#include "dimensionalisablefield.h"
#include "monitor.h"
#include "rhs.h"
#include "cJSON.h"
#include "version.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscReal dtmacro, PetscReal dtmacro_years, PetscInt step, PetscReal time, Vec sol, void *ptr, MonitorCtx *mctx)
{
  PetscErrorCode ierr;
  Ctx            *ctx = (Ctx*)ptr;
  PetscBool      test_view = PETSC_FALSE;
  Vec rhs;
  /* it remains convenient to be able to plot both short and long times together
     on the same plot, and since these are output by different settings in main.c,
     it is easiest if the output files reference the actual time in years, rather
     than timestep.  But the consequence of this approach is that small output
     times (less than a year) will overwrite each other.  For debugging this could
     be annoying, but in terms of plotting and analysis it is a non-issue since nothing
     is happening at timescales less than a year */
  Parameters const *P = &ctx->parameters;
  long long        nstep = (long long) dtmacro_years * (long long )step;

  PetscFunctionBeginUser;

  ierr = PetscOptionsGetBool(NULL,NULL,"-test_view",&test_view,NULL);CHKERRQ(ierr);

  /* Info for stdout */
  {
    double        walltime;
    long long int elapsedSeconds;
    int           days,hours,minutes,seconds;
    Vec           *subVecs;
    PetscScalar   *maxVals,*minVals;
    PetscInt      i;

    /* Compute min/max for each field */
    ierr = PetscMalloc3(ctx->numFields,&subVecs,ctx->numFields,&maxVals,ctx->numFields,&minVals);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(ctx->dm_sol,sol,ctx->numFields,NULL,subVecs);CHKERRQ(ierr);

    for (i=0; i<ctx->numFields; ++i) {
      ierr = VecMin(subVecs[i],NULL,&minVals[i]);CHKERRQ(ierr);
      ierr = VecMax(subVecs[i],NULL,&maxVals[i]);CHKERRQ(ierr);
    }
    ierr = DMCompositeRestoreAccessArray(ctx->dm_sol,sol,ctx->numFields,NULL,subVecs);CHKERRQ(ierr);

    /* Compute timing information */
    walltime = MPI_Wtime();
    mctx->walltimeprev = walltime;
    elapsedSeconds = (long long int) (walltime - mctx->walltime0);
    seconds = elapsedSeconds % 60;
    minutes = ((elapsedSeconds - seconds)/60) % 60;
    hours = ((elapsedSeconds - 60*minutes - seconds)/(60*60)) % 24;
    days = (elapsedSeconds - 24*60*hours - 60*minutes - seconds)/(24*60*60);

    /* Print */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"***  Writing output at macro Step %D, t=%f. [%lld:%02d:%02d:%02d]\n",step,(double)time,days,hours,minutes,seconds);CHKERRQ(ierr);
    for (i=0; i<ctx->numFields; ++i) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Field %D Min/Max %f/%f\n",i,(double)minVals[i],(double)maxVals[i]);CHKERRQ(ierr);
    }

    /* Cleanup */
    ierr = PetscFree3(subVecs,maxVals,minVals);CHKERRQ(ierr);
  }

  /* Reevaluate the RHS */
  ierr = VecDuplicate(sol,&rhs);CHKERRQ(ierr);
  ierr = RHSFunction(ts,time,sol,rhs,ctx);CHKERRQ(ierr);

  /* Ensure that the output directory exists */
  if (!mctx->outputDirectoryExistenceConfirmed) {
    PetscBool exists;
    ierr = PetscTestDirectory(P->outputDirectory,'w',&exists);CHKERRQ(ierr);
    if (!exists) {
      ierr = PetscMkdir(P->outputDirectory);CHKERRQ(ierr);
    }
    mctx->outputDirectoryExistenceConfirmed = PETSC_TRUE;
  }

  /* Dump the solution to a file named for the timestep */
  {
    PetscViewer viewer;
    char filename[PETSC_MAX_PATH_LEN],vecname[PETSC_MAX_PATH_LEN];
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/sol_%lld.m",P->outputDirectory,nstep);CHKERRQ(ierr);
    ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"sol",nstep);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr); //Annoyingly, PETSc wants you to use binary output so badly that this is the easiest way to get full-precision ASCII..
    ierr = PetscObjectSetName((PetscObject)sol,vecname);CHKERRQ(ierr);
    ierr = VecView(sol,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  if (test_view) {
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"--- Printing sol for testing ---\n",time);CHKERRQ(ierr);
    ierr = VecView(sol,viewer);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"--- Printing rhs for testing ---\n",time);CHKERRQ(ierr);
    ierr = VecView(rhs,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  ierr = VecDestroy(&rhs);CHKERRQ(ierr);

  // TODO: once JSON output stable, remove the above old output (leaving the stdout output!)
  /* Dump a JSON file for this timestep */
  {
    PetscMPIInt rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

    if (rank==0) {
      cJSON       *json,*data;
      char        *outputString;
      char        filename[PETSC_MAX_PATH_LEN];
      PetscInt    i;

      ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/%lld.json",P->outputDirectory,nstep);CHKERRQ(ierr);

      /* Create a new JSON object */
      json = cJSON_CreateObject();

      /* Add some header info */
      cJSON_AddItemToObject(json,"SPIDER major version",cJSON_CreateString(SPIDER_MAJOR_VERSION));
      cJSON_AddItemToObject(json,"SPIDER minor version",cJSON_CreateString(SPIDER_MINOR_VERSION));
      cJSON_AddItemToObject(json,"SPIDER patch version",cJSON_CreateString(SPIDER_PATCH_VERSION));
      // TODO: some of these variable names lack clear meaning (e.g., nstep)
      cJSON_AddItemToObject(json,"step",cJSON_CreateNumber(step));
      cJSON_AddItemToObject(json,"dtmacro_years",cJSON_CreateNumber(dtmacro_years));
      cJSON_AddItemToObject(json,"nstep",cJSON_CreateNumber(nstep));
      // TODO : add other stuff we might like in the header

      /* Add solution to file */
      {
        cJSON *solJSON;
        ierr = DimensionalisableFieldToJSON(ctx->solDF,&solJSON);CHKERRQ(ierr);
        cJSON_AddItemToObject(json,"solution",solJSON);
      }

      /* Re-evaluate rhs and add to file */
      {
        DimensionalisableField rhsDF;
        Vec                    rhs;
        cJSON                  *rhsJSON;
        ierr = DimensionalisableFieldDuplicate(ctx->solDF,&rhsDF);CHKERRQ(ierr);
        ierr = DimensionalisableFieldSetName(rhsDF,"SPIDER rhs");CHKERRQ(ierr);
        ierr = DimensionalisableFieldGetGlobalVec(rhsDF,&rhs);CHKERRQ(ierr);
        ierr = RHSFunction(ts,time,sol,rhs,ctx);CHKERRQ(ierr);
        ierr = DimensionalisableFieldToJSON(rhsDF,&rhsJSON);CHKERRQ(ierr);
        cJSON_AddItemToObject(json,"rhs",rhsJSON);CHKERRQ(ierr);
        ierr = DimensionalisableFieldDestroy(&rhsDF);CHKERRQ(ierr);
      }

      /* Add data of interest */
      // Note: this is duplicative, but it's too painful to flatten out the Ctx
      data = cJSON_CreateArray();

      for (i=0; i<NUMMESHVECS_S; ++i) {
        cJSON *item;
        DimensionalisableField curr = ctx->mesh.meshFields_s[i];
        ierr = DimensionalisableFieldToJSON(curr,&item);CHKERRQ(ierr);
        cJSON_AddItemToArray(data,item);
      }
      for (i=0; i<NUMMESHVECS_B; ++i) {
        cJSON *item;
        DimensionalisableField curr = ctx->mesh.meshFields_b[i];
        ierr = DimensionalisableFieldToJSON(curr,&item);CHKERRQ(ierr);
        cJSON_AddItemToArray(data,item);
      }
      for (i=0; i<NUMSOLUTIONVECS_S; ++i) {
        cJSON *item;
        DimensionalisableField curr = ctx->solution.solutionFields_s[i];
        ierr = DimensionalisableFieldToJSON(curr,&item);CHKERRQ(ierr);
        cJSON_AddItemToArray(data,item);
      }
      for (i=0; i<NUMSOLUTIONVECS_B; ++i) {
        cJSON *item;
        DimensionalisableField curr = ctx->solution.solutionFields_b[i];
        ierr = DimensionalisableFieldToJSON(curr,&item);CHKERRQ(ierr);
        cJSON_AddItemToArray(data,item);
      }

      /* atmosphere */
      {
        Atmosphere const *A = &ctx->atmosphere;
        Parameters           const *P = &ctx->parameters;
        Constants            const *C  = &P->constants;
        AtmosphereParameters const *Ap = &P->atmosphere_parameters;
        VolatileParameters   const *CO2 = &Ap->CO2_volatile_parameters;
        VolatileParameters   const *H2O = &Ap->H2O_volatile_parameters;
        Mesh                 const *M = &ctx->mesh;
        const PetscInt             ind0 = 0;
        PetscScalar                x0,x1;
        Vec                        *subVecs;

        ierr = PetscMalloc1(ctx->numFields,&subVecs);CHKERRQ(ierr);
        ierr = DMCompositeGetAccessArray(ctx->dm_sol,sol,ctx->numFields,NULL,subVecs);CHKERRQ(ierr);

        /* CO2 content of magma ocean (liquid phase) */
        ierr = VecGetValues(subVecs[ctx->solutionSlots[SPIDER_SOLUTION_FIELD_MO_CO2]],1,&ind0,&x0);CHKERRQ(ierr);

        /* H2O content of magma ocean (liquid phase) */
        ierr = VecGetValues(subVecs[ctx->solutionSlots[SPIDER_SOLUTION_FIELD_MO_H2O]],1,&ind0,&x1);CHKERRQ(ierr);

        ierr = DMCompositeRestoreAccessArray(ctx->dm_sol,sol,ctx->numFields,NULL,subVecs);CHKERRQ(ierr);
        ierr = PetscFree(subVecs);CHKERRQ(ierr);

        /* TODO: scalings are incorporated into the outputs below, but
           keeping these scalings to document at some point */
        //MASS = 4.0 * PETSC_PI * C->MASS; // includes 4*PI for spherical geometry
        //FAC = C->VOLATILE / 1.0E6;
        //Msol = A->Msol * MASS;
        //Mliq = A->Mliq * MASS;
        // CO2
        //sol0 = FAC * x0 * CO2->kdist * Msol; // solid
        //liq0 = FAC * x0 * Mliq; // liquid
        //atm0 = A->m0 * MASS; // atmosphere
        //tot0 = FAC * CO2->initial * M->mantle_mass * MASS; // total
        // H2O
        //sol1 = FAC * x1 * H2O->kdist * Msol; // solid
        //liq1 = FAC * x1 * Mliq; // liquid
        //atm1 = A->m1 * MASS; // atmosphere
        //tot1 = FAC * H2O->initial * M->mantle_mass * MASS; // total

        /* total liquid mass of mantle, kg */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = 4.0 * PETSC_PI * C->MASS; // includes 4*PI for spherical geometry
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"mass_liquid");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->Mliq,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* total solid mass of mantle, kg */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = 4.0 * PETSC_PI * C->MASS; // includes 4*PI for spherical geometry
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"mass_solid");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->Msol,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* total mass of mantle, kg (for sanity check) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = 4.0 * PETSC_PI * C->MASS; // includes 4*PI for spherical geometry
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"mass_mantle");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,M->mantle_mass,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* surface temperature, K */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = C->TEMP;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"temperature_surface");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"K");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->tsurf,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* optical depth, non-dimensional */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = 1.0;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"optical_depth");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"None");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->tau,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* (effective) emissivity, non-dimensional */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = 1.0;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"emissivity");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"None");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->emissivity,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* initial CO2 volatile (ppm) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = C->VOLATILE;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"CO2_initial_ppm");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"ppm");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,CO2->initial,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* initial CO2 volatile (kg) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = (C->VOLATILE/1.0E6) * 4.0 * PETSC_PI * C->MASS;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"CO2_initial_kg");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,CO2->initial*M->mantle_mass,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* CO2 volatile in liquid mantle (ppm) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = C->VOLATILE;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"CO2_liquid_ppm");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"ppm");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,x0,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* CO2 volatile in liquid mantle (kg) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = (C->VOLATILE/1.0E6) * 4.0 * PETSC_PI * C->MASS;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"CO2_liquid_kg");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,x0*A->Mliq,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* CO2 volatile in solid mantle (ppm) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = C->VOLATILE;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"CO2_solid");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"ppm");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,x0*CO2->kdist,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* CO2 volatile in solid mantle (kg) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = (C->VOLATILE/1.0E6) * 4.0 * PETSC_PI * C->MASS;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"CO2_solid_kg");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,x0*CO2->kdist*A->Msol,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* CO2 volatile in atmosphere (bar) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = C->PRESSURE / 1.0E5; /* for bars */
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"CO2_atmosphere_bar");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"bar");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->p0,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* CO2 volatile in atmosphere (kg) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = (C->VOLATILE/1.0E6) * 4.0 * PETSC_PI * C->MASS;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"CO2_atmosphere_kg");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->m0,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* CO2 optical depth (non-dimensional) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = 1.0;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"CO2_optical_depth");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"None");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->tau0,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* initial H2O volatile (ppm) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = C->VOLATILE;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"H2O_initial_ppm");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"ppm");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,H2O->initial,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* initial H2O volatile (kg) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = (C->VOLATILE/1.0E6) * 4.0 * PETSC_PI * C->MASS;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"H2O_initial_kg");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,H2O->initial*M->mantle_mass,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* H2O volatile in liquid mantle (ppm) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = C->VOLATILE;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"H2O_liquid_ppm");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"ppm");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,x1,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* H2O volatile in liquid mantle (kg) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = (C->VOLATILE/1.0E6) * 4.0 * PETSC_PI * C->MASS;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"H2O_liquid_kg");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,x1*A->Mliq,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* H2O volatile in solid mantle (ppm) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = C->VOLATILE;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"H2O_solid_ppm");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"ppm");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,x1*H2O->kdist,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* H2O volatile in solid mantle (kg) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = (C->VOLATILE/1.0E6) * 4.0 * PETSC_PI * C->MASS;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"H2O_solid_kg");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,x1*H2O->kdist*A->Msol,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* H2O volatile in atmosphere (bar) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = C->PRESSURE / 1.0E5; /* for bar */
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"H2O_atmosphere_bar");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"bar");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->p1,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* H2O volatile in atmosphere (kg) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = (C->VOLATILE/1.0E6) * 4.0 * PETSC_PI * C->MASS;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"H2O_atmosphere_kg");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->m1,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

        /* H2O optical depth (non-dimensional) */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = 1.0;
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"H2O_optical_depth");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"None");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,A->tau1,INSERT_VALUES);CHKERRQ(ierr);
          ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
          ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
          cJSON_AddItemToArray(data,item);
          ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);
        }

      }

      /* now add all the data array to the output */
      cJSON_AddItemToObject(json,"data",data);

      /* Print to a string */
      outputString = cJSON_Print(json);

      /* Print to file (the dumb way since PETSc's ASCII viewer seems buggy) */
      {
        FILE *fp;
        fp = fopen(filename,"w");
        fprintf(fp,"%s\n",outputString);
        fclose(fp);
      }

      /* Free the string */
      free(outputString);

      /* Delete the JSON object and viewer (close file)*/
      cJSON_Delete(json);
    }
  }

  /* At this point, we should have processed all command line options, so we
     check to see if any are stray (this usually would indicate a typo) */
  ierr = PetscOptionsLeft(NULL);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode TSMonitorWalltimed(TS ts,PetscInt steps,PetscReal time,Vec x,void *mctx)
{
  PetscErrorCode ierr;
  MonitorCtx     *ctx = mctx;
  double         walltime;
  const double   split = 10.0; /* seconds between outputs */
  long long int  elapsedSeconds;
  int            hours,days,minutes,seconds;
  PetscReal      dt;

  PetscFunctionBeginUser;
  walltime = MPI_Wtime();
  if (walltime - ctx->walltimeprev > split){
    ctx->walltimeprev = walltime;
    elapsedSeconds = (int) (walltime - ctx->walltime0);
    seconds = elapsedSeconds % 60;
    minutes = ((elapsedSeconds - seconds)/60) % 60;
    hours = ((elapsedSeconds - 60*minutes - seconds)/(60*60)) % 24;
    days = (elapsedSeconds - 24*60*hours - 60*minutes - seconds)/(24*60*60);
    ierr = TSGetTimeStep(ts,&dt);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  micro step %D, t=%.3f, dt=%.3g [%d:%02d:%02d:%02d]\n",steps,(double)time,(double)dt,days,hours,minutes,seconds,elapsedSeconds);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
