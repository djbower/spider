#include "dimensionalisablefield.h"
#include "monitor.h"
#include "output.h"
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

  /* Globally, output double, if we are computing with quad precision */
#if (defined PETSC_USE_REAL___FLOAT128)
  ierr = PetscOptionsSetValue(NULL,"-binary_write_double","");CHKERRQ(ierr);
#endif

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
  {

  ierr = VecDuplicate(sol,&rhs);CHKERRQ(ierr);
  ierr = RHSFunction(ts,time,sol,rhs,ctx);CHKERRQ(ierr);
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

/* below can be removed */
#if 0
  /* Dump several PETSc binary vectors to a file named for the timestep */
  {
    /* Set up a binary viewer */
    PetscViewer          viewer;
    char                 filename[PETSC_MAX_PATH_LEN];
    Parameters           *P = &ctx->parameters;
    Constants            *C = &P->constants;
    Mesh                 *M = &ctx->mesh;
    Solution             *S = &ctx->solution;

    /* for debugging, switch 'nstep' below to 'step' to over overwriting output data files */
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/petscbin.%lld",P->outputDirectory,nstep);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);

    /* the following order of export must match exactly the order of the fields in plot_figure.py
       see get_vector_index in plot_figure.py */

    /* Add a data vector */
    {
      Vec data;
      PetscMPIInt rank;
      //char vecname[PETSC_MAX_PATH_LEN];
      const int nData = 2;
      ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      if (!rank) {
        ierr = VecCreate(PETSC_COMM_SELF,&data);CHKERRQ(ierr);
        ierr = VecSetType(data,VECSEQ);CHKERRQ(ierr);
        ierr = VecSetSizes(data,nData,nData);CHKERRQ(ierr);
        ierr = VecSetValue(data,0,(PetscScalar)step,INSERT_VALUES);CHKERRQ(ierr);
        ierr = VecSetValue(data,1,(PetscScalar)dtmacro,INSERT_VALUES);CHKERRQ(ierr);
        ierr = VecAssemblyBegin(data);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(data);CHKERRQ(ierr);
        //ierr = PetscObjectSetName((PetscObject)data,vecname);CHKERRQ(ierr);
        ierr = VecView(data,viewer);CHKERRQ(ierr);
        ierr = VecDestroy(&data);CHKERRQ(ierr);
      }
    }

    /* Add constants */
    {
      Vec data;
      PetscMPIInt rank;
      //char vecname[PETSC_MAX_PATH_LEN];
      /* FIXME: it's really easy to update output.c and forget to change the size
         of the array in the next line */
      const int nData = 28;
      ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      if (!rank) {
        ierr = VecCreate(PETSC_COMM_SELF,&data);CHKERRQ(ierr);
        ierr = VecSetType(data,VECSEQ);CHKERRQ(ierr);
        ierr = VecSetSizes(data,nData,nData);CHKERRQ(ierr);
        ierr = constants_struct_to_vec( C, data ); CHKERRQ(ierr);
        //ierr = PetscObjectSetName((PetscObject)data,vecname);CHKERRQ(ierr);
        ierr = VecView(data,viewer);CHKERRQ(ierr);
        ierr = VecDestroy(&data);CHKERRQ(ierr);
      }
    }

    /* Add atmosphere */
    {
      Vec data;
      PetscMPIInt rank;
      //char vecname[PETSC_MAX_PATH_LEN];
      /* FIXME: it's really easy to update output.c and forget to change the size
         of the array in the next line */
      const int nData = 20; //34;
      ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      if (!rank) {
        ierr = VecCreate(PETSC_COMM_SELF,&data);CHKERRQ(ierr);
        ierr = VecSetType(data,VECSEQ);CHKERRQ(ierr);
        ierr = VecSetSizes(data,nData,nData);CHKERRQ(ierr);
        ierr = atmosphere_structs_to_vec( sol, ctx, data ); CHKERRQ(ierr);
        //ierr = PetscObjectSetName((PetscObject)data,vecname);CHKERRQ(ierr);
        ierr = VecView(data,viewer);CHKERRQ(ierr);
        ierr = VecDestroy(&data);CHKERRQ(ierr);
      }
    }

    /* Add the solution vector */
    {
      char vecname[PETSC_MAX_PATH_LEN];
      ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"sol_%d",step);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)sol,vecname);CHKERRQ(ierr);
      ierr = VecView(sol,viewer);CHKERRQ(ierr);
    }

    /* write mesh information for basic node quantities */
    ierr = scale_vectors_and_output( M->meshFields_b, NUMMESHVECS_B, viewer );CHKERRQ(ierr);

    /* write mesh information for staggered node quantities */
    ierr = scale_vectors_and_output( M->meshFields_s, NUMMESHVECS_S, viewer );CHKERRQ(ierr);

    /* write all vectors associated with basic node quantities */
    ierr = scale_vectors_and_output( S->solutionFields_b, NUMSOLUTIONVECS_B, viewer );CHKERRQ(ierr);

    /* write all vectors associated with staggered node quantities */
    ierr = scale_vectors_and_output( S->solutionFields_s, NUMSOLUTIONVECS_S, viewer );CHKERRQ(ierr);

    /* Close the viewer */
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

#endif

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
        PetscScalar                Msol,Mliq;
        PetscScalar                sol0,liq0,atm0,tot0,sol1,liq1,atm1,tot1;
        PetscScalar                x0,x1;
        PetscScalar                FAC, MASS;
        Vec                        *subVecs;

        ierr = PetscMalloc1(ctx->numFields,&subVecs);CHKERRQ(ierr);
        ierr = DMCompositeGetAccessArray(ctx->dm_sol,sol,ctx->numFields,NULL,subVecs);CHKERRQ(ierr);

        /* CO2 content of magma ocean (liquid phase) */
        ierr = VecGetValues(subVecs[ctx->solutionSlots[SPIDER_SOLUTION_FIELD_MO_CO2]],1,&ind0,&x0);CHKERRQ(ierr);

        /* H2O content of magma ocean (liquid phase) */
        ierr = VecGetValues(subVecs[ctx->solutionSlots[SPIDER_SOLUTION_FIELD_MO_H2O]],1,&ind0,&x1);CHKERRQ(ierr);

        ierr = DMCompositeRestoreAccessArray(ctx->dm_sol,sol,ctx->numFields,NULL,subVecs);CHKERRQ(ierr);
        ierr = PetscFree(subVecs);CHKERRQ(ierr);

        /* scalings */
        MASS = 4.0 * PETSC_PI * C->MASS; // includes 4*PI for spherical geometry
        FAC = C->VOLATILE / 1.0E6;

        // TODO: this was previous, might break output below if commented!
        //Msol = A->Msol * MASS;
        //Mliq = A->Mliq * MASS;

        // CO2
        sol0 = FAC * x0 * CO2->kdist * Msol; // solid
        liq0 = FAC * x0 * Mliq; // liquid
        atm0 = A->m0 * MASS; // atmosphere
        tot0 = FAC * CO2->initial * M->mantle_mass * MASS; // total

        // H2O
        sol1 = FAC * x1 * H2O->kdist * Msol; // solid
        liq1 = FAC * x1 * Mliq; // liquid
        atm1 = A->m1 * MASS; // atmosphere
        tot1 = FAC * H2O->initial * M->mantle_mass * MASS; // total

        /* PDS: the next block is what I will duplicate for all other atmosphere outputs
           hence wanting to get it right before I duplicate! */
        {
          cJSON *item;
          DimensionalisableField dfield;
          PetscScalar scaling = 4.0 * PETSC_PI * C->MASS; // includes 4*PI for spherical geometry
          ierr = DimensionalisableFieldCreate(&dfield,ctx->da_point,&scaling,PETSC_FALSE);CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetName(dfield,"Mliq");CHKERRQ(ierr);
          ierr = DimensionalisableFieldSetUnits(dfield,"kg");CHKERRQ(ierr);
          ierr = VecSetValue(dfield->vecGlobal,0,Mliq,INSERT_VALUES);CHKERRQ(ierr);
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
