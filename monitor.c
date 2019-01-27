#include "dimensionalisablefield.h"
#include "monitor.h"
#include "rhs.h"
#include "cJSON.h"
#include "version.h"
#include "rheologicalfront.h"
#include "atmosphere.h"
#include "twophase.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscReal dtmacro, PetscReal dtmacro_years, PetscInt step, PetscReal time, Vec sol, void *ptr, MonitorCtx *mctx)
{
  PetscErrorCode ierr;
  Ctx            *ctx = (Ctx*)ptr;
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

  /* Ensure that the output directory exists */
  if (!mctx->outputDirectoryExistenceConfirmed) {
    PetscBool exists;
    ierr = PetscTestDirectory(P->outputDirectory,'w',&exists);CHKERRQ(ierr);
    if (!exists) {
      ierr = PetscMkdir(P->outputDirectory);CHKERRQ(ierr);
    }
    mctx->outputDirectoryExistenceConfirmed = PETSC_TRUE;
  }

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
      {
        data = cJSON_CreateObject();

        for (i=0; i<NUMMESHVECS_S; ++i) {
          cJSON *item;
          DimensionalisableField curr = ctx->mesh.meshFields_s[i];
          ierr = DimensionalisableFieldToJSON(curr,&item);CHKERRQ(ierr);
          cJSON_AddItemToObject(data,curr->name,item);
        } 
        for (i=0; i<NUMMESHVECS_B; ++i) {
          cJSON *item;
          DimensionalisableField curr = ctx->mesh.meshFields_b[i];
          ierr = DimensionalisableFieldToJSON(curr,&item);CHKERRQ(ierr);
          cJSON_AddItemToObject(data,curr->name,item);
        }
        for (i=0; i<NUMSOLUTIONVECS_S; ++i) {
          cJSON *item;
          DimensionalisableField curr = ctx->solution.solutionFields_s[i];
          ierr = DimensionalisableFieldToJSON(curr,&item);CHKERRQ(ierr);
          cJSON_AddItemToObject(data,curr->name,item);
        }
        for (i=0; i<NUMSOLUTIONVECS_B; ++i) {
          cJSON *item;
          DimensionalisableField curr = ctx->solution.solutionFields_b[i];
          ierr = DimensionalisableFieldToJSON(curr,&item);CHKERRQ(ierr);
          cJSON_AddItemToObject(data,curr->name,item);
        }

        /* now add all the data array to the output */
        cJSON_AddItemToObject(json,"data",data);
      }

      /* rheological front */
      ierr = JSON_add_rheological_front( ctx->da_point, &ctx->parameters.constants, &ctx->rheological_front_phi, "rheological_front_phi", json); CHKERRQ(ierr);

      /* atmosphere */
      ierr = JSON_add_atmosphere( ctx->da_point, &ctx->parameters, &ctx->atmosphere, "atmosphere", json); CHKERRQ(ierr);

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
