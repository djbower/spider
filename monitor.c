#include "dimensionalisablefield.h"
#include "eos_output.h"
#include "eos_composite.h"
#include "monitor.h"
#include "rhs.h"
#include "cJSON.h"
#include "version.h"
#include "rheologicalfront.h"
#include "atmosphere.h"
#include "twophase.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscInt step, PetscReal time, Vec sol, void *ptr, MonitorCtx *mctx)
{
  PetscErrorCode ierr;
  Ctx            *ctx = (Ctx*)ptr;
  Parameters const P = ctx->parameters;
  ScalingConstants  const SC = P->scaling_constants;

  PetscFunctionBeginUser;

  /* Info for stdout */
  {
    double        walltime;
    long long int elapsedSeconds, time_years_actual, time_years_desired;
    int           days,hours,minutes,seconds;
    Vec           *subVecs;
    PetscBool     *fieldActive;
    PetscScalar   *maxVals,*minVals;
    PetscInt      i,len;

    /* Compute min/max for each field */
    ierr = PetscMalloc4(ctx->numFields,&fieldActive,ctx->numFields,&subVecs,ctx->numFields,&maxVals,ctx->numFields,&minVals);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(ctx->dm_sol,sol,ctx->numFields,NULL,subVecs);CHKERRQ(ierr);

    for (i=0; i<ctx->numFields; ++i) {
      ierr = VecGetSize(subVecs[i],&len);CHKERRQ(ierr);
      fieldActive[i] = len > 0;
      if (fieldActive[i]) {
        ierr = VecMin(subVecs[i],NULL,&minVals[i]);CHKERRQ(ierr);
        ierr = VecMax(subVecs[i],NULL,&maxVals[i]);CHKERRQ(ierr);
      }
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

    /* actual and desired time in integer number of years */
    time_years_actual = llround( time * SC->TIMEYRS );
    time_years_desired = llround( (P->t0 + ((step-P->stepmacro) * P->dtmacro)) * SC->TIMEYRS );

    /* Print */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"***  Writing output at macro Step %D, t=%f. [%lld:%02d:%02d:%02d]\n",step,(double)time,days,hours,minutes,seconds);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"     Actual time: %lld years, Desired time: %lld years\n",time_years_actual,time_years_desired);CHKERRQ(ierr);
    for (i=0; i<ctx->numFields; ++i) {
      if (fieldActive[i]) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Field %D Min/Max %f/%f\n",i,(double)minVals[i],(double)maxVals[i]);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Field %D Min/Max N/A\n",i);CHKERRQ(ierr);
      }
    }

    /* Cleanup */
    ierr = PetscFree4(fieldActive,subVecs,maxVals,minVals);CHKERRQ(ierr);
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
      cJSON         *json,*data;
      char          *outputString;
      char          filename[PETSC_MAX_PATH_LEN];
      long long int time_years_actual;
      PetscInt      i;

      /* current age in integer number of years */
      time_years_actual = llround( time * SC->TIMEYRS );

      ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/%lld.json",P->outputDirectory,time_years_actual);CHKERRQ(ierr);

      /* Create a new JSON object */
      json = cJSON_CreateObject();

      /* Add some header info */
      cJSON_AddItemToObject(json,"SPIDER major version",cJSON_CreateString(SPIDER_MAJOR_VERSION));
      cJSON_AddItemToObject(json,"SPIDER minor version",cJSON_CreateString(SPIDER_MINOR_VERSION));
      cJSON_AddItemToObject(json,"SPIDER patch version",cJSON_CreateString(SPIDER_PATCH_VERSION));
      cJSON_AddItemToObject(json,"step",cJSON_CreateNumber(step));
      cJSON_AddItemToObject(json,"dtmacro",cJSON_CreateNumber(P->dtmacro));
      cJSON_AddItemToObject(json,"dtmacro_years",cJSON_CreateNumber(P->dtmacro*SC->TIMEYRS));
      cJSON_AddItemToObject(json,"time",cJSON_CreateNumber(time));
      cJSON_AddItemToObject(json,"time_years",cJSON_CreateNumber(time*SC->TIMEYRS));
      // add other stuff we might like in the header

      /* Add solution to file */
      {
        cJSON *solJSON;
        ierr = DimensionalisableFieldToJSON(ctx->solDF,&solJSON);CHKERRQ(ierr);
        cJSON_AddItemToObject(json,"solution",solJSON);
      }

      /* Re-evaluate rhs and add to file */
      {
        Vec                    rhs;
        cJSON                  *rhsJSON;
        ierr = DimensionalisableFieldGetGlobalVec(ctx->rhsDF,&rhs);CHKERRQ(ierr);
        ierr = RHSFunction(ts,time,sol,rhs,ctx);CHKERRQ(ierr);
        ierr = DimensionalisableFieldToJSON(ctx->rhsDF,&rhsJSON);CHKERRQ(ierr);
        cJSON_AddItemToObject(json,"rhs",rhsJSON);CHKERRQ(ierr);
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

        /* we compute the phase boundaries here for the basic nodes, but still store them
           in the solution container */
        if( ctx->parameters->n_phases > 1 ){
          /* phase boundary */
          EOS      *sub_eos;
          PetscInt should_be_two;

          ierr = EOSCompositeGetSubEOS(ctx->parameters->eos, &sub_eos, &should_be_two);CHKERRQ(ierr);
          if (should_be_two!=2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Expecting two sub-EOSs");
          ierr = JSON_add_phase_boundary( ctx, sub_eos[0], "liquidus", data ); CHKERRQ(ierr);
          ierr = JSON_add_phase_boundary( ctx, sub_eos[1], "solidus", data ); CHKERRQ(ierr);
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

      /* rheological front phi */
      ierr = JSON_add_rheological_front( ctx->da_point, ctx->parameters->scaling_constants, &ctx->rheological_front_phi, "rheological_front_phi", json); CHKERRQ(ierr);

      /* rheological front dynamic */
      ierr = JSON_add_rheological_front( ctx->da_point, ctx->parameters->scaling_constants, &ctx->rheological_front_dynamic, "rheological_front_dynamic", json); CHKERRQ(ierr);

      /* atmosphere */
      ierr = JSON_add_atmosphere( ctx->da_point, ctx->parameters, &ctx->atmosphere, "atmosphere", json); CHKERRQ(ierr);

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
  (void) x; // mark explicitly as unused
  (void) mctx; // mark explicitly as unused
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

/* A custom, verbose SNES monitor */
PetscErrorCode SNESMonitorVerbose(SNES snes, PetscInt its, PetscReal norm, void *mctx)
{
  PetscErrorCode ierr;
  Vec            x,r;
  PetscErrorCode (*func)(SNES,Vec,Vec,void*);
  void           *ctx;

  PetscFunctionBeginUser;
  (void) mctx; // mark explicitly as unused
  ierr = SNESGetSolution(snes,&x);CHKERRQ(ierr);
  ierr = SNESGetFunction(snes,&r,&func,&ctx);CHKERRQ(ierr); /* ctx should be the same as mctx */
  ierr = func(snes,x,r,ctx);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\e[36m### Iteration %D\e[0m\n",its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\e[32mresidual function norm: %g\nx:\e[0m\n",(double)norm);CHKERRQ(ierr);
  ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\e[32mresidual function:\e[0m\n");CHKERRQ(ierr);
  ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
