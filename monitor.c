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
  ierr = VecDuplicate(sol,&rhs);CHKERRQ(ierr);
  ierr = RHSFunction(ts,time,sol,rhs,ctx);CHKERRQ(ierr);

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

    // PDS TODO: replace all this with calls to write individual output files (dimensionalized) with DimensionalisableField
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
      cJSON       *json,*str,*number,*data;
      char        *outputString;
      PetscViewer viewer;
      char        filename[PETSC_MAX_PATH_LEN];

      ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"%s/%lld.json",P->outputDirectory,nstep);CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);

      /* Create a new JSON object */
      json = cJSON_CreateObject();

      /* Add some header info */
      str = cJSON_CreateString(SPIDER_MAJOR_VERSION);
      cJSON_AddItemToObject(json,"SPIDER major version",str);
      str = cJSON_CreateString(SPIDER_MINOR_VERSION);
      cJSON_AddItemToObject(json,"SPIDER minor version",str);
      str = cJSON_CreateString(SPIDER_PATCH_VERSION);
      cJSON_AddItemToObject(json,"SPIDER patch version",str);
      number = cJSON_CreateNumber(nstep);
      cJSON_AddItemToObject(json,"step",number);
      // TODO : add other stuff we might like in the header

      /* Add data of interest */
      data = cJSON_CreateArray();
      {
        cJSON *item;
        ierr = DimensionalisableFieldToJSON(ctx->solution.solutionFields_b[0],&item);CHKERRQ(ierr);
        cJSON_AddItemToArray(data,item);
      }
      // TODO : add all desired data

      cJSON_AddItemToObject(json,"data",data);

      /* Print to a string */
      outputString = cJSON_Print(json);

      // DEBUG (TODO remove)
      ierr = PetscPrintf(PETSC_COMM_WORLD,outputString);CHKERRQ(ierr);

      /* Print to file */
      ierr = PetscViewerASCIIPrintf(viewer,outputString);CHKERRQ(ierr);

      /* Delete the JSON object and viewer (close file)*/
      cJSON_Delete(json);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
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
