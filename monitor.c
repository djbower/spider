/* Here, we define a monitor to produce output for debugging */

#include "monitor.h"
#include "output.h"
#include "rhs.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscReal dtmacro, PetscInt step, PetscReal time, PetscReal time0, PetscReal timeprev, Vec x_aug, void * ptr, double walltime0, double *walltimeprev)
{
  PetscErrorCode ierr;
  Ctx            *ctx = (Ctx*)ptr;
  PetscBool      test_view = PETSC_FALSE;
  /* it remains convenient to be able to plot both short and long times together
     on the same plot, and since these are output by different settings in main.c,
     it is easiest if the output files reference the actual time in years, rather
     than timestep.  But the consequence of this approach is that small output
     times (less than a year) will overwrite each other.  For debugging this could
     be annoying, but in terms of plotting and analysis it is a non-issue since nothing
     is happening at timescales less than a year */  
  long long      nstep = (long long) dtmacro * (long long) step;

  PetscFunctionBeginUser;

  ierr = PetscOptionsGetBool(NULL,NULL,"-test_view",&test_view,NULL);CHKERRQ(ierr);

  {
    PetscReal minval,maxval;
    ierr = VecMin(x_aug,NULL,&minval);CHKERRQ(ierr); /* Note that this includes the extra point now, so might not be as meaningful! */
    ierr = VecMax(x_aug,NULL,&maxval);CHKERRQ(ierr);
    double walltime = MPI_Wtime();
    double walltimetotal = walltime - walltime0;
    double walltimestep = walltime - *walltimeprev;
    double time_per_walltime_hour_total = 3600.0 * (time-time0) / walltimetotal;
    double time_per_walltime_hour_step = 3600.0 * (time-timeprev) / walltimestep;
    ierr = PetscPrintf(PETSC_COMM_WORLD,
        "***  Writing output at macro Step %D, t=%f. Min/Max %f/%f [%.2f s : %.2f time units / hour (%.2f total)]\n",
        step,(double)time,(double)minval,(double)maxval,walltime-walltime0,time_per_walltime_hour_step,time_per_walltime_hour_total);CHKERRQ(ierr);
    *walltimeprev = walltime;
  }

  /* Dump the solution to a file named for the timestep */
  {
    PetscViewer viewer;
    char filename[PETSC_MAX_PATH_LEN],vecname[PETSC_MAX_PATH_LEN];
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"output/dSdr_b_aug.%lld.m",nstep);CHKERRQ(ierr);
    ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"dSdr_b_aug_%lld",nstep);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr); //Annoyingly, PETSc wants you to use binary output so badly that this is the easiest way to get full-precision ASCII..
    ierr = PetscObjectSetName((PetscObject)x_aug,vecname);CHKERRQ(ierr);
    ierr = VecView(x_aug,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Dump several PETSc binary vectors to a file named for the timestep */
  {
    /* Set up a binary viewer */
    PetscViewer viewer;
    char filename[PETSC_MAX_PATH_LEN];
    Mesh     *M = &ctx->mesh;
    Solution *S = &ctx->solution;
    PetscInt i;

    /* for debugging, switch 'nstep' below to 'step' to over overwriting output data files */
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"output/petscbin.%lld",nstep);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);

    /* the following order of export must match exactly the order of the fields in plot_figure.py
       see get_vector_index in plot_figure.py */

    /* Add a data vector */
    {
      Vec data;
      PetscMPIInt rank;
      char vecname[PETSC_MAX_PATH_LEN];
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
        ierr = PetscObjectSetName((PetscObject)data,vecname);CHKERRQ(ierr);
        ierr = VecView(data,viewer);CHKERRQ(ierr);
        ierr = VecDestroy(&data);CHKERRQ(ierr);
      }
    }

    /* Add the solution vector */
    {
      char vecname[PETSC_MAX_PATH_LEN];
      ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"dSdr_b_aug_%d",step);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)x_aug,vecname);CHKERRQ(ierr);
      ierr = VecView(x_aug,viewer);CHKERRQ(ierr);
    }

    /* write mesh information for basic node quantities */
    for (i=0;i<NUMMESHVECS_B;++i){
        add_vector_to_binary_output( M->meshVecs_b[i], viewer ); CHKERRQ(ierr);
    }

    /* write mesh information for staggered node quantities */
    for (i=0;i<NUMMESHVECS_S;++i){
        add_vector_to_binary_output( M->meshVecs_s[i], viewer ); CHKERRQ(ierr);
    }

    /* write all vectors associated with basic node quantities */
    for (i=0;i<NUMSOLUTIONVECS_B;++i){
        add_vector_to_binary_output( S->solutionVecs_b[i], viewer ); CHKERRQ(ierr);
    }

    /* write all vectors associated with staggered node quantities */
    for (i=0;i<NUMSOLUTIONVECS_S;++i){
        add_vector_to_binary_output( S->solutionVecs_s[i], viewer ); CHKERRQ(ierr);
    }

    /* Close the viewer */
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  if (test_view) {
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"--- Printing dSdr_b_aug for testing ---\n",time);CHKERRQ(ierr);
    ierr = VecView(x_aug,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

#if 0
  /* Recompute the rhs to a file named for the timestep */
  {
    Vec rhs_b_aug;

    ierr = VecDuplicate(x_aug,&rhs_b_aug);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)rhs_b_aug,"rhs_b_aug");CHKERRQ(ierr);
    ierr = RHSFunction(ts,time,x_aug,rhs_b_aug,ctx);CHKERRQ(ierr);
    {
      PetscViewer viewer;
      char filename[PETSC_MAX_PATH_LEN],vecname[PETSC_MAX_PATH_LEN];

      ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"output/rhs_b.%D.m",step);CHKERRQ(ierr);
      ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"rhs_step_%D",step);CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);
      ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)x,vecname);CHKERRQ(ierr);
      ierr = VecView(rhs_b,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }

    if (test_view) {
      PetscViewer viewer;

      ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
      ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
      ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"--- Printing rhs_b_aug for testing ---\n",time);CHKERRQ(ierr);
      ierr = VecView(rhs_b_aug,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }
  }
#endif

  PetscFunctionReturn(0);
}
