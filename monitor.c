/* Here, we define a monitor to produce output for debugging */

#include "monitor.h"
#include "rhs.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscInt step, PetscReal time, PetscReal time0, PetscReal timeprev, Vec x_aug, void * ptr, double walltime0, double *walltimeprev)
{
  PetscErrorCode ierr;
  Ctx            *ctx = (Ctx*)ptr;
  PetscBool      test_view = PETSC_FALSE;

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

    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"output/dSdr_b_aug_%D.m",step);CHKERRQ(ierr);
    ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"dSdr_b_aug_step_%D",step);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr); //Annoyingly, PETSc wants you to use binary output so badly that this is the easiest way to get full-precision ASCII..
    ierr = PetscObjectSetName((PetscObject)x_aug,vecname);CHKERRQ(ierr);
    ierr = VecView(x_aug,viewer);CHKERRQ(ierr);
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

  /* Recompute the rhs to a file named for the timestep */
  {
    Vec rhs_b_aug; 

    ierr = VecDuplicate(x_aug,&rhs_b_aug);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)rhs_b_aug,"rhs_b_aug");CHKERRQ(ierr);
    ierr = RHSFunction(ts,time,x_aug,rhs_b_aug,ctx);CHKERRQ(ierr);
  // NOTE: we turn off dumping of the RHS for now, but retain this code as it may become useful later
#if 0
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
#endif
    
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

  PetscFunctionReturn(0);
}
