/* Here, we define a monitor to produce output for debugging */

#include "monitor.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscInt step, PetscReal time, Vec x, void * ptr)
{
  PetscErrorCode ierr;
  Ctx            *ctx = (Ctx*)ptr;

  PetscFunctionBeginUser;

  {
    PetscReal minval,maxval;
    ierr = VecMin(x,NULL,&minval);CHKERRQ(ierr);
    ierr = VecMax(x,NULL,&maxval);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Custom monitor: Step %D, t=%f. Min/Max %f/%f\n",step,time,minval,maxval);CHKERRQ(ierr);
    
  }

  /* Dump the solution to a file named for the timestep */
  {
    PetscViewer viewer;
    char filename[PETSC_MAX_PATH_LEN],vecname[PETSC_MAX_PATH_LEN];


    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"output/TIMESTEPPER/S_s/S_s_%D.m",step);CHKERRQ(ierr);
    ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"S_s_step_%D",step);CHKERRQ(ierr);

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)x,vecname);CHKERRQ(ierr);
    ierr = VecView(x,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Dump the rhs stored in the ctx to a file named for the timestep */
  {
    PetscViewer viewer;
    char filename[PETSC_MAX_PATH_LEN],vecname[PETSC_MAX_PATH_LEN];
    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"output/TIMESTEPPER/rhs/rhs_%D.m",step);CHKERRQ(ierr);
    ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"rhs_step_%D",step);CHKERRQ(ierr);

    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)x,vecname);CHKERRQ(ierr);
    ierr = VecView(ctx->solution.rhs_s,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
