/* Here, we define a monitor to produce output for debugging */

#include "monitor.h"
#include "rhs.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscInt step, PetscReal time, Vec x, void * ptr)
{
  PetscErrorCode ierr;
  Ctx            *ctx = (Ctx*)ptr;
  PetscBool      test_view = PETSC_FALSE;

  PetscFunctionBeginUser;

  ierr = PetscOptionsGetBool(NULL,NULL,"-test_view",&test_view,NULL);CHKERRQ(ierr);

  {
    PetscReal minval,maxval;
    ierr = VecMin(x,NULL,&minval);CHKERRQ(ierr);
    ierr = VecMax(x,NULL,&maxval);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Writing output at macro Step %D, t=%f. Min/Max %f/%f\n",step,time,minval,maxval);CHKERRQ(ierr);
  }

  // TODO PDS: discuss with Dan what output format is actually best here. Currently it dumps out MATLAB-style files for S_s and rhs_s at each timestep

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

  if (test_view) {
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"--- Printing S_s for testing ---\n",time);CHKERRQ(ierr);
    ierr = VecView(x,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Recompute the rhs stored in the ctx to a file named for the timestep */
  {
    Vec rhs_s; 

    ierr = VecDuplicate(x,&rhs_s);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)rhs_s,"rhs_s");CHKERRQ(ierr);
    ierr = RHSFunction(ts,time,x,rhs_s,ctx);CHKERRQ(ierr);
    {
      PetscViewer viewer;
      char filename[PETSC_MAX_PATH_LEN],vecname[PETSC_MAX_PATH_LEN];

      ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"output/TIMESTEPPER/rhs/rhs_%D.m",step);CHKERRQ(ierr);
      ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"rhs_step_%D",step);CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);
      ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)x,vecname);CHKERRQ(ierr);
      ierr = VecView(rhs_s,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }
    
    if (test_view) {
      PetscViewer viewer;

      ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
      ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
      ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"--- Printing rhs_s for testing ---\n",time);CHKERRQ(ierr);
      ierr = VecView(rhs_s,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }
    ierr = VecDestroy(&rhs_s);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
