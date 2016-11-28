/* Here, we define a monitor to produce output for debugging */

#include "monitor.h"
#include "rhs.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscInt step, PetscReal time, PetscReal time0, Vec x, void * ptr, double walltime0)
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
    double walltime = MPI_Wtime() - walltime0;
    double time_per_walltime_hour = 3600.0 * (time-time0) / walltime;
    ierr = PetscPrintf(PETSC_COMM_WORLD,
        "*** [%.2f s : %.2f time units / hour] Writing output at macro Step %D, t=%f. Min/Max %f/%f\n",
        walltime,time_per_walltime_hour,step,(double)time,(double)minval,(double)maxval);CHKERRQ(ierr);
  }

  /* Dump the solution to a file named for the timestep */
  {
    PetscViewer viewer;
    char filename[PETSC_MAX_PATH_LEN],vecname[PETSC_MAX_PATH_LEN];

    ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"output/S_s_%D.m",step);CHKERRQ(ierr);
    ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"S_s_step_%D",step);CHKERRQ(ierr);
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr); //Annoyingly, PETSc wants you to use binary output so badly that this is the easiest way to get full-precision ASCII..
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
  // NOTE: we turn off dumping of the RHS for now, but retain this code as it may become useful later
#if 0
    {
      PetscViewer viewer;
      char filename[PETSC_MAX_PATH_LEN],vecname[PETSC_MAX_PATH_LEN];

      ierr = PetscSNPrintf(filename,PETSC_MAX_PATH_LEN,"output/rhs_s.%D.m",step);CHKERRQ(ierr);
      ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"rhs_step_%D",step);CHKERRQ(ierr);
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);CHKERRQ(ierr);
      ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject)x,vecname);CHKERRQ(ierr);
      ierr = VecView(rhs_s,viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    }
#endif
    
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
