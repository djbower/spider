static char help[] ="Parallel magma ocean timestepper";

#include "petsc.h"
#include "ctx.h" 

PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);

int main(int argc, char ** argv)
{
  PetscErrorCode ierr;
  TS             ts;      /* ODE solver object */
  Vec            S_s;     /* Solution Vector */
  Ctx            ctx;     /* Solver context */

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);

  /* Perform all initialization for our problem, allocating data */
  ierr = setup_ctx(&ctx);CHKERRQ(ierr);

  /* We will use this solution vector as our data object for timestepping */
  S_s = ctx.solution.S_s;

  {
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," *** Viewing S_s initial condition ***\n");CHKERRQ(ierr);
    ierr = VecView(ctx.solution.S_s,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  
  /* Set up the Jacobian function (omitted for now) */

  /* Set up timestepper */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,S_s);CHKERRQ(ierr);
 // ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  // TODO: More CVODE setup ...

  /* Set up the RHS Function */
  ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&ctx);CHKERRQ(ierr);

  /* Set up the integration period */
  ierr = TSSetDuration(ts,1,1.e12);CHKERRQ(ierr); /* One time step, huge final time */
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);

  /* Accept command line options */
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  /* Set up some kind of monitor or viewer (omitted for now) */

  /* Solve */
  ierr = TSSolve(ts,S_s);CHKERRQ(ierr);

  /* For debugging, view some things stored in the context.
      Note that there are other viewer implementations, including 
      those that will write to a file which we should be able
      to open with python.
  */
  {
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," *** Viewing S_s ***\n");CHKERRQ(ierr);
    ierr = VecView(ctx.solution.S_s,viewer);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," *** Viewing rhs_s ***\n");CHKERRQ(ierr);
    ierr = VecView(ctx.solution.rhs_s,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Free allocated data in the context */
  ierr = destroy_ctx(&ctx);CHKERRQ(ierr);

  /* Cleanup and finalize */
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}
