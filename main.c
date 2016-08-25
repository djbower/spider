static char help[] ="Parallel magma ocean timestepper\n\
                     -n : specify the number of staggered points\n\
                     -sinit : specify an entropy value to base the initial condition upon\n\
                     -monitor : use custom monitor to dump output\n\
                     -nstepsmacro : specify the number of macro (output dump) steps\n";

/* Note: if you would like more verbose output, see the preprocessor define
         in global_defs.h */

#include "petsc.h"
#include "ic.h" 
#include "ctx.h"
#include "monitor.h"
#include "rhs.h"

int main(int argc, char ** argv)
{
  PetscErrorCode ierr;
  TS             ts;                   /* ODE solver object */
  Vec            S_s;                  /* Solution Vector */
  Ctx            ctx;                  /* Solver context */
  PetscBool      test_view;            /* View vectors for testing purposes */
  PetscBool      monitor;              /* Macros step custom monitor (monitor.c) */
  const PetscReal t0 = 0;              /* Initial time */
  const PetscInt maxsteps = 1000000000;/* Max internal steps */
  PetscInt       nstepsmacro = 100;    /* Max macros steps */
  PetscReal      dtmacro = 1.0;        /* Macros step size (see stepover note) */

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);

  /* Obtain a command-line argument for testing */
  test_view = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-test_view",&test_view,NULL);CHKERRQ(ierr);
  monitor = PETSC_TRUE;
  ierr = PetscOptionsGetBool(NULL,NULL,"-monitor",&monitor,NULL);CHKERRQ(ierr);
  nstepsmacro = 100;
  ierr = PetscOptionsGetInt(NULL,NULL,"-nstepsmacro",&nstepsmacro,NULL);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD," *** Will perform %D macro steps of length %f\n",
      nstepsmacro,dtmacro);CHKERRQ(ierr);

  // Note: it might make for a less-confusing code if all command-line 
  //       processing was here, instead of hidden in the ctx setup

  /* Perform all initialization for our problem, allocating data 
     Note that this checks for a command line option -n */
  ierr = setup_ctx(&ctx);CHKERRQ(ierr);

  /* We will use this solution vector as our data object for timestepping */
  ierr = DMCreateGlobalVector( ctx.da_s, &S_s );CHKERRQ(ierr);

  /* must call this after setup_ctx */
  set_initial_condition(&ctx,S_s);CHKERRQ(ierr);

  if (test_view) {
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," *** Viewing S_s initial condition ***\n");CHKERRQ(ierr);
    ierr = VecView(S_s,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  
  /* Set up the Jacobian function (omitted for now) */

  /* Set up timestepper */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,S_s);CHKERRQ(ierr);

#if (defined QUAD)
  ierr = TSSetType(ts,TSBEULER);
  ierr = TSSetTolerances( ts, 1e-10, NULL, 1e-10, NULL );
#else
  ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  ierr = TSSundialsSetTolerance(ts,1e-10,1e-10);CHKERRQ(ierr);
#endif

  /* Set up the RHS Function */
  ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&ctx);CHKERRQ(ierr);

  /* Set up the integration period for first macro step */
  ierr = TSSetDuration(ts,maxsteps,dtmacro);CHKERRQ(ierr); 
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr); /* Note: will not compute at precisely the requested time (TODO: see if interpolate can be made to work) */

  /* Accept command line options */
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  /* Solve macro steps (could also use a monitor which keeps some state) */
  {
    PetscReal time = t0;
    PetscInt stepmacro=0;
    if (monitor) {
      ierr = TSCustomMonitor(ts,stepmacro,time,S_s,&ctx);CHKERRQ(ierr);
    }
    for (stepmacro=1;stepmacro<=nstepsmacro;++stepmacro){
      const PetscInt time_next = (stepmacro + 1) * dtmacro;
      ierr = TSSolve(ts,S_s);CHKERRQ(ierr);
      ierr = TSGetTime(ts,&time);CHKERRQ(ierr);
      if (monitor) {
        ierr = TSCustomMonitor(ts,stepmacro,time,S_s,&ctx);CHKERRQ(ierr);
      }
      ierr = TSSetDuration(ts,maxsteps,time_next);CHKERRQ(ierr); 
    }
  }

  /* For testing, view some things stored in the context.
      Note that there are other viewer implementations, including 
      those that will write to a file which we should be able
      to open with python.
  */
  // TODO: update this to view after every macro step so we can automatically test 
  //       that as well
  if (test_view) {
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," *** Viewing S_s ***\n");CHKERRQ(ierr);
    ierr = VecView(S_s,viewer);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," *** Viewing rhs_s ***\n");CHKERRQ(ierr);
    // TODO: don't store the rhs like this. Just recompute it here
    ierr = VecView(ctx.solution.rhs_s,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Free allocated data in the context */
  ierr = destroy_ctx(&ctx);CHKERRQ(ierr);
  ierr = VecDestroy(&S_s);CHKERRQ(ierr);

  /* Cleanup and finalize */
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}
