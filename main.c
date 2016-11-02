static char help[] =
"Parallel magma ocean timestepper\n\
-n : specify the number of staggered points\n\
-sinit : specify an entropy value to base the initial condition upon\n\
-monitor : use custom monitor to dump output\n\
-test_view : with -monitor, also print out solution and rhs for testing\n\
-nstepsmacro : specify the number of macro (output dump) steps\n\
-dtmacro : specify the macro (output dump) time step (will not be exact!)\n\
";

/* Note: if you would like more verbose output, see the preprocessor defines
         in global_defs.h */

#include "petsc.h"
#include "ic.h" 
#include "ctx.h"
#include "monitor.h"
#include "rhs.h"

int main(int argc, char ** argv)
{
  PetscErrorCode ierr;
  TS             ts;                    /* ODE solver object */
  Vec            S_s;                   /* Solution Vector */
  Ctx            ctx;                   /* Solver context */
  PetscBool      monitor = PETSC_TRUE;  /* Macros step custom monitor (monitor.c) */
  const PetscReal t0 = 0;               /* Initial time */
  const PetscInt maxsteps = 1000000000; /* Max internal steps */
  PetscInt       nstepsmacro = 100;     /* Max macros steps */
  PetscReal      dtmacro = 1.0;         /* Macros step size (see stepover note) */

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);

  // Emit some warnings as a placeholder to ensure that we don't forget the QUAD flag if its needed.
#if ((defined QUAD) && !(defined (PETSC_USE_REAL___FLOAT128)))
#error Configuration error. You appear to be using QUAD without a quad-precision build of PETSc
#endif
#if (!(defined QUAD) && (defined (PETSC_USE_REAL___FLOAT128)))
#error Configuration error. You appear to be using a quad-precision build of PETSc without definiing QUAD here
#endif

  /* Obtain a command-line argument for testing */
  ierr = PetscOptionsGetBool(NULL,NULL,"-monitor",&monitor,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-nstepsmacro",&nstepsmacro,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dtmacro",&dtmacro,NULL);CHKERRQ(ierr);

  /* This code proceeds by performing multiple solves with a TS object,
     pausing to optionally produce output before updating the "final" time
     and proceeding again */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Will perform %D macro (output) steps of length %f\n",
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
  
  /* Set up the Jacobian function (omitted for now) */

  /* Set up timestepper */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,S_s);CHKERRQ(ierr);

#if (defined QUAD)
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr); //PDS TODO: replace this with TSARKIMEX or something else that's more sophisticated and native
  ierr = TSSetTolerances( ts, 1e-10, NULL, 1e-10, NULL );CHKERRQ(ierr);
#else
  ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  ierr = TSSundialsSetTolerance(ts,1e-10,1e-10);CHKERRQ(ierr);
#endif

  /* Set up the RHS Function */
  ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&ctx);CHKERRQ(ierr);

  /* Set up the integration period for first macro step */
  ierr = TSSetDuration(ts,maxsteps,dtmacro);CHKERRQ(ierr); 
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr); /* Note: will not compute at precisely the requested time (PDS TODO: see if interpolate or matchstep can be made to work) */

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

  /* Free allocated data in the context */
  ierr = destroy_ctx(&ctx);CHKERRQ(ierr);

  /* Destroy solution vector */
  ierr = VecDestroy(&S_s);CHKERRQ(ierr);

  /* Cleanup and finalize */
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}
