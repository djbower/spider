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
  PetscErrorCode  ierr;
  TS              ts;                        /* ODE solver object */
  /* DJB BELOW NEEDS TO BE LENGTH NUM_BASIC+1 */
  /* AND RENAME TO SOMETHING MORE DESCRIPTIVE */
  Vec             S_s;                       /* Solution Vector */
  Ctx             ctx;                       /* Solver context */
  PetscBool       monitor = PETSC_TRUE;      /* Macro step custom monitor (monitor.c) */
  const PetscReal t0 = 0;                    /* Initial time */
  const PetscInt  maxsteps =    1000000000;  /* Unlimited Max internal steps */ /* TODO: figure out if PETSc actually has a way to specify unlimited steps */
  PetscInt        nstepsmacro = 1000000000;  /* Max macros steps */
  PetscReal       dtmacro = 1000.0;          /* Macro step size */
  PetscMPIInt     size;

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);

  /* We don't want to take the time to debug things in MPI (though the
     problems are likely minor), so don't allow multi-rank runs */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");

  /* Obtain a command-line argument for testing */
  ierr = PetscOptionsGetBool(NULL,NULL,"-monitor",&monitor,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-nstepsmacro",&nstepsmacro,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dtmacro",&dtmacro,NULL);CHKERRQ(ierr);

  /* This code proceeds by performing multiple solves with a TS object,
     pausing to optionally produce output before updating the "final" time
     and proceeding again */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Will perform %D macro (output) steps of length %f\n",
      nstepsmacro,(double) dtmacro);CHKERRQ(ierr);

  /* Note: it might make for a less-confusing code if all command-line 
           processing was here, instead of hidden in the ctx setup */

  /* Perform all initialization for our problem, allocating data 
     Note that this checks for a command line option -n */
  ierr = setup_ctx(&ctx);CHKERRQ(ierr);

  /* We will use this solution vector as our data object for timestepping */
  /* DJB THIS GLOBAL VECTOR NEEDS TO HAVE LENGTH NUM_BASIC+1 AND RENAME */
  ierr = DMCreateGlobalVector( ctx.da_s, &S_s );CHKERRQ(ierr);

  /* must call this after setup_ctx */
  set_initial_condition(&ctx,S_s);CHKERRQ(ierr);
  
  /* Set up the Jacobian function (omitted for now) */

  /* Set up timestepper */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,S_s);CHKERRQ(ierr);

#if (defined PETSC_USE_REAL___FLOAT128)
  /* PDS TODO: pick a solver that we can test quad precision for.
     Perhaps current arkimex is OK, but it seems to run slow for
     double precision, compared to BDF methods in SUNDIALS */
  ierr = TSSetType(ts,TSARKIMEX);CHKERRQ(ierr); 
  ierr = TSARKIMEXSetType(ts, TSARKIMEX1BEE ); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue(NULL,"-snes_mf","1");CHKERRQ(ierr);
  ierr = TSSetTolerances(ts, 1e-11, NULL, 1e-11, NULL);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
#else
  /* SUNDIALS */
  ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  ierr = TSSundialsSetTolerance(ts,1e-11,1e-11);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
#endif

  /* Set up the RHS Function */
  ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&ctx);CHKERRQ(ierr);

  /* Set up the integration period for first macro step */
  ierr = TSSetDuration(ts,maxsteps,dtmacro);CHKERRQ(ierr); 

  /* Set a very small initial timestep to prevent problems with 
     challenging initial conditions and adaptive steppers */
  ierr = TSSetInitialTimeStep(ts,0.0,1e-10);CHKERRQ(ierr);

  /* Accept command line options */
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  /* Solve macro steps (could also use a monitor which keeps some state) */
  {
    PetscReal time = t0;
    double walltime0 = MPI_Wtime();
    PetscInt stepmacro=0;
    if (monitor) {
      ierr = TSCustomMonitor(ts,stepmacro,time,t0,S_s,&ctx,walltime0);CHKERRQ(ierr);
    }
    for (stepmacro=1;stepmacro<=nstepsmacro;++stepmacro){
      ierr = TSSolve(ts,S_s);CHKERRQ(ierr);
      ierr = TSGetTime(ts,&time);CHKERRQ(ierr);
      if (monitor) {
        ierr = TSCustomMonitor(ts,stepmacro,time,t0,S_s,&ctx,walltime0);CHKERRQ(ierr);
      }
      ierr = TSSetDuration(ts,maxsteps,(stepmacro + 1) * dtmacro);CHKERRQ(ierr); 
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
