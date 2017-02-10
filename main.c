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
#include "aug.h"

int main(int argc, char ** argv)
{
  PetscErrorCode  ierr;
  TS              ts;                        /* ODE solver object */
  Vec             dSdr_b;
  Vec             dSdr_b_aug;                /* Solution Vector (basic points, plus an extra point) */
  Ctx             ctx;                       /* Solver context */
  PetscBool       monitor = PETSC_TRUE;      /* Macro step custom monitor (monitor.c) */
  const PetscReal t0 = 0;                    /* Initial time */
  const PetscInt  maxsteps =    1000000000;  /* Unlimited Max internal steps */ /* TODO: figure out if PETSc actually has a way to specify unlimited steps */
  PetscInt        nstepsmacro = 1000000000;  /* Max macros steps */
  PetscReal       dtmacro = 100000.0;          /* Macro step size */
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

  ierr = DMCreateGlobalVector( ctx.da_b, &dSdr_b );CHKERRQ(ierr);

  /* must call this after setup_ctx */
  set_initial_condition(&ctx,dSdr_b);CHKERRQ(ierr);
  
  /* Set up the Jacobian function (omitted for now) */

  /* Create a special vector with an extra point and transfer IC */
  ierr = CreateAug(dSdr_b,&dSdr_b_aug);CHKERRQ(ierr);
  ierr = ToAug(dSdr_b,dSdr_b_aug);CHKERRQ(ierr);

  /* include initial entropy at staggered node */
  /* TODO: is this over-rideable from the command line? */
  ierr = VecSetValue(dSdr_b_aug,0,SINIT_DEFAULT,INSERT_VALUES);CHKERRQ(ierr);

  /* Set up timestepper */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,dSdr_b_aug);CHKERRQ(ierr);

#if (defined PETSC_USE_REAL___FLOAT128)
  /* PDS TODO: pick a solver that we can test quad precision for.
     Perhaps current arkimex is OK, but it seems to run slow for
     double precision, compared to BDF methods in SUNDIALS */
  ierr = TSSetType(ts,TSARKIMEX);CHKERRQ(ierr); 
  ierr = TSARKIMEXSetType(ts, TSARKIMEX5 ); CHKERRQ(ierr);
  ierr = PetscOptionsSetValue(NULL,"-snes_mf","1");CHKERRQ(ierr);
  ierr = TSSetTolerances(ts, 1e-11, NULL, 1e-11, NULL);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
#else
  /* SUNDIALS */
  ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  ierr = TSSundialsSetTolerance(ts,1e-11,1e-11);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  //ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_INTERPOLATE);CHKERRQ(ierr);
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
    PetscReal time = t0,timeprev;
    double walltime0 = MPI_Wtime(),walltimeprev;
    timeprev = t0;
    walltimeprev = walltime0;
    PetscInt stepmacro=0;
    if (monitor) {
      ierr = TSCustomMonitor(ts,stepmacro,time,t0,timeprev,dSdr_b_aug,&ctx,walltime0,&walltimeprev);CHKERRQ(ierr);
    }
    for (stepmacro=1;stepmacro<=nstepsmacro;++stepmacro){
      ierr = TSSolve(ts,dSdr_b_aug);CHKERRQ(ierr);
      timeprev = time;
      ierr = TSGetTime(ts,&time);CHKERRQ(ierr);
      if (monitor) {
        ierr = TSCustomMonitor(ts,stepmacro,time,t0,timeprev,dSdr_b_aug,&ctx,walltime0,&walltimeprev);CHKERRQ(ierr);
      }
      ierr = TSSetDuration(ts,maxsteps,(stepmacro + 1) * dtmacro);CHKERRQ(ierr); 
    }
  }

  /* Free allocated data in the context */
  ierr = destroy_ctx(&ctx);CHKERRQ(ierr);

  /* Destroy solution vector */
  ierr = VecDestroy(&dSdr_b);CHKERRQ(ierr);
  ierr = VecDestroy(&dSdr_b_aug);CHKERRQ(ierr);

  /* Cleanup and finalize */
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}
