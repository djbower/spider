static char help[] =
"Magma ocean timestepper\n\
-n : specify the number of staggered points\n\
-sinit : specify an entropy value to base the initial condition upon\n\
-monitor : use custom monitor to dump output\n\
-test_view : with -monitor, also print out solution and rhs for testing\n\
-nstepsmacro : specify the number of macro (output dump) steps\n\
-dtmacro : specify the macro (output dump) time step, in nondimensionalised time (will not be exact!)\n\
-dtmacro_years : specify the macro (output dump) time step, in years\n\
-curves : specify a set of data files to load. Currently accepts 'andrault2011' (default) and 'stixrude2009'\n\
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
  Vec             dSdr_b_aug;                /* Solution Vector (basic points, plus extra points) */
  Ctx             ctx;                       /* Solver context */
  PetscBool       monitor = PETSC_TRUE;      /* Macro step custom monitor (monitor.c) */
  const PetscReal t0 = 0;                    /* Initial time */
  const PetscInt  maxsteps =    1000000000;  /* Unlimited Max internal steps */ /* TODO: figure out if PETSc actually has a way to specify unlimited steps */
  PetscReal       nexttime;                  /* next time to integrate to */

  /* FIXME: at the moment we uncomment two lines below to switch between
     early, middle, and late evolution */
  /* early evolution to about 10 kyr */
  PetscInt        nstepsmacro = 1000;  /* Max macros steps */
  PetscReal       dtmacro = 1000000;        /* Macro step size (in nondimensional units. About a year) */
  /* middle evolution to about 100 Myr */
  //PetscInt        nstepsmacro = 10000;  /* Max macros steps */
  //PetscReal       dtmacro = 1000000000; /* Macro step size (roughly years) */
  /* late evolution to 4.55 Byr */
  //PetscInt        nstepsmacro = 455;  /* Max macros steps */
  //PetscReal       dtmacro = 1000000000000; /* Macro step size (roughly years) */

  PetscMPIInt     size;

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);

  /* We don't want to take the time to debug things in MPI (though the
     problems are likely minor), so don't allow multi-rank runs */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");

  /* 
   **All** parameters/settings should be inside ctx->parameters . Here, we
     1. Initialize these to default values
     2. Obtain any inputs from the user at the command line and update
     3. Print out all paramters
     After this point, **no** command line options should be processed,
     and ctx->paramters should **never** change 
    */
  // TODO : move everything into Parameters
  // TODO: initialize Parameters here
  // TODO: move all commadn line processing into a function to update parameters
  // TODO : print out all params with scaled value / dim value / scaling  formula and value

  /* Obtain command-line options for simulation time frame and monitoring */
  {
    PetscReal dtmacro_years = dtmacro * TIME0;
    PetscBool dtmacro_set = PETSC_FALSE, dtmacro_years_set = PETSC_FALSE;
    ierr = PetscOptionsGetBool(NULL,NULL,"-monitor",&monitor,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,NULL,"-nstepsmacro",&nstepsmacro,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-dtmacro",&dtmacro,&dtmacro_set);CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL,NULL,"-dtmacro_years",&dtmacro_years,&dtmacro_years_set);CHKERRQ(ierr);
    if (dtmacro_set) {
        if (dtmacro_years_set) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: both -dtmacro and -dtmacro_years provided. Using -dtmacro\n");CHKERRQ(ierr);
        }
        dtmacro_years = dtmacro * TIME0; /* FIXME : TIME0 currently wrong */
    } else if (dtmacro_years_set) {
        dtmacro = dtmacro_years / TIME0; /* FIXME : TIME0 currently wrong */
    }

    /* This code proceeds by performing multiple solves with a TS object,
       pausing to optionally produce output before updating the "final" time
       and proceeding again */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Will perform %D macro (output) steps of length %f = %f years\n",
        nstepsmacro,(double) dtmacro, (double) dtmacro_years);CHKERRQ(ierr);
  }

  /* Note: it might make for a less-confusing code if all command-line 
           processing was here, instead of hidden in the ctx setup */

  /* Perform all initialization for our problem, allocating data 
     Note that this checks for a command line option -n */
  ierr = setup_ctx(&ctx);CHKERRQ(ierr);

  ///////////////////////
  /* initial condition */
  ///////////////////////
  /* must call after setup_ctx */

  ierr = DMCreateGlobalVector( ctx.da_b, &dSdr_b );CHKERRQ(ierr);

  set_ic_dSdr( &ctx, dSdr_b ); CHKERRQ(ierr);

  /* Create a special vector with extra points and transfer IC */
  ierr = CreateAug(dSdr_b,&dSdr_b_aug);CHKERRQ(ierr);
  ierr = ToAug(dSdr_b,dSdr_b_aug);CHKERRQ(ierr);

  /* add initial values of other quantities to the augmented vector */
  set_ic_aug( &ctx, dSdr_b_aug ); CHKERRQ(ierr);

  ///////////////////////////
  /* end initial condition */
  ///////////////////////////

  /* Set up timestepper */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,dSdr_b_aug);CHKERRQ(ierr);

  /* always use SUNDIALS CVODE (BDF) */
  /* must use direct solver, so requires Patrick's hacks */
  ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  // next seems to be set by default
  ierr = TSSundialsSetType(ts,SUNDIALS_BDF);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  //ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_INTERPOLATE);CHKERRQ(ierr);

  /* can tighten tolerances for quad */
#if (defined PETSC_USE_REAL___FLOAT128)
  // 1.0E-12 crashes at time step 5
  ierr = TSSundialsSetTolerance(ts, 1.0e-14, 1.0e-14 );CHKERRQ(ierr);
#else
  // TODO: 1.0E-4 works for simple bottom up case
  ierr = TSSundialsSetTolerance(ts, 1.0e-4, 1.0e-4 );CHKERRQ(ierr);
#endif

  /* Set up the RHS Function */
  ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&ctx);CHKERRQ(ierr);

  /* Set up the integration period for first macro step */
  nexttime = dtmacro; // non-dim time
  ierr = TSSetDuration(ts,maxsteps,nexttime);CHKERRQ(ierr);

  /* Set a very small initial timestep to prevent problems with 
     challenging initial conditions and adaptive steppers */
  /* DJB with revised code this probably is not necessary anymore */
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
      ierr = TSCustomMonitor(ts,dtmacro,stepmacro,time,t0,timeprev,dSdr_b_aug,&ctx,walltime0,&walltimeprev);CHKERRQ(ierr);
    }
    for (stepmacro=1;stepmacro<=nstepsmacro;++stepmacro){
      ierr = TSSolve(ts,dSdr_b_aug);CHKERRQ(ierr);
      timeprev = time;
      ierr = TSGetTime(ts,&time);CHKERRQ(ierr);
      if (monitor) {
        ierr = TSCustomMonitor(ts,dtmacro,stepmacro,time,t0,timeprev,dSdr_b_aug,&ctx,walltime0,&walltimeprev);CHKERRQ(ierr);
      }
      nexttime = (stepmacro + 1) * dtmacro; // non-dim
      ierr = TSSetDuration(ts,maxsteps,nexttime);CHKERRQ(ierr); 
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
