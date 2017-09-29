static char help[] =
"Magma ocean timestepper\n\
-n : specify the number of staggered points\n\
-sinit : specify an entropy value to base the initial condition upon\n\
-monitor : use custom monitor to dump output\n\
-nstepsmacro : specify the number of macro (output dump) steps\n\
-dtmacro : specify the macro (output dump) time step, in nondimensionalised time (will not be exact!)\n\
-dtmacro_years : specify the macro (output dump) time step, in years\n\
-curves : specify a set of data files to load. Currently accepts 'andrault2011' (default) and 'stixrude2009'\n\
-early, -middle, -late: shortcuts to set -dtmacro and -nstepsmacro (overrides these)\n\
";

#include "petsc.h"
#include "ctx.h"
#include "ic.h"
#include "monitor.h"
#include "rhs.h"
#include "aug.h"

int main(int argc, char ** argv)
{
  PetscErrorCode   ierr;
  TS               ts;                        /* ODE solver object */
  Vec              dSdr_b;
  Vec              dSdr_b_aug;                /* Solution Vector (basic points, plus extra points) */
  Ctx              ctx;                       /* Solver context */
  Parameters const *P;
  Constants const  *C;

  PetscMPIInt     size;

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);

  /* We don't want to take the time to debug things in MPI (though the
     problems are likely minor), so don't allow multi-rank runs */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");


  /* Perform all initialization for our problem, allocating data
     Note that this checks all command-line options. */
  ierr = SetupCtx(&ctx);CHKERRQ(ierr);
  P = &ctx.parameters;
  C = &P->constants;

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
  /* TODO: 1.0E-4 works for simple bottom up case with grey-body atmosphere
         with constant emissivity */
  ierr = TSSundialsSetTolerance(ts, 1.0e-6, 1.0e-6 );CHKERRQ(ierr);
#endif

  /* Set up the RHS Function */
  ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&ctx);CHKERRQ(ierr);


  /* Set a very small initial timestep to prevent problems with
     challenging initial conditions and adaptive steppers */
  /* DJB with revised code this probably is not necessary anymore */
  //ierr = TSSetInitialTimeStep(ts,0.0,1e-10);CHKERRQ(ierr);

  /* Accept command line options */
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  /* Solve macro steps (could also use a monitor which keeps some state) */
  {
    PetscReal time,nexttime;
    double    walltime0,walltimeprev;

    time = P->t0;
    walltime0 = MPI_Wtime();
    walltimeprev = walltime0;
    PetscInt stepmacro=0;
    /* This code proceeds by performing multiple solves with a TS object,
       pausing to optionally produce output before updating the "final" time
       and proceeding again */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Will perform %D macro (output) steps of length %f = %f years\n",
        P->nstepsmacro,(double) P->dtmacro, (double) (P->dtmacro * C->TIMEYRS));CHKERRQ(ierr);
    nexttime = P->dtmacro; // non-dim time
    ierr = TSSetDuration(ts,P->maxsteps,nexttime);CHKERRQ(ierr);
    if (P->monitor) {
      ierr = TSCustomMonitor(ts,P->dtmacro,stepmacro,time,dSdr_b_aug,&ctx,walltime0,&walltimeprev);CHKERRQ(ierr);
    }
    for (stepmacro=1; stepmacro<=P->nstepsmacro; ++stepmacro){
      ierr = TSSolve(ts,dSdr_b_aug);CHKERRQ(ierr);
      ierr = TSGetTime(ts,&time);CHKERRQ(ierr);
      if (P->monitor) {
        ierr = TSCustomMonitor(ts,P->dtmacro,stepmacro,time,dSdr_b_aug,&ctx,walltime0,&walltimeprev);CHKERRQ(ierr);
      }
      nexttime = (stepmacro + 1) * P->dtmacro; // non-dim
      ierr = TSSetDuration(ts,P->maxsteps,nexttime);CHKERRQ(ierr);
    }
  }

  /* Free allocated data in the context */
  ierr = DestroyCtx(&ctx);CHKERRQ(ierr);

  /* Destroy solution vector */
  ierr = VecDestroy(&dSdr_b);CHKERRQ(ierr);
  ierr = VecDestroy(&dSdr_b_aug);CHKERRQ(ierr);

  /* Cleanup and finalize */
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}
