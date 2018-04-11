static char help[] =
"SPIDER. Magma ocean timestepper\n\
-n : specify the number of staggered points\n\
-monitor : use custom monitor to dump output\n\
-nstepsmacro : specify the number of macro (output dump) steps\n\
-dtmacro : specify the macro (output dump) time step, in nondimensionalised time (will not be exact!)\n\
-dtmacro_years : specify the macro (output dump) time step, in years\n\
-early, -middle, -late : shortcuts to set -dtmacro and -nstepsmacro (overrides these)\n\
See example input files in examples/ for many more available options\n\
";

#include "petsc.h"
#include "ctx.h"
#include "ic.h"
#include "monitor.h"
#include "rhs.h"
#include "version.h"

static PetscErrorCode PrintSPIDERHeader();
static PetscErrorCode PrintFields(Ctx*);

int main(int argc, char ** argv)
{
  PetscErrorCode   ierr;
  TS               ts;                 /* ODE solver object */
  Vec              sol;                /* Solution Vector (packed vector from several DMs) */
  Ctx              ctx;                /* Solver context */
  Parameters const *P=&ctx.parameters;

  PetscMPIInt     size;

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);

  ierr = PrintSPIDERHeader();CHKERRQ(ierr);

  /* We don't want to take the time to debug things in MPI (though the
     problems are likely minor), so don't allow multi-rank runs */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");


  /* Perform all initialization for our problem, allocating data
     Note that this checks all command-line options. */
  ierr = SetupCtx(&ctx);CHKERRQ(ierr);

  /* Print out the quantities we're solving for */
  ierr = PrintFields(&ctx);CHKERRQ(ierr);

  /* Create Vector for solution provided to timestepper */
  ierr = DimensionalisableFieldGetGlobalVec(ctx.solDF,&sol);CHKERRQ(ierr);

  /* see ic.c : must call after setup_ctx */
  ierr = set_initial_condition( &ctx, sol ); CHKERRQ(ierr);

  /* Set up timestepper */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,sol);CHKERRQ(ierr);

  /* always use SUNDIALS CVODE (BDF) */
  /* must use direct solver, so requires Patrick's hacks */
  ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  // next seems to be set by default
  ierr = TSSundialsSetType(ts,SUNDIALS_BDF);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  //ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_INTERPOLATE);CHKERRQ(ierr);

  /* can tighten tolerances for quad */
#if (defined PETSC_USE_REAL___FLOAT128)
  ierr = TSSundialsSetTolerance(ts, 1.0e-14, 1.0e-14 );CHKERRQ(ierr);
#else
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

  /* Solve macro steps. We have our own outer loop which repeatedly calls TSSolve
     and then calls a custom monitor. We also provide a special monitor to the
     TS object, to output periodically within our macro steps, useful if a
     macro step takes a longer amount of time. */
  {
    PetscReal  time,nexttime;
    MonitorCtx mctx;

    time = P->t0;
    mctx.walltime0 = MPI_Wtime();
    mctx.walltimeprev = mctx.walltime0;
    PetscInt stepmacro=0;
    /* This code proceeds by performing multiple solves with a TS object,
       pausing to optionally produce output before updating the "final" time
       and proceeding again */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Will perform %D macro (output) steps of length %f = %f years\n",
        P->nstepsmacro,(double) P->dtmacro, (double) P->dtmacro_years);CHKERRQ(ierr);
    ierr = TSMonitorSet(ts,TSMonitorWalltimed,&mctx,NULL);CHKERRQ(ierr);
    nexttime = P->dtmacro; // non-dim time
    ierr = TSSetDuration(ts,P->maxsteps,nexttime);CHKERRQ(ierr);
    if (P->monitor) {
      ierr = TSCustomMonitor(ts,P->dtmacro,P->dtmacro_years,stepmacro,time,sol,&ctx,&mctx);CHKERRQ(ierr);
    }
    for (stepmacro=1; stepmacro<=P->nstepsmacro; ++stepmacro){
      ierr = TSSolve(ts,sol);CHKERRQ(ierr);
      if (stepmacro == 1) {
        ierr = TSView(ts,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      }
      ierr = TSGetTime(ts,&time);CHKERRQ(ierr);
      if (P->monitor) {
        ierr = TSCustomMonitor(ts,P->dtmacro,P->dtmacro_years,stepmacro,time,sol,&ctx,&mctx);CHKERRQ(ierr);
      }
      nexttime = (stepmacro + 1) * P->dtmacro; // non-dim
      ierr = TSSetDuration(ts,P->maxsteps,nexttime);CHKERRQ(ierr);
    }
  }

  /* Free allocated data and clean up */
  ierr = DestroyCtx(&ctx);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = VecDestroy(&sol);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PrintSPIDERHeader()
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(
      PETSC_COMM_WORLD,"::::::::::::::::: SPIDER version %s.%s.%s :::::::::::::::::\n\n",
      SPIDER_MAJOR_VERSION,SPIDER_MINOR_VERSION,SPIDER_PATCH_VERSION);CHKERRQ(ierr);
#if defined(PETSC_USE_REAL___FLOAT128)
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using quadruple (128-bit,__float128) precision scalars\n\n");CHKERRQ(ierr);
#elif defined(PETSC_USE_REAL_DOUBLE)
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using double (64-bit) precision scalars\n\n");CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode PrintFields(Ctx *ctx)
{
  PetscErrorCode ierr;
  PetscInt f;
  DM       *sub_dms;

  PetscFunctionBeginUser;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n*** Timestepper will solve for coupled fields:\n");CHKERRQ(ierr);
  ierr = PetscMalloc1(ctx->numFields,&sub_dms);CHKERRQ(ierr);
  ierr = DMCompositeGetEntriesArray(ctx->dm_sol,sub_dms);CHKERRQ(ierr);
  for (f=0; f<ctx->numFields; ++f) {
    Vec vecTemp;
    PetscInt nEntries;
    const char * entriesString;
    ierr = DMGetGlobalVector(sub_dms[f],&vecTemp);CHKERRQ(ierr);
    ierr = VecGetSize(vecTemp,&nEntries);CHKERRQ(ierr);
    entriesString = nEntries > 1 ? "entries" : "entry";
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  %D: %s (%D %s)\n",f,SpiderSolutionFieldDescriptions[ctx->solutionFieldIDs[f]],nEntries,entriesString);CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(sub_dms[f],&vecTemp);CHKERRQ(ierr);
  }
  ierr = PetscFree(sub_dms);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
