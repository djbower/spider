static char help[] =
"SPIDER. Magma ocean timestepper\n\
-n : specify the number of staggered points\n\
-monitor : use custom monitor to dump output\n\
-nstepsmacro : specify the number of macro (output dump) steps\n\
-dtmacro : specify the macro (output dump) time step, in years\n\
See example input files in examples/ for many more available options\n\
";

#include "petsc.h"
#include "ctx.h"
#include "ic.h"
#include "monitor.h"
#include "rhs.h"
#include "version.h"
#include "rollback.h"
#include "poststep.h"
#include "util.h"

static PetscErrorCode PrintSPIDERHeader();
static PetscErrorCode PrintFields(Ctx*);

int main(int argc, char ** argv)
{
  PetscErrorCode   ierr;
  TS               ts;                 /* ODE solver object */
  Vec              sol;                /* Solution Vector (packed vector from several DMs) */
  Ctx              ctx;                /* Solver context */

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);

  ierr = PrintSPIDERHeader();CHKERRQ(ierr);

  /* Print options file(s) being used. If no options file is specified, use a default */
  {
    PetscBool options_file_provided = PETSC_FALSE;
    char      default_options_filename[PETSC_MAX_PATH_LEN] = "tests/opts/blackbody50.opts";

    for (int i=0; i<argc; ++i) {
      PetscBool detected;

      ierr = PetscStrcmp(argv[i],"-options_file",&detected);CHKERRQ(ierr);
      if (detected && i < argc-1) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Using options file: %s\n",argv[i+1]);CHKERRQ(ierr);
      }
      options_file_provided = options_file_provided || detected;
    }
    if (!options_file_provided) {
      if (argc > 1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"If you supply any options, you must supply an options file explicitly. For the default, use \n   -options_file %s",default_options_filename);
      ierr = MakeRelativeToSourcePathAbsolute(default_options_filename);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Using default options file: %s\n",default_options_filename);CHKERRQ(ierr);
      ierr = PetscOptionsInsertFile(PETSC_COMM_WORLD,NULL,default_options_filename,PETSC_TRUE);CHKERRQ(ierr);
    }
  }

  /* We haven't debugged in parallel, and usually use a serial, dense
     solver within the timestepper, so don't allow multi-rank runs */
  {
    PetscMPIInt size;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");
  }

  /* Perform all initialization for our problem, allocating data
     Note that this checks all command-line options. */
  ierr = SetupCtx(&ctx);CHKERRQ(ierr);

  // important this is here
  Parameters const P=ctx.parameters;

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

  /* always use SUNDIALS CVODE (BDF) with a dense (serial) inner linear solver */
  ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  ierr = TSSundialsSetUseDense(ts,PETSC_TRUE);CHKERRQ(ierr);
  ierr = TSSundialsSetType(ts,SUNDIALS_BDF);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_INTERPOLATE);CHKERRQ(ierr);

  /* can tighten tolerances for quad */
#if (defined PETSC_USE_REAL___FLOAT128)
  ierr = TSSundialsSetTolerance(ts, 1.0e-14, 1.0e-14 );CHKERRQ(ierr);
#else
  ierr = TSSundialsSetTolerance(ts, 1.0e-6, 1.0e-6 );CHKERRQ(ierr);
#endif

  /* Set up the RHS Function */
  ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&ctx);CHKERRQ(ierr);

  /* Accept command line options */
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);


  /* Solve macro steps. We have our own outer loop which repeatedly calls TSSolve
     and then calls a custom monitor. We also provide a special monitor to the
     TS object, to output periodically within our macro steps, useful if a
     macro step takes a longer amount of time. */
  {
    PetscReal  time,nexttime;
    PetscInt   stepmacro;
    MonitorCtx mctx;

    /* Monitor to output periodically, based on wall time */
    mctx.walltime0 = MPI_Wtime();
    mctx.walltimeprev = mctx.walltime0;
    mctx.outputDirectoryExistenceConfirmed = PETSC_FALSE;
    ierr = TSMonitorSet(ts,TSMonitorWalltimed,&mctx,NULL);CHKERRQ(ierr);

    /* Proceed by performing multiple solves with a TS object,
       pausing to optionally produce output before updating the "final" time
       and proceeding again */

    ierr = PetscPrintf(PETSC_COMM_WORLD,"*** Starting at t0 = %f, Will perform %D macro (output) steps of length %f = %f years\n",
        P->t0,P->nstepsmacro,(double) P->dtmacro, (double) (P->dtmacro*(*P->scaling_constants).TIMEYRS) );CHKERRQ(ierr);
    time = P->t0;
    ierr = TSSetTime(ts,time);CHKERRQ(ierr);
    nexttime = P->t0 + P->dtmacro;
    stepmacro = P->stepmacro; /* start macrostep */
    ierr = TSSetMaxSteps(ts,P->maxsteps);CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts,nexttime);CHKERRQ(ierr);

    /* Final setup logic (needs to be here because some things, like TSSetTime(),
       won't work properly for the SUNDIALS implementation if called after TSSetUp()) */
    ierr = TSSetUp(ts);CHKERRQ(ierr);

    /* output first time step (initial condition) */
    /* need this here to update quantities to make poststep data */
    if (P->monitor) {
      ierr = TSCustomMonitor(ts,stepmacro,time,sol,&ctx,&mctx);CHKERRQ(ierr);
    }

    /* Activate rollback capability and PostStep logic (needs to happen post-SetUp()) */
    if (P->rollBackActive) {
      ierr = TSRollBackGenericActivate(ts);CHKERRQ(ierr);
    }
    if (P->postStepActive) {
      ierr = TSSetApplicationContext(ts,&ctx);CHKERRQ(ierr);
      ierr = PostStepDataInitialize(&ctx);CHKERRQ(ierr); /* use initial condition here */
      ierr = TSSetPostStep(ts,PostStep);CHKERRQ(ierr);
    }

    for (stepmacro=P->stepmacro+1; stepmacro<=P->nstepsmacro; ++stepmacro){
      ierr = TSSolve(ts,sol);CHKERRQ(ierr);
      /* view the timestepper parameters once */
      if (stepmacro == P->stepmacro+1) {
        ierr = TSView(ts,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      }
      ierr = TSGetSolveTime(ts,&time);CHKERRQ(ierr);
      if (P->monitor) {
        ierr = TSCustomMonitor(ts,stepmacro,time,sol,&ctx,&mctx);CHKERRQ(ierr);
      }
      if (ctx.stopEarly) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Stopping macro timestepping loop early!\n");CHKERRQ(ierr);
        ierr = TSCustomMonitor(ts,stepmacro,time,sol,&ctx,&mctx);CHKERRQ(ierr);
        break;
      }
      nexttime = P->t0 + ((stepmacro - P->stepmacro + 1) * P->dtmacro);
      ierr = TSSetMaxSteps(ts,P->maxsteps);CHKERRQ(ierr);
      ierr = TSSetMaxTime(ts,nexttime);CHKERRQ(ierr);
    }
  }

  /* Free allocated data and clean up */
  ierr = DestroyCtx(&ctx);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
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
