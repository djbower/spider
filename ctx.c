#include "ctx.h"
#include "mesh.h"
#include "twophase.h"
#include "util.h"

static PetscErrorCode CtxCreateFields(Ctx* ctx);

/* Set up the Context */
PetscErrorCode SetupCtx(Ctx* ctx)
{
  PetscErrorCode             ierr;

  PetscFunctionBeginUser;

  /*
   **All** parameters/settings should be inside ctx->parameters . Here, we
     1. Initialize these to default values
     2. Obtain any inputs from the user at the command line and update
     3. Print out all parameters
     After this point, no custom command line options should be processed,
     and ctx->parameters should never change (meaning that we can and should use
     pointers-to-constants to refer to ctx->parameters).

   Generally, the above is adhered to, although when we set the initial condition
   we must update the parameters to be self-consistent (e.g., the initial volatile
   abundances cannot be known a priori if chemical reactions are turned on).
  */

  ierr = ParametersCreate(&ctx->parameters);CHKERRQ(ierr);

  Parameters const           P = ctx->parameters;
  AtmosphereParameters const Ap = P->atmosphere_parameters;

  ierr = PrintParameters(P);CHKERRQ(ierr);

  /* Set up a parallel structured grid as DMComposite with several included DMDAs
     This is used to define vectors which hold the solution. The included
     DMDAs are used to create vectors which hold various properties
     living on primal and staggered nodes as well as on auxiliary domains.
  */

  const PetscInt stencilWidth = 1;
  const PetscInt dof = 1;
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,P->numpts_b,      dof,stencilWidth,NULL,&ctx->da_b        );CHKERRQ(ierr);
  ierr = DMSetUp(ctx->da_b);CHKERRQ(ierr);

  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,P->numpts_s,      dof,stencilWidth,NULL,&ctx->da_s        );CHKERRQ(ierr);
  ierr = DMSetUp(ctx->da_s);CHKERRQ(ierr);

  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,1,                dof,0,           NULL,&ctx->da_point    );CHKERRQ(ierr);
  ierr = DMSetUp(ctx->da_point);CHKERRQ(ierr);
  if (Ap->n_volatiles > 0) {
    ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,Ap->n_volatiles,dof,0           ,NULL,&ctx->da_volatiles);CHKERRQ(ierr); /* A collection of 0-dimensional bulk quantities */
    ierr = DMSetUp(ctx->da_volatiles);CHKERRQ(ierr);
  } else {
    ctx->da_volatiles = NULL;
  }
  if (Ap->n_reactions > 0) {
    ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,Ap->n_reactions,dof,0           ,NULL,&ctx->da_reactions);CHKERRQ(ierr); /* A collection of 0-dimensional bulk quantities */
    ierr = DMSetUp(ctx->da_reactions);CHKERRQ(ierr);
  } else {
    ctx->da_reactions = NULL;
  }

  /* Create a composite DM of the basic nodes plus additional quantities.
    This allows us to create a vector to solve for all of these quantities.  */
  ierr = DMCompositeCreate(PETSC_COMM_WORLD,&ctx->dm_sol);CHKERRQ(ierr);

  /* Define some more-descriptive names for the solution fields, for output */

  /* Note: the order in which DMCompositeAddDM() is called matters!
     This should be the ONLY place this function is called in the code.
     To try to avoid confusion, we maintain a list of strings describing each
     entry (indexed in the same way as the sub-DMs).
     This is intended to catch errors when we change the number of fields we are
     simultaneously solving for. */
  {
    ScalingConstants const SC = ctx->parameters->scaling_constants;
    PetscInt i,f;
    PetscScalar *sol_scalings, *rhs_scalings;

    /* We (over)allocate assuming a maximum of SPIDER_NUM_FIELD_IDS fields, that is that
       you will never have more than one copy of the same field.
       Note that ctx->solutionSlots has an extra entry to account for the extra,
       "undefined" entry in SpiderSolutionFieldID */
    ierr = PetscMalloc2(SPIDER_NUM_FIELD_IDS,&ctx->solutionFieldIDs,SPIDER_NUM_FIELD_IDS+1,&ctx->solutionSlots);CHKERRQ(ierr);
    ierr = PetscMalloc1(SPIDER_NUM_FIELD_IDS,&sol_scalings);CHKERRQ(ierr);
    ierr = PetscMalloc1(SPIDER_NUM_FIELD_IDS,&rhs_scalings);CHKERRQ(ierr);
    for (i=0; i<SPIDER_NUM_FIELD_IDS; ++i) ctx->solutionSlots[i] = -1;
    f = 0;

    ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_b);CHKERRQ(ierr);
    ctx->solutionFieldIDs[f] = SPIDER_SOLUTION_FIELD_DSDXI_B;
    ctx->solutionSlots[ctx->solutionFieldIDs[f]] = f;
    sol_scalings[f] = SC->DSDR;
    rhs_scalings[f] = SC->DSDR / SC->TIME;
    ++f;

    ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_point);CHKERRQ(ierr);
    ctx->solutionFieldIDs[f] = SPIDER_SOLUTION_FIELD_S0;
    ctx->solutionSlots[ctx->solutionFieldIDs[f]] = f;
    sol_scalings[f] = SC->ENTROPY;
    rhs_scalings[f] = SC->ENTROPY / SC->TIME;
    ++f;

    if (ctx->da_volatiles) {
      ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_volatiles);CHKERRQ(ierr);
      ctx->solutionFieldIDs[f] = SPIDER_SOLUTION_FIELD_MO_VOLATILES;
      ctx->solutionSlots[ctx->solutionFieldIDs[f]] = f;
      sol_scalings[f] = SC->PRESSURE;
      rhs_scalings[f] = SC->PRESSURE / SC->TIME;
      ++f;
    }

    if (ctx->da_reactions) {
      ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_reactions);CHKERRQ(ierr);
      ctx->solutionFieldIDs[f] = SPIDER_SOLUTION_FIELD_MO_REACTIONS;
      ctx->solutionSlots[ctx->solutionFieldIDs[f]] = f;
      sol_scalings[f] = SC->VOLATILE * 4.0 * PETSC_PI * SC->MASS; // physical volatile reservoir mass scaling
      rhs_scalings[f] = SC->VOLATILE * 4.0 * PETSC_PI * SC->MASS / SC->TIME;
      ++f;
    }

    ierr = DMCompositeGetNumberDM(ctx->dm_sol,&ctx->numFields);CHKERRQ(ierr); /* For convenience */

    /* Create a DimensionalisableField, referring to this DMComposite,
       for keeping track of the "solution" as seen by the solver. We will
       extra the vector from this to use as our solution vector */
    ierr = DimensionalisableFieldCreate(&ctx->solDF,ctx->dm_sol,sol_scalings,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetName(ctx->solDF,"SPIDER solution");CHKERRQ(ierr);
    for (f=0; f<ctx->numFields; ++f) {
      ierr = DimensionalisableFieldSetSubdomainName(ctx->solDF,f,SpiderSolutionFieldDescriptions[ctx->solutionFieldIDs[f]]);CHKERRQ(ierr);
      ierr = DimensionalisableFieldSetSubdomainUnits(ctx->solDF,f,SpiderSolutionFieldUnits[ctx->solutionFieldIDs[f]]);CHKERRQ(ierr);
    }

    /* for rhs */
    ierr = DimensionalisableFieldCreate(&ctx->rhsDF,ctx->dm_sol,rhs_scalings,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetName(ctx->rhsDF,"SPIDER rhs");CHKERRQ(ierr);
    for (f=0; f<ctx->numFields; ++f) {
      ierr = DimensionalisableFieldSetSubdomainName(ctx->rhsDF,f,SpiderSolutionFieldDescriptions[ctx->solutionFieldIDs[f]]);CHKERRQ(ierr);
      ierr = DimensionalisableFieldSetSubdomainUnits(ctx->rhsDF,f,SpiderRhsFieldUnits[ctx->solutionFieldIDs[f]]);CHKERRQ(ierr);
    }

    ierr = PetscFree(sol_scalings);CHKERRQ(ierr);
    ierr = PetscFree(rhs_scalings);CHKERRQ(ierr);
  }

  /* Continue to initialize context with distributed data */
  ierr = CtxCreateFields(ctx);CHKERRQ(ierr);

  /* Create a work vector */
  ierr = DMCreateLocalVector(ctx->da_b,&ctx->work_local_b);CHKERRQ(ierr);

  /* Populate vectors (initial condition is set in main.c) */
  set_mesh(ctx);

  ierr = initialise_atmosphere( &ctx->atmosphere, ctx->parameters->atmosphere_parameters, ctx->parameters->scaling_constants );CHKERRQ(ierr);

  /* Initialize PostStep reference data */
  ctx->postStepData = NULL;

  /* Initialize control flags */
  ctx->stopEarly = PETSC_FALSE;

  PetscFunctionReturn(0);
}

PetscErrorCode DestroyCtx(Ctx* ctx)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBeginUser;

  /* Destroy data allocated in Ctx */
  ierr = DimensionalisableFieldDestroy(&ctx->solDF);CHKERRQ(ierr);
  ierr = DimensionalisableFieldDestroy(&ctx->rhsDF);CHKERRQ(ierr);
  for (i=0;i<NUMMESHVECS_B;++i){
    ierr = DimensionalisableFieldDestroy(&ctx->mesh.meshFields_b[i]);CHKERRQ(ierr);
  }
  for (i=0;i<NUMMESHVECS_S;++i){
    ierr = DimensionalisableFieldDestroy(&ctx->mesh.meshFields_s[i]);CHKERRQ(ierr);
  }
  for (i=0;i<NUMSOLUTIONVECS_B;++i){
    ierr = DimensionalisableFieldDestroy(&ctx->solution.solutionFields_b[i]);CHKERRQ(ierr);
  }
  for (i=0;i<NUMSOLUTIONVECS_S;++i){
    ierr = DimensionalisableFieldDestroy(&ctx->solution.solutionFields_s[i]);CHKERRQ(ierr);
  }

  ierr = destroy_atmosphere( &ctx->atmosphere );CHKERRQ(ierr);

  /* Destroy parameter data */
  ierr = ParametersDestroy(&ctx->parameters);CHKERRQ(ierr);

  ierr = PetscFree2(ctx->solutionFieldIDs,ctx->solutionSlots);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->work_local_b);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_s);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_b);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_point);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_volatiles);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_reactions);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->dm_sol);CHKERRQ(ierr);

  /* Destroy PostStep data */
  if (ctx->postStepData) {
    ierr = PetscFree(ctx->postStepData);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

static PetscErrorCode CtxCreateFields(Ctx* ctx)
{
  PetscErrorCode ierr;
  ScalingConstants const SC = ctx->parameters->scaling_constants;

  /* note: the scalings include the factors of 4 pi associated with spherical
     geometry, since these are required for the output routines */

  PetscFunctionBeginUser;
  /* basic nodes */
  { // alpha
    PetscScalar scaling = 1.0 / SC->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[0],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[0],&ctx->solution.alpha); // Just for convenience - can always get this vecotr out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[0],"alpha_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[0],"K$^{-1}$");CHKERRQ(ierr);
  }
  { // cond
    PetscScalar scaling = SC->COND;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[1],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[1],&ctx->solution.cond); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[1],"cond_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[1],"W m$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  { // cp
    PetscScalar scaling = SC->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[2],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[2],&ctx->solution.cp); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[2],"cp_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[2],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  { 
    PetscScalar scaling = SC->DSDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[3],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[3],&ctx->solution.dSdxi); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[3],"dSdxi_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[3],"J kg$^{-1}$ K$^{-1}$ m$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->DTDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[4],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[4],&ctx->solution.dTdxis); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[4],"dTdxis_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[4],"K m$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->POWER * 4.0 * PETSC_PI; // total energy flow over spherical surface
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[5],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[5],&ctx->solution.Etot); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[5],"Etot_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[5],"W");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->GSUPER;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[8],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[8],&ctx->solution.gsuper); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[8],"gsuper_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[8],"K s$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[9],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[9],&ctx->solution.Jcond); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[9],"Jcond_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[9],"W m$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[10],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[10],&ctx->solution.Jconv); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[10],"Jconv_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[10],"W m$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[11],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[11],&ctx->solution.Jgrav); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[11],"Jgrav_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[11],"W m$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[12],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[12],&ctx->solution.Jmix); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[12],"Jmix_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[12],"W m$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[13],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[13],&ctx->solution.Jtot); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[13],"Jtot_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[13],"W m$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->KAPPA;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[14],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[14],&ctx->solution.kappah); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[14],"kappah_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[14],"m$^2$ s$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->KAPPA;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[15],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[15],&ctx->solution.kappac); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[15],"kappah_c");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[15],"m$^2$ s$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->KAPPA;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[16],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[16],&ctx->solution.nu); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[16],"nu_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[16],"m s$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // melt fraction is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[17],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[17],&ctx->solution.phi); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[17],"phi_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[17],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // dynamic regime (subadiabatic,inviscid,viscous)
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[18],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[18],&ctx->solution.regime); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[18],"regime_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[18],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[19],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[19],&ctx->solution.rho); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[19],"rho_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[19],"kg m$^{-3}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[20],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[20],&ctx->solution.S); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[20],"S_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[20],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[6],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[6],&ctx->solution.temp); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[6],"temp_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[6],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->VISC;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[7],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[7],&ctx->solution.visc); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[7],"visc_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[7],"Pa s");CHKERRQ(ierr);
  }
  /* staggered nodes */
  {
    PetscScalar scaling = SC->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[0],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[0],&ctx->solution.cp_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[0],"cp_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[0],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->ENTROPY / SC->TIME;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[1],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[1],&ctx->solution.dSdt_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[1],"dSdt_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[1],"Various units");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->SENERGY / SC->TIME; // W/kg
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[4],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[4],&ctx->solution.Hradio_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[4],"Hradio_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[4],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->SENERGY / SC->TIME; // W/kg
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[5],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[5],&ctx->solution.Htidal_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[5],"Htidal_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[5],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->SENERGY / SC->TIME; // W/kg
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[6],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[6],&ctx->solution.Htot_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[6],"Htot_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[6],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->MASS * SC->TEMP * 4.0 * PETSC_PI; // kg K
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[7],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[7],&ctx->solution.capacitance_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[7],"capacitance_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[7],"kg K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // melt fraction is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[9],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[9],&ctx->solution.phi_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[9],"phi_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[9],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[2],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[2],&ctx->solution.rho_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[2],"rho_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[2],"kg m$^{-3}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[8],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[8],&ctx->solution.S_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[8],"S_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[8],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[3],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[3],&ctx->solution.temp_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[3],"temp_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[3],"K");CHKERRQ(ierr);
  }
  /* mesh basic nodes */
  {
    PetscScalar scaling = SC->AREA * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[0],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[0],&ctx->mesh.area_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[0],"area_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[0],"m$^2$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = SC->DPDR;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[1],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[1],&ctx->mesh.dPdr_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[1],"dPdr_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[1],"Pa m$^{-1}$");CHKERRQ(ierr);
  }
  { // pressure_b
    PetscScalar scaling = SC->PRESSURE;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[2],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[2],&ctx->mesh.pressure_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[2],"pressure_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[2],"Pa");CHKERRQ(ierr);
  }
  { // radius_b
    PetscScalar scaling = SC->RADIUS;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[3],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[3],&ctx->mesh.radius_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[3],"radius_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[3],"m");CHKERRQ(ierr);
  }
  { // layer_b
    PetscScalar scaling = 1.0;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[4],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[4],&ctx->mesh.layer_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[4],"layer_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[4],"None");CHKERRQ(ierr);
  }
  { // xi_b (mass coordinate)
    PetscScalar scaling = SC->RADIUS;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[5],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[5],&ctx->mesh.xi_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[5],"xi_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[5],"m");CHKERRQ(ierr);
  }
  { // dxidr_b (derivative of mass coordinate)
    PetscScalar scaling = 1.0;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[6],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[6],&ctx->mesh.dxidr_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[6],"dxidr_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[6],"None");CHKERRQ(ierr);
  }
  /* mesh staggered nodes */
  { // pressure_s
    PetscScalar scaling = SC->PRESSURE;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[0],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[0],&ctx->mesh.pressure_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[0],"pressure_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[0],"Pa");CHKERRQ(ierr);
  }
  { // radius_s
    PetscScalar scaling = SC->RADIUS;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[1],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[1],&ctx->mesh.radius_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[1],"radius_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[1],"m");CHKERRQ(ierr);
  }
  { // volume_s
    PetscScalar scaling = SC->VOLUME * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[2],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[2],&ctx->mesh.volume_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[2],"volume_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[2],"m$^3$");CHKERRQ(ierr);
  }
  { // dPdr_s
    PetscScalar scaling = SC->DPDR;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[3],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[3],&ctx->mesh.dPdr_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[3],"dPdr_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[3],"Pa m$^{-1}$");CHKERRQ(ierr);
  }
  { // area_s
    PetscScalar scaling = SC->AREA * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[4],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[4],&ctx->mesh.area_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[4],"area_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[4],"m$^2$");CHKERRQ(ierr);
  }
  { // mass_s
    PetscScalar scaling = SC->MASS * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[5],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[5],&ctx->mesh.mass_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[5],"mass_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[5],"kg");CHKERRQ(ierr);
  }
  { // xi_s (mass coordinate)
    PetscScalar scaling = SC->RADIUS;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[6],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[6],&ctx->mesh.xi_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[6],"xi_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[6],"m");CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
