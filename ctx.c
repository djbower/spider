#include "ctx.h"
#include "eos.h"
#include "mesh.h"
#include "twophase.h"
#include "util.h"
// FIXME
//#include "composition.h"

static PetscErrorCode CtxCreateFields(Ctx* ctx);

/* Set up the Context */
PetscErrorCode SetupCtx(Ctx* ctx)
{
  PetscErrorCode             ierr;
  Parameters const           *P = &ctx->parameters;
  AtmosphereParameters const *Ap = &P->atmosphere_parameters;

  PetscFunctionBeginUser;

  /*
   **All** parameters/settings should be inside ctx->parameters . Here, we
     1. Initialize these to default values
     2. Obtain any inputs from the user at the command line and update
     3. Print out all parameters
     After this point, no custom command line options should be processed,
     and ctx->parameters should never change (meaning that we can and should use
     pointers-to-constants to refer to ctx->parameters).
    */
  ierr = InitializeParametersAndSetFromOptions(&ctx->parameters);CHKERRQ(ierr); /* Note we use ctx->parameters, not P! */
  ierr = PrintParameters(P);CHKERRQ(ierr);

  /* Set up a parallel structured grid as DMComposite with several included DMDAs
     This is used to define vectors which hold the solution. The included
     DMDAs are used to create vectors which hold various properties
     living on primal and staggered nodes as well as on auxiliary domains.
  */
  const PetscInt stencilWidth = 1;
  const PetscInt dof = 1;
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,P->numpts_b,    dof,stencilWidth,NULL,&ctx->da_b        );CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,P->numpts_s,    dof,stencilWidth,NULL,&ctx->da_s        );CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,1,              dof,0,           NULL,&ctx->da_point    );CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,Ap->n_volatiles,dof,0           ,NULL,&ctx->da_volatiles);CHKERRQ(ierr); /* A collection of 0-dimensional bulk quantities */
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,Ap->n_reactions,dof,0           ,NULL,&ctx->da_reactions);CHKERRQ(ierr); /* A collection of 0-dimensional bulk quantities */

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
    Constants const *C = &ctx->parameters.constants;
    PetscInt i,f;
    PetscScalar *sol_scalings;

    /* We (over)allocate assuming a maximum of SPIDER_NUM_FIELD_IDS fields, that is that
       you will never have more than one copy of the same field.
       Note that ctx->solutionSlots has an extra entry to account for the extra,
       "undefined" entry in SpiderSolutionFieldID */
    ierr = PetscMalloc2(SPIDER_NUM_FIELD_IDS,&ctx->solutionFieldIDs,SPIDER_NUM_FIELD_IDS+1,&ctx->solutionSlots);CHKERRQ(ierr);
    ierr = PetscMalloc1(SPIDER_NUM_FIELD_IDS,&sol_scalings);CHKERRQ(ierr);
    for (i=0; i<SPIDER_NUM_FIELD_IDS; ++i) ctx->solutionSlots[i] = -1;
    f = 0;

    ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_b);CHKERRQ(ierr);
    ctx->solutionFieldIDs[f] = SPIDER_SOLUTION_FIELD_DSDR_B;
    ctx->solutionSlots[ctx->solutionFieldIDs[f]] = f;
    sol_scalings[f] = C->DSDR;
    ++f;

    ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_point);CHKERRQ(ierr);
    ctx->solutionFieldIDs[f] = SPIDER_SOLUTION_FIELD_S0;
    ctx->solutionSlots[ctx->solutionFieldIDs[f]] = f;
    sol_scalings[f] = C->ENTROPY;
    ++f;

    ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_volatiles);CHKERRQ(ierr);
    ctx->solutionFieldIDs[f] = SPIDER_SOLUTION_FIELD_MO_VOLATILES;
    ctx->solutionSlots[ctx->solutionFieldIDs[f]] = f;
    sol_scalings[f] = C->PRESSURE;
    ++f;

    ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_reactions);CHKERRQ(ierr);
    ctx->solutionFieldIDs[f] = SPIDER_SOLUTION_FIELD_MO_REACTIONS;
    ctx->solutionSlots[ctx->solutionFieldIDs[f]] = f;
    sol_scalings[f] = C->VOLATILE * 4.0 * PETSC_PI * C->MASS; // physical reservoir mass scaling
    ++f;

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
    ierr = PetscFree(sol_scalings);CHKERRQ(ierr);
  }

  /* Continue to initialize context with distributed data */
  ierr = CtxCreateFields(ctx);CHKERRQ(ierr);

  /* Create a work vector */
  ierr = DMCreateLocalVector(ctx->da_b,&ctx->work_local_b);CHKERRQ(ierr);

  /* Populate vectors (initial condition is set in main.c) */
  set_mesh(ctx);
  set_d_dr( ctx );
  set_twophase(ctx);

  ierr = initialise_atmosphere( &ctx->atmosphere, &ctx->parameters.atmosphere_parameters, &ctx->parameters.constants );CHKERRQ(ierr);

  // FIXME
  //if(P->COMPOSITION){
  //    initialise_composition(ctx);
 // }

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
  Parameters *P = &ctx->parameters;

  PetscFunctionBeginUser;

  /* destroy malloc arrays in interp objects */
  /* interp1d */
  /* TODO: if we have different options for the melting curves (e.g. analytic)
     then we may not allocate memory (and thus should not try to deallocate
     here) */
  Interp1dDestroy( &P->eos2_parameters.lookup.liquidus );
  Interp1dDestroy( &P->eos2_parameters.lookup.solidus );
  Interp1dDestroy( &P->eos1_parameters.lookup.liquidus );
  Interp1dDestroy( &P->eos1_parameters.lookup.solidus );

  /* interp2d if allocated */
  if( P->SOLID_EOS == 1 ){
      EosParametersInterp2dDestroy( &P->eos2_parameters );
  }

  /* interp2d if allocated */
  if( P->MELT_EOS == 1 ){
      EosParametersInterp2dDestroy( &P->eos1_parameters );
  }

  /* Destroy data allocated in Ctx */
  ierr = DimensionalisableFieldDestroy(&ctx->solDF);CHKERRQ(ierr);
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
  ierr = ParametersDestroy(P);CHKERRQ(ierr);

  ierr = PetscFree2(ctx->solutionFieldIDs,ctx->solutionSlots);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->work_local_b);CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->qty_at_b);CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->ddr_at_b);CHKERRQ(ierr);
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
  Constants const *C = &ctx->parameters.constants;

  PetscFunctionBeginUser;
  /* basic nodes */
  { // alpha
    PetscScalar scaling = 1.0 / C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[0],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[0],&ctx->solution.alpha); // Just for convenience - can always get this vecotr out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[0],"alpha_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[0],"K$^{-1}$");CHKERRQ(ierr);
  }
  { // alpha_mix
    PetscScalar scaling = 1.0 / C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[1],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[1],&ctx->solution.alpha_mix); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[1],"alpha_mix_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[1],"K$^{-1}$");CHKERRQ(ierr);
  }
  { // cond
    PetscScalar scaling = C->COND;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[2],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[2],&ctx->solution.cond); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[2],"cond_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[2],"W m$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  { // cp
    PetscScalar scaling =  C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[3],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[3],&ctx->solution.cp); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[3],"cp_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[3],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  { // cp_mix
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[4],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[4],&ctx->solution.cp_mix); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[4],"cp_mix_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[4],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  { // dfusdr
    PetscScalar scaling = C->DSDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[5],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[5],&ctx->solution.dfusdr); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[5],"dfusdr_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[5],"J kg$^{-1}$ K$^{-1}$ m$^{-1}$");CHKERRQ(ierr);
  }
  { // dfusdr_temp
    PetscScalar scaling = C->DTDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[6],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[6],&ctx->solution.dfusdr_temp); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[6],"dfusdr_temp_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[6],"K m$^{-1}$");CHKERRQ(ierr);
  }
  { // dsdr
    PetscScalar scaling = C->DSDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[7],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[7],&ctx->solution.dSdr); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[7],"dSdr_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[7],"J kg$^{-1}$ K$^{-1}$ m$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DSDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[8],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[8],&ctx->solution.dSliqdr); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[8],"dSliqdr_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[8],"J kg$^{-1}$ K$^{-1}$ m$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DSDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[9],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[9],&ctx->solution.dSsoldr); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[9],"dSsoldr_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[9],"J kg$^{-1}$ K$^{-1}$ m$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DTDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[10],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[10],&ctx->solution.dTdrs); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[10],"dTdrs_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[10],"K m$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DTDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[11],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[11],&ctx->solution.dTdrs_mix); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[11],"dTdrs_mix_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[11],"K m$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->POWER * 4.0 * PETSC_PI; // total energy flow over spherical surface
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[12],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[12],&ctx->solution.Etot); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[12],"Etot_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[12],"W$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[13],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[13],&ctx->solution.fusion); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[13],"fusion_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[13],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[14],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[14],&ctx->solution.fusion_curve); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[14],"fusion_curve_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[14],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[15],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[15],&ctx->solution.fusion_curve_temp); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[15],"fusion_curve_temp_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[15],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[16],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[16],&ctx->solution.fusion_rho); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[16],"fusion_rho_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[16],"kg m$^{-3}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[17],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[17],&ctx->solution.fusion_temp); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[17],"fusion_temp_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[17],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // weight is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[18],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[18],&ctx->solution.fwtl); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[18],"fwtl_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[18],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // weight is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[19],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[19],&ctx->solution.fwts); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[19],"fwts_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[19],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // (generalised) melt fraction is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[20],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[20],&ctx->solution.gphi); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[20],"gphi_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[20],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->GSUPER;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[21],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[21],&ctx->solution.gsuper); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[21],"gsuper_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[21],"K s$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[22],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[22],&ctx->solution.Jcond); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[22],"Jcond_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[22],"W m$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[23],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[23],&ctx->solution.Jconv); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[23],"Jconv_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[23],"W m$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[24],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[24],&ctx->solution.Jgrav); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[24],"Jgrav_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[24],"W m$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[25],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[25],&ctx->solution.Jmix); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[25],"Jmix_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[25],"W m$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[26],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[26],&ctx->solution.Jtot); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[26],"Jtot_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[26],"W m$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->KAPPA;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[27],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[27],&ctx->solution.kappah); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[27],"kappah_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[27],"m$^2$ s$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->KAPPA;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[28],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[28],&ctx->solution.kappac); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[28],"kappah_c");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[28],"m$^2$ s$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[29],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[29],&ctx->solution.liquidus); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[29],"liquidus_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[29],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[30],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[30],&ctx->solution.liquidus_rho); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[30],"liquidus_rho_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[30],"kg m$^{-3}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[31],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[31],&ctx->solution.liquidus_temp); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[31],"liquidus_temp_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[31],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->KAPPA;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[32],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[32],&ctx->solution.nu); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[32],"nu_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[32],"m s$^{-2}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // melt fraction is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[33],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[33],&ctx->solution.phi); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[33],"phi_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[33],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[34],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[34],&ctx->solution.Ra); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[34],"Ra_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[34],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // dynamic regime (subadiabatic,inviscid,viscous)
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[35],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[35],&ctx->solution.regime); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[35],"regime_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[35],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[36],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[36],&ctx->solution.rho); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[36],"rho_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[36],"kg m$^{-3}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[37],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[37],&ctx->solution.S); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[37],"S_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[37],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[38],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[38],&ctx->solution.solidus); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[38],"solidus_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[38],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[39],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[39],&ctx->solution.solidus_rho); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[39],"solidus_rho_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[39],"kg m$^{-3}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[40],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[40],&ctx->solution.solidus_temp); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[40],"solidus_temp_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[40],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[41],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[41],&ctx->solution.temp); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[41],"temp_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[41],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->VISC;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[42],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[42],&ctx->solution.visc); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_b[42],"visc_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_b[42],"Pa s");CHKERRQ(ierr);
  }
  /* staggered nodes */
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[0],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[0],&ctx->solution.cp_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[0],"cp_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[0],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[1],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[1],&ctx->solution.cp_mix_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[1],"cp_mix_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[1],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->RHS; // note: C->RHS is 1.0 (see parameters.c)
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[2],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[2],&ctx->solution.dSdt_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[2],"dSdt_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[2],"Various units");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[3],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[3],&ctx->solution.fusion_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[3],"fusion_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[3],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[4],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[4],&ctx->solution.fusion_curve_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[4],"fusion_curve_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[4],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[5],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[5],&ctx->solution.fusion_curve_temp_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[5],"fusion_curve_temp_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[5],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[6],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[6],&ctx->solution.fusion_temp_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[6],"fusion_temp_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[6],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // non-dimensional weight
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[7],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[7],&ctx->solution.fwtl_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[7],"fwtl_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[7],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // non-dimenisonal weight
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[8],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[8],&ctx->solution.fwts_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[8],"fwts_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[8],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // generalised melt fraction is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[9],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[9],&ctx->solution.gphi_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[9],"gphi_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[9],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->SENERGY / C->TIME; // W/kg
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[10],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[10],&ctx->solution.Hradio_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[10],"Hradio_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[10],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->SENERGY / C->TIME; // W/kg
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[11],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[11],&ctx->solution.Htidal_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[11],"Htidal_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[11],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->SENERGY / C->TIME; // W/kg
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[12],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[12],&ctx->solution.Htot_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[12],"Htot_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[12],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->LHS * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[13],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[13],&ctx->solution.lhs_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[13],"lhs_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[13],"kg K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[14],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[14],&ctx->solution.liquidus_rho_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[14],"liquidus_rho_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[14],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[15],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[15],&ctx->solution.liquidus_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[15],"liquidus_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[15],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[16],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[16],&ctx->solution.liquidus_temp_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[16],"liquidus_temp_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[16],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = 1.0; // melt fraction is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[17],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[17],&ctx->solution.phi_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[17],"phi_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[17],"None");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[18],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[18],&ctx->solution.rho_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[18],"rho_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[18],"kg m$^{-3}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[19],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[19],&ctx->solution.S_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[19],"S_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[19],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[20],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[20],&ctx->solution.solidus_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[20],"solidus_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[20],"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[21],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[21],&ctx->solution.solidus_rho_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[21],"solidus_rho_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[21],"kg m$^{-3}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[22],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[22],&ctx->solution.solidus_temp_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[22],"solidus_temp_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[22],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[23],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[23],&ctx->solution.temp_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[23],"temp_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[23],"K");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->HEATGEN;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[24],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[24],&ctx->solution.Hal26_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[24],"Hal26_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[24],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->HEATGEN;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[25],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[25],&ctx->solution.Hk40_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[25],"Hk40_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[25],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->HEATGEN;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[26],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[26],&ctx->solution.Hfe60_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[26],"Hfe60_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[26],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->HEATGEN;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[27],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[27],&ctx->solution.Hth232_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[27],"Hth232_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[27],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->HEATGEN;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[28],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[28],&ctx->solution.Hu235_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[28],"Hu235_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[28],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->HEATGEN;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[29],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[29],&ctx->solution.Hu238_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->solution.solutionFields_s[29],"Hu238_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->solution.solutionFields_s[29],"W kg$^{-1}$");CHKERRQ(ierr);
  }
  /* mesh basic nodes */
  {
    PetscScalar scaling = C->AREA * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[0],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[0],&ctx->mesh.area_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[0],"area_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[0],"m$^2$");CHKERRQ(ierr);
  }
  {
    PetscScalar scaling = C->DPDR;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[1],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[1],&ctx->mesh.dPdr_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[1],"dPdr_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[1],"Pa m$^{-1}$");CHKERRQ(ierr);
  }
  { // pressure_b
    PetscScalar scaling = C->PRESSURE;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[2],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[2],&ctx->mesh.pressure_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[2],"pressure_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[2],"Pa");CHKERRQ(ierr);
  }
  { // radius_b
    PetscScalar scaling = C->RADIUS;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[3],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[3],&ctx->mesh.radius_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[3],"radius_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[3],"m");CHKERRQ(ierr);
  }
  { // mix_b
    PetscScalar scaling = C->RADIUS;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[4],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[4],&ctx->mesh.mix_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[4],"mix_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[4],"m");CHKERRQ(ierr);
  }
  { // layer_b
    PetscScalar scaling = 1.0;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[5],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[5],&ctx->mesh.layer_b); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_b[5],"layer_b");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_b[5],"None");CHKERRQ(ierr);
  }
  /* mesh staggered nodes */
  { // pressure_s
    PetscScalar scaling = C->PRESSURE;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[0],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[0],&ctx->mesh.pressure_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[0],"pressure_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[0],"Pa");CHKERRQ(ierr);
  }
  { // radius_s
    PetscScalar scaling = C->RADIUS;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[1],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[1],&ctx->mesh.radius_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[1],"radius_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[1],"m");CHKERRQ(ierr);
  }
  { // volume_s
    PetscScalar scaling = C->VOLUME * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[2],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[2],&ctx->mesh.volume_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[2],"volume_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[2],"m$^3$");CHKERRQ(ierr);
  }
  { // dPdr_s
    PetscScalar scaling = C->DPDR;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[3],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[3],&ctx->mesh.dPdr_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[3],"dPdr_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[3],"Pa m$^{-1}$");CHKERRQ(ierr);
  }
  { // area_s
    PetscScalar scaling = C->AREA * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[4],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[4],&ctx->mesh.area_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[4],"area_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[4],"m$^2$");CHKERRQ(ierr);
  }
  { // rho_s
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[5],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[5],&ctx->mesh.rho_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[5],"rho_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[5],"kg m$^{-3}$");CHKERRQ(ierr);
  }
  { // mass_s
    PetscScalar scaling = C->MASS * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[6],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[6],&ctx->mesh.mass_s); // Just for convenience - can always get this vector out when you need it
    ierr = DimensionalisableFieldSetName(ctx->mesh.meshFields_s[6],"mass_s");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(ctx->mesh.meshFields_s[6],"kg");CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
