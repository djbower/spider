#include "ctx.h"
#include "lookup.h"
#include "mesh.h"
#include "twophase.h"
#include "util.h"
#include "output.h"

static PetscErrorCode CtxCreateFields(Ctx* ctx);

/* Set up the Context */
PetscErrorCode SetupCtx(Ctx* ctx)
{
  PetscErrorCode   ierr;
  Parameters const *P = &ctx->parameters;

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

  /* Set up a parallel structured grid as DMComposite with two included DMDAs
     This is used to define vectors which hold the solution. The included
     DMDAs are used to create vectors which hold various properties
     living on the same primal and staggered nodes.
  */
  const PetscInt stencilWidth = 1;
  const PetscInt dof = 1;
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,P->numpts_b,dof,stencilWidth,NULL,&ctx->da_b      );CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,P->numpts_s,dof,stencilWidth,NULL,&ctx->da_s      );CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,1          ,dof,0           ,NULL,&ctx->da_surface);CHKERRQ(ierr); // single-point DMDA
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,1          ,dof,0           ,NULL,&ctx->da_mo_co2 );CHKERRQ(ierr); // single-point DMDA
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,1          ,dof,0           ,NULL,&ctx->da_mo_h2o );CHKERRQ(ierr); // single-point DMDA

  /* Create a composite DM of the basic nodes plus additional quantities.
    This allows us to create a vector to solve for all of these quantites.

     */
  ierr = DMCompositeCreate(PETSC_COMM_WORLD,&ctx->dm_sol);CHKERRQ(ierr);

  /* Define some more-descriptive names for the solution fields, for output */

  /* Note: the order in which DMCompositeAddDM() is called matters!
     This should be the ONLY place this function is called in the code.
     To try to avoid confusion, we maintain a list of strings describing each
     entry (indexed in the same way as the sub-DMs).
     This is intended to catch errors when we change the number of fields we are
     simultaneously solving for. */
  {
    PetscInt f = 0;

    ctx->numFields = 4;
    ierr = PetscMalloc1(4,&ctx->solutionFieldDescription);CHKERRQ(ierr);

    ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_b);CHKERRQ(ierr);
    ctx->solutionFieldDescription[f] = "dS/dr (basic nodes)";
    ++f;

    ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_surface);CHKERRQ(ierr);
    ctx->solutionFieldDescription[f] = "S at surface";
    ++f;

    ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_mo_co2);CHKERRQ(ierr);
    ctx->solutionFieldDescription[f] = "Magma ocean CO2 content";
    ++f;

    ierr = DMCompositeAddDM(ctx->dm_sol,(DM)ctx->da_mo_h2o);CHKERRQ(ierr);
    ctx->solutionFieldDescription[f] = "Magma ocean H20 content";
    ++f;
  }

  {
    PetscInt numFieldsCheck;
    ierr = DMCompositeGetNumberDM(ctx->dm_sol,&numFieldsCheck);CHKERRQ(ierr);
    if (numFieldsCheck != ctx->numFields) SETERRQ2(PetscObjectComm((PetscObject)ctx->dm_sol),PETSC_ERR_ARG_INCOMP,"The number of sub-DMs in ctx->dm_sol (%D) should be equal to %D",numFieldsCheck,ctx->numFields);
  }

  /* Continue to initialize context with distributed data */
  ierr = CtxCreateFields(ctx);

  /* Create a work vector */
  ierr = DMCreateLocalVector(ctx->da_b,&ctx->work_local_b);CHKERRQ(ierr);

  /* Populate vectors (initial condition is set in main.c) */
  set_mesh(ctx);
  set_d_dr( ctx );
  set_twophase(ctx);

  PetscFunctionReturn(0);
}

PetscErrorCode DestroyCtx(Ctx* ctx)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBeginUser;

  /* Destroy data allocated in Ctx */
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
  ierr = PetscFree(ctx->solutionFieldDescription);CHKERRQ(ierr);
  ierr = VecDestroy(&ctx->work_local_b);CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->qty_at_b);CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->ddr_at_b);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_s);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_b);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_surface);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_mo_co2);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_mo_h2o);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->dm_sol);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

static PetscErrorCode CtxCreateFields(Ctx* ctx)
{
  PetscErrorCode ierr;
  Constants const *C = &ctx->parameters.constants;

  PetscFunctionBeginUser;
  { // alpha
    PetscScalar scaling = 1.0 / C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[0],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[0],&ctx->solution.alpha); // Just for convenience - can always get this vecotr out when you need it
  }
  { // alpha_mix
    PetscScalar scaling = 1.0 / C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[1],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[1],&ctx->solution.alpha_mix); // Just for convenience - can always get this vector out when you need it
  }
  { // cond
    PetscScalar scaling = C->COND;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[2],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[2],&ctx->solution.cond); // Just for convenience - can always get this vector out when you need it
  }
  { // cp
    PetscScalar scaling =  C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[3],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[3],&ctx->solution.cp); // Just for convenience - can always get this vector out when you need it
  }
  { // cp_mix
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[4],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[4],&ctx->solution.cp_mix); // Just for convenience - can always get this vector out when you need it
  }
  { // dfusdr
    PetscScalar scaling = C->DSDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[5],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[5],&ctx->solution.dfusdr); // Just for convenience - can always get this vector out when you need it
  }
  { // dfustdr_temp
    PetscScalar scaling = C->DTDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[6],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[6],&ctx->solution.dfusdr_temp); // Just for convenience - can always get this vector out when you need it
  }
  { // dsdr
    PetscScalar scaling = C->DSDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[7],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[7],&ctx->solution.dSdr); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DSDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[8],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[8],&ctx->solution.dSliqdr); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DSDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[9],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[9],&ctx->solution.dSsoldr); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DTDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[10],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[10],&ctx->solution.dTdrs); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DTDR;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[11],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[11],&ctx->solution.dTdrs_mix); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->POWER * 4.0 * PETSC_PI; // total energy flow over spherical surface
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[12],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[12],&ctx->solution.Etot); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[13],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[13],&ctx->solution.fusion); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[14],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[14],&ctx->solution.fusion_curve); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[15],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[15],&ctx->solution.fusion_curve_temp); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[16],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[16],&ctx->solution.fusion_rho); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[17],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[17],&ctx->solution.fusion_temp); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = 1.0; // weight is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[18],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[18],&ctx->solution.fwtl); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = 1.0; // weight is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[19],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[19],&ctx->solution.fwts); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = 1.0; // (generalised) melt fraction is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[20],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[20],&ctx->solution.gphi); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->GSUPER;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[21],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[21],&ctx->solution.gsuper); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[22],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[22],&ctx->solution.Jcond); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[23],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[23],&ctx->solution.Jconv); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[24],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[24],&ctx->solution.Jgrav); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[25],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[25],&ctx->solution.Jmix); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->FLUX;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[26],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[26],&ctx->solution.Jtot); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->KAPPA;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[27],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[27],&ctx->solution.kappah); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[28],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[28],&ctx->solution.liquidus); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[29],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[29],&ctx->solution.liquidus_rho); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[30],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[30],&ctx->solution.liquidus_temp); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->KAPPA;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[31],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[31],&ctx->solution.nu); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = 1.0; // melt fraction is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[32],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[32],&ctx->solution.phi); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = 1.0; // dynamic regime (subadiabatic,inviscid,viscous)
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[33],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[33],&ctx->solution.regime); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[34],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[34],&ctx->solution.rho); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[35],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[35],&ctx->solution.S); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[36],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[36],&ctx->solution.solidus); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[37],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[37],&ctx->solution.solidus_rho); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[38],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[38],&ctx->solution.solidus_temp); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[39],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[39],&ctx->solution.temp); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->VISC;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_b[40],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_b[40],&ctx->solution.visc); // Just for convenience - can always get this vector out when you need it
  }

  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[0],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[0],&ctx->solution.cp_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[1],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[1],&ctx->solution.cp_mix_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->RHS; // note: C->RHS is 1.0 (see parameters.c)
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[2],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[2],&ctx->solution.dSdt_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[3],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[3],&ctx->solution.fusion_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[4],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[4],&ctx->solution.fusion_curve_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[5],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[5],&ctx->solution.fusion_curve_temp_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[6],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[6],&ctx->solution.fusion_temp_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = 1.0; // non-dimensional weight
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[7],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[7],&ctx->solution.fwtl_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = 1.0; // non-dimenisonal weight
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[8],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[8],&ctx->solution.fwts_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = 1.0; // generalised melt fraction is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[9],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[9],&ctx->solution.gphi_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->SENERGY / C->TIME; // W/kg
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[10],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[10],&ctx->solution.Hradio_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->SENERGY / C->TIME; // W/kg
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[11],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[11],&ctx->solution.Htidal_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->SENERGY / C->TIME; // W/kg
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[12],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[12],&ctx->solution.Htot_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->LHS * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[13],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[13],&ctx->solution.lhs_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[14],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[14],&ctx->solution.liquidus_rho_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[15],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[15],&ctx->solution.liquidus_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[16],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[16],&ctx->solution.liquidus_temp_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = 1.0; // melt fraction is non-dimensional
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[17],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[17],&ctx->solution.phi_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[18],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[18],&ctx->solution.rho_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[19],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[19],&ctx->solution.S_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->ENTROPY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[20],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[20],&ctx->solution.solidus_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[21],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[21],&ctx->solution.solidus_rho_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[22],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[22],&ctx->solution.solidus_temp_s); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&ctx->solution.solutionFields_s[23],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->solution.solutionFields_s[23],&ctx->solution.temp_s); // Just for convenience - can always get this vector out when you need it
  }

  {
    PetscScalar scaling = C->AREA * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[0],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[0],&ctx->mesh.area_b); // Just for convenience - can always get this vector out when you need it
  }
  {
    PetscScalar scaling = C->DPDR;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[1],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[1],&ctx->mesh.dPdr_b); // Just for convenience - can always get this vector out when you need it
  }
  { // pressure_b
    PetscScalar scaling = C->PRESSURE;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[2],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[2],&ctx->mesh.pressure_b); // Just for convenience - can always get this vector out when you need it
  }
  { // radius_b
    PetscScalar scaling = C->RADIUS;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[3],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[3],&ctx->mesh.radius_b); // Just for convenience - can always get this vector out when you need it
  }
  { // mix_b
    PetscScalar scaling = C->RADIUS;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_b[4],ctx->da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_b[4],&ctx->mesh.mix_b); // Just for convenience - can always get this vector out when you need it
  }

  { // pressure_s
    PetscScalar scaling = C->PRESSURE;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[0],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[0],&ctx->mesh.pressure_s); // Just for convenience - can always get this vector out when you need it
  }
  { // radius_s
    PetscScalar scaling = C->RADIUS;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[1],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[1],&ctx->mesh.radius_s); // Just for convenience - can always get this vector out when you need it
  }
  { // volume_s
    PetscScalar scaling = C->VOLUME * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[2],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[2],&ctx->mesh.volume_s); // Just for convenience - can always get this vector out when you need it
  }
  { // dPdr_s
    PetscScalar scaling = C->DPDR;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[3],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[3],&ctx->mesh.dPdr_s); // Just for convenience - can always get this vector out when you need it
  }
  { // area_s
    PetscScalar scaling = C->AREA * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[4],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[4],&ctx->mesh.area_s); // Just for convenience - can always get this vector out when you need it
  }
  { // rho_s
    PetscScalar scaling = C->DENSITY;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[5],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[5],&ctx->mesh.rho_s); // Just for convenience - can always get this vector out when you need it
  }
  { // mass_s
    PetscScalar scaling = C->MASS * 4.0 * PETSC_PI;
    ierr = DimensionalisableFieldCreate(&ctx->mesh.meshFields_s[6],ctx->da_s,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(ctx->mesh.meshFields_s[6],&ctx->mesh.mass_s); // Just for convenience - can always get this vector out when you need it
  }
  PetscFunctionReturn(0);
}
