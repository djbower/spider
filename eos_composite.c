#include "eos_composite.h"
#include "util.h"

/* Prototypes for helpers called from interface functions */
static PetscErrorCode EOSEval_Composite_TwoPhase(EOS,PetscScalar,PetscScalar,EOSEvalData*);
static PetscErrorCode EOSCompositeGetTwoPhaseLiquidus(EOS,PetscScalar,PetscScalar*);
static PetscErrorCode EOSCompositeGetTwoPhaseSolidus(EOS,PetscScalar,PetscScalar*);

/* EOS Interface functions */
static PetscErrorCode EOSEval_Composite(EOS eos, PetscScalar P, PetscScalar S, EOSEvalData *eval)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  /* currently only simple two phase implemented */
  ierr = EOSEval_Composite_TwoPhase(eos, P, S, eval); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* This destroy function destroys the sub-EOSs */
static PetscErrorCode EOSDestroy_Composite(EOS eos)
{
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  /* Note that we do NOT destroy the sub-EOSs. That must be done 
     by whoever supplied them */
  ierr = PetscFree(eos->impl_data);CHKERRQ(ierr);
  eos->impl_data = NULL;
  PetscFunctionReturn(0);
}

static PetscErrorCode EOSSetUpFromOptions_Composite(EOS eos, const char *prefix, const FundamentalConstants FC, const ScalingConstants SC)
{
  PetscErrorCode  ierr;
  data_EOSComposite *composite = (data_EOSComposite*) eos->impl_data;

  PetscFunctionBegin;
  (void) prefix; // unused
  (void) FC; // unused
  (void) SC; // unused
  if (!composite->eos) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONGSTATE,"Sub-EOSs must be added to composite before setting up");
  ierr = PetscOptionsGetScalar(NULL,NULL,"-matprop_smooth_width",&composite->matprop_smooth_width,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-phi_critical",&composite->phi_critical,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(NULL,NULL,"-phi_width",&composite->phi_width,NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Creation */
PetscErrorCode EOSCreate_Composite(EOS eos) {
  PetscErrorCode    ierr;

  PetscFunctionBeginUser;
  eos->eval = EOSEval_Composite;
  eos->destroy = EOSDestroy_Composite;
  eos->setupfromoptions = EOSSetUpFromOptions_Composite;

  ierr = PetscMalloc1(1, (data_EOSComposite**) (&eos->impl_data));CHKERRQ(ierr);
  {
    data_EOSComposite *composite = (data_EOSComposite*) eos->impl_data;

    composite->matprop_smooth_width = 0.0;
    composite->phi_critical = 0.4;
    composite->phi_width = 0.15;
    composite->eos = NULL;
    composite->n_eos = 0;

    composite->melt_slot = 0;
    composite->liquidus_slot = 0;
    composite->solid_slot= 1;
    composite->solidus_slot = 1;
  }
  PetscFunctionReturn(0);
}

/* EOSComposite interface functions */
PetscErrorCode EOSCompositeGetMatpropSmoothWidth(EOS eos, PetscScalar *matprop_smooth_width)
{
  data_EOSComposite *composite = (data_EOSComposite*) eos->impl_data;

  PetscFunctionBeginUser;
  *matprop_smooth_width = composite->matprop_smooth_width;
  PetscFunctionReturn(0);
}

PetscErrorCode EOSCompositeGetTwoPhasePhaseFractionNoTruncation(EOS eos, PetscScalar P, PetscScalar S, PetscScalar *phase_fraction)
{
  PetscErrorCode     ierr;
  PetscScalar        solidus, liquidus;

  PetscFunctionBeginUser;
  ierr = EOSCompositeGetTwoPhaseSolidus(eos, P, &solidus ); CHKERRQ(ierr);
  ierr = EOSCompositeGetTwoPhaseLiquidus(eos, P, &liquidus ); CHKERRQ(ierr);
  *phase_fraction = ( S - solidus ) / (liquidus-solidus);
  PetscFunctionReturn(0);
}

PetscErrorCode EOSCompositeGetSubEOS(EOS eos, EOS **sub_eos, PetscInt *n_sub_eos)
{
  data_EOSComposite *composite = (data_EOSComposite*) eos->impl_data;

  PetscFunctionBeginUser;
  if (!composite->eos) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONGSTATE,"No sub-EOS to get");
  *sub_eos = composite->eos;
  *n_sub_eos = composite->n_eos;
  PetscFunctionReturn(0);
}

PetscErrorCode EOSCompositeSetSubEOS(EOS eos, EOS *sub_eos, PetscInt n_sub_eos)
{
  data_EOSComposite *composite = (data_EOSComposite*) eos->impl_data;

  PetscFunctionBeginUser;
  if (eos->is_setup) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONGSTATE,"Can only set sub-EOSs before setup");
  if (composite->eos) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONGSTATE,"Can only set sub-EOSs once");
  composite->eos = sub_eos;
  composite->n_eos = n_sub_eos;
  PetscFunctionReturn(0);
}

/* Helper Functions */
static PetscErrorCode EOSCompositeGetTwoPhaseLiquidus(EOS eos, PetscScalar P, PetscScalar *liquidus )
{
  PetscErrorCode     ierr;
  data_EOSComposite *composite = (data_EOSComposite*) eos->impl_data;

  PetscFunctionBeginUser;
  ierr = EOSGetPhaseBoundary(composite->eos[composite->liquidus_slot], P, liquidus, NULL );CHKERRQ(ierr); /* liquidus entropy */
  PetscFunctionReturn(0);
}

static PetscErrorCode EOSCompositeGetTwoPhaseSolidus(EOS eos, PetscScalar P, PetscScalar *solidus )
{
  PetscErrorCode     ierr;
  data_EOSComposite *composite = (data_EOSComposite*) eos->impl_data;

  PetscFunctionBeginUser;
  ierr = EOSGetPhaseBoundary(composite->eos[composite->solidus_slot], P, solidus, NULL );CHKERRQ(ierr); /* solidus entropy */
  PetscFunctionReturn(0);
}

static PetscErrorCode EOSEval_Composite_TwoPhase(EOS eos, PetscScalar P, PetscScalar S, EOSEvalData *eval) 
{
  PetscErrorCode     ierr;
  data_EOSComposite *composite = (data_EOSComposite*) eos->impl_data;
  EOSEvalData        eval2; // pure phase for blending across phase boundary
  EOSEvalData        eval_melt, eval_solid;
  PetscScalar        gphi, smth, liquidus, solidus, fwt;

  /* this function is called alot, if we have two phases.  Is it therefore better to store all 
     the EOSEvalData in Ctx?  Or isn't this really a speed issue? (prob not in comparison to the
     re-evaluation of functions as described below */

  PetscFunctionBeginUser;

  /* TODO: the functions below are clean, in the sense that each updates one material property.  This might
     be an advantage for setting up function pointers or the like.  However, a consequence of this approach
     is that the code is not optimised for speed, since many of these functions call the same functions
     to evaluate properties.  So an obvious speed enhancement is to perhaps scrap these individual functions
     and just update everything together in this function (still populating the eval struct, so the end
     result is the same). */

  eval->P = P;
  eval->S = S;

  /* these are strictly only valid for the mixed phase region, and not for general P and S
     conditions */
  /* FIXME: unsure what the best approach is here.  The following functions are highly modular,
     but I think it slows the code down a lot since many of the functions repeat the same lookups
     It would reduce the modularity, but for speed the better option would be to have an
     aggregate function that only evaluates things once.  This would be trivial to implement,
     but leaving as is for the time being until PS formalises the EosParameters and EosComposite structs */

  ierr = EOSCompositeGetTwoPhaseLiquidus(eos, P, &liquidus ); CHKERRQ(ierr);
  ierr = EOSCompositeGetTwoPhaseSolidus(eos, P, &solidus ); CHKERRQ(ierr);
  eval->fusion = liquidus - solidus;
  gphi = ( S - solidus ) / eval->fusion;
  eval->phase_fraction = gphi;

  /* truncation */
  if( eval->phase_fraction > 1.0 ){
    eval->phase_fraction = 1.0;
  }
  if( eval->phase_fraction < 0.0 ){
    eval->phase_fraction = 0.0;
  }

  /* properties along melting curves */
  ierr = EOSEval( composite->eos[composite->liquidus_slot], P, liquidus, &eval_melt );CHKERRQ(ierr);
  ierr = EOSEval( composite->eos[composite->solidus_slot], P, solidus, &eval_solid );CHKERRQ(ierr);

  /* linear temperature between liquidus and solidus */
  eval->T = eval->phase_fraction * eval_melt.T;
  eval->T += (1-eval->phase_fraction) * eval_solid.T;

  /* Cp */
  eval->Cp = eval_melt.S - eval_solid.S;
  eval->Cp /= eval_melt.T - eval_solid.T;
  eval->Cp *= eval_solid.T + 0.5 * (eval_melt.T - eval_solid.T);

  /* Rho */
  eval->rho = eval->phase_fraction * (1.0/eval_melt.rho) + (1-eval->phase_fraction) * (1.0/eval_solid.rho);
  eval->rho = 1.0/(eval->rho);

  /* Alpha */
  /* FIXME: positive for MgSiO3 since solid rho > melt rho.  But need to adjust for compositional
     effects */
  eval->alpha = (eval_solid.rho - eval_melt.rho) / (eval_melt.T - eval_solid.T) / eval->rho;

  /* dTdPs */
  eval->dTdPs = eval->alpha * eval->T / ( eval->rho * eval->Cp );

  /* Conductivity */
  eval->cond = eval->phase_fraction * eval_melt.cond;
  eval->cond += (1.0-eval->phase_fraction) * eval_solid.cond;

  /* Viscosity */
  ierr = EOSEvalSetViscosity(composite->eos[composite->liquidus_slot], &eval_melt);CHKERRQ(ierr);
  ierr = EOSEvalSetViscosity(composite->eos[composite->solidus_slot], &eval_solid);CHKERRQ(ierr);
  fwt = tanh_weight( eval->phase_fraction, composite->phi_critical, composite->phi_width );
  eval->log10visc = fwt * eval_melt.log10visc + (1.0-fwt) * eval_solid.log10visc;

  /* lookup does not know about these quantities, since they are not used by
     SPIDER, but for completeness zero them here */
  eval->Cv = 0.0;
  eval->V = 0.0;

  smth = get_smoothing( composite->matprop_smooth_width, gphi);

  /* now blend mixed phase EOS with single phase EOS across the phase boundary */
  if( gphi > 0.5 ){
    /* melt only properties */
    ierr = EOSEval( composite->eos[composite->melt_slot], P, S, &eval2 );CHKERRQ(ierr);
  }
  else{
    /* solid only properties */
    ierr = EOSEval( composite->eos[composite->solid_slot], P, S, &eval2 );CHKERRQ(ierr);
  }

  /* blend mixed phase with single phase, across phase boundary */
  eval->alpha = combine_matprop( smth, eval->alpha, eval2.alpha );
  eval->rho = combine_matprop( smth, eval->rho, eval2.rho );
  eval->T = combine_matprop( smth, eval->T, eval2.T );
  eval->Cp = combine_matprop( smth, eval->Cp, eval2.Cp );
  eval->dTdPs = combine_matprop( smth, eval->dTdPs, eval2.dTdPs );
  eval->cond = combine_matprop( smth, eval->cond, eval2.cond );
  eval->log10visc = combine_matprop( smth, eval->log10visc, eval2.log10visc );

  PetscFunctionReturn(0);
}
