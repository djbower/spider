#include "eos_composite.h"
#include "util.h"

/* Prototypes for helpers called from interface functions */
static PetscErrorCode EOSEval_Composite_TwoPhase(EOS,PetscScalar,PetscScalar,EosEval*);
static PetscErrorCode EOSCompositeGetTwoPhaseLiquidus(EOS,PetscScalar,PetscScalar*);
static PetscErrorCode EOSCompositeGetTwoPhaseSolidus(EOS,PetscScalar,PetscScalar*);

/* EOS Interface functions */
static PetscErrorCode EOSEval_Composite(EOS eos, PetscScalar P, PetscScalar S, EosEval *eval)
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
  }
  PetscFunctionReturn(0);
}

/* EOSComposite interface functions */
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
  if (composite->eos) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONGSTATE,"No sub-EOS to get");
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
  // TODO bad - hard-coded slot
  ierr = EOSGetPhaseBoundary(composite->eos[0], P, liquidus, NULL );CHKERRQ(ierr); /* liquidus entropy */
  PetscFunctionReturn(0);
}

static PetscErrorCode EOSCompositeGetTwoPhaseSolidus(EOS eos, PetscScalar P, PetscScalar *solidus )
{
  PetscErrorCode     ierr;
  data_EOSComposite *composite = (data_EOSComposite*) eos->impl_data;

  PetscFunctionBeginUser;
  // TODO bad - hard-coded slot
  ierr = EOSGetPhaseBoundary(composite->eos[1], P, solidus, NULL );CHKERRQ(ierr); /* solidus entropy */
  PetscFunctionReturn(0);
}

static PetscErrorCode EOSEval_Composite_TwoPhase(EOS eos, PetscScalar P, PetscScalar S, EosEval *eval) 
{
  PetscErrorCode     ierr;
  data_EOSComposite *composite = (data_EOSComposite*) eos->impl_data;
  EosEval            eval2; // pure phase for blending across phase boundary
  EosEval            eval_melt, eval_solid;
  PetscScalar        gphi, smth, liquidus, solidus, fwt;

  /* this function is called alot, if we have two phases.  Is it therefore better to store all 
     the EosEval in Ctx?  Or isn't this really a speed issue? (prob not in comparison to the
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
  // TODO bad - hard-coded slots
  ierr = EOSEval( composite->eos[0], P, liquidus, &eval_melt );CHKERRQ(ierr);
  ierr = EOSEval( composite->eos[1], P, solidus, &eval_solid );CHKERRQ(ierr);

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
  // TODO bad - hard-coded slots
  ierr = EOSEvalSetViscosity(composite->eos[0], &eval_melt);CHKERRQ(ierr);
  ierr = EOSEvalSetViscosity(composite->eos[1], &eval_solid);CHKERRQ(ierr);
  fwt = tanh_weight( eval->phase_fraction, composite->phi_critical, composite->phi_width );
  eval->log10visc = fwt * eval_melt.log10visc + (1.0-fwt) * eval_solid.log10visc;

  /* lookup does not know about these quantities, since they are not used by
     SPIDER, but for completeness zero them here */
  eval->Cv = 0.0;
  eval->V = 0.0;

  smth = get_smoothing( composite->matprop_smooth_width, gphi);

  /* now blend mixed phase EOS with single phase EOS across the phase boundary */
  // TODO bad - hard-coded slots
  if( gphi > 0.5 ){
    /* melt only properties */
    ierr = EOSEval( composite->eos[0], P, S, &eval2 );CHKERRQ(ierr);
  }
  else{
    /* solid only properties */
    ierr = EOSEval( composite->eos[1], P, S, &eval2 );CHKERRQ(ierr);
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


// TODO here, be especially careful that all code below, and relevant code in parameters.c relating to setting from options, actually made it into the new class

// TODO ---- remove the old functions below, once new EOS class is ready ----------

PetscErrorCode EosCompositeCreateTwoPhase( EosComposite *eos_composite_ptr, const EosParameters eos_parameters[], PetscInt n_phases )
{
    PetscErrorCode ierr;
    PetscInt i,j;
    char *composite_phase_names[SPIDER_MAX_PHASES];
    PetscInt n_composite_phases = SPIDER_MAX_PHASES;
    PetscBool set,flg;
    EosComposite eos_composite;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(1,eos_composite_ptr);CHKERRQ(ierr);
    eos_composite = *eos_composite_ptr;
    eos_composite->prefix = "twophase";

    ierr = PetscOptionsGetStringArray(NULL,NULL,"-eos_composite_two_phase_names",composite_phase_names,&n_composite_phases,&set);CHKERRQ(ierr);

    /* must only be two phases selected */
    if (n_composite_phases!=2) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"-eos_composite_two_phase_names only supports 2 phases (currently %d)",n_composite_phases);

    for(j=0; j<n_composite_phases; ++j){
        for(i=0; i<n_phases; ++i) {
            ierr = PetscStrcmp(eos_parameters[i]->prefix,composite_phase_names[j],&flg);CHKERRQ(ierr);
            if(flg){
                /* FIXME: this relies on the user specifying the liquidus phase first in the list, and the
                   solidus phase second.  The input file notes this should be the case, but still a possibility
                   for bugs to be introduced here. */
                eos_composite->eos_parameters[j] = eos_parameters[i];
                break;
            }
        }
        ierr = PetscFree(composite_phase_names[j]);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode EosCompositeDestroy( EosComposite *eos_composite_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscFree(*eos_composite_ptr);CHKERRQ(ierr);
    *eos_composite_ptr = NULL;
    PetscFunctionReturn(0);
}

static PetscErrorCode SetTwoPhaseLiquidus( const EosComposite eos_composite, PetscScalar P, PetscScalar *liquidus )
{
    PetscFunctionBeginUser;
    SetPhaseBoundary( eos_composite->eos_parameters[0], P, liquidus, NULL ); /* liquidus entropy */
    PetscFunctionReturn(0);
}

static PetscErrorCode SetTwoPhaseSolidus( const EosComposite eos_composite, PetscScalar P, PetscScalar *solidus )
{
    PetscFunctionBeginUser;
    SetPhaseBoundary( eos_composite->eos_parameters[1], P, solidus, NULL ); /* solidus entropy */
    PetscFunctionReturn(0);
}

PetscErrorCode SetTwoPhasePhaseFractionNoTruncation( const EosComposite eos_composite, PetscScalar P, PetscScalar S, PetscScalar *phase_fraction )
{
    PetscErrorCode ierr;
    PetscScalar solidus, liquidus;

    PetscFunctionBeginUser;

    ierr = SetTwoPhaseSolidus( eos_composite, P, &solidus ); CHKERRQ(ierr);
    ierr = SetTwoPhaseLiquidus( eos_composite, P, &liquidus ); CHKERRQ(ierr);

    *phase_fraction = ( S - solidus ) / (liquidus-solidus);

    PetscFunctionReturn(0);
}

static PetscErrorCode SetEosCompositeEvalFromTwoPhase( const EosComposite eos_composite, PetscScalar P, PetscScalar S, EosEval *eos_eval)
{
    PetscErrorCode ierr;
    /* this function is called alot, if we have two phases.  Is it therefore better to store all 
       the EosEval in Ctx?  Or isn't this really a speed issue? (prob not in comparison to the
       re-evaluation of functions as described below */
    EosEval eos_eval2; // pure phase for blending across phase boundary
    PetscScalar gphi, smth, liquidus, solidus, fwt;
    EosEval eos_eval_melt, eos_eval_solid;

    PetscFunctionBeginUser;

    /* TODO: the functions below are clean, in the sense that each updates one material property.  This might
       be an advantage for setting up function pointers or the like.  However, a consequence of this approach
       is that the code is not optimised for speed, since many of these functions call the same functions
       to evaluate properties.  So an obvious speed enhancement is to perhaps scrap these individual functions
       and just update everything together in this function (still populating the eos_eval struct, so the end
       result is the same). */

    eos_eval->P = P;
    eos_eval->S = S;

    /* these are strictly only valid for the mixed phase region, and not for general P and S
       conditions */
    /* FIXME: unsure what the best approach is here.  The following functions are highly modular,
       but I think it slows the code down a lot since many of the functions repeat the same lookups
       It would reduce the modularity, but for speed the better option would be to have an
       aggregate function that only evaluates things once.  This would be trivial to implement,
       but leaving as is for the time being until PS formalises the EosParameters and EosComposite structs */

    ierr = SetTwoPhaseLiquidus( eos_composite, P, &liquidus ); CHKERRQ(ierr);
    ierr = SetTwoPhaseSolidus( eos_composite, P, &solidus ); CHKERRQ(ierr);
    eos_eval->fusion = liquidus - solidus;
    gphi = ( S - solidus ) / eos_eval->fusion;
    eos_eval->phase_fraction = gphi;

    /* truncation */
    if( eos_eval->phase_fraction > 1.0 ){
        eos_eval->phase_fraction = 1.0;
    }
    if( eos_eval->phase_fraction < 0.0 ){
        eos_eval->phase_fraction = 0.0;
    }

    /* properties along melting curves */
    ierr = SetEosEval( eos_composite->eos_parameters[0], P, liquidus, &eos_eval_melt );CHKERRQ(ierr);
    ierr = SetEosEval( eos_composite->eos_parameters[1], P, solidus, &eos_eval_solid );CHKERRQ(ierr);

    /* linear temperature between liquidus and solidus */
    eos_eval->T = eos_eval->phase_fraction * eos_eval_melt.T;
    eos_eval->T += (1-eos_eval->phase_fraction) * eos_eval_solid.T;

    /* Cp */
    eos_eval->Cp = eos_eval_melt.S - eos_eval_solid.S;
    eos_eval->Cp /= eos_eval_melt.T - eos_eval_solid.T;
    eos_eval->Cp *= eos_eval_solid.T + 0.5 * (eos_eval_melt.T - eos_eval_solid.T);

    /* Rho */
    eos_eval->rho = eos_eval->phase_fraction * (1.0/eos_eval_melt.rho) + (1-eos_eval->phase_fraction) * (1.0/eos_eval_solid.rho);
    eos_eval->rho = 1.0/(eos_eval->rho);

    /* Alpha */
    /* FIXME: positive for MgSiO3 since solid rho > melt rho.  But need to adjust for compositional
       effects */
    eos_eval->alpha = (eos_eval_solid.rho - eos_eval_melt.rho) / (eos_eval_melt.T - eos_eval_solid.T) / eos_eval->rho;

    /* dTdPs */
    eos_eval->dTdPs = eos_eval->alpha * eos_eval->T / ( eos_eval->rho * eos_eval->Cp );

    /* Conductivity */
    eos_eval->cond = eos_eval->phase_fraction * eos_eval_melt.cond;
    eos_eval->cond += (1.0-eos_eval->phase_fraction) * eos_eval_solid.cond;

    /* Viscosity */
    ierr = SetEosEvalViscosity( eos_composite->eos_parameters[0], &eos_eval_melt );CHKERRQ(ierr);
    ierr = SetEosEvalViscosity( eos_composite->eos_parameters[1], &eos_eval_solid );CHKERRQ(ierr);
    fwt = tanh_weight( eos_eval->phase_fraction, eos_composite->phi_critical, eos_composite->phi_width );
    eos_eval->log10visc = fwt * eos_eval_melt.log10visc + (1.0-fwt) * eos_eval_solid.log10visc;

    /* lookup does not know about these quantities, since they are not used by
       SPIDER, but for completeness zero them here */
    eos_eval->Cv = 0.0;
    eos_eval->V = 0.0;

    smth = get_smoothing( eos_composite->matprop_smooth_width, gphi);

    /* now blend mixed phase EOS with single phase EOS across the phase boundary */
    if( gphi > 0.5 ){
        /* melt only properties */
        ierr = SetEosEval( eos_composite->eos_parameters[0], P, S, &eos_eval2 );CHKERRQ(ierr);
    }
    else{
        /* solid only properties */
        ierr = SetEosEval( eos_composite->eos_parameters[1], P, S, &eos_eval2 );CHKERRQ(ierr);
    }

    /* blend mixed phase with single phase, across phase boundary */
    eos_eval->alpha = combine_matprop( smth, eos_eval->alpha, eos_eval2.alpha );
    eos_eval->rho = combine_matprop( smth, eos_eval->rho, eos_eval2.rho );
    eos_eval->T = combine_matprop( smth, eos_eval->T, eos_eval2.T );
    eos_eval->Cp = combine_matprop( smth, eos_eval->Cp, eos_eval2.Cp );
    eos_eval->dTdPs = combine_matprop( smth, eos_eval->dTdPs, eos_eval2.dTdPs );
    eos_eval->cond = combine_matprop( smth, eos_eval->cond, eos_eval2.cond );
    eos_eval->log10visc = combine_matprop( smth, eos_eval->log10visc, eos_eval2.log10visc );

    PetscFunctionReturn(0);

}

PetscErrorCode SetEosCompositeEval( const EosComposite eos_composite, PetscScalar P, PetscScalar S, EosEval *eos_eval )
{

  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  /* currently only simple two phase implemented */
  ierr = SetEosCompositeEvalFromTwoPhase( eos_composite, P, S, eos_eval ); CHKERRQ(ierr);

  PetscFunctionReturn(0);

}
