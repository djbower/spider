#include "eos.h"
#include "parameters.h"
#include "matprop.h"
#include "twophase.h"
#include "util.h"

static PetscErrorCode set_matprop_staggered( Ctx * );
static PetscScalar GetModifiedMixingLength( PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar GetConstantMixingLength( PetscScalar outer_radius, PetscScalar inner_radius );
static PetscScalar GetMixingLength( const Parameters, PetscScalar );
static PetscErrorCode apply_log10visc_cutoff( Parameters const, PetscScalar * );
static PetscErrorCode GetEddyDiffusivity( const EOSEvalData, const Parameters, PetscScalar, PetscScalar, PetscScalar, PetscScalar *, PetscScalar *, PetscScalar * );

PetscErrorCode set_capacitance_staggered( Ctx *E )
{
    PetscErrorCode    ierr;
    Mesh              *M = &E->mesh;
    Solution          *S = &E->solution;

    PetscFunctionBeginUser;

    ierr = set_matprop_staggered( E ); CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->capacitance_s,S->temp_s,S->rho_s); CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->capacitance_s,S->capacitance_s,M->volume_s); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_phase_fraction_staggered( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    DM                da_s=E->da_s;
    Mesh              *M = &E->mesh;
    Parameters        P = E->parameters;
    Solution          *S = &E->solution;
    PetscScalar       *arr_phi, *arr_S, *arr_pres, PP, SS;
    EOSEvalData       eos_eval;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = DMDAVecGetArrayRead(da_s,M->pressure_s,&arr_pres);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->S_s,&arr_S);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->phi_s,&arr_phi);CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){
        PP = arr_pres[i];
        SS = arr_S[i];
        ierr = EOSEval(P->eos, PP, SS, &eos_eval );CHKERRQ(ierr);
        arr_phi[i] = eos_eval.phase_fraction;
    }

    ierr = DMDAVecRestoreArrayRead(da_s,M->pressure_s,&arr_pres);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->S_s,&arr_S);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->phi_s,&arr_phi);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_matprop_staggered( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    DM                da_s=E->da_s;
    Mesh              *M = &E->mesh;
    Parameters const  P = E->parameters;
    Solution          *S = &E->solution;
    Vec               pres_s = M->pressure_s;
    PetscScalar       *arr_rho_s, *arr_temp_s, *arr_cp_s;
    const PetscScalar *arr_pres_s, *arr_S_s;
    EOSEvalData       eos_eval;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->cp_s,&arr_cp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){
        ierr = EOSEval( P->eos, arr_pres_s[i], arr_S_s[i], &eos_eval );CHKERRQ(ierr);
        arr_rho_s[i] = eos_eval.rho;
        arr_temp_s[i] = eos_eval.T;
        arr_cp_s[i] = eos_eval.Cp;
    }

    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->cp_s,&arr_cp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_matprop_basic( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_b,ihi_b,w_b,numpts_b;
    DM                da_b=E->da_b;
    PetscScalar       *arr_Ra, *arr_phi, *arr_nu, *arr_gsuper, *arr_kappac, *arr_kappah, *arr_dTdxis, *arr_alpha, *arr_temp, *arr_cp, *arr_cond, *arr_visc, *arr_regime, *arr_rho;
    const PetscScalar *arr_dSdxi, *arr_S_b, *arr_pres, *arr_dPdr_b, *arr_radius_b, *arr_dxidr_b;
    const PetscInt    *arr_layer_b;
    Mesh              *M = &E->mesh;
    Parameters const  P = E->parameters;
    Solution          *S = &E->solution;
    EOSEvalData       eos_eval;

    PetscFunctionBeginUser;

    /* loop over all basic nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    ierr = DMDAVecGetArrayRead(da_b,S->dSdxi,&arr_dSdxi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->S,&arr_S_b); CHKERRQ(ierr);
    /* mesh quantities */
    ierr = DMDAVecGetArrayRead(da_b,M->dPdr_b,&arr_dPdr_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->pressure_b,&arr_pres); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->layer_b,&arr_layer_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->radius_b,&arr_radius_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->dxidr_b,&arr_dxidr_b); CHKERRQ(ierr);
    /* material properties */
    ierr = DMDAVecGetArray(    da_b,S->alpha,&arr_alpha); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cond,&arr_cond); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cp,&arr_cp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->dTdxis,&arr_dTdxis); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->gsuper,&arr_gsuper); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->kappac,&arr_kappac); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->kappah,&arr_kappah); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->nu,&arr_nu); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->phi,&arr_phi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Ra,&arr_Ra); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->rho,&arr_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->temp,&arr_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->visc,&arr_visc); CHKERRQ(ierr);
    /* regime: not convecting (0), inviscid (1), viscous (2) */
    ierr = DMDAVecGetArray(    da_b,S->regime,&arr_regime); CHKERRQ(ierr);

    /* loop over all basic nodes */ 
    for(i=ilo_b; i<ihi_b; ++i){
      ierr = EOSEval( P->eos, arr_pres[i], arr_S_b[i], &eos_eval );CHKERRQ(ierr);
      arr_phi[i] = eos_eval.phase_fraction;
      arr_rho[i] = eos_eval.rho;
      arr_dTdxis[i] = arr_dPdr_b[i] * eos_eval.dTdPs / arr_dxidr_b[i];
      arr_cp[i] = eos_eval.Cp;
      arr_temp[i] = eos_eval.T;
      arr_alpha[i] = eos_eval.alpha;
      arr_cond[i] = eos_eval.cond;
      arr_visc[i] = eos_eval.log10visc;

      /* apply viscosity cutoff */
      ierr = apply_log10visc_cutoff( P, &arr_visc[i] );
      arr_visc[i] = PetscPowScalar( 10.0, arr_visc[i] );

      /* below are computed in GetEddyDiffusivity, but not yet stored to the arrays */
      /* kinematic viscosity */
      arr_nu[i] = arr_visc[i] / arr_rho[i];
      /* gravity * super-adiabatic temperature gradient */
      arr_gsuper[i] = P->gravity * arr_temp[i] / arr_cp[i] * arr_dSdxi[i] * arr_dxidr_b[i];

      ierr = GetEddyDiffusivity( eos_eval, P, arr_radius_b[i], arr_dSdxi[i], arr_dxidr_b[i], &arr_kappah[i], &arr_kappac[i], &arr_regime[i] );CHKERRQ(ierr);

      /* FIXME: below */
#if 0
      /* Rayleigh number */
      /* FIXME: should use domain size not mixing length */
      arr_Ra[i] = arr_gsuper[i];
      arr_Ra[i] *= arr_alpha[i] * PetscPowScalar(arr_mix_b[i],4) * arr_rho[i] * arr_cp[i];
      arr_Ra[i] /= arr_nu[i] * arr_cond[i];
#endif
      arr_Ra[i] = 0.0; // placeholder, to update

    }

    ierr = DMDAVecRestoreArrayRead(da_b,S->dSdxi,&arr_dSdxi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->S,&arr_S_b); CHKERRQ(ierr);
    /* mesh quantities */
    ierr = DMDAVecRestoreArrayRead(da_b,M->dPdr_b,&arr_dPdr_b); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->pressure_b,&arr_pres); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->layer_b,&arr_layer_b); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->dxidr_b,&arr_dxidr_b); CHKERRQ(ierr);
    /* material properties */
    ierr = DMDAVecRestoreArray(    da_b,S->alpha,&arr_alpha); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cond,&arr_cond); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cp,&arr_cp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->dTdxis,&arr_dTdxis); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->gsuper,&arr_gsuper); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->kappac,&arr_kappac); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->kappah,&arr_kappah); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->nu,&arr_nu); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->phi,&arr_phi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Ra,&arr_Ra); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->rho,&arr_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->temp,&arr_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->visc,&arr_visc); CHKERRQ(ierr);
    /* regime */
    ierr = DMDAVecRestoreArray(    da_b,S->regime,&arr_regime); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode GetEddyDiffusivity( const EOSEvalData eos_eval, const Parameters P, PetscScalar radius, PetscScalar dSdxi, PetscScalar dxidr, PetscScalar *kappah_ptr, PetscScalar *kappac_ptr, PetscScalar *regime_ptr )
{
    PetscErrorCode ierr;
    PetscScalar    visc, kvisc, gsuper, kh, crit, mix, kappah, kappac, regime;

    PetscFunctionBeginUser;

    visc = eos_eval.log10visc;
    ierr = apply_log10visc_cutoff( P, &visc );CHKERRQ(ierr);
    visc = PetscPowScalar( 10.0, visc );
    kvisc = visc / eos_eval.rho; // kinematic viscosity
    gsuper = P->gravity * eos_eval.T / eos_eval.Cp * dSdxi * dxidr; // g * super adiabatic gradient

    crit = 81.0 * PetscPowScalar(kvisc,2);
    mix = GetMixingLength( P, radius);
    crit /= 4.0 * eos_eval.alpha * PetscPowScalar(mix,4);

    if( gsuper <= 0.0 ){
      /* no convection, subadiabatic */
      kh = 0.0;
      regime = 0.0;
    } else if( gsuper > crit ){
      /* inviscid scaling from Vitense (1953) */
      kh = 0.25 * PetscPowScalar(mix,2) * PetscSqrtScalar(eos_eval.alpha * gsuper);
      regime = 2.0;
    } else{
      /* viscous scaling */
      kh = eos_eval.alpha * gsuper * PetscPowScalar(mix,4) / ( 18.0 * kvisc );
      regime = 1.0;
    }

    /* thermal eddy diffusivity */
    if (P->eddy_diffusivity_thermal > 0.0){
      /* scale */
      kappah = P->eddy_diffusivity_thermal * kh;
    }
    else{
      /* else set (and negate to account for sign flag) */
      kappah = -P->eddy_diffusivity_thermal;
    }

    /* chemical eddy diffusivity */
    if (P->eddy_diffusivity_chemical > 0.0){
      /* scale */
      kappac = P->eddy_diffusivity_chemical * kh;
    }
    else{
      /* else set (and negate to account for sign flag) */
      kappac = -P->eddy_diffusivity_chemical;
    }

    /* update pointers with data */
    if(kappah_ptr != NULL){
        *kappah_ptr = kappah;
    }
    if(kappac_ptr != NULL){
        *kappac_ptr = kappac;
    }
    if(regime_ptr != NULL){
        *regime_ptr = regime;
    }

    PetscFunctionReturn(0);
}

static PetscScalar GetModifiedMixingLength( PetscScalar a, PetscScalar b, PetscScalar outer_radius, PetscScalar inner_radius, PetscScalar radius )
{
    /* See Kamata, 2018, JGR */
    /* conventional mixing length theory has a = b = 0.5 */
    /* a is location of peak in depth/radius space,
       b is amplitude of the peak */

    PetscScalar mix_length1, mix_length2, mix_length;

    mix_length1 = (radius - inner_radius) * b / (1.0 - a);
    mix_length2 = (outer_radius - radius) * b / a;

    mix_length = PetscMin( mix_length1, mix_length2 );

    return mix_length;
}

static PetscScalar GetConstantMixingLength( PetscScalar outer_radius, PetscScalar inner_radius )
{
    PetscScalar mix_length;

    mix_length = 0.25 * (outer_radius - inner_radius );

    return mix_length;
}

static PetscScalar GetMixingLength( const Parameters P, PetscScalar radius )
{
    PetscScalar outer_radius, inner_radius;
    PetscScalar mix_length = 0.0;
    PetscScalar eps = 1.0E-10;

    /* for a single layer, P->layer_interface_radius = P->coresize (parameters.c),
       enabling this single expression to work for both a single and double
       layered mantle */
    /* due to floating points, sometimes the radius at the innermost boundary is slightly
       less than P->layer_interface_radius, so account for this with a small offset (eps) */
    /* if we try to resolve the ultra-thin boundary layer, eps might not be smaller enough? */

    if( radius >= P->radius * P->layer_interface_radius - eps ){
        outer_radius = P->radius;
        inner_radius = P->radius * P->layer_interface_radius;
    }
    else{
        outer_radius = P->radius * P->layer_interface_radius;
        inner_radius = P->radius * P->coresize;
    }

    if( P->mixing_length == 1){
        mix_length = GetModifiedMixingLength( P->mixing_length_a, P->mixing_length_b, outer_radius, inner_radius, radius );
    }
    else if( P->mixing_length == 2){
        mix_length = GetConstantMixingLength( outer_radius, inner_radius );
    }

    /* parameters.c ensures that P->mixing_length must be 1 or 2, so we cannot
       fall outside this if statement */

    return mix_length;
}

static PetscErrorCode apply_log10visc_cutoff( Parameters const P, PetscScalar *viscosity )
{
    PetscFunctionBeginUser;

    if(P->log10visc_min > 0.0){
        if(*viscosity < P->log10visc_min){
            *viscosity = P->log10visc_min;
        }
    }
    if(P->log10visc_max > 0.0){
        if(*viscosity > P->log10visc_max){
            *viscosity = P->log10visc_max;
        }
    }

    PetscFunctionReturn(0);
}
