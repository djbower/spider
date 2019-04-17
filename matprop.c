#include "matprop.h"
#include "util.h"
#include "lookup.h"
#include "twophase.h"
// FIXME
//#include "composition.h"

static PetscErrorCode set_matprop_staggered( Ctx * );
static PetscScalar get_melt_fraction_truncated( PetscScalar );
static PetscScalar get_log10_viscosity_solid( PetscScalar, PetscScalar, PetscInt, PetscScalar, Parameters const *);
static PetscScalar add_compositional_viscosity( PetscScalar, PetscScalar );
static PetscScalar get_log10_viscosity_melt( PetscScalar, PetscScalar, PetscInt, Parameters const *);
static PetscScalar get_log10_viscosity_mix( PetscScalar, PetscScalar, PetscScalar, Parameters const * );
static PetscScalar get_log10_viscosity_cutoff( PetscScalar, Parameters const * );
static PetscScalar get_viscosity_mix_no_skew( PetscScalar, Parameters const * );

PetscErrorCode set_capacitance_staggered( Ctx *E )
{
    PetscErrorCode    ierr;
    Mesh              *M = &E->mesh;
    //Parameters        *P = &E->parameters;
    Solution          *S = &E->solution;

    PetscFunctionBeginUser;

    /* useful to passively compute crystal and Bridgmanite fraction,
       even if they do not feedback into the density calculation */
    /* determine values to implement compositional differentiation */
    // FIXME
    //if(P->COMPOSITION){
    //ierr = set_composition( E ); CHKERRQ(ierr);
    //}

    ierr = set_matprop_staggered( E ); CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->lhs_s,S->temp_s,S->rho_s); CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->lhs_s,S->lhs_s,M->volume_s); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscScalar get_melt_fraction_truncated( PetscScalar phi )
{
    /* truncate [0,1] */
    if (phi > 1.0){
      /* superliquidus */
      return 1.0;
    }
    else if (phi < 0.0){
      /* subsolidus */
      return 0.0;
    }
    else
      return phi;
}

PetscErrorCode set_melt_fraction_staggered( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    DM                da_s=E->da_s;
    Parameters        *P = &E->parameters;
    Solution          *S = &E->solution;
    PetscScalar       *arr_phi_s;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    // compute melt fraction
    if(P->SOLID_CONVECTION_ONLY){
        ierr = VecSet( S->phi_s, 0.0 );CHKERRQ(ierr); // by definition
    }
    else if(P->LIQUID_CONVECTION_ONLY){
        ierr = VecSet( S->phi_s, 1.0 );CHKERRQ(ierr); // by definition
    }
    else{
        ierr = VecWAXPY(S->phi_s,-1.0,S->solidus_s,S->S_s);CHKERRQ(ierr);
        ierr = VecPointwiseDivide(S->phi_s,S->phi_s,S->fusion_s);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
        /* TODO: can we remove this loop and truncate using Petsc Vec
           operations instead? */
        for(i=ilo_s; i<ihi_s; ++i){
            /* truncate melt fraction */
            arr_phi_s[i] = get_melt_fraction_truncated( arr_phi_s[i] );
        }
        ierr = DMDAVecRestoreArray(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

static PetscErrorCode set_matprop_staggered( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    DM                da_s=E->da_s;
    Lookup const      *L;
    Mesh              *M = &E->mesh;
    Parameters const  *P = &E->parameters;
    // FIXME
    //CompositionalParameters const *Comp = &P->compositional_parameters;
    Solution          *S = &E->solution;
    Vec               pres_s = M->pressure_s;
    // material properties that are updated here
    PetscScalar       *arr_rho_s, *arr_temp_s, *arr_cp_s;
    // material properties used to update above
    const PetscScalar *arr_pres_s, *arr_liquidus_rho_s, *arr_solidus_rho_s, *arr_liquidus_temp_s, *arr_solidus_temp_s, *arr_S_s, *arr_liquidus_s, *arr_solidus_s, *arr_fusion_s, *arr_cp_mix_s, *arr_phi_s;
    // for smoothing properties across liquidus and solidus
    const PetscScalar *arr_fwtl_s, *arr_fwts_s;
    PetscScalar       fwtl, fwts;
    PetscScalar       rho_sol, temp_sol, cp_sol;
    PetscScalar       rho_mel, temp_mel, cp_mel;
    PetscScalar       rho_mix, temp_mix, cp_mix;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = DMDAVecGetArrayRead(da_s,S->cp_mix_s,&arr_cp_mix_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fusion_s,&arr_fusion_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fwtl_s,&arr_fwtl_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fwts_s,&arr_fwts_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->cp_s,&arr_cp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){

        /* solid phase */
        L = &P->solid_prop;
        rho_sol = get_val2d( &L->rho, arr_pres_s[i], arr_S_s[i] );
        temp_sol = get_val2d( &L->temp, arr_pres_s[i], arr_S_s[i] );
        cp_sol = get_val2d( &L->cp, arr_pres_s[i], arr_S_s[i] );

        /* melt phase */
        L = &P->melt_prop;
        rho_mel = get_val2d( &L->rho, arr_pres_s[i], arr_S_s[i] );
        temp_mel = get_val2d( &L->temp, arr_pres_s[i], arr_S_s[i] );
        cp_mel = get_val2d( &L->cp, arr_pres_s[i], arr_S_s[i] );

        /* mixed phase */
        /* volume additivity, excluding temperature (since phase effect dominant) */
        rho_mix = combine_matprop( arr_phi_s[i], 1.0/arr_liquidus_rho_s[i], 1.0/arr_solidus_rho_s[i] );
        rho_mix = 1.0 / rho_mix;

        /* FIXME: run this past Aaron */
        /*if(P->COMPOSITION){
            rho_mel *= Comp->mass_ratio_liquidus;
            rho_mix = arr_liquidus_rho_s[i];
            if(i < Comp->rheological_front_index){
                rho_mix *= Comp->mo_mass_ratio;
            }
            else{
                rho_mix *= Comp->mass_ratio_liquidus;
            }
        }
        else{
            rho_mix = combine_matprop( arr_phi_s[i], 1.0/arr_liquidus_rho_s[i], 1.0/arr_solidus_rho_s[i] );
            rho_mix = 1.0 / rho_mix;
        }*/

        temp_mix = combine_matprop( arr_phi_s[i], arr_liquidus_temp_s[i], arr_solidus_temp_s[i] );
        cp_mix = arr_cp_mix_s[i];

        if(P->SOLID_CONVECTION_ONLY){
            arr_rho_s[i] = rho_sol;
            arr_temp_s[i] = temp_sol;
            arr_cp_s[i] = cp_sol;
        }
        else if(P->LIQUID_CONVECTION_ONLY){
            arr_rho_s[i] = rho_mel;
            arr_temp_s[i] = temp_mel;
            arr_cp_s[i] = cp_mel;
        }
        else{
            if (arr_phi_s[i] > 0.5){
                fwtl = arr_fwtl_s[i]; // for smoothing
                arr_rho_s[i]   = combine_matprop( fwtl, rho_mel, rho_mix );
                arr_temp_s[i]  = combine_matprop( fwtl, temp_mel, temp_mix );
                arr_cp_s[i]    = combine_matprop( fwtl, cp_mel, cp_mix );
            }
            else{ 
                fwts = arr_fwts_s[i]; // for smoothing
                arr_rho_s[i]   = combine_matprop( fwts, rho_mix, rho_sol );
                arr_temp_s[i]  = combine_matprop( fwts, temp_mix, temp_sol );
                arr_cp_s[i]    = combine_matprop( fwts, cp_mix, cp_sol );
            }
        }

    }

    ierr = DMDAVecRestoreArrayRead(da_s,S->cp_mix_s,&arr_cp_mix_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fusion_s,&arr_fusion_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fwtl_s,&arr_fwtl_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fwts_s,&arr_fwts_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->cp_s,&arr_cp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_matprop_basic( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_b,ihi_b,w_b,ilo,ihi,numpts_b;
    DM                da_b=E->da_b;
    // material properties that are updated here
    PetscScalar       *arr_Ra, *arr_phi, *arr_nu, *arr_gsuper, *arr_kappac, *arr_kappah, *arr_dTdrs, *arr_alpha, *arr_temp, *arr_cp, *arr_cond, *arr_visc, *arr_regime, *arr_rho;
    // material properties used to update above
    const PetscScalar *arr_dSdr, *arr_S_b, *arr_solidus, *arr_fusion, *arr_pres, *arr_dPdr_b, *arr_liquidus, *arr_liquidus_rho, *arr_solidus_rho, *arr_cp_mix, *arr_dTdrs_mix, *arr_liquidus_temp, *arr_solidus_temp, *arr_fusion_rho, *arr_fusion_temp, *arr_mix_b, *arr_radius_b;
    const PetscInt *arr_layer_b;
    // for smoothing properties across liquidus and solidus
    const PetscScalar *arr_fwtl, *arr_fwts;
    PetscScalar       fwtl, fwts;
    PetscScalar       rho_sol, dTdrs_sol, cp_sol, temp_sol, alpha_sol, cond_sol, log10visc_sol;
    PetscScalar       rho_mel, dTdrs_mel, cp_mel, temp_mel, alpha_mel, cond_mel, log10visc_mel;
    PetscScalar       rho_mix, dTdrs_mix, cp_mix, temp_mix, alpha_mix, cond_mix, log10visc_mel_mix, log10visc_sol_mix, log10visc_mix;
    Lookup const      *L;
    Mesh              *M = &E->mesh;
    Parameters const  *P = &E->parameters;
    // FIXME
    //CompositionalParameters const *Comp = &P->compositional_parameters;
    Solution          *S = &E->solution;

    PetscFunctionBeginUser;

    /* melt fraction, not truncated (this happens below) */
    ierr = VecWAXPY(S->phi,-1.0,S->solidus,S->S);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(S->phi,S->phi,S->fusion);CHKERRQ(ierr);

    /* loop over all basic nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    /* loop over all basic nodes */
    ilo = ilo_b;
    ihi = ihi_b;
    /* restrict to basic nodes by uncommenting below */
    //ilo = ilo_b == 0        ? 1            : ilo_b;
    //ihi = ihi_b == numpts_b ? numpts_b - 1 : ihi_b;

    ierr = DMDAVecGetArrayRead(da_b,S->dSdr,&arr_dSdr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->S,&arr_S_b); CHKERRQ(ierr);
    /* mesh quantities */
    ierr = DMDAVecGetArrayRead(da_b,M->dPdr_b,&arr_dPdr_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->mix_b,&arr_mix_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->pressure_b,&arr_pres); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->layer_b,&arr_layer_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->radius_b,&arr_radius_b); CHKERRQ(ierr);
    /* material properties */
    ierr = DMDAVecGetArray(    da_b,S->alpha,&arr_alpha); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cond,&arr_cond); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cp,&arr_cp); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->cp_mix,&arr_cp_mix); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->dTdrs,&arr_dTdrs); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dTdrs_mix,&arr_dTdrs_mix); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion,&arr_fusion); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion_rho,&arr_fusion_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion_temp,&arr_fusion_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fwtl,&arr_fwtl); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fwts,&arr_fwts); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->gsuper,&arr_gsuper); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->kappac,&arr_kappac); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->kappah,&arr_kappah); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus,&arr_liquidus); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus_rho,&arr_liquidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus_temp,&arr_liquidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->nu,&arr_nu); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->phi,&arr_phi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Ra,&arr_Ra); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->rho,&arr_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus,&arr_solidus); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus_rho,&arr_solidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus_temp,&arr_solidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->temp,&arr_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->visc,&arr_visc); CHKERRQ(ierr);
    /* regime: not convecting (0), inviscid (1), viscous (2) */
    ierr = DMDAVecGetArray(    da_b,S->regime,&arr_regime); CHKERRQ(ierr);

    for(i=ilo; i<ihi; ++i){

      /* truncate melt fraction */
      arr_phi[i] = get_melt_fraction_truncated( arr_phi[i] );

      /* solid phase */
      L = &P->solid_prop;
      rho_sol = get_val2d( &L->rho, arr_pres[i], arr_S_b[i] );
      dTdrs_sol = arr_dPdr_b[i] * get_val2d( &L->dTdPs, arr_pres[i], arr_S_b[i] );
      cp_sol = get_val2d( &L->cp, arr_pres[i], arr_S_b[i] );
      temp_sol = get_val2d( &L->temp, arr_pres[i], arr_S_b[i] );
      alpha_sol = get_val2d( &L->alpha, arr_pres[i], arr_S_b[i] );
      cond_sol = P->cond_sol;
      log10visc_sol = get_log10_viscosity_solid( temp_sol, arr_pres[i], arr_layer_b[i], arr_radius_b[i], P );

      /* melt phase */
      L = &P->melt_prop;
      rho_mel = get_val2d( &L->rho, arr_pres[i], arr_S_b[i] );
      dTdrs_mel = arr_dPdr_b[i] * get_val2d( &L->dTdPs, arr_pres[i], arr_S_b[i] );
      cp_mel = get_val2d( &L->cp, arr_pres[i], arr_S_b[i] );
      temp_mel = get_val2d( &L->temp, arr_pres[i], arr_S_b[i] );
      alpha_mel = get_val2d( &L->alpha, arr_pres[i], arr_S_b[i] );
      cond_mel = P->cond_mel;
      log10visc_mel = get_log10_viscosity_melt( temp_mel, arr_pres[i], arr_layer_b[i], P );

      /* mixed phase */
      rho_mix = combine_matprop( arr_phi[i], 1.0/arr_liquidus_rho[i], 1.0/arr_solidus_rho[i] );
      rho_mix = 1.0 / rho_mix;

      /* FIXME: run past Aaron */
      /*if(P->COMPOSITION){
          rho_mel *= Comp->mass_ratio_liquidus;
          rho_mix = arr_liquidus_rho[i];
          if(i <= Comp->rheological_front_index){
              rho_mix *= Comp->mo_mass_ratio;
          }
          else{
              rho_mix *= Comp->mass_ratio_liquidus;
          }
      }
      else{
          rho_mix = combine_matprop( arr_phi[i], 1.0/arr_liquidus_rho[i], 1.0/arr_solidus_rho[i] );
          rho_mix = 1.0 / rho_mix;
      }*/

      dTdrs_mix = arr_dTdrs_mix[i];
      cp_mix = arr_cp_mix[i];
      temp_mix = combine_matprop( arr_phi[i], arr_liquidus_temp[i], arr_solidus_temp[i] );
      alpha_mix = -arr_fusion_rho[i] / arr_fusion_temp[i] / rho_mix;
      cond_mix = combine_matprop( arr_phi[i], P->cond_mel, P->cond_sol );
      /* need to get viscosity of melt and solid phases at the liquidus and solidus temperature,
         since this is consistent with the notion of ignoring temperature effects in the mixed
         phase region (e.g., for density) */
      log10visc_mel_mix = get_log10_viscosity_melt( arr_liquidus_temp[i], arr_pres[i], arr_layer_b[i], P );
      log10visc_sol_mix = get_log10_viscosity_solid( arr_solidus_temp[i], arr_pres[i], arr_layer_b[i], arr_radius_b[i], P );
      log10visc_mix = get_log10_viscosity_mix( arr_phi[i], log10visc_mel_mix, log10visc_sol_mix, P );

      if(P->SOLID_CONVECTION_ONLY){
          arr_phi[i] = 0.0; // by definition
          arr_rho[i] = rho_sol;
          arr_dTdrs[i] = dTdrs_sol;
          arr_cp[i] = cp_sol;
          arr_temp[i] = temp_sol;
          arr_alpha[i] = alpha_sol;
          arr_cond[i] = cond_sol;
          arr_visc[i] = log10visc_sol;
      }
      else if(P->LIQUID_CONVECTION_ONLY){
          arr_phi[i] = 1.0; // by definition
          arr_rho[i] = rho_mel;
          arr_dTdrs[i] = dTdrs_mel;
          arr_cp[i] = cp_mel;
          arr_temp[i] = temp_mel;
          arr_alpha[i] = alpha_mel;
          arr_cond[i] = cond_mel;
          arr_visc[i] = log10visc_mel;
      }
      else{
          if(arr_phi[i] > 0.5){
              fwtl = arr_fwtl[i]; // for smoothing
              arr_rho[i] = combine_matprop( fwtl, rho_mel, rho_mix );
              arr_dTdrs[i] = combine_matprop( fwtl, dTdrs_mel, dTdrs_mix );
              arr_cp[i] = combine_matprop( fwtl, cp_mel, cp_mix );
              arr_temp[i] = combine_matprop( fwtl, temp_mel, temp_mix );
              arr_alpha[i] = combine_matprop( fwtl, alpha_mel, alpha_mix );
              arr_cond[i] = combine_matprop( fwtl, cond_mel, cond_mix );
              arr_visc[i] = combine_matprop( fwtl, log10visc_mel, log10visc_mix );
          }
          else{
              fwts = arr_fwts[i]; // for smoothing
              arr_rho[i] = combine_matprop( fwts, rho_mix, rho_sol );
              arr_dTdrs[i] = combine_matprop( fwts, dTdrs_mix, dTdrs_sol );
              arr_cp[i] = combine_matprop( fwts, cp_mix, cp_sol );
              arr_temp[i] = combine_matprop( fwts, temp_mix, temp_sol );
              arr_alpha[i] = combine_matprop( fwts, alpha_mix, alpha_sol );
              arr_cond[i] = combine_matprop( fwts, cond_mix, cond_sol );
              arr_visc[i] = combine_matprop( fwts, log10visc_mix, log10visc_sol );
          }
      }
 
      /* compute viscosity */
      arr_visc[i] = PetscPowScalar( 10.0, arr_visc[i] );

      /* other useful material properties */
      /* kinematic viscosity */
      arr_nu[i] = arr_visc[i] / arr_rho[i];

      /* gravity * super-adiabatic temperature gradient */
      arr_gsuper[i] = P->gravity * arr_temp[i] / arr_cp[i] * arr_dSdr[i];

      /* Rayleigh number */
      /* FIXME: should use domain size not mixing length */
      arr_Ra[i] = arr_gsuper[i];
      arr_Ra[i] *= arr_alpha[i] * PetscPowScalar(arr_mix_b[i],4) * arr_rho[i] * arr_cp[i];
      arr_Ra[i] /= arr_nu[i] * arr_cond[i];

      /* eddy diffusivity */
      {
        /* always compute based on force balance and then select below */
        PetscScalar kh, crit;
        crit = 81.0 * PetscPowScalar(arr_nu[i],2);
        crit /= 4.0 * arr_alpha[i] * PetscPowScalar(arr_mix_b[i],4);

        if( arr_gsuper[i] <= 0.0 ){
          /* no convection, subadiabatic */
          kh = 0.0;
          arr_regime[i] = 0.0;
        } else if( arr_gsuper[i] > crit ){
          /* inviscid scaling from Vitense (1953) */
          kh = 0.25 * PetscPowScalar(arr_mix_b[i],2) * PetscSqrtScalar(arr_alpha[i]*arr_gsuper[i]);
          arr_regime[i] = 2.0;
        } else{
          /* viscous scaling */
          kh = arr_alpha[i] * arr_gsuper[i] * PetscPowScalar(arr_mix_b[i],4) / (18.0*arr_nu[i]);
          arr_regime[i] = 1.0;
        }

        /* thermal eddy diffusivity */
        if (P->eddy_diffusivity_thermal > 0.0){
          /* scale */
          arr_kappah[i] = P->eddy_diffusivity_thermal * kh;
        }
        else{
          /* else set (and negate to account for sign flag) */
          arr_kappah[i] = -P->eddy_diffusivity_thermal;
        }

        /* chemical eddy diffusivity */
        if (P->eddy_diffusivity_chemical > 0.0){
          /* scale */
          arr_kappac[i] = P->eddy_diffusivity_chemical * kh;
        }
        else{
          /* else set (and negate to account for sign flag) */
          arr_kappac[i] = -P->eddy_diffusivity_chemical;
        }
      }
    }

    ierr = DMDAVecRestoreArrayRead(da_b,S->dSdr,&arr_dSdr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->S,&arr_S_b); CHKERRQ(ierr);
    /* mesh quantities */
    ierr = DMDAVecRestoreArrayRead(da_b,M->dPdr_b,&arr_dPdr_b); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->mix_b,&arr_mix_b); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->pressure_b,&arr_pres); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->layer_b,&arr_layer_b); CHKERRQ(ierr);
    /* material properties */
    ierr = DMDAVecRestoreArray(    da_b,S->alpha,&arr_alpha); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cond,&arr_cond); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cp,&arr_cp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->cp_mix,&arr_cp_mix); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->dTdrs,&arr_dTdrs); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dTdrs_mix,&arr_dTdrs_mix); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion,&arr_fusion); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion_rho,&arr_fusion_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion_temp,&arr_fusion_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fwtl,&arr_fwtl); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fwts,&arr_fwts); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->gsuper,&arr_gsuper); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->kappac,&arr_kappac); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->kappah,&arr_kappah); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus,&arr_liquidus); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus_rho,&arr_liquidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus_temp,&arr_liquidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->nu,&arr_nu); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->phi,&arr_phi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Ra,&arr_Ra); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->rho,&arr_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus,&arr_solidus); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus_rho,&arr_solidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus_temp,&arr_solidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->temp,&arr_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->visc,&arr_visc); CHKERRQ(ierr);
    /* regime */
    ierr = DMDAVecRestoreArray(    da_b,S->regime,&arr_regime); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscScalar get_log10_viscosity_cutoff( PetscScalar in_visc, Parameters const *P )
{

    PetscScalar out_visc;

    out_visc = in_visc;

    if(P->log10visc_min > 0.0){
        if(in_visc < P->log10visc_min){
            out_visc = P->log10visc_min;
        }
    }
    if(P->log10visc_max > 0.0){
        if(in_visc > P->log10visc_max){
            out_visc = P->log10visc_max;
        }
    }

    return out_visc;

}

static PetscScalar get_log10_viscosity_solid( PetscScalar temperature, PetscScalar pressure, PetscInt layer, PetscScalar radius, Parameters const *P )
{

    PetscScalar Ea = P->activation_energy_sol; // activation energy (non-dimensional)
    PetscScalar Va = P->activation_volume_sol; // activation volume (non-dimensional)
    PetscScalar Mg_Si0 = P->Mg_Si0; // layer 0 (default) Mg/Si ratio
    PetscScalar Mg_Si1 = P->Mg_Si1; // layer 1 (basal layer) Mg/Si ratio
    PetscInt    VISCOUS_LID = P->VISCOUS_LID;
    PetscScalar lid_log10visc = P->lid_log10visc;
    PetscScalar lid_thickness = P->lid_thickness;

    PetscScalar A, lvisc;
    

    /* reference viscosity */
    lvisc = P->log10visc_sol; // i.e., log10(eta_0)

    /* temperature and pressure contribution
    A(T,P) = (E_a + V_a P) / RT
    eta = eta_0 * exp(A)
    log10(eta) = log10(eta0) + log10(exp(A))
    log10(eta) = P->log10visc_sol + A/ln(10) */
    A = 0.0;
    if(Ea>0.0)
        A += Ea;
    if(Va>0.0)
        A += Va*pressure;
    A *= 1.0 / temperature;
    lvisc += A / PetscLogReal(10);

    /* compositional contribution (based on Mg/Si ratio) */
    /* regular (default) layer */
    if(layer == 0){
        if(Mg_Si0 > 0.0)
            lvisc = add_compositional_viscosity( lvisc, Mg_Si0 );
    }
    /* optional, compositionally distinct layer */
    if(layer == 1){
        if(Mg_Si1 > 0.0)
            lvisc = add_compositional_viscosity( lvisc, Mg_Si1 );
    }

    /* viscous lid added by Rob Spaargaren */
    if(VISCOUS_LID){
        /* TODO: make this tanh to be smoother? */
        if(radius > P->radius - lid_thickness){
            lvisc += lid_log10visc;
        }
     }
    
    lvisc = get_log10_viscosity_cutoff( lvisc, P );

    return lvisc;

}

static PetscScalar add_compositional_viscosity( PetscScalar lvisc, PetscScalar Mg_Si ){

    /* These expressions were worked out by Rob Spaargaren as part
       of his MSc thesis (2018) */

    if(Mg_Si <= 0.5)
        lvisc += 2;
    else if (Mg_Si <= 0.7)
        lvisc += 2 - 1.4815 * (Mg_Si - 0.5)/0.2; // 1.4815 = 2 - log10(3.3)
    else if (Mg_Si <= 1.0)
        lvisc += 0.5185 * (1 - Mg_Si)/0.3; // 0.5185 = log10(3.3)
    else if (Mg_Si <= 1.25)
        lvisc += -1.4815 * (Mg_Si - 1)/0.25; // -1.4815 = log10(0.033)
    else if (Mg_Si <= 1.5)
        lvisc += -2 + (0.5185) * (1.5 - Mg_Si)/0.25; // 0.5185 = log10(0.033) - -2
    else
        lvisc += -2;

    return lvisc;
}

static PetscScalar get_log10_viscosity_melt( PetscScalar temperature, PetscScalar pressure, PetscInt layer, Parameters const *P )
{

    /* melt viscosity is currently a constant, but this retains symmetry with the function used
       to compute solid viscosity */

    PetscScalar lvisc;

    lvisc = P->log10visc_mel;

    lvisc = get_log10_viscosity_cutoff( lvisc, P );

    return lvisc;

}

static PetscScalar get_log10_viscosity_mix( PetscScalar meltf, PetscScalar log10visc_mel, PetscScalar log10visc_sol, Parameters const *P )
{
    PetscScalar fwt, lvisc;

    /* below needs revising to use critical melt fraction and not
       critical solid fraction */
    //fwt = viscosity_mix_skew( meltf );

    fwt = get_viscosity_mix_no_skew( meltf, P );
    lvisc = fwt * log10visc_mel + (1.0 - fwt) * log10visc_sol;

    return lvisc;
}

static PetscScalar get_viscosity_mix_no_skew( PetscScalar meltf, Parameters const *P )
{
    /* viscosity in mixed phase region with no skew */

    PetscScalar fwt;

    fwt = tanh_weight( meltf, P->phi_critical, P->phi_width );

    return fwt;

}

/* below is for the skewed viscosity formulation worked out by ASW */

#if 0
static PetscScalar viscosity_mix_skew( PetscScalar meltf )
{
    /* skewed viscosity in mixed phase region */

    PetscScalar fsol, fwt, fwt2, fwt3, z;

    fsol = 1.0 - meltf;
    z = (fsol-F_THRESHOLD) / DF_TRANSITION;
    fwt = zmap( z );

    z = (1.0-F_THRESHOLD) / DF_TRANSITION;
    fwt2 = zmap( z );

    z = -F_THRESHOLD / DF_TRANSITION;
    fwt3 = zmap( z );

    fwt = -(fwt-fwt3)/(fwt3-fwt2);

    return fwt;
}

static PetscScalar zmap( PetscScalar z )
{
    /* for skewed viscosity profile */

    PetscScalar fac, fwt, zmap, shp;

    shp = PHI_SKEW;

    fac = PetscSqrtScalar(PetscSqr(shp*z)+1.0);

    zmap = shp*PetscSqr(z)+z*fac+1.0/shp \
        *PetscLogScalar(shp*z+fac);
    zmap *= 0.5;

    fwt = tanh_weight( zmap, 0.0, 1.0 );

    return fwt;
}
#endif
