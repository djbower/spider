#include "parameters.h"
#include "matprop.h"
#include "util.h"
#include "twophase.h"
#include "eos.h"

static PetscErrorCode set_matprop_staggered( Ctx * );
static PetscScalar get_melt_fraction_truncated( PetscScalar );
static PetscErrorCode apply_log10visc_cutoff( Parameters const, PetscScalar * );
static PetscScalar GetModifiedMixingLength( PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar GetConstantMixingLength( PetscScalar outer_radius, PetscScalar inner_radius );
static PetscScalar GetMixingLength( const Parameters, PetscScalar );


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
    Parameters        P = E->parameters;
    Solution          *S = &E->solution;
    PetscScalar       *arr_phi_s;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    /* FIXME: update this with EOS two phase */

    // compute melt fraction
    /* FIXME: implicitly assumes solid state convection, but actually the melt
       fraction shouldn't do or mean anything for single phase systems */
    if(P->n_phases == 1){
        ierr = VecSet( S->phi_s, 0.0 );CHKERRQ(ierr); // by definition
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
    Mesh              *M = &E->mesh;
    Parameters const  P = E->parameters;
    Solution          *S = &E->solution;
    Vec               pres_s = M->pressure_s;
    // material properties that are updated here
    PetscScalar       *arr_rho_s, *arr_temp_s, *arr_cp_s;
    // for smoothing properties across liquidus and solidus
    const PetscScalar *arr_pres_s, *arr_S_s, *arr_phi_s, *arr_fwtl_s, *arr_fwts_s;
    PetscScalar       fwtl, fwts;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = DMDAVecGetArrayRead(da_s,S->fwtl_s,&arr_fwtl_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fwts_s,&arr_fwts_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->cp_s,&arr_cp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){

        /* single phase */
        if( P->n_phases==1 ){
            SetEosEval( P->eos_parameters[0], arr_pres_s[i], arr_S_s[i], &E->eos_evals[0] );
            arr_rho_s[i] = E->eos_evals[0].rho;
            arr_temp_s[i] = E->eos_evals[0].T;
            arr_cp_s[i] = E->eos_evals[0].Cp;
        }   

        else{
            /* mixed phase */
            SetEosCompositeEval( P->eos_composites[0], arr_pres_s[i], arr_S_s[i], &E->eos_evals[2] );
            if (arr_phi_s[i] > 0.5){
                /* blend melt (typically) and mixed */
                SetEosEval( P->eos_parameters[0], arr_pres_s[i], arr_S_s[i], &E->eos_evals[0] );
                fwtl = arr_fwtl_s[i]; // for smoothing
                arr_rho_s[i] = combine_matprop( fwtl, E->eos_evals[0].rho, E->eos_evals[2].rho );
                arr_temp_s[i] = combine_matprop( fwtl, E->eos_evals[0].T, E->eos_evals[2].T );
                arr_cp_s[i] = combine_matprop( fwtl, E->eos_evals[0].Cp, E->eos_evals[2].Cp );     
            }
            else{ 
                /* blend solid (typically) and mixed */
                SetEosEval( P->eos_parameters[1], arr_pres_s[i], arr_S_s[i], &E->eos_evals[1] );
                fwts = arr_fwts_s[i]; // for smoothing
                arr_rho_s[i] = combine_matprop( fwts, E->eos_evals[2].rho, E->eos_evals[1].rho );
                arr_temp_s[i] = combine_matprop( fwts, E->eos_evals[2].T, E->eos_evals[1].T );
                arr_cp_s[i] = combine_matprop( fwts, E->eos_evals[2].Cp, E->eos_evals[1].Cp );
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(da_s,S->fwtl_s,&arr_fwtl_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fwts_s,&arr_fwts_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
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
    PetscInt          i,ilo_b,ihi_b,w_b,ilo,ihi,numpts_b;
    DM                da_b=E->da_b;
    // material properties that are updated here
    PetscScalar       *arr_Ra, *arr_phi, *arr_nu, *arr_gsuper, *arr_kappac, *arr_kappah, *arr_dTdrs, *arr_alpha, *arr_temp, *arr_cp, *arr_cond, *arr_visc, *arr_regime, *arr_rho;
    // material properties used to update above
    const PetscScalar *arr_dSdr, *arr_S_b, *arr_solidus, *arr_fusion, *arr_pres, *arr_dPdr_b, *arr_liquidus, *arr_radius_b;
    const PetscInt *arr_layer_b;
    // for smoothing properties across liquidus and solidus
    const PetscScalar *arr_fwtl, *arr_fwts;
    PetscScalar       fwtl, fwts, mix;
    Mesh              *M = &E->mesh;
    Parameters const  P = E->parameters;
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
    ierr = DMDAVecGetArrayRead(da_b,M->pressure_b,&arr_pres); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->layer_b,&arr_layer_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->radius_b,&arr_radius_b); CHKERRQ(ierr);
    /* material properties */
    ierr = DMDAVecGetArray(    da_b,S->alpha,&arr_alpha); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cond,&arr_cond); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cp,&arr_cp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->dTdrs,&arr_dTdrs); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion,&arr_fusion); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fwtl,&arr_fwtl); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fwts,&arr_fwts); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->gsuper,&arr_gsuper); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->kappac,&arr_kappac); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->kappah,&arr_kappah); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus,&arr_liquidus); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->nu,&arr_nu); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->phi,&arr_phi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Ra,&arr_Ra); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->rho,&arr_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus,&arr_solidus); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->temp,&arr_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->visc,&arr_visc); CHKERRQ(ierr);
    /* regime: not convecting (0), inviscid (1), viscous (2) */
    ierr = DMDAVecGetArray(    da_b,S->regime,&arr_regime); CHKERRQ(ierr);

    for(i=ilo; i<ihi; ++i){

      /* single phase */
      if( P->n_phases==1 ){
          SetEosEval( P->eos_parameters[0], arr_pres[i], arr_S_b[i], &E->eos_evals[0] );
          arr_phi[i] = 0.0; // by definition
          arr_rho[i] = E->eos_evals[0].rho;
          arr_dTdrs[i] = arr_dPdr_b[i] * E->eos_evals[0].dTdPs;
          arr_cp[i] = E->eos_evals[0].Cp;
          arr_temp[i] = E->eos_evals[0].T;
          arr_alpha[i] = E->eos_evals[0].alpha;
          arr_cond[i] = E->eos_evals[0].cond;
          arr_visc[i] = E->eos_evals[0].log10visc;
      }

      else{
          /* truncate melt fraction */
          arr_phi[i] = get_melt_fraction_truncated( arr_phi[i] );
          /* mixed phase */
          SetEosCompositeEval( P->eos_composites[0], arr_pres[i], arr_S_b[i], &E->eos_evals[2] );
    
          if(arr_phi[i] > 0.5){
              /* blend melt (typically) and mixed */
              SetEosEval( P->eos_parameters[0], arr_pres[i], arr_S_b[i], &E->eos_evals[0] );
              fwtl = arr_fwtl[i]; // for smoothing
              arr_rho[i] = combine_matprop( fwtl, E->eos_evals[0].rho, E->eos_evals[2].rho );
#if 0
              /* this is the original formulation */
              arr_dTdrs[i] = combine_matprop( fwtl, E->eos_evals[0].dTdPs * arr_dPdr_b[i], arr_dTdrs_mix[i] );
#endif
#if 1
              /* this uses a different formulation for computing the mixed phase alpha */
              arr_dTdrs[i] = combine_matprop( fwtl, E->eos_evals[0].dTdPs * arr_dPdr_b[i], E->eos_evals[2].dTdPs * arr_dPdr_b[i] );
#endif
              arr_cp[i] = combine_matprop( fwtl, E->eos_evals[0].Cp, E->eos_evals[2].Cp );
              arr_temp[i] = combine_matprop( fwtl, E->eos_evals[0].T, E->eos_evals[2].T );
              arr_alpha[i] = combine_matprop( fwtl, E->eos_evals[0].alpha, E->eos_evals[2].alpha );
              arr_cond[i] = combine_matprop( fwtl, E->eos_evals[0].cond, E->eos_evals[2].cond );
              arr_visc[i] = combine_matprop( fwtl, E->eos_evals[0].log10visc, E->eos_evals[2].log10visc );
          }
          else{
              /* blend solid (typically) and mixed */
              SetEosEval( P->eos_parameters[1], arr_pres[i], arr_S_b[i], &E->eos_evals[1] );
              fwts = arr_fwts[i]; // for smoothing
              arr_rho[i] = combine_matprop( fwts, E->eos_evals[2].rho, E->eos_evals[1].rho );
#if 0
              /* this is the original formulation */
              arr_dTdrs[i] = combine_matprop( fwts, arr_dTdrs_mix[i], E->eos_evals[1].dTdPs * arr_dPdr_b[i] );
#endif
#if 1
              /* this uses a different formulation for computing the mixed phase alpha */
              arr_dTdrs[i] = combine_matprop( fwts, E->eos_evals[2].dTdPs * arr_dPdr_b[i], E->eos_evals[1].dTdPs * arr_dPdr_b[i] );
#endif
              arr_cp[i] = combine_matprop( fwts, E->eos_evals[2].Cp, E->eos_evals[1].Cp );
              arr_temp[i] = combine_matprop( fwts, E->eos_evals[2].T, E->eos_evals[1].T );
              arr_alpha[i] = combine_matprop( fwts, E->eos_evals[2].alpha, E->eos_evals[1].alpha );
              arr_cond[i] = combine_matprop( fwts, E->eos_evals[2].cond, E->eos_evals[1].cond );
              arr_visc[i] = combine_matprop( fwts, E->eos_evals[2].log10visc, E->eos_evals[1].log10visc );
          }
      }
 
      /* compute viscosity */
      /* note that prior versions of the code applied a cutoff to each individual
         phase, rather than to the aggregate.  But I think it makes the most sense to
         apply the cutoff once */
      ierr = apply_log10visc_cutoff( P, &arr_visc[i] );
      arr_visc[i] = PetscPowScalar( 10.0, arr_visc[i] );

      /* other useful material properties */
      /* kinematic viscosity */
      arr_nu[i] = arr_visc[i] / arr_rho[i];

      /* gravity * super-adiabatic temperature gradient */
      arr_gsuper[i] = P->gravity * arr_temp[i] / arr_cp[i] * arr_dSdr[i];

      /* FIXME: below */
#if 0
      /* Rayleigh number */
      /* FIXME: should use domain size not mixing length */
      arr_Ra[i] = arr_gsuper[i];
      arr_Ra[i] *= arr_alpha[i] * PetscPowScalar(arr_mix_b[i],4) * arr_rho[i] * arr_cp[i];
      arr_Ra[i] /= arr_nu[i] * arr_cond[i];
#endif
      arr_Ra[i] = 0.0; // HACK FOR RA TO AVOID UNINITIALISED VALUES

      /* eddy diffusivity */
      {
        /* always compute based on force balance and then select below */
        PetscScalar kh, crit;
        crit = 81.0 * PetscPowScalar(arr_nu[i],2);
        mix = GetMixingLength( P, arr_radius_b[i]);
        crit /= 4.0 * arr_alpha[i] * PetscPowScalar(mix,4);

        if( arr_gsuper[i] <= 0.0 ){
          /* no convection, subadiabatic */
          kh = 0.0;
          arr_regime[i] = 0.0;
        } else if( arr_gsuper[i] > crit ){
          /* inviscid scaling from Vitense (1953) */
          kh = 0.25 * PetscPowScalar(mix,2) * PetscSqrtScalar(arr_alpha[i]*arr_gsuper[i]);
          arr_regime[i] = 2.0;
        } else{
          /* viscous scaling */
          kh = arr_alpha[i] * arr_gsuper[i] * PetscPowScalar(mix,4) / (18.0*arr_nu[i]);
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
    ierr = DMDAVecRestoreArrayRead(da_b,M->pressure_b,&arr_pres); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->layer_b,&arr_layer_b); CHKERRQ(ierr);
    /* material properties */
    ierr = DMDAVecRestoreArray(    da_b,S->alpha,&arr_alpha); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cond,&arr_cond); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cp,&arr_cp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->dTdrs,&arr_dTdrs); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion,&arr_fusion); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fwtl,&arr_fwtl); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fwts,&arr_fwts); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->gsuper,&arr_gsuper); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->kappac,&arr_kappac); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->kappah,&arr_kappah); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus,&arr_liquidus); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->nu,&arr_nu); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->phi,&arr_phi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Ra,&arr_Ra); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->rho,&arr_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus,&arr_solidus); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->temp,&arr_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->visc,&arr_visc); CHKERRQ(ierr);
    /* regime */
    ierr = DMDAVecRestoreArray(    da_b,S->regime,&arr_regime); CHKERRQ(ierr);

    PetscFunctionReturn(0);
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

    /* for a single layer, P->layer_interface_radius = P->coresize (parameters.c),
       enabling this single expression to work for both a single and double
       layered mantle */

    if( radius > P->radius * P->layer_interface_radius ){
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
