#include "matprop.h"

static PetscErrorCode set_matprop_staggered( Ctx * );
static PetscScalar log_viscosity_mix( PetscScalar );
static PetscScalar viscosity_mix_no_skew( PetscScalar );

PetscErrorCode set_capacitance_staggered( Ctx *E )
{
    PetscErrorCode    ierr;
    Mesh              *M;
    Solution          *S;

    PetscFunctionBeginUser;

    M = &E->mesh;
    S = &E->solution;

    ierr = set_matprop_staggered( E ); CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->lhs_s,S->temp_s,S->rho_s); CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->lhs_s,S->lhs_s,M->volume_s); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_matprop_staggered( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    DM                da_s=E->da_s;
    Lookup            *L; 
    Mesh              *M; 
    Solution          *S; 
    Vec               pres_s;
    // material properties that are updated here
    PetscScalar       *arr_phi_s, *arr_rho_s, *arr_temp_s, *arr_cp_s;
    // material properties used to update above
    const PetscScalar *arr_pres_s, *arr_liquidus_rho_s, *arr_solidus_rho_s, *arr_liquidus_temp_s, *arr_solidus_temp_s, *arr_S_s, *arr_liquidus_s, *arr_solidus_s, *arr_fusion_s, *arr_cp_mix_s;
    // for smoothing properties across liquidus and solidus
    const PetscScalar *arr_gphi_s, *arr_fwtl_s, *arr_fwts_s;
    PetscScalar       gphi, fwtl, fwts;

    PetscFunctionBeginUser;

    M = &E->mesh;
    S = &E->solution;
    pres_s = M->pressure_s;

    ierr = VecWAXPY(S->phi_s,-1.0,S->solidus_s,S->S_s);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(S->phi_s,S->phi_s,S->fusion_s);CHKERRQ(ierr);

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = DMDAVecGetArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fusion_s,&arr_fusion_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->gphi_s,&arr_gphi_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fwtl_s,&arr_fwtl_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fwts_s,&arr_fwts_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->cp_s,&arr_cp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->cp_mix_s,&arr_cp_mix_s);CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){

      // for smoothing
      gphi = arr_gphi_s[i];
      fwtl = arr_fwtl_s[i];
      fwts = arr_fwts_s[i];

      /* melt fraction, also truncated here [0,1] */
      if (arr_phi_s[i] > 1.0){
        /* superliquidus */
        arr_phi_s[i] = 1.0;
      }
      if (arr_phi_s[i] < 0.0){
        /* subsolidus */
        arr_phi_s[i] = 0.0;
      }

      /////////////////
      /* mixed phase */
      /////////////////
      /* density */
      arr_rho_s[i] = combine_matprop( arr_phi_s[i], 1.0/arr_liquidus_rho_s[i], 1.0/arr_solidus_rho_s[i] );
      arr_rho_s[i] = 1.0 / arr_rho_s[i];
      /* temperature */
      arr_temp_s[i] = combine_matprop( arr_phi_s[i], arr_liquidus_temp_s[i], arr_solidus_temp_s[i] );
      /* heat capacity */
      arr_cp_s[i] = arr_cp_mix_s[i];

      ////////////////
      /* melt phase */
      ////////////////
      if (gphi > 0.5){

          /* get melt properties */
          L = &E->melt_prop;
          /* density */
          arr_rho_s[i]   *= 1.0 - fwtl;
          arr_rho_s[i]   += fwtl * get_val2d( &L->rho, arr_pres_s[i], arr_S_s[i] );
          /* temperature */
          arr_temp_s[i]  *= 1.0 - fwtl;
          arr_temp_s[i]  += fwtl * get_val2d( &L->temp, arr_pres_s[i], arr_S_s[i] );
          /* heat capacity */
          arr_cp_s[i]    *= 1.0 - fwtl;
          arr_cp_s[i]    += fwtl * get_val2d( &L->cp, arr_pres_s[i], arr_S_s[i] );
      }

      /////////////////
      /* solid phase */
      /////////////////
      else if (gphi<=0.5){

          /* get solid properties */
          L = &E->solid_prop;
          /* density */
          arr_rho_s[i]   *= fwts;
          arr_rho_s[i]   += ( 1.0 - fwts ) * get_val2d( &L->rho, arr_pres_s[i], arr_S_s[i] );
          /* temperature */
          arr_temp_s[i]  *= fwts;
          arr_temp_s[i]  += ( 1.0 - fwts ) * get_val2d( &L->temp, arr_pres_s[i], arr_S_s[i] );
          /* heat capacity */
          arr_cp_s[i]    *= fwts;
          arr_cp_s[i]    += ( 1.0 - fwts ) * get_val2d( &L->cp, arr_pres_s[i], arr_S_s[i] );
      }

    }

    ierr = DMDAVecRestoreArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fusion_s,&arr_fusion_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->cp_mix_s,&arr_cp_mix_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->gphi_s,&arr_gphi_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fwtl_s,&arr_fwtl_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fwts_s,&arr_fwts_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->cp_s,&arr_cp_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_matprop_basic( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_b,ihi_b,w_b,ilo,ihi,numpts_b;
    DM                da_b=E->da_b;
    // material properties that are updated here
    PetscScalar       *arr_phi, *arr_nu, *arr_gsuper, *arr_kappah, *arr_dTdrs, *arr_alpha, *arr_temp, *arr_cp, *arr_cond, *arr_visc, *arr_rho;
    // material properties used to update above
    const PetscScalar *arr_dSdr, *arr_S_b, *arr_dSliqdr, *arr_dSsoldr, *arr_solidus, *arr_fusion, *arr_pres, *arr_dPdr_b, *arr_liquidus, *arr_liquidus_rho, *arr_solidus_rho, *arr_cp_mix, *arr_dTdrs_mix, *arr_liquidus_temp, *arr_solidus_temp, *arr_fusion_rho, *arr_fusion_temp, *arr_mix_b;
    // for smoothing properties across liquidus and solidus
    const PetscScalar *arr_gphi, *arr_fwtl, *arr_fwts;
    PetscScalar       gphi, fwtl, fwts;
    Lookup            *L; 
    Mesh              *M; 
    Solution          *S; 

    PetscFunctionBeginUser;

    M = &E->mesh;
    S = &E->solution;

    ierr = VecWAXPY(S->phi,-1.0,S->solidus,S->S);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(S->phi,S->phi,S->fusion);CHKERRQ(ierr);

    /* loop over all basic internal nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ilo = ilo_b == 0        ? 1            : ilo_b;
    ihi = ihi_b == numpts_b ? numpts_b - 1 : ihi_b;

    ierr = DMDAVecGetArrayRead(da_b,S->dSdr,&arr_dSdr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->S,&arr_S_b); CHKERRQ(ierr);
    /* mesh quantities */
    ierr = DMDAVecGetArrayRead(da_b,M->dPdr_b,&arr_dPdr_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->mix_b,&arr_mix_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->pressure_b,&arr_pres); CHKERRQ(ierr);
    /* material properties */
    ierr = DMDAVecGetArray(    da_b,S->alpha,&arr_alpha); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cond,&arr_cond); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cp,&arr_cp); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->cp_mix,&arr_cp_mix); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dSliqdr,&arr_dSliqdr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dSsoldr,&arr_dSsoldr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->dTdrs,&arr_dTdrs); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dTdrs_mix,&arr_dTdrs_mix); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion,&arr_fusion); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion_rho,&arr_fusion_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion_temp,&arr_fusion_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fwtl,&arr_fwtl); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fwts,&arr_fwts); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->gphi,&arr_gphi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->gsuper,&arr_gsuper); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->kappah,&arr_kappah); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus,&arr_liquidus); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus_rho,&arr_liquidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus_temp,&arr_liquidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->nu,&arr_nu); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->phi,&arr_phi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->rho,&arr_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus,&arr_solidus); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus_rho,&arr_solidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus_temp,&arr_solidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->temp,&arr_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->visc,&arr_visc); CHKERRQ(ierr);

    for(i=ilo; i<ihi; ++i){

      /* for smoothing */
      gphi = arr_gphi[i];
      fwtl = arr_fwtl[i];
      fwts = arr_fwts[i];

      /* truncate melt fraction */
      if (arr_phi[i] > 1.0){
        /* superliquidus */
        arr_phi[i] = 1.0;
      }
      if (arr_phi[i] < 0.0){
        /* subsolidus */
        arr_phi[i] = 0.0;
      }
      /////////////////
      /* mixed phase */
      /////////////////

      /* TODO: confirm that this cannot break for any value of i.
         by inspection, I do not think it can */
      /* density */
      arr_rho[i] = combine_matprop( arr_phi[i], 1.0/arr_liquidus_rho[i], 1.0/arr_solidus_rho[i] );
      arr_rho[i] = 1.0 / arr_rho[i];
      /* adiabatic temperature gradient */
      arr_dTdrs[i] = arr_dTdrs_mix[i];
      /* heat capacity */
      arr_cp[i] = arr_cp_mix[i];
      /* temperature */
      arr_temp[i] = combine_matprop( arr_phi[i], arr_liquidus_temp[i], arr_solidus_temp[i] );
      /* thermal expansion coefficient */
      arr_alpha[i] = -arr_fusion_rho[i] / arr_fusion_temp[i] / arr_rho[i];
      /* thermal conductivity */
      arr_cond[i] = combine_matprop( arr_phi[i], COND_MEL, COND_SOL );
      /* log viscosity */
      arr_visc[i] = log_viscosity_mix( arr_phi[i] );

      ////////////////
      /* melt phase */
      ////////////////
      if (gphi > 0.5){

        /* get melt properties */
        L = &E->melt_prop;
        /* density */
        arr_rho[i] *= 1.0 - fwtl;
        arr_rho[i] += fwtl * get_val2d( &L->rho, arr_pres[i], arr_S_b[i] );
        /* adiabatic temperature gradient */
        arr_dTdrs[i] *= 1.0 - fwtl;
        arr_dTdrs[i] += fwtl * arr_dPdr_b[i] * get_val2d( &L->dTdPs, arr_pres[i], arr_S_b[i] );
        /* heat capacity */
        arr_cp[i] *= 1.0 - fwtl;
        arr_cp[i] += fwtl * get_val2d( &L->cp, arr_pres[i], arr_S_b[i] );
        /* temperature */
        arr_temp[i] *= 1.0 - fwtl;
        arr_temp[i] += fwtl * get_val2d( &L->temp, arr_pres[i], arr_S_b[i] );
        /* thermal expansion coefficient */
        arr_alpha[i] *= 1.0 - fwtl;
        arr_alpha[i] += fwtl * get_val2d( &L->alpha, arr_pres[i], arr_S_b[i] );
        /* thermal conductivity */
        arr_cond[i] *= 1.0 - fwtl;
        arr_cond[i] += fwtl * COND_MEL;
        /* viscosity */
        arr_visc[i] *= 1.0 - fwtl;
        arr_visc[i] += fwtl * LOG10VISC_MEL;
      }

      else if (gphi <= 0.5){

        /* get solid properties */
        L = &E->solid_prop;
        /* density */
        arr_rho[i] *= fwts;
        arr_rho[i] += (1.0-fwts) * get_val2d( &L->rho, arr_pres[i], arr_S_b[i] );
        /* adiabatic temperature gradient */
        arr_dTdrs[i] *= fwts;
        arr_dTdrs[i] += (1.0-fwts) * arr_dPdr_b[i] * get_val2d( &L->dTdPs, arr_pres[i], arr_S_b[i] );
        /* heat capacity */
        arr_cp[i] *= fwts;
        arr_cp[i] += (1.0-fwts) * get_val2d( &L->cp, arr_pres[i], arr_S_b[i] );
        /* temperature */
        arr_temp[i] *= fwts;
        arr_temp[i] += (1.0-fwts) * get_val2d( &L->temp, arr_pres[i], arr_S_b[i] );
        /* thermal expansion coefficient */
        arr_alpha[i] *= fwts;
        arr_alpha[i] += (1.0-fwts) * get_val2d( &L->alpha, arr_pres[i], arr_S_b[i] );
        /* thermal conductivity */
        arr_cond[i] *= fwts;
        arr_cond[i] += (1.0-fwts) * COND_SOL;
        /* viscosity */
        arr_visc[i] *= fwts;
        arr_visc[i] += (1.0-fwts) * LOG10VISC_SOL;
      }

      /* compute viscosity */
      arr_visc[i] = PetscPowScalar( 10.0, arr_visc[i] );

      /* other useful material properties */
      /* kinematic viscosity */
      arr_nu[i] = arr_visc[i] / arr_rho[i];

      /* gravity * super-adiabatic temperature gradient */
      arr_gsuper[i] = GRAVITY * arr_temp[i] / arr_cp[i] * arr_dSdr[i];

      /* eddy diffusivity */
      PetscScalar kh, crit;
      crit = 81.0 * PetscPowScalar(arr_nu[i],2);
      crit /= 4.0 * arr_alpha[i] * PetscPowScalar(arr_mix_b[i],4);

      if( arr_gsuper[i] <= 0.0 ){
        /* no convection, subadiabatic */
        kh = 0.0;
      } else if( arr_gsuper[i] > crit ){
        /* inviscid scaling from Vitense (1953) */
        kh = 0.25 * PetscPowScalar(arr_mix_b[i],2) * PetscSqrtScalar(arr_alpha[i]*arr_gsuper[i]);
      } else{
        /* viscous scaling */
        kh = arr_alpha[i] * arr_gsuper[i] * PetscPowScalar(arr_mix_b[i],4) / (18.0*arr_nu[i]);
      }
      arr_kappah[i] = kh;

    }

    ierr = DMDAVecRestoreArrayRead(da_b,S->dSdr,&arr_dSdr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->S,&arr_S_b); CHKERRQ(ierr);
    /* mesh quantities */
    ierr = DMDAVecRestoreArrayRead(da_b,M->dPdr_b,&arr_dPdr_b); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->mix_b,&arr_mix_b); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->pressure_b,&arr_pres); CHKERRQ(ierr);
    /* material properties */
    ierr = DMDAVecRestoreArray(    da_b,S->alpha,&arr_alpha); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cond,&arr_cond); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cp,&arr_cp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->cp_mix,&arr_cp_mix); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dSliqdr,&arr_dSliqdr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dSsoldr,&arr_dSsoldr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->dTdrs,&arr_dTdrs); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dTdrs_mix,&arr_dTdrs_mix); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion,&arr_fusion); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion_rho,&arr_fusion_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion_temp,&arr_fusion_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fwtl,&arr_fwtl); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fwts,&arr_fwts); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->gphi,&arr_gphi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->gsuper,&arr_gsuper); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->kappah,&arr_kappah); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus,&arr_liquidus); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus_rho,&arr_liquidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus_temp,&arr_liquidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->nu,&arr_nu); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->phi,&arr_phi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->rho,&arr_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus,&arr_solidus); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus_rho,&arr_solidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus_temp,&arr_solidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->temp,&arr_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->visc,&arr_visc); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscScalar log_viscosity_mix( PetscScalar meltf )
{
    PetscScalar fwt, lvisc;

    /* below needs revising to use critical melt fraction and not
       critical solid fraction */
    //fwt = viscosity_mix_skew( meltf );

    fwt = viscosity_mix_no_skew( meltf );
    lvisc = fwt * LOG10VISC_MEL + (1.0 - fwt) * LOG10VISC_SOL;

    return lvisc;
}

static PetscScalar viscosity_mix_no_skew( PetscScalar meltf )
{
    /* viscosity in mixed phase region with no skew */

    PetscScalar fwt;

    fwt = tanh_weight( meltf, PHI_CRITICAL, PHI_WIDTH );

    return fwt;

}
