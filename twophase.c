#include "twophase.h"
#include "util.h"
#include "lookup.h"

static PetscErrorCode set_liquidus( Ctx * );
static PetscErrorCode set_solidus( Ctx * );
static PetscErrorCode set_fusion( Ctx * );
static PetscErrorCode set_fusion_curve( Ctx * );
static PetscErrorCode set_mixed_phase( Ctx * );

PetscErrorCode set_twophase( Ctx *E )
{
    PetscFunctionBeginUser;

    set_liquidus( E );
    set_solidus( E );

    /* these all need the liquidus and solidus to be set */
    set_fusion( E );
    set_fusion_curve( E );
    set_mixed_phase( E );

    PetscFunctionReturn(0);
}

PetscErrorCode set_gphi_smooth( Ctx *E )
{
    /* smoothing at each radial coordinate as a function
       of generalised melt fraction */

    PetscErrorCode ierr;
    DM             da_s=E->da_s, da_b=E->da_b;
    Parameters     *P = &E->parameters;
    Solution       *S = &E->solution;
    PetscInt       i, ilo_s, ihi_s, w_s, ilo, ihi, w;
    PetscScalar    *arr_gphi_s, *arr_fwtl_s, *arr_fwts_s, *arr_gphi, *arr_fwtl, *arr_fwts;

    PetscFunctionBeginUser;

    /* basic nodes */
    ierr = VecWAXPY(S->gphi,-1.0,S->solidus,S->S);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(S->gphi,S->gphi,S->fusion);CHKERRQ(ierr);
    ierr = DMDAGetCorners(da_b,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;

    ierr = DMDAVecGetArray(da_b,S->fwtl,&arr_fwtl);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->fwts,&arr_fwts);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->gphi,&arr_gphi);CHKERRQ(ierr);

    for(i=ilo; i<ihi; ++i){
        /* fwtl -> 1.0 for gphi > 1.0 */
        /* fwtl -> 0.0 for gphi < 1.0 */
        arr_fwtl[i] = tanh_weight( arr_gphi[i], 1.0, P->swidth );
        /* fwts -> 1.0 for gphi > 0.0 */
        /* fwts -> 0.0 for gphi < 0.0 */
        arr_fwts[i] = tanh_weight( arr_gphi[i], 0.0, P->swidth );
    }

    ierr = DMDAVecRestoreArray(da_b,S->fwtl,&arr_fwtl);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->fwts,&arr_fwts);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->gphi,&arr_gphi);CHKERRQ(ierr);

    /* staggered nodes */
    ierr = VecWAXPY(S->gphi_s,-1.0,S->solidus_s,S->S_s);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(S->gphi_s,S->gphi_s,S->fusion_s);CHKERRQ(ierr);
    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = DMDAVecGetArray(da_s,S->fwtl_s,&arr_fwtl_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->fwts_s,&arr_fwts_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->gphi_s,&arr_gphi_s);CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){
        arr_fwtl_s[i] = tanh_weight( arr_gphi_s[i], 1.0, P->swidth );
        arr_fwts_s[i] = tanh_weight( arr_gphi_s[i], 0.0, P->swidth );
    }

    ierr = DMDAVecRestoreArray(da_s,S->fwtl_s,&arr_fwtl_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->fwts_s,&arr_fwts_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->gphi_s,&arr_gphi_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_liquidus( Ctx *E )
{
    /* liquidus */

    PetscErrorCode ierr;
    PetscInt       i,ilo_b,ihi_b,ilo_s,ihi_s,w_s,w_b;
    DM             da_s=E->da_s, da_b=E->da_b;
    Vec            pres_b,pres_s;
    PetscScalar    z,*arr_liquidus,*arr_liquidus_rho,*arr_liquidus_temp,*arr_liquidus_s,*arr_liquidus_rho_s,*arr_liquidus_temp_s;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Interp1d       *I; 
    Solution       *S; 
    Interp2d       *IR, *IT;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_liquidus:\n");CHKERRQ(ierr);
#endif

    S = &E->solution;
    I = &E->melt_prop.liquidus;
    IR = &E->melt_prop.rho;
    IT = &E->melt_prop.temp;

    pres_b = E->mesh.pressure_b;
    pres_s = E->mesh.pressure_s;


    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    /* basic nodes */
    ierr = DMDAVecGetArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->liquidus,&arr_liquidus);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->liquidus_temp,&arr_liquidus_temp);CHKERRQ(ierr);
    for(i=ilo_b;i<ihi_b;++i){
        z = get_val1d( I, arr_pres_b[i] );
        arr_liquidus[i] = z;
        arr_liquidus_rho[i] = get_val2d( IR, arr_pres_b[i], z );
        arr_liquidus_temp[i] = get_val2d( IT, arr_pres_b[i], z );
    }   
    ierr = DMDAVecRestoreArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->liquidus,&arr_liquidus);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->liquidus_temp,&arr_liquidus_temp);CHKERRQ(ierr);


    /* staggered nodes */
    /* need in order to compute dfus/dr at basic internal nodes */ 
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        z = get_val1d( I, arr_pres_s[i] );
        arr_liquidus_s[i] = z;
        arr_liquidus_rho_s[i] = get_val2d( IR, arr_pres_s[i], z );
        arr_liquidus_temp_s[i] = get_val2d( IT, arr_pres_s[i], z );
    }
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_solidus( Ctx *E )
{
    /* solidus */

    PetscErrorCode ierr;
    PetscInt       i,ilo_b,ihi_b,w_b,ilo_s,ihi_s,w_s;
    DM             da_s=E->da_s,da_b=E->da_b;
    Vec            pres_b,pres_s;
    PetscScalar    z,*arr_solidus,*arr_solidus_rho,*arr_solidus_temp,*arr_solidus_s,*arr_solidus_rho_s,*arr_solidus_temp_s;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Interp1d       *I;
    Solution       *S;
    Interp2d       *IR, *IT;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_solidus:\n");CHKERRQ(ierr);
#endif
    S = &E->solution;
    I = &E->solid_prop.solidus;
    IR = &E->solid_prop.rho;
    IT = &E->solid_prop.temp;

    pres_b = E->mesh.pressure_b;
    pres_s = E->mesh.pressure_s;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    /* basic nodes */
    ierr = DMDAVecGetArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->solidus_temp,&arr_solidus_temp);CHKERRQ(ierr);
    for(i=ilo_b;i<ihi_b;++i){
        z = get_val1d( I, arr_pres_b[i] );
        arr_solidus[i] = z;
        arr_solidus_rho[i] = get_val2d( IR, arr_pres_b[i], z );
        arr_solidus_temp[i] = get_val2d( IT, arr_pres_b[i], z );
    }
    ierr = DMDAVecRestoreArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->solidus_temp,&arr_solidus_temp);CHKERRQ(ierr);


    /* staggered nodes */
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        z = get_val1d( I, arr_pres_s[i] );
        arr_solidus_s[i] = z;
        arr_solidus_rho_s[i] = get_val2d( IR, arr_pres_s[i], z );
        arr_solidus_temp_s[i] = get_val2d( IT, arr_pres_s[i], z );
    }
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_fusion( Ctx *E )
{
    /* entropy of fusion */
    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_fusion:\n");CHKERRQ(ierr);
#endif
    S = &E->solution;

    /* fusion = liquidus - solidus */
    /* basic nodes */
    ierr = VecWAXPY(S->fusion,-1.0,S->solidus,S->liquidus);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_rho,-1.0,S->solidus_rho,S->liquidus_rho);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_temp,-1.0,S->solidus_temp,S->liquidus_temp);CHKERRQ(ierr);

    /* staggered nodes */
    ierr = VecWAXPY(S->fusion_s,-1.0,S->solidus_s,S->liquidus_s);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_temp_s,-1.0,S->solidus_temp_s,S->liquidus_temp_s);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode set_fusion_curve( Ctx *E )
{
    /* fusion curve is defined by the 50% melt fraction contour */

    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_fusion_curve:\n");CHKERRQ(ierr);
#endif
    S = &E->solution;

    /* basic nodes fusion_curve = solidus + 0.5*fusion*/
    ierr = VecWAXPY(S->fusion_curve,0.5,S->fusion,S->solidus);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_curve_temp,0.5,S->fusion_temp,S->solidus_temp);CHKERRQ(ierr);

    /* staggered nodes */
    ierr = VecWAXPY(S->fusion_curve_s,0.5,S->fusion_s,S->solidus_s);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_curve_temp_s,0.5,S->fusion_temp_s,S->solidus_temp_s);CHKERRQ(ierr);

    /* basic nodes */
    MatMult( E->ddr_at_b, S->fusion_curve_s, S->dfusdr ); CHKERRQ(ierr);
    MatMult( E->ddr_at_b, S->fusion_curve_temp_s, S->dfusdr_temp ); CHKERRQ(ierr);
    MatMult( E->ddr_at_b, S->liquidus_s, S->dSliqdr ); CHKERRQ(ierr);
    MatMult( E->ddr_at_b, S->solidus_s, S->dSsoldr ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_mixed_phase( Ctx *E )
{
    /* dTdrs and heat capacity can be precomputed for the mixed phase
       because they are only pressure-dependent */

    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_mixed_phase\n");CHKERRQ(ierr);
#endif

    S = &E->solution;

    /* dTdrs_mix = -dfusdr * (fusion_temp/fusion) + dfusdr_temp */
    ierr = VecPointwiseDivide(S->dTdrs_mix,S->fusion_temp,S->fusion);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->dTdrs_mix,S->dTdrs_mix,S->dfusdr);CHKERRQ(ierr);
    ierr = VecScale( S->dTdrs_mix, -1.0 );
    ierr = VecAXPY(S->dTdrs_mix,1.0,S->dfusdr_temp);CHKERRQ(ierr);

    /* cp_mix */
    /* basic nodes */
    ierr = VecPointwiseDivide(S->cp_mix,S->fusion,S->fusion_temp);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->cp_mix,S->cp_mix,S->fusion_curve_temp);CHKERRQ(ierr);

    /* staggered nodes */
    ierr = VecPointwiseDivide(S->cp_mix_s,S->fusion_s,S->fusion_temp_s);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->cp_mix_s,S->cp_mix_s,S->fusion_curve_temp_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
