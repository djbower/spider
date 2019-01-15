#include "dimensionalisablefield.h"
#include "twophase.h"
#include "util.h"
#include "lookup.h"

static PetscErrorCode set_liquidus( Ctx * );
static PetscErrorCode set_solidus( Ctx * );
static PetscErrorCode set_fusion( Ctx * );
static PetscErrorCode set_fusion_curve( Ctx * );
static PetscErrorCode set_mixed_phase( Ctx * );

static PetscErrorCode set_rheological_front_mask_phi( Ctx *, PetscInt *, Vec );
static PetscErrorCode set_rheological_front_mantle_properties( Ctx *, RheologicalFront *, PetscInt const, Vec const );

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

        /* no smoothing */
        if( P->matprop_smooth_width == 0.0 ){
            /* by default, mixed properties only */
            arr_fwtl[i] = 0.0;
            arr_fwts[i] = 1.0;
            // liquid properties only
            if( arr_gphi[i] > 1.0 ){
                arr_fwtl[i] = 1.0;
                //arr_fwts[i] = 0.0; // not used
            }
            // solid properties only
            if( arr_gphi[i] < 0.0 ){
                //arr_fwtl[i] = 0.0; // not used
                arr_fwts[i] = 0.0;
            }
        }

        /* tanh smoothing */
        else{
            /* fwtl -> 1.0 for gphi > 1.0 */
            /* fwtl -> 0.0 for gphi < 1.0 */
            arr_fwtl[i] = tanh_weight( arr_gphi[i], 1.0, P->matprop_smooth_width );
            /* fwts -> 1.0 for gphi > 0.0 */
            /* fwts -> 0.0 for gphi < 0.0 */
            arr_fwts[i] = tanh_weight( arr_gphi[i], 0.0, P->matprop_smooth_width );
        }
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
        arr_fwtl_s[i] = tanh_weight( arr_gphi_s[i], 1.0, P->matprop_smooth_width );
        arr_fwts_s[i] = tanh_weight( arr_gphi_s[i], 0.0, P->matprop_smooth_width );
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
    PetscInt          i,ilo_b,ihi_b,ilo_s,ihi_s,w_s,w_b;
    DM                da_s=E->da_s, da_b=E->da_b;
    Vec               pres_b,pres_s;
    PetscScalar       z,*arr_liquidus,*arr_liquidus_rho,*arr_liquidus_temp,*arr_liquidus_s,*arr_liquidus_rho_s,*arr_liquidus_temp_s;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Interp1d const    *interp;
    Interp2d const    *interpR, *interpT;
    Solution          *S;
    Parameters const  *P = &E->parameters;

    PetscFunctionBeginUser;

    S = &E->solution;
    interp = &P->melt_prop.liquidus;
    interpR = &P->melt_prop.rho;
    interpT = &P->melt_prop.temp;

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
        z = get_val1d( interp, arr_pres_b[i] );
        arr_liquidus[i] = z;
        arr_liquidus_rho[i] = get_val2d( interpR, arr_pres_b[i], z );
        arr_liquidus_temp[i] = get_val2d( interpT, arr_pres_b[i], z );
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
        z = get_val1d( interp, arr_pres_s[i] );
        arr_liquidus_s[i] = z;
        arr_liquidus_rho_s[i] = get_val2d( interpR, arr_pres_s[i], z );
        arr_liquidus_temp_s[i] = get_val2d( interpT, arr_pres_s[i], z );
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
    PetscInt          i,ilo_b,ihi_b,w_b,ilo_s,ihi_s,w_s;
    DM                da_s=E->da_s,da_b=E->da_b;
    Vec               pres_b,pres_s;
    PetscScalar       z,*arr_solidus,*arr_solidus_rho,*arr_solidus_temp,*arr_solidus_s,*arr_solidus_rho_s,*arr_solidus_temp_s;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Interp1d const    *interp;
    Interp2d const    *interpR, *interpT;
    Solution          *S;
    Parameters const  *P = &E->parameters;

    PetscFunctionBeginUser;
    S = &E->solution;
    interp = &P->solid_prop.solidus;
    interpR = &P->solid_prop.rho;
    interpT = &P->solid_prop.temp;

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
        z = get_val1d( interp, arr_pres_b[i] );
        arr_solidus[i] = z;
        arr_solidus_rho[i] = get_val2d( interpR, arr_pres_b[i], z );
        arr_solidus_temp[i] = get_val2d( interpT, arr_pres_b[i], z );
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
        z = get_val1d( interp, arr_pres_s[i] );
        arr_solidus_s[i] = z;
        arr_solidus_rho_s[i] = get_val2d( interpR, arr_pres_s[i], z );
        arr_solidus_temp_s[i] = get_val2d( interpT, arr_pres_s[i], z );
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

PetscErrorCode set_Mliq( Ctx *E )
{
    PetscErrorCode ierr;
    Atmosphere     *A = &E->atmosphere;
    Solution       const *S = &E->solution;
    Mesh           const *M = &E->mesh;
    Vec            mass_s;

    PetscFunctionBeginUser;

    // Mliq = sum[ phi*dm ]
    ierr = VecDuplicate( S->phi_s, &mass_s ); CHKERRQ(ierr);
    ierr = VecCopy( S->phi_s, mass_s ); CHKERRQ(ierr);
    ierr = VecPointwiseMult( mass_s, mass_s, M->mass_s ); CHKERRQ(ierr);
    ierr = VecSum( mass_s, &A->Mliq );

    VecDestroy( &mass_s );

    PetscFunctionReturn(0);
}

PetscErrorCode set_Msol( Ctx *E )
{
    PetscErrorCode ierr;
    Atmosphere     *A = &E->atmosphere;
    Solution       const *S = &E->solution;
    Mesh           const *M = &E->mesh;
    Vec            mass_s;

    PetscFunctionBeginUser;

    // Msol = sum[ (1-phi)*dm ]
    ierr = VecDuplicate( S->phi_s, &mass_s ); CHKERRQ(ierr);
    ierr = VecCopy( S->phi_s, mass_s ); CHKERRQ(ierr);
    ierr = VecScale( mass_s, -1.0 ); CHKERRQ(ierr);
    ierr = VecShift( mass_s, 1.0 ); CHKERRQ(ierr);
    ierr = VecPointwiseMult( mass_s, mass_s, M->mass_s ); CHKERRQ(ierr);
    ierr = VecSum( mass_s, &A->Msol );

    VecDestroy( &mass_s );

    PetscFunctionReturn(0);

}

PetscErrorCode set_dMliqdt( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    DM                da_s = E->da_s;
    Atmosphere        *A = &E->atmosphere;
    Solution          *S = &E->solution;
    Mesh              *M = &E->mesh;
    Vec               result_s;
    PetscScalar       *arr_result_s;
    const PetscScalar *arr_dSdt_s, *arr_fusion_s, *arr_fwtl_s, *arr_fwts_s, *arr_phi_s, *arr_mass_s;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = VecDuplicate( S->dSdt_s, &result_s); CHKERRQ(ierr);
    ierr = VecCopy( S->dSdt_s, result_s ); CHKERRQ(ierr);

    ierr = DMDAVecGetArrayRead(da_s,M->mass_s,&arr_mass_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->dSdt_s,&arr_dSdt_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fusion_s,&arr_fusion_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fwtl_s,&arr_fwtl_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fwts_s,&arr_fwts_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->phi_s,&arr_phi_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,result_s,&arr_result_s); CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){
        arr_result_s[i] = arr_dSdt_s[i] * arr_mass_s[i];
        arr_result_s[i] /= arr_fusion_s[i];

        /* with smoothing */
#if 0
        if (arr_phi_s[i] > 0.5){
            arr_result_s[i] *= ( 1.0 - arr_fwtl_s[i] );
        }

        else if (arr_phi_s[i] <=0.5){
            arr_result_s[i] *= arr_fwts_s[i];
        }
#endif

// generally had more luck without smoothing
#if 1
        /* no smoothing approach */
        if (arr_phi_s[i] <= 0.0){
            arr_result_s[i] = 0.0;
        }   
        else if (arr_phi_s[i] >=1.0){
            arr_result_s[i] = 0.0;
        }   
#endif

    }   

    ierr = DMDAVecRestoreArrayRead(da_s,M->mass_s,&arr_mass_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->dSdt_s, &arr_dSdt_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fusion_s,&arr_fusion_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fwtl_s,&arr_fwtl_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fwts_s,&arr_fwts_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->phi_s,&arr_phi_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,result_s, &arr_result_s); CHKERRQ(ierr);

    ierr = VecSum( result_s, &A->dMliqdt ); CHKERRQ(ierr);

    VecDestroy( &result_s );

    PetscFunctionReturn(0);

}

/* ----------------- */
/* rheological front */
/* ----------------- */

static PetscErrorCode set_rheological_front_mask_phi( Ctx *E, PetscInt *index, Vec mask_s )
{
    PetscErrorCode    ierr;
    DM                da_s = E->da_s;
    PetscInt          i,ilo_s,ihi_s,w_s;
    Parameters        *P = &E->parameters;
    Solution          *S = &E->solution;
    const PetscScalar *arr_phi_s;
    PetscScalar phi;
    const PetscScalar one=1.0, zero=0.0;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    /* this simple algorithm counts up from the surface towards the
       CMB until the melt fraction at a staggered node is larger than
       the melt fraction value (rheological_front_phi) used to define
       the base of the magma ocean.  Once this value is reached, the 
       loop is exited and the index stored for later use.  There are 
       many instances when this algorithm could return nonsense 
       values, but as long as the magma ocean is generally 
       crystalising from the bottom-up, it should be OK */

    /* end-member cases:
           i = 0 if surface is below rheological transition
           i = ihi_s-1 if all mantle is above rheological transition */

    /* this loop should always return a meaningful value if the cooling
       sequence can be adequately modelled as bottom-up */

    ierr = DMDAVecGetArrayRead(da_s,S->phi_s,&arr_phi_s); CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){
        phi = arr_phi_s[i];
        ierr = VecSetValues( mask_s, 1, &i, &one, INSERT_VALUES ); CHKERRQ(ierr);
        if( phi < P->phi_critical ){
            ierr = VecSetValues( mask_s, 1, &i, &zero, INSERT_VALUES ); CHKERRQ(ierr);
            break;
        }   
    }   

    ierr = VecAssemblyBegin(mask_s);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(mask_s);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);

    *index = i;

    PetscFunctionReturn(0);

}

PetscErrorCode set_rheological_front_phi( Ctx *E )
{
    PetscErrorCode    ierr;
    DM                da_s = E->da_s;
    RheologicalFront  *Rf = &E->rheological_front_phi;
    PetscInt          numpts_s, index;
    Vec               mask_s;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    index = numpts_s;

    /* create mask vector */
    ierr = VecCreate( PETSC_COMM_WORLD, &mask_s ); CHKERRQ(ierr);
    ierr = VecSetSizes( mask_s, PETSC_DECIDE, numpts_s ); CHKERRQ(ierr);
    ierr = VecSetFromOptions( mask_s ); CHKERRQ(ierr);
    ierr = VecSetUp( mask_s ); CHKERRQ(ierr);
    ierr = VecSet( mask_s, 0.0 ); CHKERRQ(ierr);

    /* TODO: this computes the rheological front mask based on the
       melt fraction, but another option is to compute based on the
       dynamic criterion instead */
    ierr = set_rheological_front_mask_phi( E, &index, mask_s ); CHKERRQ(ierr);

    ierr = set_rheological_front_mantle_properties( E, Rf, index, mask_s );

    ierr = VecDestroy( &mask_s ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_rheological_front_mantle_properties( Ctx *E, RheologicalFront *Rf, PetscInt const index, Vec mask_s )
{
    PetscErrorCode ierr;
    DM             da_s = E->da_s;
    // FIXME: make some of these const?
    Mesh           *M = &E->mesh;
    Parameters     *P = &E->parameters;
    Solution       *S = &E->solution;
    PetscScalar    phi, radius, pressure, temperature;
    PetscInt       numpts_s, i_above_avg, i_below_avg;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* rheological front coordinates */
    Rf->mesh_index = index;
    ierr = VecGetValues(M->radius_b,1,&index,&radius);CHKERRQ(ierr);
    Rf->depth = P->radius - radius;
    ierr = VecGetValues(M->pressure_b,1,&index,&pressure);CHKERRQ(ierr);
    Rf->pressure = pressure;

    /* mantle properties in magma ocean (above rheological front) */
    /* middle of layer */
    i_above_avg = index/2;
    ierr = VecGetValues(S->phi,1,&i_above_avg,&phi);CHKERRQ(ierr);
    Rf->above_middle.phi = phi;
    ierr = VecGetValues(M->radius_b,1,&i_above_avg,&radius);CHKERRQ(ierr);
    Rf->above_middle.depth = P->radius - radius;
    ierr = VecGetValues(M->pressure_b,1,&i_above_avg,&pressure);CHKERRQ(ierr);
    Rf->above_middle.pressure = pressure;
    ierr = VecGetValues(S->temp,1,&i_above_avg,&temperature);CHKERRQ(ierr);
    Rf->above_middle.temperature = temperature;
    /* average by mass */
    ierr = average_by_mass_staggered( E, S->phi_s, mask_s, &Rf->above_mass_avg.phi); CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, M->radius_s, mask_s, &Rf->above_mass_avg.depth); CHKERRQ(ierr);
    Rf->above_mass_avg.depth = P->radius - Rf->above_mass_avg.depth;
    ierr = average_by_mass_staggered( E, M->pressure_s, mask_s, &Rf->above_mass_avg.pressure); CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, S->temp_s, mask_s, &Rf->above_mass_avg.temperature); CHKERRQ(ierr);

    /* only compute properties in the solid layer once the rheological front
       begins advancing through the mantle */
    if( index < numpts_s){
        /* mantle properties in he solid layer (below rheological front) */
        /* middle of layer */
        i_below_avg = (numpts_s - index)/2 + index;
        ierr = VecGetValues(S->phi,1,&i_below_avg,&phi);CHKERRQ(ierr);
        Rf->below_middle.phi = phi;
        ierr = VecGetValues(M->radius_b,1,&i_below_avg,&radius);CHKERRQ(ierr);
        Rf->below_middle.depth = P->radius - radius;
        ierr = VecGetValues(M->pressure_b,1,&i_below_avg,&pressure);CHKERRQ(ierr);
        Rf->below_middle.pressure = pressure;
        ierr = VecGetValues(S->temp,1,&i_below_avg,&temperature);CHKERRQ(ierr);
        Rf->below_middle.temperature = temperature;
        /* average by mass */
        ierr = invert_vec_mask( mask_s ); CHKERRQ(ierr);
        /* mask is now one below the rheological front */
        ierr = average_by_mass_staggered( E, S->phi_s, mask_s, &Rf->below_mass_avg.phi); CHKERRQ(ierr);
        ierr = average_by_mass_staggered( E, M->radius_s, mask_s, &Rf->below_mass_avg.depth); CHKERRQ(ierr);
        Rf->below_mass_avg.depth = P->radius - Rf->below_mass_avg.depth;
        ierr = average_by_mass_staggered( E, M->pressure_s, mask_s, &Rf->below_mass_avg.pressure); CHKERRQ(ierr);
        ierr = average_by_mass_staggered( E, S->temp_s, mask_s, &Rf->below_mass_avg.temperature); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);

}

PetscErrorCode add_rheological_front_to_cJSON( Ctx const *E, RheologicalFront const *Rf, cJSON *json )
{
    PetscErrorCode        ierr;
    cJSON                *data;
    cJSON             *subdata;
    Constants          const *C = &E->parameters.constants;

    PetscFunctionBeginUser;

    data = cJSON_CreateObject();

    ierr = AddSingleValueToJSONArray(E->da_point, 1.0, "mesh_index", "None", Rf->mesh_index, data);CHKERRQ(ierr);
    ierr = AddSingleValueToJSONArray(E->da_point, C->RADIUS, "depth", "m", Rf->depth, data);CHKERRQ(ierr);
    ierr = AddSingleValueToJSONArray(E->da_point, C->PRESSURE, "pressure", "Pa", Rf->pressure, data);CHKERRQ(ierr);

    /* above rheological front */
    /* middle */
    subdata = cJSON_CreateArray();
    ierr = AddSingleValueToJSONArray(E->da_point, 1.0, "phi", "None", Rf->above_middle.phi, subdata);CHKERRQ(ierr);
    ierr = AddSingleValueToJSONArray(E->da_point, 1.0, "depth", "m", Rf->above_middle.depth, subdata);CHKERRQ(ierr);
    ierr = AddSingleValueToJSONArray(E->da_point, 1.0, "pressure", "Pa", Rf->above_middle.pressure, subdata);CHKERRQ(ierr);
    ierr = AddSingleValueToJSONArray(E->da_point, 1.0, "temperature", "K", Rf->above_middle.temperature, subdata);CHKERRQ(ierr);
    cJSON_AddItemToObject(data,"above_middle",subdata);

    /* mass averaged */
    /*subdata = cJSON_CreateArray();
    ierr = AddSingleValueToJSON(E->da_point, 1.0, "phi", "None", Rf->above_mass_avg.phi, subdata);CHKERRQ(ierr);
    ierr = AddSingleValueToJSON(E->da_point, 1.0, "depth", "m", Rf->above_mass_avg.depth, subdata);CHKERRQ(ierr);
    ierr = AddSingleValueToJSON(E->da_point, 1.0, "pressure", "Pa", Rf->above_mass_avg.pressure, subdata);CHKERRQ(ierr);
    ierr = AddSingleValueToJSON(E->da_point, 1.0, "temperature", "K", Rf->above_mass_avg.temperature, subdata);CHKERRQ(ierr);
    cJSON_AddItemToObject(data,"above_mass_avg",subdata);
    */

    cJSON_AddItemToObject(json,"rheological_front",data);


    PetscFunctionReturn(0);

}

