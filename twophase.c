#include "parameters.h"
#include "dimensionalisablefield.h"
#include "twophase.h"
#include "util.h"
#include "eos.h"
#include "rheologicalfront.h"

static PetscErrorCode set_liquidus( Ctx *, PetscInt index );
static PetscErrorCode set_solidus( Ctx *, PetscInt index );
static PetscErrorCode set_fusion( Ctx * );

static PetscErrorCode set_rheological_front_mantle_properties( Ctx *, RheologicalFront *, PetscInt, Vec * );

PetscErrorCode set_twophase( Ctx *E )
{
    PetscFunctionBeginUser;

    Parameters const P = E->parameters;

    /* FIXME: below a bit hacky, but basically the slot number for the melt and solid EOSs changes
       depending on the number of phases in the system */
    if( P->n_phases==1 ){
        set_liquidus( E, 0 );
        set_solidus( E, 0 );
    }
    else{
        set_liquidus( E, 0);
        set_solidus( E, 1 );
    }

    /* these all need the liquidus and solidus to be set */
    set_fusion( E );

    PetscFunctionReturn(0);
}

PetscErrorCode set_gphi_smooth( Ctx *E )
{
    /* smoothing at each radial coordinate as a function
       of generalised melt fraction */

    PetscErrorCode ierr;
    DM             da_s=E->da_s, da_b=E->da_b;
    Parameters const P = E->parameters;
    Solution       *S = &E->solution;
    PetscInt       i, ilo_s, ihi_s, w_s, ilo, ihi, w;
    PetscScalar    *arr_gphi_s, *arr_gphi;

    PetscFunctionBeginUser;

    /* basic nodes */
    ierr = VecWAXPY(S->gphi,-1.0,S->solidus,S->S);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(S->gphi,S->gphi,S->fusion);CHKERRQ(ierr);

    /* staggered nodes */
    ierr = VecWAXPY(S->gphi_s,-1.0,S->solidus_s,S->S_s);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(S->gphi_s,S->gphi_s,S->fusion_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_liquidus( Ctx *E, PetscInt index )
{
    /* liquidus */

    PetscErrorCode ierr;
    PetscInt          i,ilo_b,ihi_b,ilo_s,ihi_s,w_s,w_b;
    DM                da_s=E->da_s, da_b=E->da_b;
    Vec               pres_b,pres_s;
    PetscScalar       z,*arr_liquidus,*arr_liquidus_s;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Interp1d          interp;
    Interp2d          interpR;
    Solution          *S;
    Parameters const  P = E->parameters;

    PetscFunctionBeginUser;

    S = &E->solution;
    interp = P->eos_parameters[index]->phase_boundary;
    interpR = P->eos_parameters[index]->lookup->rho;

    pres_b = E->mesh.pressure_b;
    pres_s = E->mesh.pressure_s;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    /* basic nodes */
    ierr = DMDAVecGetArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->liquidus,&arr_liquidus);CHKERRQ(ierr);
    for(i=ilo_b;i<ihi_b;++i){
        ierr = SetInterp1dValue( interp, arr_pres_b[i], &z, NULL );CHKERRQ(ierr);
        arr_liquidus[i] = z;
    }
    ierr = DMDAVecRestoreArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->liquidus,&arr_liquidus);CHKERRQ(ierr);

    /* staggered nodes */
    /* need in order to compute dfus/dr at basic internal nodes */
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        ierr = SetInterp1dValue( interp, arr_pres_s[i], &z, NULL );CHKERRQ(ierr);
        arr_liquidus_s[i] = z;
    }
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_solidus( Ctx *E, PetscInt index )
{
    /* solidus */

    PetscErrorCode ierr;
    PetscInt          i,ilo_b,ihi_b,w_b,ilo_s,ihi_s,w_s;
    DM                da_s=E->da_s,da_b=E->da_b;
    Vec               pres_b,pres_s;
    PetscScalar       z,*arr_solidus,*arr_solidus_s;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Interp1d          interp;
    Interp2d          interpR;
    Solution          *S;
    Parameters const  P = E->parameters;

    PetscFunctionBeginUser;
    S = &E->solution;
    interp = P->eos_parameters[index]->phase_boundary;
    interpR = P->eos_parameters[index]->lookup->rho;

    pres_b = E->mesh.pressure_b;
    pres_s = E->mesh.pressure_s;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    /* basic nodes */
    ierr = DMDAVecGetArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->solidus,&arr_solidus);CHKERRQ(ierr);
    for(i=ilo_b;i<ihi_b;++i){
        ierr = SetInterp1dValue( interp, arr_pres_b[i], &z, NULL );CHKERRQ(ierr);
        arr_solidus[i] = z;
    }
    ierr = DMDAVecRestoreArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->solidus,&arr_solidus);CHKERRQ(ierr);

    /* staggered nodes */
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        ierr = SetInterp1dValue( interp, arr_pres_s[i], &z, NULL );CHKERRQ(ierr);
        arr_solidus_s[i] = z;
    }
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);

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
    ierr = VecWAXPY(S->fusion,-1.0,S->solidus,S->liquidus);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_s,-1.0,S->solidus_s,S->liquidus_s);CHKERRQ(ierr);

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
    Parameters        P = E->parameters;
    Vec               result_s;
    PetscScalar       *arr_result_s;
    const PetscScalar *arr_dSdt_s, *arr_phi_s, *arr_mass_s, *arr_pres, *arr_S;
    EosEval           eos_eval;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = VecDuplicate( S->dSdt_s, &result_s); CHKERRQ(ierr);
    ierr = VecCopy( S->dSdt_s, result_s ); CHKERRQ(ierr);

    ierr = DMDAVecGetArrayRead(da_s,M->mass_s,&arr_mass_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,M->pressure_s,&arr_pres);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->dSdt_s,&arr_dSdt_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,result_s,&arr_result_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->S_s,&arr_S);CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){
        ierr = SetEosCompositeEval( P->eos_composites[0], arr_pres[i], arr_S[i], &eos_eval );CHKERRQ(ierr);
        arr_result_s[i] = arr_dSdt_s[i] * arr_mass_s[i];
        arr_result_s[i] /= eos_eval.fusion;

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
    ierr = DMDAVecRestoreArrayRead(da_s,S->phi_s,&arr_phi_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,result_s, &arr_result_s); CHKERRQ(ierr);

    ierr = VecSum( result_s, &A->dMliqdt ); CHKERRQ(ierr);

    VecDestroy( &result_s );

    PetscFunctionReturn(0);

}

/* FIXME: the next two functions relating to the rheological front
   require Ctx, but this is outside of the scope of rheologicalfront.h.
   Since ctx.h must know about rheologicalfront.h, to avoid a circular 
   dependency these two functions live here.  Ask PS about the best
   coding practice for dealing with this situation */

PetscErrorCode set_rheological_front( Ctx *E )
{
    PetscErrorCode    ierr;
    DM                da_s = E->da_s;
    DM                da_b = E->da_b;
    Parameters        P = E->parameters;
    Solution          *S = &E->solution;
    RheologicalFront  *Rfp = &E->rheological_front_phi;
    RheologicalFront  *Rfd = &E->rheological_front_dynamic;
    PetscInt          numpts_s, index;
    Vec               mask_s;
    PetscScalar       phi_global;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* critical melt fraction */
    index = get_crossover_index( da_s, S->phi_s, P->phi_critical, 1 );
    ierr = make_vec_mask( da_s, index, &mask_s );CHKERRQ(ierr);
    ierr = set_rheological_front_mantle_properties( E, Rfp, index, &mask_s );CHKERRQ(ierr);
    ierr = VecDestroy( &mask_s );CHKERRQ(ierr);

    /* inviscid to viscous regime crossover */
    index = get_crossover_index( da_b, S->regime, 1.5, 0 );
    ierr = make_vec_mask( da_s, index, &mask_s );CHKERRQ(ierr);
    ierr = set_rheological_front_mantle_properties( E, Rfd, index, &mask_s );CHKERRQ(ierr);
    ierr = VecDestroy( &mask_s );CHKERRQ(ierr);

    /* global melt fraction */
    ierr = make_vec_mask( da_s, numpts_s, &mask_s );CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, S->phi_s, &mask_s, &phi_global );CHKERRQ(ierr);
    /* store to both rheological structs */
    Rfp->phi_global = phi_global;
    Rfd->phi_global = phi_global;
    ierr = VecDestroy( &mask_s );CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_rheological_front_mantle_properties( Ctx *E, RheologicalFront *Rf, PetscInt index, Vec * mask_ptr_s )
{
    PetscErrorCode   ierr;
    const DM         da_s = E->da_s;
    const Mesh       *M = &E->mesh;
    const Parameters P = E->parameters;
    const Solution   *S = &E->solution;
    PetscScalar      phi, radius, pressure, temperature;
    PetscInt         numpts_s, index_above, index_below;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* rheological front coordinates */
    Rf->mesh_index = index;
    ierr = VecGetValues(M->radius_b,1,&index,&radius);CHKERRQ(ierr);
    Rf->depth = P->radius - radius;
    ierr = VecGetValues(M->pressure_b,1,&index,&pressure);CHKERRQ(ierr);
    Rf->pressure = pressure;
    ierr = VecGetValues(S->temp,1,&index,&temperature);CHKERRQ(ierr);
    Rf->temperature = temperature;

    /* mantle properties in magma ocean (above rheological front) */
    /* middle of layer */
    index_above = index/2; /* integer algebra */
    ierr = VecGetValues(S->phi,1,&index_above,&phi);CHKERRQ(ierr);
    Rf->above_middle.phi = phi;
    ierr = VecGetValues(M->radius_b,1,&index_above,&radius);CHKERRQ(ierr);
    Rf->above_middle.depth = P->radius - radius;
    ierr = VecGetValues(M->pressure_b,1,&index_above,&pressure);CHKERRQ(ierr);
    Rf->above_middle.pressure = pressure;
    ierr = VecGetValues(S->temp,1,&index_above,&temperature);CHKERRQ(ierr);
    Rf->above_middle.temperature = temperature;
    /* average by mass */
    ierr = average_by_mass_staggered( E, S->phi_s, mask_ptr_s, &Rf->above_mass_avg.phi); CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, M->radius_s, mask_ptr_s, &Rf->above_mass_avg.depth); CHKERRQ(ierr);
    Rf->above_mass_avg.depth = P->radius - Rf->above_mass_avg.depth;
    ierr = average_by_mass_staggered( E, M->pressure_s, mask_ptr_s, &Rf->above_mass_avg.pressure); CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, S->temp_s, mask_ptr_s, &Rf->above_mass_avg.temperature); CHKERRQ(ierr);

    // TODO: need to initialise values?
    /* middle of layer */
    index_below = (numpts_s - index)/2 + index;
    ierr = VecGetValues(S->phi,1,&index_below,&phi);CHKERRQ(ierr);
    Rf->below_middle.phi = phi;
    ierr = VecGetValues(M->radius_b,1,&index_below,&radius);CHKERRQ(ierr);
    Rf->below_middle.depth = P->radius - radius;
    ierr = VecGetValues(M->pressure_b,1,&index_below,&pressure);CHKERRQ(ierr);
    Rf->below_middle.pressure = pressure;
    ierr = VecGetValues(S->temp,1,&index_below,&temperature);CHKERRQ(ierr);
    Rf->below_middle.temperature = temperature;
    /* average by mass */
    ierr = invert_vec_mask( mask_ptr_s ); CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, S->phi_s, mask_ptr_s, &Rf->below_mass_avg.phi); CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, M->radius_s, mask_ptr_s, &Rf->below_mass_avg.depth); CHKERRQ(ierr);
    Rf->below_mass_avg.depth = P->radius - Rf->below_mass_avg.depth;
    ierr = average_by_mass_staggered( E, M->pressure_s, mask_ptr_s, &Rf->below_mass_avg.pressure); CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, S->temp_s, mask_ptr_s, &Rf->below_mass_avg.temperature); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
