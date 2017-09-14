/*
Atmosphere-related functions which interact with the general context
*/

#include "set_dxdt.h"
#include "atmosphere.h"

static PetscErrorCode set_Mliq( Ctx *E )
{
    PetscErrorCode ierr;
    Atmosphere     *A = &E->atmosphere;
    Solution       *S = &E->solution;
    Mesh           *M = &E->mesh;
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

static PetscErrorCode set_Msol( Ctx *E )
{
    PetscErrorCode ierr;
    Atmosphere     *A = &E->atmosphere;
    Solution       *S = &E->solution;
    Mesh           *M = &E->mesh;
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

static PetscErrorCode set_dMliqdt( Ctx *E )
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

PetscErrorCode set_dxdt( Ctx *E )
{
    PetscErrorCode ierr;
    Atmosphere           *A = &E->atmosphere;
    Parameters           *P = &E->parameters;
    AtmosphereParameters *Ap = &P->atmosphere_parameters;

    PetscFunctionBeginUser;

    // current stage of cooling
    ierr = set_Mliq( E ); CHKERRQ(ierr);
    ierr = set_Msol( E ); CHKERRQ(ierr);
    ierr = set_dMliqdt( E ); CHKERRQ(ierr);

    // CO2
    ierr = set_dx0dt(A,Ap); CHKERRQ(ierr);

    // H20
    ierr = set_dx1dt(A,Ap); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
