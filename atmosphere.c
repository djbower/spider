#include "atmosphere.h"

static PetscErrorCode set_Mliq( Ctx * );
static PetscErrorCode set_Msol( Ctx * );
static PetscErrorCode set_dMliqdt( Ctx * );
static PetscScalar get_atmosphere_mass( PetscScalar );
static PetscScalar get_optical_depth( PetscScalar, PetscScalar );
static PetscErrorCode set_pCO2( Atmosphere * );
static PetscErrorCode set_dpCO2dx( Atmosphere * );
static PetscErrorCode set_pH2O( Atmosphere * );
static PetscErrorCode set_dpH2Odx( Atmosphere * );
static PetscScalar get_partial_pressure_volatile( PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_partial_pressure_derivative_volatile( PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_initial_volatile( Ctx *, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar f( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar f_prim( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar newton( PetscScalar, PetscScalar, PetscScalar, PetscScalar );

/* general functions */

PetscErrorCode set_emissivity( Ctx *E )
{
    PetscErrorCode ierr;
    Atmosphere     *A = &E->atmosphere;

    PetscFunctionBeginUser;

    /* CO2 */
    ierr = set_pCO2( A ); CHKERRQ(ierr);
    A->m0 = get_atmosphere_mass( A->p0 );
    A->tau0 = get_optical_depth( A->m0, CO2_KABS );

    /* H2O */
    ierr = set_pH2O( A ); CHKERRQ(ierr);
    A->m1 = get_atmosphere_mass( A->p1 );
    A->tau1 = get_optical_depth( A->m1, H2O_KABS );

    /* total */
    A->tau = A->tau0 + A->tau1;
    A->emissivity = 2.0 / (A->tau + 2.0);

    /* for the case with no coupled atmospheric growth (H2O_INITIAL
       and CO2_INITIAL both zero), we can scale by the EMISSIVITY
       to give a grey-body model with constant emissivity.  Since
       above result will give emissivity=1.0 with no atmospheric
       growth (but see below, now over-written to avoid negative
       emissivity which could arise if the user enters a negative
       initial volatile content, presumably in error!) */
    /* this is really a check, since it makes no sense to have a
       negative initial volatile content, and therefore always
       switch to a constant emissivity grey-body if the volatile
       content is either zero or negative */
    if (CO2_INITIAL <= 0.0 && H2O_INITIAL <= 0.0){
        A->emissivity = EMISSIVITY;
    }

    PetscFunctionReturn(0);

}

static PetscScalar get_optical_depth( PetscScalar mass_atm, PetscScalar kabs )
{
    PetscScalar tau;

    tau = 3.0 * mass_atm / (8.0*PETSC_PI*PetscSqr(RADIUS));
    // note negative gravity!
    tau *= PetscSqrtScalar( kabs*-GRAVITY/(3.0*P0) );

    return tau;
}

static PetscScalar get_atmosphere_mass( PetscScalar p )
{
    PetscScalar mass_atm;

    mass_atm = 4.0*PETSC_PI*PetscSqr(RADIUS) * p / -GRAVITY;

    return mass_atm; // kg

}

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

        if (arr_phi_s[i] > 0.5){
            arr_result_s[i] *= ( 1.0 - arr_fwtl_s[i] );
        }

        else if (arr_phi_s[i] <=0.5){
            arr_result_s[i] *= arr_fwts_s[i];
        }
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

PetscErrorCode set_dx0dt( Ctx *E )
{
   /* update for dissolved CO2 content in the magma ocean */
    PetscErrorCode ierr;
    Atmosphere     *A = &E->atmosphere;
    Mesh           *M = &E->mesh;
    PetscScalar    num, den;

    PetscFunctionBeginUser;

    ierr = set_dpCO2dx( A ); CHKERRQ(ierr);
    ierr = set_Mliq( E ); CHKERRQ(ierr);
    ierr = set_Msol( E ); CHKERRQ(ierr);
    ierr = set_dMliqdt( E ); CHKERRQ(ierr);

    num = A->x0 * (CO2_KDIST-1.0) * A->dMliqdt;
    den = CO2_KDIST * M->mass0 + (1.0-CO2_KDIST) * A->Mliq;
    den += ( 4.0*PETSC_PI*PetscSqr(RADIUS) / -GRAVITY ) * A->dp0dx;

    A->dx0dt = num / den;

    PetscFunctionReturn(0);
}

/////////////////////////////////////
/* generic functions for volatiles */
/////////////////////////////////////

static PetscScalar get_partial_pressure_volatile( PetscScalar x, PetscScalar henry, PetscScalar henry_pow)
{

    /* partial pressure of volatile */

    PetscScalar p;

    // in x is wt %
    p = (x / 100.0) / henry; // numerator must be mass fraction, not wt %
    p = PetscPowScalar( p, henry_pow );

    return p; // Pa

}

static PetscScalar get_partial_pressure_derivative_volatile( PetscScalar x, PetscScalar henry, PetscScalar henry_pow )
{

    /* derivative of partial pressure with respect to x */

    PetscScalar dpdx;

    // in x is wt %
    dpdx = (x / 100.0) / henry;
    dpdx = PetscPowScalar( dpdx, henry_pow-1.0 );
    dpdx *= henry_pow / henry;

    return dpdx; // Pa per mass fraction (NOT wt %)

}


static PetscScalar get_initial_volatile( Ctx *E, PetscScalar xinit, PetscScalar henry, PetscScalar henry_pow )
{

    /* initial wt. % of volatiles in the aqueous phase */

    PetscScalar alpha, beta, gamma;
    PetscScalar M0 = E->mesh.mass0;
    PetscScalar x;

    x = xinit; // initial guess (wt. %)
    alpha = 4.0 * PETSC_PI * PetscSqr(RADIUS);
    alpha /= -GRAVITY * PetscPowScalar( henry, henry_pow ) * M0;
    beta = henry_pow;
    gamma = xinit;

    x = newton( x, alpha, beta, gamma );

    return x; // wt %

}


////////////////////////////
/* CO2 specific functions */
////////////////////////////

PetscScalar get_initial_xCO2( Ctx *E )
{
    Atmosphere *A = &E->atmosphere;

    A->x0 = get_initial_volatile( E, CO2_INITIAL, CO2_HENRY, CO2_HENRY_POW );
    A->x0init = CO2_INITIAL;

    // need to return ot add to augmented vector
    return A->x0;

}

static PetscErrorCode set_pCO2( Atmosphere *A )
{
    /* partial pressure of CO2 */
    /* see magma_ocean_notes.tex */
    /* Lebrun et al. (2013) eqn. 17
       Massol et al. (2016) eqn. 18
       Salvador et al. (2017) eqn. C4 */

    PetscFunctionBeginUser;

    A->p0 = get_partial_pressure_volatile( A->x0, CO2_HENRY, CO2_HENRY_POW );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_dpCO2dx( Atmosphere *A )
{
    /* dpCO2/dx. i.e., derivative of partial pressure (Pa) with respect
       to concentration */

    PetscFunctionBeginUser;

    A->dp0dx = get_partial_pressure_derivative_volatile( A->x0, CO2_HENRY, CO2_HENRY_POW );

    PetscFunctionReturn(0);

}

////////////////////////////
/* H2O specific functions */
////////////////////////////

PetscScalar get_initial_xH2O( Ctx *E )
{
    Atmosphere *A = &E->atmosphere;

    A->x1 = get_initial_volatile( E, H2O_INITIAL, H2O_HENRY, H2O_HENRY_POW );
    A->x1init = H2O_INITIAL;

    // need to return to add to augmented vector
    return A->x1;
}

static PetscErrorCode set_pH2O( Atmosphere *A )
{
    /* partial pressure of H2O */
    /* Lebrun et al. (2013) eqn. 16 */

    PetscFunctionBeginUser;

    A->p1 = get_partial_pressure_volatile( A->x1, H2O_HENRY, H2O_HENRY_POW );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_dpH2Odx( Atmosphere *A )
{
    /* dpH2O/dx. i.e., derivative of partial pressure (Pa) with respect
       to concentration */

    PetscFunctionBeginUser;

    A->dp1dx = get_partial_pressure_derivative_volatile( A->x1, H2O_HENRY, H2O_HENRY_POW );

    PetscFunctionReturn(0);

}

/////////////////////
/* Newton's method */
/////////////////////

/* for determining the initial mass fraction of volatiles in the
   melt.  The initial condition can be expressed as:

       x + alpha * x ** beta = gamma */

static PetscScalar f( PetscScalar x, PetscScalar alpha, PetscScalar beta, PetscScalar gamma )
{
    PetscScalar result;

    result = x + alpha * PetscPowScalar( x, beta ) - gamma;

    return result;

}

static PetscScalar f_prim( PetscScalar x, PetscScalar alpha, PetscScalar beta, PetscScalar gamma )
{
    PetscScalar result;

    result = 1.0 + alpha*beta*PetscPowScalar( x, beta-1.0 );

    return result;

}

static PetscScalar newton( PetscScalar x0, PetscScalar alpha, PetscScalar beta, PetscScalar gamma )
{
    PetscInt i=0;
    PetscScalar x;
    x = x0; // initial guess
    while(i < 10){
        x = x - f( x, alpha, beta, gamma ) / f_prim( x, alpha, beta, gamma );
        i++;
    }
    return x;

}
