#include "atmosphere.h"

static PetscErrorCode set_Mliq( Ctx * );
static PetscErrorCode set_Msol( Ctx * );
static PetscErrorCode set_dMliqdt( Ctx * );

static PetscScalar get_atmosphere_mass( PetscScalar );
static PetscScalar get_optical_depth( PetscScalar, PetscScalar );
static PetscScalar get_dxdt( Atmosphere *, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_initial_volatile( Atmosphere *, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_partial_pressure_volatile( PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_partial_pressure_derivative_volatile( PetscScalar, PetscScalar, PetscScalar );

static PetscErrorCode set_dx0dt( Atmosphere * );
static PetscErrorCode set_pCO2( Atmosphere * );
static PetscErrorCode set_dpCO2dx( Atmosphere * );
static PetscErrorCode set_dx1dt( Atmosphere * );
static PetscErrorCode set_pH2O( Atmosphere * );
static PetscErrorCode set_dpH2Odx( Atmosphere * );

/* to solve for the initial volatile content of the magma ocean (liquid)
   we use Newton's method */
static PetscScalar f( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar f_prim( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar newton( PetscScalar, PetscScalar, PetscScalar );

/* general functions */

PetscErrorCode set_emissivity( Atmosphere *A )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    if (CO2_INITIAL <= 0.0 && H2O_INITIAL <= 0.0){
      A->emissivity = EMISSIVITY;
    }
    else{
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
    }

    PetscFunctionReturn(0);

}

static PetscScalar get_optical_depth( PetscScalar mass_atm, PetscScalar kabs )
{
    PetscScalar tau;

    tau = 3.0 * mass_atm / (8.0*PETSC_PI*PetscSqr(RADIUS));
    // note negative gravity!
    tau *= PetscSqrtScalar( kabs*-GRAVITY/(3.0*P0) );

    return tau; // dimensionless (by definition)
}

static PetscScalar get_atmosphere_mass( PetscScalar p )
{
    PetscScalar mass_atm;

    // p is partial pressure in Pa
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

    PetscFunctionBeginUser;

    // current stage of cooling
    ierr = set_Mliq( E ); CHKERRQ(ierr);
    ierr = set_Msol( E ); CHKERRQ(ierr);
    ierr = set_dMliqdt( E ); CHKERRQ(ierr);

    // CO2
    ierr = set_dx0dt( &E->atmosphere ); CHKERRQ(ierr);

    // H20
    ierr = set_dx1dt( &E->atmosphere ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/////////////////////////////////////
/* generic functions for volatiles */
/////////////////////////////////////
static PetscScalar get_dxdt( Atmosphere *A, PetscScalar x, PetscScalar kdist, PetscScalar dpdx )
{
    PetscScalar dxdt;
    PetscScalar num, den;
    PetscScalar M0 = A->M0;

    num = x * (kdist-1.0) * A->dMliqdt;
    den = kdist * M0 + (1.0-kdist) * A->Mliq;
    den += (4.0*PETSC_PI*PetscSqr(RADIUS) / -GRAVITY) * VOLSCALE * dpdx;

    dxdt = num / den;

    return dxdt;

}

static PetscScalar get_partial_pressure_volatile( PetscScalar x, PetscScalar henry, PetscScalar henry_pow)
{

    /* partial pressure of volatile where x is wt % */

    PetscScalar p;

    p = (x / VOLSCALE) / henry;
    p = PetscPowScalar( p, henry_pow );

    return p; // Pa

}

static PetscScalar get_partial_pressure_derivative_volatile( PetscScalar x, PetscScalar henry, PetscScalar henry_pow )
{

    /* derivative of partial pressure wrt x where x is wt % */

    PetscScalar dpdx;

    dpdx = 1.0 / PetscPowScalar( VOLSCALE*henry, henry_pow );
    dpdx *= henry_pow * PetscPowScalar( x, henry_pow-1.0);

    return dpdx; // Pa per wt %

}


static PetscScalar get_initial_volatile( Atmosphere *A, PetscScalar xinit, PetscScalar henry, PetscScalar henry_pow )
{

    /* initial wt. % of volatile in the aqueous phase */

    PetscScalar beta, gamma;
    PetscScalar x;

    gamma = 4.0 * PETSC_PI * PetscSqr(RADIUS);
    gamma /= -GRAVITY * A->M0;
    gamma *= PetscPowScalar( VOLSCALE, 1.0-henry_pow );
    gamma /= PetscPowScalar( henry, henry_pow );
    beta = henry_pow;

    x = newton( gamma, beta, xinit );

    return x; // wt %

}


////////////////////////////
/* CO2 specific functions */
////////////////////////////

PetscErrorCode set_dx0dt( Atmosphere *A )
{
   /* update for dissolved CO2 content in the magma ocean */
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = set_dpCO2dx( A ); CHKERRQ(ierr);
    A->dx0dt = get_dxdt( A, A->x0, CO2_KDIST, A->dp0dx );

    PetscFunctionReturn(0);
}

PetscErrorCode set_initial_xCO2( Atmosphere *A )
{
    PetscFunctionBeginUser;

    A->x0 = get_initial_volatile( A, CO2_INITIAL, CO2_HENRY, CO2_HENRY_POW );
    A->x0init = CO2_INITIAL;

    PetscFunctionReturn(0);

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

PetscErrorCode set_dx1dt( Atmosphere *A )
{
   /* update for dissolved H2O content in the magma ocean */
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = set_dpH2Odx( A ); CHKERRQ(ierr);
    A->dx1dt = get_dxdt( A, A->x1, H2O_KDIST, A->dp1dx );

    PetscFunctionReturn(0);
}

PetscErrorCode set_initial_xH2O( Atmosphere *A )
{
    PetscFunctionBeginUser;

    A->x1 = get_initial_volatile( A, H2O_INITIAL, H2O_HENRY, H2O_HENRY_POW );
    A->x1init = H2O_INITIAL;

    PetscFunctionReturn(0);

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

       x + gamma * x ** beta = xinit */

static PetscScalar f( PetscScalar x, PetscScalar A, PetscScalar B, PetscScalar C )
{
    PetscScalar result;

    result = x + A * PetscPowScalar( x, B ) - C;

    return result;

}

static PetscScalar f_prim( PetscScalar x, PetscScalar A, PetscScalar B, PetscScalar C )
{
    PetscScalar result;

    result = 1.0 + A*B*PetscPowScalar( x, B-1.0 );

    return result;

}

static PetscScalar newton( PetscScalar A, PetscScalar B, PetscScalar xinit )
{
    PetscInt i=0;
    PetscScalar x;
    x = xinit;
    while(i < 50){
        x = x - f( x, A, B, xinit ) / f_prim( x, A, B, xinit );
        i++;
    }
    return x;

}
