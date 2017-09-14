#include "atmosphere.h"
#include "ctx.h"

static PetscScalar get_atmosphere_mass( AtmosphereParameters const *, PetscScalar );
static PetscScalar get_optical_depth( AtmosphereParameters const *, PetscScalar, PetscScalar );
static PetscScalar get_dxdt( Atmosphere const *, AtmosphereParameters const *,PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_partial_pressure_volatile( AtmosphereParameters const *, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_partial_pressure_derivative_volatile( AtmosphereParameters const *, PetscScalar, PetscScalar, PetscScalar );
static PetscErrorCode set_pCO2( Atmosphere *, AtmosphereParameters const * );
static PetscErrorCode set_dpCO2dx( Atmosphere *, AtmosphereParameters const *  );
static PetscErrorCode set_pH2O( Atmosphere *, AtmosphereParameters const * );
static PetscErrorCode set_dpH2Odx( Atmosphere *, AtmosphereParameters const *  );

/* to solve for the initial volatile content of the magma ocean (liquid)
   we use Newton's method */
static PetscScalar f( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar f_prim( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar newton( PetscScalar, PetscScalar, PetscScalar );

/* general functions */

PetscErrorCode set_emissivity_abe_matsui( Atmosphere *A, AtmosphereParameters const *Ap )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    /* CO2 */
    ierr = set_pCO2( A, Ap ); CHKERRQ(ierr);
    A->m0 = get_atmosphere_mass( Ap, A->p0 );
    A->tau0 = get_optical_depth( Ap, A->m0, Ap->CO2_KABS );

    /* H2O */
    ierr = set_pH2O( A, Ap ); CHKERRQ(ierr);
    A->m1 = get_atmosphere_mass( Ap, A->p1 );
    A->tau1 = get_optical_depth( Ap, A->m1, Ap->H2O_KABS );

    /* total */
    A->tau = A->tau0 + A->tau1;
    A->emissivity = 2.0 / (A->tau + 2.0);

    PetscFunctionReturn(0);

}

static PetscScalar get_optical_depth( AtmosphereParameters const *Ap, PetscScalar mass_atm, PetscScalar kabs )
{
    PetscScalar tau;

    tau = 3.0 * mass_atm / (8.0*PETSC_PI*PetscSqr(Ap->RADIUS));
    // note negative gravity!
    tau *= PetscSqrtScalar( kabs*-Ap->GRAVITY/(3.0*Ap->P0) );

    return tau; // dimensionless (by definition)
}

static PetscScalar get_atmosphere_mass( AtmosphereParameters const *Ap, PetscScalar p )
{
    PetscScalar mass_atm;

    // p is partial pressure in Pa
    mass_atm = 4.0*PETSC_PI*PetscSqr(Ap->RADIUS) * p / -Ap->GRAVITY;

    return mass_atm; // kg

}

/////////////////////////////////////
/* generic functions for volatiles */
/////////////////////////////////////
static PetscScalar get_dxdt( Atmosphere const *A, AtmosphereParameters const *Ap, PetscScalar x, PetscScalar kdist, PetscScalar dpdx )
{
    PetscScalar       dxdt;
    PetscScalar       num, den;
    const PetscScalar M0 = A->M0;

    num = x * (kdist-1.0) * A->dMliqdt;
    den = kdist * M0 + (1.0-kdist) * A->Mliq;
    den += (4.0*PETSC_PI*PetscSqr(Ap->RADIUS) / -Ap->GRAVITY) * Ap->VOLSCALE * dpdx;

    dxdt = num / den;

    return dxdt;

}

static PetscScalar get_partial_pressure_volatile( AtmosphereParameters const *Ap, PetscScalar x, PetscScalar henry, PetscScalar henry_pow)
{

    /* partial pressure of volatile where x is wt % */

    PetscScalar p;

    p = (x / Ap->VOLSCALE) / henry;
    p = PetscPowScalar( p, henry_pow );

    return p; // Pa

}

static PetscScalar get_partial_pressure_derivative_volatile( AtmosphereParameters const *Ap, PetscScalar x, PetscScalar henry, PetscScalar henry_pow )
{

    /* derivative of partial pressure wrt x where x is wt % */

    PetscScalar dpdx;

    dpdx = 1.0 / PetscPowScalar( Ap->VOLSCALE*henry, henry_pow );
    dpdx *= henry_pow * PetscPowScalar( x, henry_pow-1.0);

    return dpdx; // Pa per wt %

}


static PetscScalar get_initial_volatile( Atmosphere const *A, AtmosphereParameters const *Ap, PetscScalar xinit, PetscScalar henry, PetscScalar henry_pow )
{

    /* initial wt. % of volatile in the aqueous phase */

    PetscScalar beta, gamma;
    PetscScalar x;

    gamma = 4.0 * PETSC_PI * PetscSqr(Ap->RADIUS);
    gamma /= -Ap->GRAVITY * A->M0;
    gamma *= PetscPowScalar( Ap->VOLSCALE, 1.0-henry_pow );
    gamma /= PetscPowScalar( henry, henry_pow );
    beta = henry_pow;

    x = newton( gamma, beta, xinit );

    return x; // wt %

}


////////////////////////////
/* CO2 specific functions */
////////////////////////////

PetscErrorCode set_dx0dt( Atmosphere *A, AtmosphereParameters const * Ap )
{
   /* update for dissolved CO2 content in the magma ocean */
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = set_dpCO2dx( A, Ap ); CHKERRQ(ierr);
    A->dx0dt = get_dxdt( A, Ap, A->x0, Ap->CO2_KDIST, A->dp0dx );

    PetscFunctionReturn(0);
}

PetscErrorCode set_initial_xCO2( Atmosphere *A, AtmosphereParameters const * Ap )
{
    PetscFunctionBeginUser;

    A->x0 = get_initial_volatile( A, Ap, Ap->CO2_INITIAL, Ap->CO2_HENRY, Ap->CO2_HENRY_POW );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_pCO2( Atmosphere *A, AtmosphereParameters const *Ap )
{
    /* partial pressure of CO2 */
    /* see magma_ocean_notes.tex */
    /* Lebrun et al. (2013) eqn. 17
       Massol et al. (2016) eqn. 18
       Salvador et al. (2017) eqn. C4 */

    PetscFunctionBeginUser;

    A->p0 = get_partial_pressure_volatile( Ap, A->x0, Ap->CO2_HENRY, Ap->CO2_HENRY_POW );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_dpCO2dx( Atmosphere *A, AtmosphereParameters const * Ap )
{
    /* dpCO2/dx. i.e., derivative of partial pressure (Pa) with respect
       to concentration */

    PetscFunctionBeginUser;

    A->dp0dx = get_partial_pressure_derivative_volatile( Ap, A->x0, Ap->CO2_HENRY, Ap->CO2_HENRY_POW );

    PetscFunctionReturn(0);

}

////////////////////////////
/* H2O specific functions */
////////////////////////////

PetscErrorCode set_dx1dt( Atmosphere *A, AtmosphereParameters const *Ap )
{
   /* update for dissolved H2O content in the magma ocean */
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = set_dpH2Odx( A, Ap ); CHKERRQ(ierr);
    A->dx1dt = get_dxdt( A, Ap, A->x1, Ap->H2O_KDIST, A->dp1dx );

    PetscFunctionReturn(0);
}

PetscErrorCode set_initial_xH2O( Atmosphere *A, AtmosphereParameters const *Ap )
{
    PetscFunctionBeginUser;

    A->x1 = get_initial_volatile( A, Ap, Ap->H2O_INITIAL, Ap->H2O_HENRY, Ap->H2O_HENRY_POW );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_pH2O( Atmosphere *A, AtmosphereParameters const *Ap )
{
    /* partial pressure of H2O */
    /* Lebrun et al. (2013) eqn. 16 */

    PetscFunctionBeginUser;

    A->p1 = get_partial_pressure_volatile( Ap, A->x1, Ap->H2O_HENRY, Ap->H2O_HENRY_POW );

    PetscFunctionReturn(0);
}

static PetscErrorCode set_dpH2Odx( Atmosphere *A, AtmosphereParameters const * Ap)
{
    /* dpH2O/dx. i.e., derivative of partial pressure (Pa) with respect
       to concentration */

    PetscFunctionBeginUser;

    A->dp1dx = get_partial_pressure_derivative_volatile( Ap, A->x1, Ap->H2O_HENRY, Ap->H2O_HENRY_POW );

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
