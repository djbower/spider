#include "ctx.h"
#include "atmosphere.h"

static PetscScalar get_atmosphere_mass( AtmosphereParameters const *, PetscScalar );
static PetscScalar get_optical_depth( AtmosphereParameters const *, PetscScalar, PetscScalar );
static PetscScalar get_dxdt( Atmosphere const *, AtmosphereParameters const *,PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_partial_pressure_volatile( AtmosphereParameters const *, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_partial_pressure_derivative_volatile( AtmosphereParameters const *, PetscScalar, PetscScalar, PetscScalar );

/* to solve for the initial volatile content of the magma ocean (liquid)
   we use Newton's method */
static PetscScalar f( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar f_prim( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar newton( PetscScalar, PetscScalar, PetscScalar );

///////////////////////
/* general functions */
///////////////////////

PetscScalar tsurf_param( PetscScalar temp, AtmosphereParameters const *Ap )
{
    PetscScalar Ts, c, fac, num, den;
    c = Ap->param_utbl_const;

    fac = 3.0*PetscPowScalar(c,3.0)*(27.0*PetscPowScalar(temp,2.0)*c+4.0);
    fac = PetscPowScalar( fac, 1.0/2.0 );
    fac += 9.0*temp*PetscPowScalar(c,2.0);
    // numerator
    num = PetscPowScalar(2.0,1.0/3)*PetscPowScalar(fac,2.0/3)-2.0*PetscPowScalar(3.0,1.0/3)*c;
    // denominator
    den = PetscPowScalar(6.0,2.0/3)*c*PetscPowScalar(fac,1.0/3);
    // surface temperature
    Ts = num / den;

    return Ts; 
}

PetscScalar grey_body( PetscScalar Tsurf, Atmosphere *A, AtmosphereParameters const *Ap )
{
    PetscScalar Fsurf;

    Fsurf = PetscPowScalar(Tsurf,4.0)-PetscPowScalar(Ap->teqm,4.0);
    Fsurf *= Ap->sigma * A->emissivity; /* Note emissivity may vary */

    return Fsurf;
}

PetscScalar steam_atmosphere_zahnle_1988( PetscScalar Tsurf, AtmosphereParameters const *Ap )
{
    PetscScalar       Fsurf;

    /* fit to Zahnle et al. (1988) from Solomatov and Stevenson (1993)
       Eqn. 40 */

    /* FIXME: will break for non-dimensional
       see commit 780b1dd to reverse this */

    Fsurf = 1.5E2 + 1.02E-5 * PetscExpScalar(0.011*Tsurf);

    return Fsurf;

}

PetscErrorCode set_emissivity_abe_matsui( Atmosphere *A, AtmosphereParameters const *Ap )
{
    VolatileParameters const *H2O = &Ap->H2O_volatile_parameters;
    VolatileParameters const *CO2 = &Ap->CO2_volatile_parameters;

    PetscFunctionBeginUser;

    /* CO2 */
    A->p0 = get_partial_pressure_volatile( Ap, A->x0, CO2->henry, CO2->henry_pow );
    A->m0 = get_atmosphere_mass( Ap, A->p0 );
    A->tau0 = get_optical_depth( Ap, A->m0, CO2->kabs );

    /* H2O */
    A->p1 = get_partial_pressure_volatile( Ap, A->x1, H2O->henry, H2O->henry_pow );
    A->m1 = get_atmosphere_mass( Ap, A->p1 );
    A->tau1 = get_optical_depth( Ap, A->m1, H2O->kabs );

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

    mass_atm = 4.0*PETSC_PI*PetscSqr(Ap->RADIUS) * p / -Ap->GRAVITY;

    return mass_atm;

}

/////////////////////////////////////
/* general functions for volatiles */
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

    /* partial pressure of volatile */

    PetscScalar p;

    p = (x / Ap->VOLSCALE) / henry;
    p = PetscPowScalar( p, henry_pow );

    return p;

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
    VolatileParameters const *CO2 = &Ap->CO2_volatile_parameters;

    PetscFunctionBeginUser;

    A->dp0dx = get_partial_pressure_derivative_volatile( Ap, A->x0, CO2->henry, CO2->henry_pow );
    A->dx0dt = get_dxdt( A, Ap, A->x0, CO2->kdist, A->dp0dx );

    PetscFunctionReturn(0);
}

PetscErrorCode set_initial_xCO2( Atmosphere *A, AtmosphereParameters const * Ap )
{
    VolatileParameters const *CO2 = &Ap->CO2_volatile_parameters;

    PetscFunctionBeginUser;

    A->x0 = get_initial_volatile( A, Ap, CO2->initial, CO2->henry, CO2->henry_pow );

    PetscFunctionReturn(0);

}

////////////////////////////
/* H2O specific functions */
////////////////////////////

PetscErrorCode set_dx1dt( Atmosphere *A, AtmosphereParameters const *Ap )
{
    /* update for dissolved H2O content in the magma ocean */
    VolatileParameters const *H2O = &Ap->H2O_volatile_parameters;

    PetscFunctionBeginUser;

    A->dp0dx = get_partial_pressure_derivative_volatile( Ap, A->x0, H2O->henry, H2O->henry_pow );
    A->dx1dt = get_dxdt( A, Ap, A->x1, H2O->kdist, A->dp1dx );

    PetscFunctionReturn(0);
}

PetscErrorCode set_initial_xH2O( Atmosphere *A, AtmosphereParameters const *Ap )
{
    VolatileParameters const *H2O = &Ap->H2O_volatile_parameters;

    PetscFunctionBeginUser;

    A->x1 = get_initial_volatile( A, Ap, H2O->initial, H2O->henry, H2O->henry_pow );

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
