#include "ctx.h"
#include "atmosphere.h"

static PetscScalar get_atmosphere_mass( AtmosphereParameters const *, PetscScalar );
static PetscScalar get_optical_depth( AtmosphereParameters const *, PetscScalar, VolatileParameters const * );
static PetscScalar get_dxdt( Atmosphere const *, AtmosphereParameters const *,PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_partial_pressure_volatile( AtmosphereParameters const *, PetscScalar, VolatileParameters const * );
static PetscScalar get_partial_pressure_derivative_volatile( AtmosphereParameters const *, PetscScalar, VolatileParameters const * );

/* to solve for the initial volatile content of the magma ocean (liquid)
   we use Newton's method */
static PetscScalar f( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar f_prim( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar newton( PetscScalar, PetscScalar, PetscScalar );

///////////////////////
/* general functions */
///////////////////////

PetscErrorCode set_atmosphere_volatile_content( Atmosphere *A, AtmosphereParameters const *Ap, PetscScalar x0, PetscScalar x1 )
{
    PetscFunctionBeginUser;

    VolatileParameters const *CO2 = &Ap->CO2_volatile_parameters;
    VolatileParameters const *H2O = &Ap->H2O_volatile_parameters;

    /* CO2 */
    A->p0 = get_partial_pressure_volatile( Ap, x0, CO2 );
    A->dp0dx = get_partial_pressure_derivative_volatile( Ap, x0, CO2 );
    A->m0 = get_atmosphere_mass( Ap, A->p0 );

    /* H2O */
    A->p1 = get_partial_pressure_volatile( Ap, x1, H2O );
    A->dp1dx = get_partial_pressure_derivative_volatile( Ap, x1, H2O );
    A->m1 = get_atmosphere_mass( Ap, A->p1 );

    PetscFunctionReturn(0);

}


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

PetscScalar grey_body( Atmosphere const *A, AtmosphereParameters const *Ap )
{
    PetscScalar Fsurf;

    Fsurf = PetscPowScalar(A->tsurf,4.0)-PetscPowScalar(Ap->teqm,4.0);
    Fsurf *= Ap->sigma * A->emissivity; /* Note emissivity may vary */

    return Fsurf;
}

PetscScalar steam_atmosphere_zahnle_1988( Atmosphere const *A, PetscScalar TEMP, PetscScalar FLUX )
{
    PetscScalar       Tsurf, Fsurf;

    /* fit to Zahnle et al. (1988) from Solomatov and Stevenson (1993)
       Eqn. 40.  Expressed dimensionally so must convert here using
       TEMP and FLUX scalings */

    Tsurf = A->tsurf * TEMP;
    Fsurf = 1.5E2 + 1.02E-5 * PetscExpScalar(0.011*Tsurf);
    Fsurf /= FLUX;

    return Fsurf;

}

PetscScalar get_emissivity_abe_matsui( Atmosphere *A, AtmosphereParameters const *Ap )
{

    PetscScalar emissivity;

    /* CO2 */
    A->tau0 = get_optical_depth( Ap, A->m0, &Ap->CO2_volatile_parameters );

    /* H2O */
    A->tau1 = get_optical_depth( Ap, A->m1, &Ap->H2O_volatile_parameters );

    /* total */
    A->tau = A->tau0 + A->tau1;
    emissivity = 2.0 / (A->tau + 2.0);

    return emissivity;

}

static PetscScalar get_optical_depth( AtmosphereParameters const *Ap, PetscScalar mass_atm, VolatileParameters const *V )
{
    PetscScalar tau;

    tau = 3.0 * mass_atm / (8.0*PETSC_PI*PetscSqr(Ap->RADIUS));
    // note negative gravity!
    tau *= PetscSqrtScalar( V->kabs*-Ap->GRAVITY/(3.0*Ap->P0) );

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
static PetscScalar get_dxdt( Atmosphere const *A, AtmosphereParameters const *Ap, PetscScalar x, PetscScalar kdist, PetscScalar dpdx, PetscScalar mantle_mass )
{
    PetscScalar       dxdt;
    PetscScalar       num, den;

    num = x * (kdist-1.0) * A->dMliqdt;
    den = kdist * mantle_mass + (1.0-kdist) * A->Mliq;
    den += (4.0*PETSC_PI*PetscSqr(Ap->RADIUS) / -Ap->GRAVITY) * Ap->volscale * dpdx;

    dxdt = num / den;

    return dxdt;

}

static PetscScalar get_partial_pressure_volatile( AtmosphereParameters const *Ap, PetscScalar x, VolatileParameters const *V)
{

    /* partial pressure of volatile */

    PetscScalar p;

    p = (x / Ap->volscale) / V->henry;
    p = PetscPowScalar( p, V->henry_pow );

    return p;

}

static PetscScalar get_partial_pressure_derivative_volatile( AtmosphereParameters const *Ap, PetscScalar x, VolatileParameters const *V )
{

    /* derivative of partial pressure wrt x where x is wt % */

    PetscScalar dpdx;

    dpdx = 1.0 / PetscPowScalar( Ap->volscale*V->henry, V->henry_pow );
    dpdx *= V->henry_pow * PetscPowScalar( x, V->henry_pow-1.0);

    return dpdx; // Pa per wt %

}


PetscScalar get_initial_volatile( AtmosphereParameters const *Ap, VolatileParameters const *V, PetscScalar mantle_mass )
{

    /* initial volatile in the aqueous phase */

    PetscScalar beta, gamma;
    PetscScalar x;

    gamma = 4.0 * PETSC_PI * PetscSqr(Ap->RADIUS);
    gamma /= -Ap->GRAVITY * mantle_mass;
    gamma *= PetscPowScalar( Ap->volscale, 1.0-V->henry_pow );
    gamma /= PetscPowScalar( V->henry, V->henry_pow );
    beta = V->henry_pow;

    x = newton( gamma, beta, V->initial );

    return x;

}

PetscScalar get_dx0dt( Atmosphere *A, AtmosphereParameters const * Ap, PetscScalar x0, PetscScalar mantle_mass )
{
    /* update for dissolved CO2 content in the magma ocean */
    VolatileParameters const *CO2 = &Ap->CO2_volatile_parameters;
    PetscScalar dx0dt;

    PetscFunctionBeginUser;

    dx0dt = get_dxdt( A, Ap, x0, CO2->kdist, A->dp0dx, mantle_mass );

    return dx0dt;
}

PetscScalar get_dx1dt( Atmosphere *A, AtmosphereParameters const *Ap, PetscScalar x1, PetscScalar mantle_mass )
{
    /* update for dissolved H2O content in the magma ocean */
    VolatileParameters const *H2O = &Ap->H2O_volatile_parameters;
    PetscScalar dx1dt;

    PetscFunctionBeginUser;

    dx1dt = get_dxdt( A, Ap, x1, H2O->kdist, A->dp1dx, mantle_mass );

    return dx1dt;
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
