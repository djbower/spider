#include "atmosphere.h"

static PetscScalar get_partial_pressure_volatile( PetscScalar, const VolatileParameters * );
static PetscScalar get_partial_pressure_derivative_volatile( PetscScalar, const VolatileParameters * );
static PetscScalar get_atmosphere_mass( const Parameters *, PetscScalar );
static PetscScalar get_optical_depth( const Parameters *, PetscScalar, const VolatileParameters * );
static PetscScalar get_newton_f( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_newton_f_prim( PetscScalar, PetscScalar, PetscScalar, PetscScalar );

static PetscScalar get_partial_pressure_volatile( PetscScalar x, const VolatileParameters *V )
{
    /* partial pressure of volatile */
    PetscScalar p;

    p = PetscPowScalar( x / V->henry, V->henry_pow);

    return p;
}

static PetscScalar get_partial_pressure_derivative_volatile( PetscScalar x, const VolatileParameters *V )
{
    /* partial pressure derivative of volatile */
    PetscScalar dpdx;

    dpdx = V->henry_pow / V->henry;
    dpdx *= PetscPowScalar( x / V->henry, V->henry_pow-1.0 );

    return dpdx;
}

static PetscScalar get_atmosphere_mass( const Parameters *P, PetscScalar p )
{
    /* mass of atmosphere */
    const Constants *C  = &P->constants;
    PetscScalar mass_atm;

    mass_atm = PetscSqr(P->radius) * p / -P->gravity;
    mass_atm *= 1.0E6 / C->VOLATILE;

    return mass_atm;
}

PetscErrorCode set_atmosphere_volatile_content( const Parameters *P, Atmosphere *A, PetscScalar x0, PetscScalar x1 )
{
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    VolatileParameters   const *CO2 = &Ap->CO2_volatile_parameters;
    VolatileParameters   const *H2O = &Ap->H2O_volatile_parameters;

    PetscFunctionBeginUser;

    /* if x0 and/or x1 are zero, the quantities below will also all
       be set to zero */

    /* CO2 */
    A->p0 = get_partial_pressure_volatile( x0, CO2 );
    A->dp0dx = get_partial_pressure_derivative_volatile( x0, CO2 );
    A->m0 = get_atmosphere_mass( P, A->p0 );

    /* H2O */
    A->p1 = get_partial_pressure_volatile( x1, H2O );
    A->dp1dx = get_partial_pressure_derivative_volatile( x1, H2O );
    A->m1 = get_atmosphere_mass( P, A->p1 );

    PetscFunctionReturn(0);
}

PetscScalar get_grey_body_flux( const Atmosphere *A, const AtmosphereParameters *Ap )
{
    PetscScalar Fsurf;

    Fsurf = PetscPowScalar(A->tsurf,4.0)-PetscPowScalar(Ap->teqm,4.0);
    Fsurf *= Ap->sigma * A->emissivity; /* Note emissivity may vary */

    return Fsurf;
}

PetscScalar get_steam_atmosphere_zahnle_1988_flux( const Atmosphere *A, const Constants *C )
{
    PetscScalar       Tsurf, Fsurf;

    /* fit to Zahnle et al. (1988) from Solomatov and Stevenson (1993)
       Eqn. 40.  Expressed dimensionally so must convert here using
       TEMP and FLUX scalings */

    Tsurf = A->tsurf * C->TEMP;
    Fsurf = 1.5E2 + 1.02E-5 * PetscExpScalar(0.011*Tsurf);
    Fsurf /= C->FLUX; // non-dimensionalise

    return Fsurf;

}
PetscScalar get_emissivity_from_flux( const Atmosphere *A, const AtmosphereParameters *Ap, PetscScalar flux )
{
    PetscScalar emissivity;

    emissivity = flux / Ap->sigma;
    emissivity /= PetscPowScalar(A->tsurf,4.0)-PetscPowScalar(Ap->teqm,4.0);

    return emissivity;

}

PetscScalar get_emissivity_abe_matsui( const Parameters *P, Atmosphere *A )
{
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;

    PetscScalar emissivity;

    /* CO2 */
    A->tau0 = get_optical_depth( P, A->m0, &Ap->CO2_volatile_parameters );

    /* H2O */
    A->tau1 = get_optical_depth( P, A->m1, &Ap->H2O_volatile_parameters );

    /* total */
    A->tau = A->tau0 + A->tau1;
    emissivity = 2.0 / (A->tau + 2.0);

    return emissivity;

}

static PetscScalar get_optical_depth( const Parameters *P, PetscScalar mass_atm, const VolatileParameters *V )
{
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    Constants            const *C  = &P->constants;
    PetscScalar tau;

    // note negative gravity
    tau = 3.0/2.0 * C->VOLATILE/1.0E6 * mass_atm / PetscSqr( P->radius );
    tau *= PetscSqrtScalar( V->kabs * -P->gravity / (3.0*Ap->P0) );

    return tau; // dimensionless (by definition)
}

/* Newton's method */
/* for determining the initial mass fraction of volatiles in the
   melt.  The initial condition can be expressed as:

       x + A * x ** B = C */

static PetscScalar get_newton_f( PetscScalar x, PetscScalar A, PetscScalar B, PetscScalar C ) 
{
    PetscScalar result;

    result = x + A * PetscPowScalar( x, B ) - C;

    return result;

}

static PetscScalar get_newton_f_prim( PetscScalar x, PetscScalar A, PetscScalar B, PetscScalar C ) 
{
    PetscScalar result;

    result = 1.0 + A*B*PetscPowScalar( x, B-1.0 );

    return result;

}

PetscScalar solve_newton_method( PetscScalar A, PetscScalar B, PetscScalar xinit )
{
    PetscInt i=0;
    PetscScalar x;
    x = xinit;
    while(i < 50){
        x = x - get_newton_f( x, A, B, xinit ) / get_newton_f_prim( x, A, B, xinit );
        i++;
    }   
    return x;

}
