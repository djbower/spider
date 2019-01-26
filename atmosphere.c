#include "atmosphere.h"
#include "dimensionalisablefield.h"

static PetscErrorCode set_partial_pressure_volatile( const VolatileParameters *, Volatile * );
static PetscErrorCode set_partial_pressure_derivative_volatile( const VolatileParameters *, Volatile * );
static PetscErrorCode set_atmosphere_mass( const AtmosphereParameters *, Volatile * );
static PetscErrorCode set_optical_depth(  const AtmosphereParameters *, const VolatileParameters *, Volatile * );
static PetscScalar solve_newton_method( PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_newton_f( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_newton_f_prim( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscErrorCode JSON_add_volatile( DM, Parameters const *, VolatileParameters const *, Volatile const *, Atmosphere const *, char const *name, cJSON * );

static PetscErrorCode set_partial_pressure_volatile( const VolatileParameters *Vp, Volatile *V )
{
    /* partial pressure of volatile */

    PetscFunctionBeginUser;

    V->p = PetscPowScalar( V->x / Vp->henry, Vp->henry_pow);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_partial_pressure_derivative_volatile( const VolatileParameters *Vp, Volatile *V )
{
    /* partial pressure derivative of volatile */

    PetscFunctionBeginUser;

    V->dpdx = Vp->henry_pow / Vp->henry;
    V->dpdx *= PetscPowScalar( V->x / Vp->henry, Vp->henry_pow-1.0 );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_atmosphere_mass( const AtmosphereParameters *Ap, Volatile *V )
{
    /* mass of atmosphere */

    PetscFunctionBeginUser;

    V->m = PetscSqr((*Ap->radius_ptr)) * V->p / -(*Ap->gravity_ptr);
    V->m *= 1.0E6 / (*Ap->VOLATILE_ptr);

    PetscFunctionReturn(0);

}

PetscErrorCode set_atmosphere_volatile_content( const AtmosphereParameters *Ap, Atmosphere *A )
{
    PetscErrorCode             ierr;
    VolatileParameters   const *CO2_parameters = &Ap->CO2_parameters;
    VolatileParameters   const *H2O_parameters = &Ap->H2O_parameters;
    Volatile                   *CO2 = &A->CO2;
    Volatile                   *H2O = &A->H2O;

    PetscFunctionBeginUser;

    /* if x0 and/or x1 are zero, the quantities below will also all
       be set to zero */

    /* CO2 */
    ierr = set_partial_pressure_volatile( CO2_parameters, CO2 );CHKERRQ(ierr);
    ierr = set_partial_pressure_derivative_volatile( CO2_parameters, CO2 );CHKERRQ(ierr);
    ierr = set_atmosphere_mass( Ap, CO2 );CHKERRQ(ierr);

    /* H2O */
    ierr = set_partial_pressure_volatile( H2O_parameters, H2O );CHKERRQ(ierr);
    ierr = set_partial_pressure_derivative_volatile( H2O_parameters, H2O );CHKERRQ(ierr);
    ierr = set_atmosphere_mass( Ap, H2O );CHKERRQ(ierr);

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

PetscScalar get_emissivity_abe_matsui( const AtmosphereParameters *Ap, Atmosphere *A )
{
    PetscErrorCode             ierr;
    VolatileParameters   const *CO2_parameters = &Ap->CO2_parameters;
    VolatileParameters   const *H2O_parameters = &Ap->H2O_parameters;
    Volatile                   *CO2 = &A->CO2;
    Volatile                   *H2O = &A->H2O;

    PetscScalar emissivity;

    /* CO2 */
    ierr = set_optical_depth( Ap, CO2_parameters, CO2 );CHKERRQ(ierr);

    /* H2O */
    ierr = set_optical_depth( Ap, H2O_parameters, H2O );CHKERRQ(ierr);

    /* total */
    A->tau = CO2->tau + CO2->tau;
    emissivity = 2.0 / (A->tau + 2.0);

    return emissivity;

}

static PetscErrorCode set_optical_depth( const AtmosphereParameters *Ap, const VolatileParameters *Vp, Volatile *V )
{
    PetscFunctionBeginUser;

    // note negative gravity
    V->tau = 3.0/2.0 * (*Ap->VOLATILE_ptr)/1.0E6 * V->m / PetscSqr( (*Ap->radius_ptr) );
    V->tau *= PetscSqrtScalar( Vp->kabs * -(*Ap->gravity_ptr) / (3.0*Ap->P0) );

    PetscFunctionReturn(0);
}

PetscErrorCode JSON_add_atmosphere( DM dm, Parameters const *P, Atmosphere const *A, const char *name, cJSON *json )
{
    PetscErrorCode ierr;
    cJSON          *data;
    PetscScalar    scaling;
    Constants      const *C = &P->constants;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;

    PetscFunctionBeginUser;

    data = cJSON_CreateObject();

    /* total liquid mass of mantle, kg */
    scaling = 4.0 * PETSC_PI * C->MASS; // includes 4*PI for spherical geometry
    ierr = JSON_add_single_value_to_object(dm, scaling, "mass_liquid", "kg", A->Mliq, data);CHKERRQ(ierr);

    /* total solid mass of mantle, kg */
    ierr = JSON_add_single_value_to_object(dm, scaling, "mass_solid", "kg", A->Msol, data);CHKERRQ(ierr);

    /* total mass of mantle, kg (for sanity check) */
    /* FIXME: will probably break due to scope of function */
    ierr = JSON_add_single_value_to_object(dm, scaling, "mass_mantle", "kg", *Ap->mantle_mass_ptr, data);CHKERRQ(ierr);

    /* surface temperature, K */
    scaling = C->TEMP;
    ierr = JSON_add_single_value_to_object(dm, scaling, "temperature_surface", "K", A->tsurf, data);CHKERRQ(ierr);

    /* optical depth, non-dimensional */
    scaling = 1.0;
    ierr = JSON_add_single_value_to_object(dm, scaling, "optical_depth", "None", A->tau, data);CHKERRQ(ierr);

    /* (effective) emissivity, non-dimensional */
    ierr = JSON_add_single_value_to_object(dm, scaling, "emissivity", "None", A->emissivity, data);CHKERRQ(ierr);

    /* CO2 */
    ierr = JSON_add_volatile(dm, P, &Ap->CO2_parameters, &A->CO2, A, "CO2", data ); CHKERRQ(ierr);

    /* H2O */
    ierr = JSON_add_volatile(dm, P, &Ap->H2O_parameters, &A->H2O, A, "H2O", data ); CHKERRQ(ierr);

    cJSON_AddItemToObject(json,name,data);

    PetscFunctionReturn(0);

}

static PetscErrorCode JSON_add_volatile( DM dm, Parameters const *P, VolatileParameters const *VP, Volatile const *V, Atmosphere const *A, char const *name, cJSON *json )
{
    PetscErrorCode  ierr;
    cJSON           *data;
    PetscScalar     scaling;
    Constants       const *C = &P->constants;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;

    PetscFunctionBeginUser;

    data = cJSON_CreateObject();

    /* parts-per-million (ppm) */
    scaling = C->VOLATILE;
    /* initial volatile (ppm) */
    ierr = JSON_add_single_value_to_object(dm, scaling, "initial_ppm", "ppm", VP->initial, data);CHKERRQ(ierr);
    /* volatile in liquid mantle (ppm) */
    ierr = JSON_add_single_value_to_object(dm, scaling, "liquid_ppm", "ppm", V->x, data);CHKERRQ(ierr);
    /* volatile in solid mantle (ppm) */
    ierr = JSON_add_single_value_to_object(dm, scaling, "solid_ppm", "ppm", V->x*VP->kdist, data);CHKERRQ(ierr);

    /* kilograms (kg) */
    scaling = (C->VOLATILE/1.0E6) * 4.0 * PETSC_PI * C->MASS;
    /* initial volatile (kg) */
    /* FIXME: mass will break here: move or copy M->mantle_mass to Constants? */
    ierr = JSON_add_single_value_to_object(dm, scaling, "initial_kg", "kg", VP->initial*(*Ap->mantle_mass_ptr), data);CHKERRQ(ierr);
    /* volatile in liquid mantle (kg) */
    ierr = JSON_add_single_value_to_object(dm, scaling, "liquid_kg", "kg", V->x*A->Mliq, data);CHKERRQ(ierr);
    /* volatile in solid mantle (kg) */
    ierr = JSON_add_single_value_to_object(dm, scaling, "solid_kg", "kg", V->x*VP->kdist*A->Msol, data);CHKERRQ(ierr);
    /* volatile in atmosphere (kg) */
    ierr = JSON_add_single_value_to_object(dm, scaling, "atmosphere_kg", "kg", V->m, data);CHKERRQ(ierr);

    /* bar (bar) */
    /* volatile in atmosphere (bar) */
    scaling = C->PRESSURE / 1.0E5; /* bar */
    ierr = JSON_add_single_value_to_object(dm, scaling, "atmosphere_bar", "bar", V->p, data);CHKERRQ(ierr);

    /* non-dimensional */
    /* optical depth (non-dimensional) */
    scaling = 1.0;
    ierr = JSON_add_single_value_to_object(dm, scaling, "optical_depth", "None", V->tau, data);CHKERRQ(ierr);

    cJSON_AddItemToObject(json,name,data);

    PetscFunctionReturn(0);

}

PetscScalar get_dxdt( const AtmosphereParameters *Ap, const Atmosphere *A, const VolatileParameters *Vp, const Volatile *V )
{
    PetscScalar dxdt;
    PetscScalar num, den;

    num = V->x * (Vp->kdist-1.0) * A->dMliqdt;
    den = Vp->kdist * (*Ap->mantle_mass_ptr) + (1.0-Vp->kdist) * A->Mliq;
    den += (1.0E6 / (*Ap->VOLATILE_ptr)) * PetscSqr( (*Ap->radius_ptr)) * V->dpdx / -(*Ap->gravity_ptr);

    dxdt = num / den;

    return dxdt;
}

PetscScalar get_initial_volatile( const AtmosphereParameters *Ap, const VolatileParameters *Vp )
{
    /* initial volatile in the aqueous phase */

    PetscScalar fac, x;

    fac = PetscSqr( (*Ap->radius_ptr) );
    fac /= -(*Ap->gravity_ptr) * (*Ap->mantle_mass_ptr);
    fac /= PetscPowScalar(Vp->henry,Vp->henry_pow);
    fac *= 1.0E6 / (*Ap->VOLATILE_ptr);

    x = solve_newton_method( fac, Vp->henry_pow, Vp->initial );

    return x;
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

static PetscScalar solve_newton_method( PetscScalar A, PetscScalar B, PetscScalar xinit )
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
