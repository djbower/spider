#include "atmosphere.h"
#include "dimensionalisablefield.h"
#include "util.h"

static PetscErrorCode set_partial_pressure_volatile( const VolatileParameters *, Volatile * );
static PetscErrorCode set_partial_pressure_derivative_volatile( const VolatileParameters *, Volatile * );
static PetscErrorCode set_atmosphere_mass( const Atmosphere *, const AtmosphereParameters *, const VolatileParameters *, Volatile * );
static PetscErrorCode set_optical_depth( const AtmosphereParameters *, const VolatileParameters *, Volatile * );
static PetscScalar get_pressure_dependent_kabs( const AtmosphereParameters *, const VolatileParameters * );
static PetscErrorCode set_mixing_ratios( Volatile *, Volatile * );
static PetscErrorCode set_jeans( const Atmosphere *, const AtmosphereParameters *, const VolatileParameters *, Volatile * );
static PetscErrorCode set_f_thermal_escape( const Atmosphere *, const AtmosphereParameters *, const VolatileParameters *, Volatile * );
// FIXME: need to calculate volatile IC using PETSc non-linear solver
//static PetscScalar solve_newton_method( PetscScalar, PetscScalar, PetscScalar );
//static PetscScalar get_newton_f( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
//static PetscScalar get_newton_f_prim( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscErrorCode JSON_add_volatile( DM, Parameters const *, VolatileParameters const *, Volatile const *, Atmosphere const *, char const *name, cJSON * );
static PetscErrorCode JSON_add_atm_struct( Atmosphere *, const AtmosphereParameters *, cJSON * );
static PetscErrorCode set_atm_struct_tau( Atmosphere * );
static PetscErrorCode set_atm_struct_temp( Atmosphere *, const AtmosphereParameters * );
static PetscErrorCode set_atm_struct_pressure( Atmosphere *, const AtmosphereParameters * );
static PetscErrorCode set_atm_struct_depth( Atmosphere *, const AtmosphereParameters * );

PetscErrorCode initialise_atmosphere( Atmosphere *A, const Constants *C )
{
    PetscErrorCode ierr;
    PetscScalar    scaling;

    PetscFunctionBeginUser;

    // pointer to da, so we can easily access it within the atmosphere structure
    const PetscInt stencilWidth = 1;
    const PetscInt dof = 1;
    const PetscInt numpts = 500; // FIXME hard-coded for atmosphere structure output
    ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,numpts,dof,stencilWidth,NULL,&A->da_atm);CHKERRQ(ierr);

    /* create dimensionalisable fields for outputting atmosphere structure */
    /* TODO: this does not really need to be a DimensionalisableField, and PS
       proposed creating a new structure call a DimenisonalisableValue instead */
    scaling = 1.0;
    ierr = DimensionalisableFieldCreate(&A->atm_struct[0],A->da_atm,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(A->atm_struct[0],&A->atm_struct_tau);
    ierr = DimensionalisableFieldSetName(A->atm_struct[0],"atm_struct_tau");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(A->atm_struct[0],"None");CHKERRQ(ierr);

    scaling = C->TEMP;
    ierr = DimensionalisableFieldCreate(&A->atm_struct[1],A->da_atm,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(A->atm_struct[1],&A->atm_struct_temp);
    ierr = DimensionalisableFieldSetName(A->atm_struct[1],"atm_struct_temp");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(A->atm_struct[1],"K");CHKERRQ(ierr);

    scaling = C->PRESSURE / 1.0E5; // bar
    ierr = DimensionalisableFieldCreate(&A->atm_struct[2],A->da_atm,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(A->atm_struct[2],&A->atm_struct_pressure);
    ierr = DimensionalisableFieldSetName(A->atm_struct[2],"atm_struct_pressure");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(A->atm_struct[2],"bar");CHKERRQ(ierr);

    scaling = C->RADIUS;
    ierr = DimensionalisableFieldCreate(&A->atm_struct[3],A->da_atm,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(A->atm_struct[3],&A->atm_struct_depth);
    ierr = DimensionalisableFieldSetName(A->atm_struct[3],"atm_struct_depth");CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(A->atm_struct[3],"m");CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

PetscErrorCode destroy_atmosphere( Atmosphere *A )
{
    PetscErrorCode ierr;
    PetscInt       i;

    PetscFunctionBeginUser;

    for (i=0;i<NUMATMSTRUCTVECS;++i){
        ierr = DimensionalisableFieldDestroy(&A->atm_struct[i]);CHKERRQ(ierr);
    }

    ierr = DMDestroy(&A->da_atm);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_atm_struct_tau( Atmosphere *A )
{

    /* builds an evenly-spaced profile of optical depth from unity at
       the top to the surface value */

    PetscErrorCode    ierr;
    PetscScalar const logtau_min = -6; // FIXME hard-coded here
    PetscScalar const logtau_max = PetscLog10Real(A->tau); // surface optical depth
    PetscScalar       tau,logdtau;
    PetscInt          i,ilo,w,ihi,numpts;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(A->da_atm,NULL,&numpts,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    ierr = DMDAGetCorners(A->da_atm,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;

    logdtau = (logtau_max - logtau_min) / (numpts-1);

    for(i=ilo; i<ihi; ++i){
        tau = PetscPowScalar(10.0,logtau_min + i*logdtau);
        ierr = VecSetValues( A->atm_struct_tau, 1, &i, &tau, INSERT_VALUES );CHKERRQ(ierr);         
    }

    ierr = VecAssemblyBegin( A->atm_struct_tau );CHKERRQ(ierr);
    ierr = VecAssemblyEnd( A->atm_struct_tau );CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_atm_struct_temp( Atmosphere *A, const AtmosphereParameters *Ap )
{
    PetscErrorCode ierr;
    PetscScalar    TO4, Teq4;

    PetscFunctionBeginUser;

    TO4 = A->Fatm / Ap->sigma;
    Teq4 = PetscPowScalar(Ap->teqm,4.0);

    ierr = VecCopy( A->atm_struct_tau, A->atm_struct_temp ); CHKERRQ(ierr);
    ierr = VecShift( A->atm_struct_temp, 1.0 ); CHKERRQ(ierr);
    ierr = VecScale( A->atm_struct_temp, 0.5 ); CHKERRQ(ierr);
    ierr = VecScale( A->atm_struct_temp, TO4 ); CHKERRQ(ierr);
    ierr = VecShift( A->atm_struct_temp, Teq4 ); CHKERRQ(ierr); 
    ierr = VecPow( A->atm_struct_temp, 1.0/4.0 ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_atm_struct_pressure( Atmosphere *A, const AtmosphereParameters *Ap )
{
    PetscErrorCode     ierr;
    PetscScalar        CO2_ratio, CO2_kabs, H2O_ratio, H2O_kabs, kabs;

    PetscFunctionBeginUser;

    /* assume atmosphere is well mixed */
    /* volume mixing ratios */
    /* TODO: store and output these quantities? */
    CO2_ratio = A->CO2.p / (A->CO2.p + A->H2O.p );
    H2O_ratio = 1.0 - CO2_ratio;

    /* effective absorption coefficient */
    CO2_kabs = get_pressure_dependent_kabs( Ap, &Ap->CO2_parameters );
    H2O_kabs = get_pressure_dependent_kabs( Ap, &Ap->H2O_parameters );
    kabs = CO2_ratio * CO2_kabs + H2O_ratio * H2O_kabs;

    ierr = VecCopy( A->atm_struct_tau, A->atm_struct_pressure ); CHKERRQ(ierr);
    ierr = VecScale( A->atm_struct_pressure, -(*Ap->gravity_ptr) ); CHKERRQ(ierr); // note negative gravity
    ierr = VecScale( A->atm_struct_pressure, 2.0 ); CHKERRQ(ierr);
    ierr = VecScale( A->atm_struct_pressure, 1.0/3.0 ); CHKERRQ(ierr);
    ierr = VecScale( A->atm_struct_pressure, 1.0/kabs); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_jeans( const Atmosphere *A, const AtmosphereParameters *Ap, const VolatileParameters *Vp, Volatile *V )
{
    /* surface Jeans parameter */

    PetscFunctionBeginUser;

    V->jeans = 4.0 * PETSC_PI * (*Ap->VOLATILE_ptr) / 1.0E6; // prefactor
    V->jeans *= -(*Ap->gravity_ptr) * (*Ap->radius_ptr) * (Vp->molar_mass/Ap->Avogadro); // note negative gravity
    V->jeans /= Ap->kB * A->tsurf;

    PetscFunctionReturn(0);

}

static PetscErrorCode set_f_thermal_escape( const Atmosphere *A, const AtmosphereParameters *Ap, const VolatileParameters *Vp, Volatile *V )
{
    /* thermal escape prefactor for atmospheric growth rate */

    PetscErrorCode ierr;
    // FIXME: placeholder
    PetscScalar const R = 1.0;

    PetscFunctionBeginUser;

    /* no thermal escape means this is simply unity */
    V->jeans = 0.0; // since not used
    V->f_thermal_escape = 1.0;

    if(Ap->THERMAL_ESCAPE){
        ierr = set_jeans( A, Ap, Vp, V );CHKERRQ(ierr);
        V->f_thermal_escape -= R * (1.0+V->jeans) * PetscExpReal(-V->jeans);
    }

    PetscFunctionReturn(0);

}

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

static PetscErrorCode set_atmosphere_mass( const Atmosphere *A, const AtmosphereParameters *Ap, const VolatileParameters *Vp, Volatile *V )
{
    /* mass of volatile in atmosphere */

    PetscFunctionBeginUser;

    V->m = PetscSqr((*Ap->radius_ptr)) * V->p / -(*Ap->gravity_ptr);
    V->m *= 1.0E6 / (*Ap->VOLATILE_ptr);
    // must weight by molar mass
    V->m *= Vp->molar_mass / A->molar_mass;

    PetscFunctionReturn(0);

}

static PetscErrorCode set_atmosphere_molar_mass( const AtmosphereParameters *Ap, Atmosphere *A )
{
    VolatileParameters   const *CO2_parameters = &Ap->CO2_parameters;
    VolatileParameters   const *H2O_parameters = &Ap->H2O_parameters;
    Volatile                   *CO2 = &A->CO2;
    Volatile                   *H2O = &A->H2O;

    PetscFunctionBeginUser;

    A->molar_mass =  CO2->mixing_ratio * CO2_parameters->molar_mass;
    A->molar_mass += H2O->mixing_ratio * H2O_parameters->molar_mass;

    PetscFunctionReturn(0);

}

static PetscErrorCode set_mixing_ratios( Volatile *CO2, Volatile *H2O )
{
    PetscFunctionBeginUser;

    CO2->mixing_ratio = CO2->p / (CO2->p + H2O->p );
    H2O->mixing_ratio = H2O->p / (CO2->p + H2O->p );

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

    /* H2O */
    ierr = set_partial_pressure_volatile( H2O_parameters, H2O );CHKERRQ(ierr);
    ierr = set_partial_pressure_derivative_volatile( H2O_parameters, H2O );CHKERRQ(ierr);

    /* mixing ratio */
    ierr = set_mixing_ratios( CO2, H2O );CHKERRQ(ierr);

    /* mean molar mass of atmosphere */
    ierr = set_atmosphere_molar_mass( Ap, A );

    /* these terms require the mean molar mass of the atmosphere */
    ierr = set_atmosphere_mass( A, Ap, CO2_parameters, CO2 );CHKERRQ(ierr);
    ierr = set_atmosphere_mass( A, Ap, H2O_parameters, H2O );CHKERRQ(ierr);

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

#if 0
// no longer used, but kept in case it is needed in the future
PetscScalar get_emissivity_from_flux( const Atmosphere *A, const AtmosphereParameters *Ap, PetscScalar flux )
{
    PetscScalar emissivity;

    emissivity = flux / Ap->sigma;
    emissivity /= PetscPowScalar(A->tsurf,4.0)-PetscPowScalar(Ap->teqm,4.0);

    return emissivity;

}
#endif

PetscErrorCode set_surface_temperature_from_flux( Atmosphere *A, const AtmosphereParameters *Ap )
{
    PetscScalar tsurf;

    PetscFunctionBeginUser;
    
    tsurf = A->Fatm / ( Ap->sigma * A->emissivity );
    tsurf += PetscPowScalar( Ap->teqm, 4.0 );
    tsurf = PetscPowScalar( tsurf, 1.0/4.0 );
    A->tsurf = tsurf;

    PetscFunctionReturn(0);

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
    A->tau = CO2->tau + H2O->tau;
    emissivity = 2.0 / (A->tau + 2.0);

    return emissivity;

}

static PetscErrorCode JSON_add_atm_struct( Atmosphere *A, const AtmosphereParameters *Ap, cJSON *data )
{
    PetscErrorCode ierr;
    PetscInt       i;

    PetscFunctionBeginUser;

    /* only compute 1-D atmosphere structure for output */
    /* TODO: if this feedsback into the equations, e.g. through
       atmospheric escape, it will need moving */
    ierr = set_atm_struct_tau( A );CHKERRQ(ierr);
    ierr = set_atm_struct_temp( A, Ap );CHKERRQ(ierr);
    ierr = set_atm_struct_pressure( A, Ap );CHKERRQ(ierr);
    ierr = set_atm_struct_depth( A, Ap ); CHKERRQ(ierr);

    /* write 1-D structure to JSON */
    for (i=0;i<NUMATMSTRUCTVECS;++i){
        cJSON *item;
        DimensionalisableField curr = A->atm_struct[i];
        ierr = DimensionalisableFieldToJSON(curr,&item);CHKERRQ(ierr);
        cJSON_AddItemToObject(data,curr->name,item);
    }

    PetscFunctionReturn(0);

}

static PetscScalar get_pressure_dependent_kabs( const AtmosphereParameters *Ap, const VolatileParameters *Vp )
{
    /* absorption coefficient in the grey atmosphere is pressure-dependent

       Abe and Matsui (1985), Eq. A21, A22, A23

    */

    PetscScalar kabs;

    kabs = PetscSqrtScalar( Vp->kabs * -(*Ap->gravity_ptr) / (3.0*Ap->P0) );

    return kabs;

}

static PetscErrorCode set_optical_depth( const AtmosphereParameters *Ap, const VolatileParameters *Vp, Volatile *V )
{
    PetscFunctionBeginUser;

    V->tau = (3.0/2.0) * V->p / -(*Ap->gravity_ptr);
    V->tau *= get_pressure_dependent_kabs( Ap, Vp );

    PetscFunctionReturn(0);
}

PetscErrorCode JSON_add_atmosphere( DM dm, Parameters const *P, Atmosphere *A, const char *name, cJSON *json )
{
    PetscErrorCode ierr;
    cJSON          *data;
    PetscScalar    scaling, val;
    Constants      const *C = &P->constants;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;

    PetscFunctionBeginUser;

    data = cJSON_CreateObject();

    /* atmosphere structure only relevant for Abe and Matsui (1985) */
    if (Ap->SURFACE_BC==3){
        ierr = JSON_add_atm_struct( A, Ap, data );CHKERRQ(ierr);
    }

    /* total liquid mass of mantle, kg */
    scaling = 4.0 * PETSC_PI * C->MASS; // includes 4*PI for spherical geometry
    ierr = JSON_add_single_value_to_object(dm, scaling, "mass_liquid", "kg", A->Mliq, data);CHKERRQ(ierr);
    /* total solid mass of mantle, kg */
    ierr = JSON_add_single_value_to_object(dm, scaling, "mass_solid", "kg", A->Msol, data);CHKERRQ(ierr);
    /* total mass of mantle, kg (for sanity check) */
    ierr = JSON_add_single_value_to_object(dm, scaling, "mass_mantle", "kg", *Ap->mantle_mass_ptr, data);CHKERRQ(ierr);

    /* kg, without 4*pi */
    val = C->MASS * A->molar_mass * 1.0E3;
    ierr = JSON_add_single_value_to_object(dm, 1.0, "molar_mass", "g/mol", val, data);CHKERRQ(ierr);

    /* surface temperature, K */
    scaling = C->TEMP;
    ierr = JSON_add_single_value_to_object(dm, scaling, "temperature_surface", "K", A->tsurf, data);CHKERRQ(ierr);

    /* optical depth, non-dimensional */
    scaling = 1.0;
    ierr = JSON_add_single_value_to_object(dm, scaling, "optical_depth", "None", A->tau, data);CHKERRQ(ierr);
    /* (effective) emissivity, non-dimensional */
    ierr = JSON_add_single_value_to_object(dm, scaling, "emissivity", "None", A->emissivity, data);CHKERRQ(ierr);

    /* net upward atmospheric flux */
    scaling = C->FLUX;
    ierr = JSON_add_single_value_to_object(dm, scaling, "Fatm", "W m$^{-2}$", A->Fatm, data);CHKERRQ(ierr);

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
    PetscScalar     scaling, val;
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
    ierr = JSON_add_single_value_to_object(dm, scaling, "initial_kg", "kg", VP->initial*(*Ap->mantle_mass_ptr), data);CHKERRQ(ierr);
    /* volatile in liquid mantle (kg) */
    ierr = JSON_add_single_value_to_object(dm, scaling, "liquid_kg", "kg", V->x*A->Mliq, data);CHKERRQ(ierr);
    /* volatile in solid mantle (kg) */
    ierr = JSON_add_single_value_to_object(dm, scaling, "solid_kg", "kg", V->x*VP->kdist*A->Msol, data);CHKERRQ(ierr);
    /* volatile in atmosphere (kg) */
    ierr = JSON_add_single_value_to_object(dm, scaling, "atmosphere_kg", "kg", V->m, data);CHKERRQ(ierr);

    /* kilograms (kg), without 4*pi */
    val = C->MASS * VP->molar_mass * 1.0E3;
    ierr = JSON_add_single_value_to_object(dm, 1.0, "molar_mass", "g/mol", val, data);CHKERRQ(ierr);

    /* bar (bar) */
    /* volatile in atmosphere (bar) */
    scaling = C->PRESSURE / 1.0E5; /* bar */
    ierr = JSON_add_single_value_to_object(dm, scaling, "atmosphere_bar", "bar", V->p, data);CHKERRQ(ierr);

    /* non-dimensional */
    scaling = 1.0;
    ierr = JSON_add_single_value_to_object(dm, scaling, "optical_depth", "None", V->tau, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, scaling, "mixing_ratio", "None", V->mixing_ratio, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, scaling, "jeans", "None", V->jeans, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, scaling, "f_thermal_escape", "None", V->f_thermal_escape, data);CHKERRQ(ierr);

    cJSON_AddItemToObject(json,name,data);

    PetscFunctionReturn(0);

}

PetscScalar get_dxdt( const AtmosphereParameters *Ap, const Atmosphere *A, const VolatileParameters *Vp, Volatile *V )
{

    PetscErrorCode ierr;
    PetscScalar    dxdt;
    PetscScalar    num, den;

    /* to update V->f_thermal_escape */
    ierr = set_f_thermal_escape( A, Ap, Vp, V ); CHKERRQ(ierr);

    num = V->x * (Vp->kdist-1.0) * A->dMliqdt;
    den = Vp->kdist * (*Ap->mantle_mass_ptr) + (1.0-Vp->kdist) * A->Mliq;
    den += (1.0E6 / (*Ap->VOLATILE_ptr)) * PetscSqr( (*Ap->radius_ptr)) * (V->dpdx / -(*Ap->gravity_ptr)) * (Vp->molar_mass / A->molar_mass) * V->f_thermal_escape;

    dxdt = num / den;

    return dxdt;
}

static PetscScalar get_dzdtau( PetscScalar tau, const AtmosphereParameters *Ap, const Atmosphere *A )
{
    PetscScalar dzdt;

    dzdt = A->Fatm / Ap->sigma; // T0**4
    dzdt *= (tau+1.0)/2.0;
    dzdt += PetscPowScalar(Ap->teqm,4.0);
    dzdt = PetscPowScalar(dzdt,1.0/4.0);
    dzdt /= tau;
    dzdt /= (*Ap->gravity_ptr) * A->molar_mass / Ap->Rgas;

    return dzdt;

}

static PetscScalar get_z_from_simpson(PetscScalar x, PetscScalar dx, const AtmosphereParameters *Ap, const Atmosphere *A)
{
    PetscScalar fa, fb, fm;

    /* Because f=dz/dtau is only a function of tau (not also z),
       RK4 reduces to Simpson's method */

    fa = get_dzdtau( x, Ap, A ); // start
    fb = get_dzdtau( x + dx, Ap, A) ; // end
    fm = get_dzdtau( x + 0.5*dx, Ap, A ); // midpoint

    return dx/6.0 * (fa + 4.0*fm + fb );

}

static PetscErrorCode set_atm_struct_depth( Atmosphere *A, const AtmosphereParameters *Ap )
{
    PetscErrorCode ierr;
    PetscScalar    *arr_tau, *arr_depth, tau, dtau, val;
    PetscInt       i, numpts;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(A->da_atm,NULL,&numpts,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    ierr = DMDAVecGetArrayRead(A->da_atm,A->atm_struct_tau,&arr_tau);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(A->da_atm,A->atm_struct_depth,&arr_depth);CHKERRQ(ierr);

    /* loop backwards, since height above the surface (z=0) occurs 
       at the largest optical depth */
    arr_depth[numpts-1] = 0.0;

    for(i=numpts-1; i>0; --i){
        tau = arr_tau[i];
        dtau = arr_tau[i-1] - arr_tau[i]; // negative
        val = get_z_from_simpson( tau, dtau, Ap, A );
        arr_depth[i-1] = arr_depth[i] + val;
    }

    ierr = DMDAVecRestoreArrayRead(A->da_atm,A->atm_struct_tau,&arr_tau);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(A->da_atm,A->atm_struct_depth,&arr_depth);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

PetscScalar get_initial_volatile_abundance( Atmosphere *A, const AtmosphereParameters *Ap, const VolatileParameters *Vp, const Volatile *V )
{
    PetscScalar out;

    out = 1.0E6 / *Ap->VOLATILE_ptr;
    out *= PetscSqr( (*Ap->radius_ptr) ) / -(*Ap->gravity_ptr);
    out /= (*Ap->mantle_mass_ptr);
    out *= (Vp->molar_mass / A->molar_mass);
    out *= V->p;
    out += V->x;
    out -= Vp->initial;

    return out;

}

#if 0
PetscScalar get_initial_volatile( const AtmosphereParameters *Ap, const VolatileParameters *Vp )
{
    /* initial volatile in the aqueous phase */

    /* FIXME: below is old, when it was easy to use Newton's method per individual volatile.
       Now the volatiles are coupled, and this must be solved as a system of non-linear
       equations instead */
    /*
    PetscScalar fac, x;

    fac = PetscSqr( (*Ap->radius_ptr) );
    fac /= -(*Ap->gravity_ptr) * (*Ap->mantle_mass_ptr);
    fac /= PetscPowScalar(Vp->henry,Vp->henry_pow);
    fac *= 1.0E6 / (*Ap->VOLATILE_ptr);

    x = solve_newton_method( fac, Vp->henry_pow, Vp->initial );

    return x;*/
}
#endif

/* FIXME: below needs replacing with PETSc non-linear solver */
#if 0
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
#endif
