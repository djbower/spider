#include "atmosphere.h"
#include "dimensionalisablefield.h"
#include "reaction.h"
#include "util.h"

static PetscErrorCode initialise_volatiles( Atmosphere *, const AtmosphereParameters * );
static PetscErrorCode set_atmosphere_pressures( Atmosphere *, const AtmosphereParameters *, const Constants * );
static PetscErrorCode set_volume_mixing_ratios( Atmosphere *, const AtmosphereParameters * );
static PetscErrorCode set_volatile_masses_in_atmosphere( Atmosphere *, const AtmosphereParameters * );
static PetscErrorCode set_escape( Atmosphere *, const AtmosphereParameters * );
static PetscScalar get_pressure_dependent_kabs( const AtmosphereParameters *, PetscInt );
static PetscErrorCode set_jeans( Atmosphere *, const AtmosphereParameters *, PetscInt );
static PetscErrorCode set_column_density_volatile( Atmosphere *, const AtmosphereParameters *, PetscInt );
static PetscErrorCode set_Knudsen_number( Atmosphere *, const AtmosphereParameters *, PetscInt );
static PetscErrorCode set_R_thermal_escape( Atmosphere *, const AtmosphereParameters *, PetscInt );
static PetscErrorCode set_f_thermal_escape( Atmosphere *, PetscInt );
static PetscErrorCode JSON_add_volatile( DM, Parameters const *, VolatileParameters const *, Volatile const *, Atmosphere const *, char const *name, cJSON * );
static PetscErrorCode JSON_add_atm_struct( Atmosphere *, const AtmosphereParameters *, cJSON * );
static PetscErrorCode set_atm_struct_tau( Atmosphere * );
static PetscErrorCode set_atm_struct_temp( Atmosphere *, const AtmosphereParameters * );
static PetscErrorCode set_atm_struct_pressure( Atmosphere *, const AtmosphereParameters * );
static PetscErrorCode set_atm_struct_depth( Atmosphere *, const AtmosphereParameters * );

PetscErrorCode initialise_atmosphere( Atmosphere *A, const AtmosphereParameters *Ap, const Constants *C )
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

    /* initialise volatiles */
    ierr = initialise_volatiles( A, Ap );

    /* other variables in struct that otherwise might not get set */
    /* below is only for Abe and Matsui atmosphere model */
    A->tau = 0.0;
    /* Ensure that A->tsurf is non-zero so that the escape-related
       functions do not return an uninitialised value warning
       (valgrind).  During the time loop, A->tsurf is updated to a
       meaningful value before it is actually used (as is A->psurf) */
    A->tsurf = 1.0;
    A->dtsurfdt = 0.0;
    A->psurf = 0.0;

    /* Initialize mass_reaction temporary field (maybe not needed) to zero */
    {
      PetscInt i;

      for (i=0; i<SPIDER_MAX_REACTIONS; ++i) A->mass_reaction[i] = 0.0;
    }

    PetscFunctionReturn(0);
}

static PetscErrorCode initialise_volatiles( Atmosphere *A, const AtmosphereParameters *Ap)
{
    PetscInt i;

    PetscFunctionBeginUser;

    for (i=0; i<Ap->n_volatiles; ++i) {
        /* these entries should match those in atmosphere.h */
        A->volatiles[i].x = 0.0;
        A->volatiles[i].p = 0.0;
        A->volatiles[i].dxdt = 0.0;
        A->volatiles[i].dpdx = 0.0;
        A->volatiles[i].m = 0.0;
        A->volatiles[i].tau = 0.0;
        A->volatiles[i].mixing_ratio = 0.0;
        A->volatiles[i].column_density = 0.0;
        A->volatiles[i].Knudsen = 0.0;
        A->volatiles[i].jeans = 0.0;
        A->volatiles[i].f_thermal_escape = 0.0;
        A->volatiles[i].R_thermal_escape = 0.0;
    }

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
    PetscInt           i;
    PetscScalar        kabs_tot,kabs;

    PetscFunctionBeginUser;

    /* effective absorption coefficient */
    /* TODO: this was formulated for 2 species, and presumably extends for n species as follows,
       but should check */
    kabs_tot = 0.0;
    for( i=0; i<Ap->n_volatiles; ++i) {
        kabs = get_pressure_dependent_kabs( Ap, i );
        kabs_tot += A->volatiles[i].mixing_ratio * kabs;
    }

    ierr = VecCopy( A->atm_struct_tau, A->atm_struct_pressure ); CHKERRQ(ierr);
    ierr = VecScale( A->atm_struct_pressure, -(*Ap->gravity_ptr) ); CHKERRQ(ierr); // note negative gravity
    ierr = VecScale( A->atm_struct_pressure, 2.0 ); CHKERRQ(ierr);
    ierr = VecScale( A->atm_struct_pressure, 1.0/3.0 ); CHKERRQ(ierr);
    ierr = VecScale( A->atm_struct_pressure, 1.0/kabs_tot); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_jeans( Atmosphere *A, const AtmosphereParameters *Ap, PetscInt i )
{
    /* surface Jeans parameter */

    Volatile                  *V = &A->volatiles[i];
    VolatileParameters const *Vp = &Ap->volatile_parameters[i];

    PetscFunctionBeginUser;

    if(Vp->jeans_value < 0.0){
        V->jeans = -(*Ap->gravity_ptr) * (*Ap->radius_ptr) * (Vp->molar_mass/Ap->Avogadro); // note negative gravity
        V->jeans /= Ap->kB * A->tsurf;
    }
    else{
        V->jeans = Vp->jeans_value;
    }

    PetscFunctionReturn(0);

}

static PetscErrorCode set_column_density_volatile( Atmosphere *A, const AtmosphereParameters *Ap, PetscInt i )
{
    /* see Johnson et al. (2015), Astrophys. J. */

    PetscFunctionBeginUser;

    Volatile                 *V = &A->volatiles[i];
    VolatileParameters const *Vp = &Ap->volatile_parameters[i];

    V->column_density = V->p;
    V->column_density /= -(*Ap->gravity_ptr) * (Vp->molar_mass/Ap->Avogadro);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_Knudsen_number( Atmosphere *A, const AtmosphereParameters *Ap, PetscInt i )
{
    /* see Johnson et al. (2015), Astrophys. J. */

    Volatile                  *V = &A->volatiles[i];
    VolatileParameters const *Vp = &Ap->volatile_parameters[i];

    PetscFunctionBeginUser;

    V->Knudsen = V->jeans * Vp->cross_section * V->column_density;
    V->Knudsen = 1.0 / V->Knudsen;

    PetscFunctionReturn(0);

}

static PetscErrorCode set_R_thermal_escape( Atmosphere *A, const AtmosphereParameters *Ap, PetscInt i )
{
    /* see Johnson et al. (2015), Astrophys. J. */

    PetscScalar              R1, R2, Rfit;
    Volatile                 *V = &A->volatiles[i];
    VolatileParameters const *Vp = &Ap->volatile_parameters[i];

    PetscFunctionBeginUser;

    if(Vp->R_thermal_escape_value < 0.0 ){
        R1 = PetscPowScalar( 1.0 / V->Knudsen, 0.09 );
        R2 = 70 * (V->Knudsen * PetscExpReal(V->jeans)) / PetscPowScalar(V->jeans,2.55);
        Rfit = 1.0 / R1 + 1.0 / R2;
        Rfit = 1.0 / Rfit;
        V->R_thermal_escape = Rfit;
    }
    else{
        V->R_thermal_escape = Vp->R_thermal_escape_value;
    }

    PetscFunctionReturn(0);

}

static PetscErrorCode set_f_thermal_escape( Atmosphere *A, PetscInt i )
{
    /* thermal escape prefactor for atmospheric growth rate */
    /* see Johnson et al. (2015), Astrophys. J. */

    Volatile    *V = &A->volatiles[i];

    PetscFunctionBeginUser;

    V->f_thermal_escape = 1.0 + V->R_thermal_escape * (1.0+V->jeans) * PetscExpReal(-V->jeans);

    /* TODO: do we need any extra checks to ensure that the asymptotic behaviour
       of escape is reasonable?  As jeans-->infty the limit looks OK, but what about
       as jeans-->0?  In reality, this would denote a switch to hydrodynamic
       escape */

    PetscFunctionReturn(0);

}

static PetscErrorCode set_atmosphere_pressures( Atmosphere *A, const AtmosphereParameters *Ap, const Constants *C )
{
    PetscInt                  i;
    Volatile                 *V;
    VolatileParameters const *Vp;

    PetscFunctionBeginUser;

    /* total surface pressure */
    A->psurf = 0.0;

    for (i=0; i<Ap->n_volatiles; ++i) {
        /* partial pressure of volatile */
        V = &A->volatiles[i];
        Vp = &Ap->volatile_parameters[i];
        V->p = PetscPowScalar( V->x / Vp->henry, Vp->henry_pow );
        /* derivative of partial pressure of volatile */
        V->dpdx = Vp->henry_pow / Vp->henry;
        V->dpdx *= PetscPowScalar( V->x / Vp->henry, Vp->henry_pow-1.0 );
        /* total pressure by Dalton's law */
        A->psurf += V->p;
    }

    /* Oxygen is often treated as a trace species from the
       perspective of computing the total atmospheric pressure, but
       it is trivial to include: */
    /* Actually, for oxidised meteorite material this is almost
       certainly required, since O2 could be a dominant species */
    /* fO2 is set in set_interior_structure_from_solution */
    A->psurf /= (1.0 - A->fO2);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_volatile_masses_in_atmosphere( Atmosphere *A, const AtmosphereParameters *Ap )
{
    /* mass of volatiles in atmosphere */

    PetscInt i;

    PetscFunctionBeginUser;

    for (i=0; i<Ap->n_volatiles; ++i) {
        A->volatiles[i].m = PetscSqr((*Ap->radius_ptr)) * A->volatiles[i].p / -(*Ap->gravity_ptr);
        A->volatiles[i].m *= 1.0E6 / (*Ap->VOLATILE_ptr);
        A->volatiles[i].m *= Ap->volatile_parameters[i].molar_mass / A->molar_mass;
    }

    PetscFunctionReturn(0);

}

static PetscErrorCode set_volume_mixing_ratios( Atmosphere *A, const AtmosphereParameters *Ap )
{
    PetscInt i;

    PetscFunctionBeginUser;

    A->molar_mass = 0.0;

    /* now compute mixing ratios */
    for (i=0; i<Ap->n_volatiles; ++i) {
        /* mixing ratio */
        A->volatiles[i].mixing_ratio = A->volatiles[i].p / A->psurf;
        /* atmosphere molar mass */
        A->molar_mass += Ap->volatile_parameters[i].molar_mass * A->volatiles[i].mixing_ratio;
    }

    PetscFunctionReturn(0);

}

PetscErrorCode set_atmosphere_volatile_content( Atmosphere *A, const AtmosphereParameters *Ap, const Constants *C )
{
    PetscErrorCode           ierr;

    PetscFunctionBeginUser;

    /* if x is zero, the quantities below will also all
       be set to zero, except the escape-related parameters that
       are useful to compute even if the feedback of escape is not
       actually included in the model */

    /* order of these functions is very important! */
    ierr = set_atmosphere_pressures( A, Ap, C );CHKERRQ(ierr);

    ierr = set_volume_mixing_ratios( A, Ap );CHKERRQ(ierr);

    ierr = set_volatile_masses_in_atmosphere( A, Ap );CHKERRQ(ierr);

    /* escape is always set, since then we can easily see the
       variables, regardless of whether the feedback is actually
       included for the volatile evolution */
    ierr = set_escape( A, Ap );CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_escape( Atmosphere *A, const AtmosphereParameters *Ap )
{
    PetscErrorCode ierr;
    PetscInt       i;

    for (i=0; i<Ap->n_volatiles; ++i) {

        ierr = set_column_density_volatile( A, Ap, i );CHKERRQ(ierr);

        ierr = set_jeans( A, Ap, i );CHKERRQ(ierr);

        ierr = set_Knudsen_number( A, Ap, i ); CHKERRQ(ierr);

        ierr = set_R_thermal_escape( A, Ap, i ); CHKERRQ(ierr);

        ierr = set_f_thermal_escape( A, i ); CHKERRQ(ierr);
    }

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

PetscScalar get_emissivity_abe_matsui( Atmosphere *A, const AtmosphereParameters *Ap )
{
    PetscInt    i;
    PetscScalar emissivity;

    A->tau = 0.0; // total surface optical depth

    /* now compute mixing ratios */
    for (i=0; i<Ap->n_volatiles; ++i) {
        Volatile                 *V = &A->volatiles[i];
        /* optical depth at surface for this volatile */
        V->tau = (3.0/2.0) * V->p / -(*Ap->gravity_ptr);
        V->tau *= get_pressure_dependent_kabs( Ap, i );
        /* total optical depth at surface */
        A->tau += V->tau;

    }

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

static PetscScalar get_pressure_dependent_kabs( const AtmosphereParameters *Ap, PetscInt i )
{
    /* absorption coefficient in the grey atmosphere is pressure-dependent

       Abe and Matsui (1985), Eq. A21, A22, A23

    */

    PetscScalar kabs;

    kabs = PetscSqrtScalar( Ap->volatile_parameters[i].kabs * -(*Ap->gravity_ptr) / (3.0*Ap->P0) );

    return kabs;

}

PetscErrorCode JSON_add_atmosphere( DM dm, Parameters const *P, Atmosphere *A, const char *name, cJSON *json )
{
    PetscErrorCode ierr;
    cJSON          *data;
    PetscScalar    scaling, val;
    Constants      const *C = &P->constants;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    PetscInt       v;

    PetscFunctionBeginUser;

    data = cJSON_CreateObject();

    /* atmosphere structure relevant for case Abe and Matsui (1985) */
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

    /* surface pressure, bar */
    scaling = C->PRESSURE / 1.0E5; /* bar */
    ierr = JSON_add_single_value_to_object(dm, scaling, "pressure_surface", "bar", A->psurf, data);CHKERRQ(ierr);

    /* oxygen fugacity */
    /* the value itself is the non-dimensional oxygen fugacity, a la Schaefer and Fegley (2017) */
    /* multiplying by the scaling gives the oxygen fugacity (partial pressure) in bar */
    scaling = A->psurf * C->PRESSURE / 1.0E5; /* bar */
    ierr = JSON_add_single_value_to_object(dm, scaling, "fO2", "bar", A->fO2, data);CHKERRQ(ierr);

    /* optical depth, non-dimensional */
    scaling = 1.0;
    ierr = JSON_add_single_value_to_object(dm, scaling, "optical_depth", "None", A->tau, data);CHKERRQ(ierr);
    /* (effective) emissivity, non-dimensional */
    ierr = JSON_add_single_value_to_object(dm, scaling, "emissivity", "None", A->emissivity, data);CHKERRQ(ierr);

    /* net upward atmospheric flux */
    scaling = C->FLUX;
    ierr = JSON_add_single_value_to_object(dm, scaling, "Fatm", "W m$^{-2}$", A->Fatm, data);CHKERRQ(ierr);

    /* Volatiles */
    for (v=0; v<Ap->n_volatiles; ++v) {
      ierr = JSON_add_volatile(dm, P, &Ap->volatile_parameters[v], &A->volatiles[v], A, Ap->volatile_parameters[v].prefix, data ); CHKERRQ(ierr);
    }

    /* Reactions */
    // (Not currently output)

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

    /* area */
    scaling = 1.0 / (C->AREA * 1.0E4); // 1/cm^2
    ierr = JSON_add_single_value_to_object(dm, scaling, "column_density", "1/cm^2", V->column_density, data);CHKERRQ(ierr);

    /* non-dimensional */
    scaling = 1.0;
    ierr = JSON_add_single_value_to_object(dm, scaling, "optical_depth", "None", V->tau, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, scaling, "mixing_ratio", "None", V->mixing_ratio, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, scaling, "jeans", "None", V->jeans, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, scaling, "Knudsen", "None", V->Knudsen, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, scaling, "R_thermal_escape", "None", V->R_thermal_escape, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, scaling, "f_thermal_escape", "None", V->f_thermal_escape, data);CHKERRQ(ierr);

    cJSON_AddItemToObject(json,name,data);

    PetscFunctionReturn(0);

}

PetscErrorCode FormFunction2( SNES snes, Vec x, Vec f, void *ptr)
{
    PetscErrorCode             ierr;
    const PetscScalar          *xx;
    PetscScalar                *ff;
    Ctx                        *E = (Ctx*) ptr;
    PetscInt                   i;
    const PetscScalar          *dmrdt;
    Atmosphere                 *A = &E->atmosphere;
    Parameters           const *P = &E->parameters;
    Constants            const *C = &P->constants;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;

    PetscFunctionBeginUser;

    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    for (i=0; i<Ap->n_volatiles; ++i) {
        A->volatiles[i].dxdt = xx[i];
    }
    dmrdt = &xx[Ap->n_volatiles];
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);

    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);
    for (i=0; i<Ap->n_volatiles; ++i) {
        ff[i] = get_dxdt( A, Ap, i, dmrdt );
    }

    /* chemical equilibrium constraints */
    for (i=0; i<Ap->n_reactions; ++i) {
        PetscScalar dQpdt, dQrdt, K, Qr, dKdT, dKdt;
        Qr = get_reaction_quotient_reactants( &Ap->reaction_parameters[i], A );
        dQpdt = get_reaction_quotient_products_time_derivative( &Ap->reaction_parameters[i], A, Ap );
        dQrdt = get_reaction_quotient_reactants_time_derivative( &Ap->reaction_parameters[i], A, Ap );
        K = get_equilibrium_constant( &Ap->reaction_parameters[i], A->tsurf, C );
        dKdT = get_equilibrium_constant_temperature_derivative( &Ap->reaction_parameters[i], A->tsurf, C );
        dKdt = dKdT * A->dtsurfdt;
        /* residual of reaction balance */
        //ff[Ap->n_volatiles + i] = 0.0; // for debugging
        ff[Ap->n_volatiles + i] = dQpdt - K * dQrdt - Qr * dKdt;
    }

#if 0
    /* Objective function */
    for (i=0; i<Ap->n_reactions; ++i) {
      PetscBool is_simple;

      ierr = PetscStrcmp(Ap->reaction_parameters[i]->type,"simple",&is_simple);CHKERRQ(ierr);
      if (is_simple) {
        const PetscInt v0 = Ap->reaction_parameters[i]->volatiles[0];
        const PetscInt v1 = Ap->reaction_parameters[i]->volatiles[1];

        ff[Ap->n_volatiles + i] = Ap->reaction_parameters[i]->epsilon[v0] * A->volatiles[v0].dxdt * A->volatiles[v0].dpdx
                                + Ap->reaction_parameters[i]->epsilon[v1] * A->volatiles[v1].dxdt * A->volatiles[v1].dpdx;
      } else {
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Reaction type %s not recognized",Ap->reaction_parameters[i]->type);
      }
    }

    /* FIXME: this turns off all reactions, but allows test cases to run */
    /* TODO: DJB needs to figure out and update the objective function for the time stepper */
    /* Objective function */
    for (i=0; i<Ap->n_reactions; ++i) {
        ff[Ap->n_volatiles + i] = 0.0;
    }
#endif


    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscScalar get_dxdt( Atmosphere *A, const AtmosphereParameters *Ap, PetscInt i, const PetscScalar *dmrdt )
{

    PetscScalar               out, out2, dpsurfdt, f_thermal_escape;
    PetscInt                  j,k;

    /* remember that to this point, V->f_thermal_escape is always
       computed but not necessarily used in the calculation */
    if(Ap->THERMAL_ESCAPE){
        f_thermal_escape = A->volatiles[i].f_thermal_escape;
    }
    else{
        f_thermal_escape = 1.0;
    }

    out2 = 0.0;

    /* dPsurf/dt */
    dpsurfdt = 0.0;
    for (k=0; k<Ap->n_volatiles; ++k) {
        dpsurfdt += A->volatiles[k].dpdx * A->volatiles[k].dxdt;
    }

    for (j=0; j<Ap->n_volatiles; ++j) {
        out = 0.0;
        out = -dpsurfdt * A->volatiles[j].p / A->psurf;
        out += A->volatiles[j].dpdx * A->volatiles[j].dxdt;
        out *= Ap->volatile_parameters[j].molar_mass;
        out2 += out;
    }

    out2 *= -A->volatiles[i].p / (A->psurf * PetscSqr(A->molar_mass));

    /* second part of atmosphere derivative */
    out2 += ( 1.0 / A->molar_mass ) * A->volatiles[i].dpdx * A->volatiles[i].dxdt;

    /* multiply by prefactors */
    out2 *= (1.0E6 / (*Ap->VOLATILE_ptr)) * PetscSqr(*Ap->radius_ptr) * Ap->volatile_parameters[i].molar_mass / -(*Ap->gravity_ptr); // note negative gravity

    /* thermal escape correction */
    out2 *= f_thermal_escape;

    /* chemical reactions. Loop through all reactions, check if they involve
      this volatile, and if so , subtract a term*/
    for (j=0; j<Ap->n_reactions; ++j) {
      for (k=0; k<Ap->reaction_parameters[j]->n_volatiles; ++k) {
        if (k==i) {
          // FIXME: STOICHIOMETRY FIELD IS WRONG
          out2 -= Ap->reaction_parameters[j]->stoichiometry[k] * dmrdt[j];
        }
      }
    }

    /* solid and liquid reservoirs */
    out2 += A->volatiles[i].dxdt * ( Ap->volatile_parameters[i].kdist * (*Ap->mantle_mass_ptr) + (1.0-Ap->volatile_parameters[i].kdist) * A->Mliq);
    out2 += A->volatiles[i].x * (1.0-Ap->volatile_parameters[i].kdist) * A->dMliqdt;

    return out2;
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

PetscScalar get_initial_volatile_abundance( Atmosphere *A, const AtmosphereParameters *Ap, const VolatileParameters *Vp, const Volatile *V)
{
    PetscScalar out;

    /* atmosphere */
    out = 1.0E6 / *Ap->VOLATILE_ptr;
    out *= PetscSqr( (*Ap->radius_ptr) ) / -(*Ap->gravity_ptr);
    out *= (Vp->molar_mass / A->molar_mass);
    out *= V->p;
    /* liquid and solid */
    out += V->x * ( Vp->kdist * A->Msol + A->Mliq );
    out -= Vp->initial * (*Ap->mantle_mass_ptr);
    /* Note: we do not include reaction masses here */
    return out;
}

PetscErrorCode set_oxygen_fugacity( Atmosphere *A, const AtmosphereParameters *Ap, const Constants *C )
{

    /* These are oxygen fugacity fits for individual meteoritic materials as
       as function of temperature from Schaefer and Fegley (2017), ApJ.
       Also, fits to IW buffer from Olson and Sharp (2019) based on Ebel and Grossman (2000) */

    PetscScalar a,b,c,d,f,func,dfuncdT;
    /* temperature must be dimensional in K */
    PetscScalar temp = A->tsurf * C->TEMP;

    PetscFunctionBeginUser;

    switch( Ap->OXYGEN_FUGACITY ){
        case 1:
            /* Schaefer and Fegley (2017), CI meteorite */
            a = 2.4976;
            b = -9.8605;
            c = -17.0701;
            d = 7.5220;
            f = -1.0404;
            break;
        case 2:
            /* Schaefer and Fegley (2017), CV meteorite */
            a = 9.0621;
            b = -31.193;
            c = 5.1092;
            d = -1.8475;
            f = 0.2000;
            break;
        case 3:
            /* Schaefer and Fegley (2017), H meteorite */
            a = 5.0743;
            b = -22.906;
            c = -5.6610;
            d = 2.0634;
            f = -0.2618;
            break;
        case 4:
            /* Schaefer and Fegley (2017), EH meteorite */
            a = 4.9495;
            b = -24.024;
            c = -4.6236;
            d = 1.7177;
            f = -0.2332;
            break;
        case 5:
            /* Schaefer and Fegley (2017), Eucrite meteorite */
            a = 5.4856;
            b = -25.127;
            c = -3.6580;
            d = 1.3014;
            f = -0.1650;
            break;
        case 6:
            /* Olson and Sharp (2019), IW buffer */
            a = -2.75E6;
            b = -1.7;
            c = 0.0;
            break;
        case 7:
            /* Olson and Sharp (2019), IW-1 */
            a = -2.75E6;
            b = -1.7;
            c = -1.0;
            break;
        case 8:
            /* Olson and Sharp (2019), IW-2 */
            a = -2.75E6;
            b = -1.7;
            c = -2.0;
            break;
        default:
            SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported OXYGEN_FUGACITY value %d provided",Ap->OXYGEN_FUGACITY);
            break;
        }

    /* Schaefer and Fegley (2017) */
    if( Ap->OXYGEN_FUGACITY <= 5 ){
        func = a + b*1E3/temp + c*1E6/PetscPowScalar(temp,2.0) + d*1E9/PetscPowScalar(temp,3.0) + f*1E12/PetscPowScalar(temp,4.0);
        dfuncdT = -b*1E3/PetscPowScalar(temp,2.0) - 2*c*1E6/PetscPowScalar(temp,3.0) - 3*d*1E9/PetscPowScalar(temp,4.0) - 4*f*1E12/PetscPowScalar(temp,5.0);
    }
    /* Olson and Sharp (2019) fits to Ebel and Grossman (2000) */
    else{
        func = a * PetscPowScalar( temp, b ) + c;
        dfuncdT = a * b * PetscPowScalar( temp, b-1 );
    }

    /* Remember that oxygen_fugacity is equivalent to a volume
       mixing ratio, and therefore does not need scaling */
    A->fO2 = PetscPowScalar(10.0,func);
    A->dfO2dT = A->fO2 * PetscLogReal(10.0) * dfuncdT;

    PetscFunctionReturn(0);

}
