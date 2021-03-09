#include "atmosphere.h"
#include "bc.h"
#include "energy.h"
#include "eos.h"
#include "ic.h"
#include "matprop.h"
#include "monitor.h"
#include "parameters.h"
#include "reaction.h"
#include "util.h"

/* interior ic */
static PetscErrorCode set_ic_interior( Ctx *, Vec );
static PetscErrorCode set_ic_interior_default( Ctx *, Vec );
static PetscErrorCode set_ic_interior_entropy( Ctx *, Vec );
static PetscErrorCode set_ic_interior_from_file( Ctx *, Vec );
static PetscErrorCode set_ic_interior_from_phase_boundary( Ctx *, Vec );
static PetscErrorCode set_start_time_from_file( Parameters , const char * );
/* atmosphere ic */
static PetscErrorCode set_ic_atmosphere( Ctx *, Vec );
static PetscErrorCode set_ic_atmosphere_from_ocean_moles( Ctx *E, Vec sol );
static PetscErrorCode set_ic_atmosphere_from_initial_total_abundance( Ctx *E, Vec sol );
static PetscErrorCode set_ic_atmosphere_from_file( Ctx *, Vec );
static PetscErrorCode set_ic_atmosphere_from_partial_pressure( Ctx *, Vec );
static PetscErrorCode solve_for_initial_partial_pressure( Ctx * );
/* general */
static PetscErrorCode set_ic_from_file( Ctx *, Vec, const char *, const PetscInt *, PetscInt );
static PetscErrorCode objective_function_initial_partial_pressure( SNES, Vec, Vec, void *);
static PetscErrorCode conform_atmosphere_parameters_to_ic( Ctx * );
static PetscErrorCode print_ocean_masses( Ctx * );

/* main function to set initial condition */
PetscErrorCode set_initial_condition( Ctx *E, Vec sol)
{
    PetscErrorCode ierr;
    Parameters const P = E->parameters;
    Atmosphere      *A = &E->atmosphere;
    AtmosphereParameters Ap = P->atmosphere_parameters;
    FundamentalConstants FC = P->fundamental_constants;
    ScalingConstants     SC = P->scaling_constants;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_initial_condition()\n");CHKERRQ(ierr);

    /* interior initial condition */
    ierr = set_ic_interior( E, sol ); CHKERRQ(ierr);

    /* need some interface quantities for setting atmosphere ic */
    ierr = set_interior_atmosphere_interface_from_surface_entropy( E ); CHKERRQ(ierr);

    /* atmosphere initial condition */
    /* must be after interior initial condition */
    ierr = set_ic_atmosphere( E, sol ); CHKERRQ(ierr);

    /* for surface bc, need an estimate of A->Fatm */
    ierr = set_atmosphere_emissivity_and_flux( A, Ap, FC, SC );

    /* P->t0 is set in parameters.c */
    /* time is needed for the radiogenic energy input */
    ierr = set_current_state_from_solution( E, P->t0, sol ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* initial condition of interior */

static PetscErrorCode set_ic_interior( Ctx *E, Vec sol)
{
    /* these functions generally work on the Vec sol directly */

    PetscErrorCode ierr;
    Parameters const P = E->parameters;
    PetscInt IC = P->IC_INTERIOR;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior()\n");CHKERRQ(ierr);

    if (IC==1){
        ierr = set_ic_interior_default( E, sol ); CHKERRQ(ierr);
    }
    else if (IC==2){
        ierr = set_ic_interior_from_file( E, sol ); CHKERRQ(ierr);
        ierr = set_start_time_from_file( P, P->ic_interior_filename );CHKERRQ(ierr);
    }
    else if (IC==3){
        ierr = set_ic_interior_from_phase_boundary( E, sol ); CHKERRQ(ierr);
    }
    else{
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported IC_INTERIOR value %d provided",P->IC_INTERIOR);
    }

    /* this copies the sol Vecs to E */
    /* the reason we set sol above, rather than Vecs in struct, is
       because the function below is called by the RHS and must take
       sol as the input.  This now sets the Vecs in E */
    ierr = set_entropy_from_solution(E, sol); CHKERRQ(ierr);

    /* set surface and core entropy for the IC (if set), but then can
       evolve (not necessarily isothermal).  Updates the Vecs in E */
    ierr = set_boundary_entropy_constant( E ); CHKERRQ(ierr);

    /* Now the Vecs in E are set and consistent, clone them
       to the sol Vec which is what is actually used by the
       time-stepper */
    ierr = set_solution_from_entropy(E, sol);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_interior_default( Ctx *E, Vec sol )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior_default()\n");CHKERRQ(ierr);

    ierr = set_ic_interior_entropy( E, sol ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_interior_from_file( Ctx *E, Vec sol )
{
    PetscErrorCode ierr;
    Parameters const P = E->parameters;
    /* subdomains to use, i.e. dS/dxi and S0 */
    PetscInt const arr[2] = {0, 1};

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior_from_file()\n");CHKERRQ(ierr);

    ierr = set_ic_from_file( E, sol, P->ic_interior_filename, arr, 2 ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_atmosphere_from_file( Ctx *E, Vec sol )
{
    PetscErrorCode ierr;
    Parameters const P = E->parameters;
    AtmosphereParameters const Ap = P->atmosphere_parameters;
    /* subdomains to use, i.e. volatile abundances */
    PetscInt const arr1[1] = {2};
    PetscInt const arr2[1] = {3};

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_atmosphere_from_file()\n");CHKERRQ(ierr);

    if(Ap->n_volatiles){
        ierr = set_ic_from_file( E, sol, Ap->ic_atmosphere_filename, arr1, 1 ); CHKERRQ(ierr);
    }

    if(Ap->n_reactions){
        ierr = set_ic_from_file( E, sol, Ap->ic_atmosphere_filename, arr2, 1 ); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_interior_entropy( Ctx *E, Vec sol )
{
    /* set initial entropy gradient to constant for all basic nodes */

    PetscErrorCode   ierr;
    PetscInt         i;
    PetscScalar      S0;
    Parameters const P = E->parameters;
    Mesh       const *M  = &E->mesh;
    Vec              dSdxi_b;
    Vec              *subVecs;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior_entropy()\n");CHKERRQ(ierr);

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    dSdxi_b = subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDXI_B]];

    /* set entropy gradient and map to mass coordinates */
    ierr = VecSet(dSdxi_b, P->ic_dsdr );CHKERRQ(ierr);
    ierr = VecPointwiseDivide(dSdxi_b,dSdxi_b,M->dxidr_b);

    /* set entropy at top of adiabat */
    S0 = P->ic_adiabat_entropy;
    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],0,S0,INSERT_VALUES);CHKERRQ(ierr);

    for (i=0; i<E->numFields; ++i) {
      ierr = VecAssemblyBegin(subVecs[i]);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(subVecs[i]);CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    ierr = PetscFree(subVecs);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_start_time_from_file( Parameters P , const char * filename )
{

    PetscErrorCode   ierr;
    cJSON            *json=NULL, *time;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_start_time_from_file()\n");CHKERRQ(ierr);

    ierr = read_JSON_file_to_JSON_object( filename, &json );

    /* time from this restart JSON must be passed to the time stepper
       to continue the integration and keep track of absolute time */
    time = cJSON_GetObjectItem(json,"time");
    P->t0 = time->valuedouble;

    cJSON_Delete( json );

    PetscFunctionReturn(0);

}

PetscErrorCode read_JSON_file_to_JSON_object( const char * filename, cJSON ** json )
{
    PetscErrorCode   ierr;
    FILE             *fp;
    long             length;
    char             *buffer = 0;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"read_JSON_file_to_JSON_object()\n");CHKERRQ(ierr);

    fp = fopen( filename, "r" );

    if(fp==NULL) {
       SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_FILE_OPEN,"Could not open file %s",filename);
    }

    /* read file to zero terminated string */
    /* https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c */
    if (fp){
        fseek(fp,0,SEEK_END);
        length = ftell(fp);
        fseek(fp,0,SEEK_SET);
        buffer = malloc(length+1);
        if (buffer){
            if(fread(buffer,1,length,fp) != length) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_READ,"fread() error");
        } else SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_MEM,"malloc error");
        fclose(fp);
        buffer[length] = '\0';
    }

    /* for future reference for returning pointers */
    /* https://stackoverflow.com/questions/2520777/pointer-as-second-argument-instead-of-returning-pointer */

    *json = cJSON_Parse( buffer );

    free( buffer );

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_from_file( Ctx *E, Vec sol, const char * filename, const PetscInt *arr, PetscInt size_arr )
{
    /* reads an initial condition from a previously output JSON file
       to enable restarting */

    PetscErrorCode   ierr;
    cJSON            *json=NULL, *solution, *subdomain, *values, *data, *item;
    char             *item_str;
    PetscInt         i, j, subdomain_num;
    PetscScalar      val = 0;
    Vec              invec, *subVecs;
#if (defined PETSC_USE_REAL___FLOAT128)
    char             val_str[30];
#endif

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_from_file()\n");CHKERRQ(ierr);

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    ierr = read_JSON_file_to_JSON_object( filename, &json );

    solution = cJSON_GetObjectItem(json,"solution");
    subdomain = cJSON_GetObjectItem(solution,"subdomain data");

    /* loop over all subdomains */
    cJSON_ArrayForEach( data, subdomain )
    {
        subdomain = cJSON_GetObjectItem( data, "subdomain" );
        subdomain_num = subdomain->valueint;
        values = cJSON_GetObjectItem( data, "values" );

        /* could break if ordering of subdomains changes */
        if (subdomain_num == 0){
            invec = subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDXI_B]];
        }
        else if (subdomain_num == 1){
            invec = subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]];
        }
        else if (subdomain_num == 2){
            invec = subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]];
        }
        else if (subdomain_num == 3){
            invec = subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]];
        }
        else {
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"unexpected number of subdomains");
        }

        /* determine if this subdomain is required for the IC by checking
           against the listed subdomains in arr */
        for( j=0; j<size_arr; j++ ){
            if(subdomain_num == arr[j]){
                for (i=0; i<cJSON_GetArraySize(values); i++ ){
                    item = cJSON_GetArrayItem( values, i );
                    item_str = item->valuestring;
#if (defined PETSC_USE_REAL___FLOAT128)
                    sscanf( item_str, "%s", val_str );
                    val = strtoflt128(val_str, NULL);
#else
                    sscanf( item_str, "%lf", &val );
#endif
                    /* add value to vec */
                    VecSetValue( invec, i, val, INSERT_VALUES );CHKERRQ(ierr);
                }
            }
        }
    }

    /* all vec assembly is done here */
    for (i=0; i<E->numFields; ++i) {
      ierr = VecAssemblyBegin(subVecs[i]);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(subVecs[i]);CHKERRQ(ierr);
    }
    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    ierr = PetscFree(subVecs);CHKERRQ(ierr);
    cJSON_Delete( json );

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_interior_from_phase_boundary( Ctx *E, Vec sol )
{
 
    /* entropy tracks a phase boundary (usually the solidus) below a cutoff value, and is
       equal to the cutoff value above */

    PetscErrorCode     ierr;
    Parameters const   P = E->parameters;
    Mesh const         *M = &E->mesh;
    PetscScalar        S0,S1,S2;
    const PetscScalar  *arr_pres_s;
    PetscScalar        *arr_dSdxi_b, *arr_xi_s;
    PetscInt           i, ihi_b, ilo_b, w_b;
    DM                 da_s = E->da_s, da_b=E->da_b;
    Vec                dSdxi_b, *subVecs;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior_from_phase_boundary()\n");CHKERRQ(ierr);

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    dSdxi_b = subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDXI_B]];

    /* for looping over basic nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    ierr = DMDAVecGetArray(da_b,dSdxi_b,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,M->xi_s,&arr_xi_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_s,M->pressure_s,&arr_pres_s);CHKERRQ(ierr);

    /* first and last points not updated by loop below */
    arr_dSdxi_b[0] = 0.0;
    arr_dSdxi_b[ihi_b-1] = 0.0;
    for(i=ilo_b+1; i<ihi_b-1; ++i){
        /* assumes the relevant phase boundary is in slot 0, which it will be for a single
           phase system (although still potential for a bug here if you want to initialise a
           two phase system from a phase boundary ic) */
        /* get staggered nodes values either side of the basic node */
        ierr = EOSGetPhaseBoundary( E->parameters->eos_phases[0], arr_pres_s[i], &S1, NULL);CHKERRQ(ierr);
        S1 = PetscMin( S1, P->ic_adiabat_entropy );
        ierr = EOSGetPhaseBoundary( E->parameters->eos_phases[0], arr_pres_s[i-1], &S2, NULL);CHKERRQ(ierr);
        S2 = PetscMin( S2, P->ic_adiabat_entropy );
        /* now compute gradient at basic node */
        arr_dSdxi_b[i] = S1 - S2;
        arr_dSdxi_b[i] /= arr_xi_s[i] - arr_xi_s[i-1];
    }

    /* for last point */
    /* TODO: update gradient for last point */

    /* set entropy at top staggered node */
    ierr = EOSGetPhaseBoundary( E->parameters->eos_phases[0], arr_pres_s[0], &S0, NULL);CHKERRQ(ierr);
    S0 = PetscMin( S0, P->ic_adiabat_entropy );
    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],0,S0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]]);CHKERRQ(ierr);

    /* restore basic node gradient array and sol Vec */
    ierr = DMDAVecRestoreArray(da_b,dSdxi_b,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    ierr = PetscFree(subVecs);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayRead(da_s,M->xi_s,&arr_xi_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,M->pressure_s,&arr_pres_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* initial condition of atmosphere */

static PetscErrorCode set_ic_atmosphere( Ctx *E, Vec sol )
{

    PetscErrorCode             ierr;
    Parameters const           P  = E->parameters;
    FundamentalConstants const FC = P->fundamental_constants;
    ScalingConstants const     SC = P->scaling_constants;
    AtmosphereParameters const Ap = P->atmosphere_parameters;
    Atmosphere                 *A = &E->atmosphere;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_atmosphere()\n");CHKERRQ(ierr);

    if(Ap->n_volatiles){

        /* set atmosphere IC from initial total abundance */
        if(Ap->IC_ATMOSPHERE==1){
            ierr = set_ic_atmosphere_from_initial_total_abundance( E, sol );CHKERRQ(ierr);
        }

        /* read atmosphere IC from file (could be a different file to
           the interior IC) */
        else if (Ap->IC_ATMOSPHERE==2){
            /* this only sets the sol Vec, and not the A->volatiles[i].x */
            ierr = set_ic_atmosphere_from_file( E, sol );CHKERRQ(ierr);
        }

        /* set IC from partial pressure in the atmosphere */
        else if (Ap->IC_ATMOSPHERE==3){
            ierr = set_ic_atmosphere_from_partial_pressure( E, sol );CHKERRQ(ierr);
        }

        /* set IC from number of ocean moles */
        else if (Ap->IC_ATMOSPHERE==4){
            ierr = set_ic_atmosphere_from_ocean_moles( E, sol );CHKERRQ(ierr);
        }

        else{
            SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported IC_ATMOSPHERE value %d provided",Ap->IC_ATMOSPHERE);
        }

        /* all of the above are guaranteed to set Vec sol, for the
           partial pressures and the mass reaction terms.  Now also 
           ensure that A->volatiles[i].p and A->mass_reaction[i]
           conform to the Vec sol, since these are subsequently 
           used to "correct" the mass reaction to zero */
        ierr = set_partial_pressures_from_solution( E, sol );CHKERRQ(ierr);

        /* ensure all atmosphere quantities are consistent with current
           solution */
        ierr = set_reservoir_volatile_content( A, Ap, FC, SC );CHKERRQ(ierr);

        ierr = conform_atmosphere_parameters_to_ic( E );CHKERRQ(ierr);

        ierr = set_solution_from_partial_pressures( E, sol );CHKERRQ(ierr);

    }

    /* with Ap->n_volatiles=0 there are no entries in
       the solution vector to initialise to zero */

    PetscFunctionReturn(0);
}

static PetscErrorCode conform_atmosphere_parameters_to_ic( Ctx *E )
{

    /* this function has authorisation to update the parameters struct
       to ensure that parameters obey constraints imposed by the
       initial condition */

    PetscErrorCode       ierr;
    PetscInt             i;
    PetscScalar          mass;
    Parameters           P = E->parameters;
    Atmosphere           *A = &E->atmosphere;
    AtmosphereParameters Ap = P->atmosphere_parameters;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"conform_atmosphere_parameters_to_ic()\n");CHKERRQ(ierr);

    /* prior to this function, A->volatiles[i].p and A->mass_reaction[i]
       are updated */

    for (i=0; i<Ap->n_volatiles; ++i){

       /* to this point, we have solved for A->volatiles[i].p.
          we must now ensure that other parameters are self-consistent
          with this pressure.  This might over-write some parameters
          that were previously used to determine p. */

        /* initial_total_abundance is used by some ICs to compute p,
           but with reactions mass can transfer between chemically
           reacting species.  Hence the initial prescribed total
           abundance is only honoured through conservation of the
           total number of moles of something (H/C/etc.).  In essence,
           A->volatiles[i].mass_reaction tells us how 'wrong" our initial
           prescribed total abundance was, and we can use this to instead
           compute the initial_total_abundance that obeys reactions */
        mass = A->volatiles[i].mass_liquid + A->volatiles[i].mass_solid + A->volatiles[i].mass_atmos + A->volatiles[i].mass_reaction;
        Ap->volatile_parameters[i]->initial_total_abundance = mass / (*Ap->mantle_mass_ptr);

        /* initial partial pressure */
        Ap->volatile_parameters[i]->initial_atmos_pressure = A->volatiles[i].p;

        /* TODO: conform initial_ocean_moles */
        /* Ap->volatile_parameters[i]->initial_ocean_moles */

    }

    /* since we updated the parameters above to adhere to the equilibrium
       state, we reset the mass_reaction to zero since our initial
       abundances are referenced to equilibrium at the initial
       interior and surface conditions (e.g., A->tsurf) */
    for(i=0; i<Ap->n_reactions; ++i){
        A->mass_reaction[i] = 0.0;
    }

    ierr = print_ocean_masses( E );CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode print_ocean_masses( Ctx *E )
{

    /* useful to report the equivalent present-day ocean mass
           1.4E21 kg of H2O (Olson and Sharp, 2018)
           1.55E20 kg of H2 (Olson and Sharp, 2018) */

    /* this is not general, and will break if other species contain H2,
       e.g., CH4.  Ideally, we should loop over all volatiles and identify
       those that include H2 */

    PetscErrorCode ierr;

    PetscInt i;
    /* for an arbitrary number of volatile species, would be better to somehow
       automate this output, rather than making it specific to reduced and oxidised
       phases of hydrogen and oxygen */
    PetscBool   FLAG_H2 = PETSC_FALSE;
    PetscBool   FLAG_H2O = PETSC_FALSE;
    PetscBool   FLAG_CO = PETSC_FALSE;
    PetscBool   FLAG_CO2 = PETSC_FALSE;
    PetscScalar mass_H2 = 0.0, molar_mass_H2 = 0.0, tmass_H2 = 0.0, p_H2 = 0.0;
    PetscScalar mass_H2O = 0.0, molar_mass_H2O = 0.0, tmass_H2O = 0.0, p_H2O = 0.0;
    PetscScalar mass_CO = 0.0, molar_mass_CO = 0.0, tmass_CO = 0.0, p_CO = 0.0;
    PetscScalar mass_CO2 = 0.0, molar_mass_CO2 = 0.0, tmass_CO2 = 0.0, p_CO2 = 0.0;

    Parameters P = E->parameters;
    AtmosphereParameters Ap = P->atmosphere_parameters;
    ScalingConstants const SC = P->scaling_constants;

    PetscScalar scaling = SC->VOLATILE * 4.0 * PETSC_PI * SC->MASS;
    PetscScalar scaling2 = (1.0/SC->VOLATILE) * PetscSqr(*Ap->radius_ptr) / -(*Ap->gravity_ptr);

    PetscFunctionBeginUser;

    /* this is ugly, but otherwise the flags get reset if they appear within
       the same loop over volatiles */
    for(i=0; i<Ap->n_volatiles; ++i){
        PetscStrcmp( Ap->volatile_parameters[i]->prefix, "H2", &FLAG_H2 );
        if ( FLAG_H2 ){
            mass_H2 = Ap->volatile_parameters[i]->initial_total_abundance;
            molar_mass_H2 = Ap->volatile_parameters[i]->molar_mass;
            break;
        }
    }

    for(i=0; i<Ap->n_volatiles; ++i){
        PetscStrcmp( Ap->volatile_parameters[i]->prefix, "H2O", &FLAG_H2O );
        if ( FLAG_H2O ){
            mass_H2O = Ap->volatile_parameters[i]->initial_total_abundance;
            molar_mass_H2O = Ap->volatile_parameters[i]->molar_mass;
            break;
        }
    }

    for(i=0; i<Ap->n_volatiles; ++i){
        PetscStrcmp( Ap->volatile_parameters[i]->prefix, "CO", &FLAG_CO );
        if ( FLAG_CO ){
            mass_CO = Ap->volatile_parameters[i]->initial_total_abundance;
            molar_mass_CO = Ap->volatile_parameters[i]->molar_mass;
            break;
        }
    }

    for(i=0; i<Ap->n_volatiles; ++i){
        PetscStrcmp( Ap->volatile_parameters[i]->prefix, "CO2", &FLAG_CO2 );
        if ( FLAG_CO2 ){
            mass_CO2 = Ap->volatile_parameters[i]->initial_total_abundance;
            molar_mass_CO2 = Ap->volatile_parameters[i]->molar_mass;
            break;
        }
    }

    if( FLAG_H2 || FLAG_H2O || FLAG_CO || FLAG_CO2 ){
         ierr = PetscPrintf(PETSC_COMM_WORLD,"\n**************** Volatile content **************\n");CHKERRQ(ierr);
    }

    /* the output below should report the same value, but it is a useful
       check that everything is working as expected */

    /* H2 */
    if( FLAG_H2 ){
        tmass_H2 = mass_H2;
        if ( FLAG_H2O ){
            /* equivalent mass of H2 in H2O */
            tmass_H2 += mass_H2O * (molar_mass_H2 / molar_mass_H2O );
        }
        tmass_H2 *= (*Ap->mantle_mass_ptr); /* total non-dimensional mass */
        p_H2 = tmass_H2 / scaling2; /* equivalent surface pressure */
        tmass_H2 *= scaling; /* total physical mass */
        tmass_H2 /= 1.55E20; /* non-dimensionalise according to ocean mass of H2 */
        p_H2 *= SC->PRESSURE / 1.0E5; /* to bar */
        ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %-15.6g\n","Equivalent present-day mass of ocean water from H2 (non-dimensional)",(double)tmass_H2);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %-15.6g\n","Equivalent atmospheric pressure of H2 (bar)",(double)p_H2);CHKERRQ(ierr);

    }

    /* H2O */
    if( FLAG_H2O ){
        tmass_H2O = mass_H2O;
        if ( FLAG_H2 ){
            /* equivalent mass of H2O in H2 */
            tmass_H2O += mass_H2 * (molar_mass_H2O / molar_mass_H2 );
        }
        tmass_H2O *= (*Ap->mantle_mass_ptr); /* total non-dimensional mass */
        p_H2O = tmass_H2O / scaling2; /* equivalent surface pressure */
        tmass_H2O *= scaling; /* total physical mass */
        //tmass_H2O /= 1.4E21; /* non-dimensionalise according to ocean mass of H2O */
        /* below I have taken the Olson and Sharp value for H2, and then scaled by the typical
           molar mass for H2 and H2O, such that by construction, the ocean mass estimates from
           both H2 and H2O are identical */
        tmass_H2O /= 1.385185824553049e+21;
        p_H2O *= SC->PRESSURE / 1.0E5; /* to bar */ 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %-15.6g\n","Equivalent present-day mass of ocean water from H2O (non-dimensional)",(double)tmass_H2O);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %-15.6g\n","Equivalent atmospheric pressure of H2O (bar)",(double)p_H2O);CHKERRQ(ierr);
    }

    /* CO */
    if( FLAG_CO ){
        tmass_CO = mass_CO;
        if ( FLAG_CO2 ){
            /* equivalent mass of CO in CO2 */
            tmass_CO += mass_CO2 * (molar_mass_CO / molar_mass_CO2 );
        }
        tmass_CO *= (*Ap->mantle_mass_ptr); /* total non-dimensional mass */
        p_CO = tmass_CO / scaling2; /* equivalent surface pressure */
        tmass_CO *= scaling; /* total physical mass */
        tmass_CO /= 2.153674821914003e+21; /* non-dimensionalise according to ocean mass of CO */
        p_CO *= SC->PRESSURE / 1.0E5; /* to bar */
        ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %-15.6g\n","Equivalent present-day mass of ocean water from CO (non-dimensional)",(double)tmass_CO);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %-15.6g\n","Equivalent atmospheric pressure of CO (bar)",(double)p_CO);CHKERRQ(ierr);

    }

    /* CO2 */
    if( FLAG_CO2 ){
        tmass_CO2 = mass_CO2;
        if ( FLAG_CO ){
            /* equivalent mass of H2O in H2 */
            tmass_CO2 += mass_CO * (molar_mass_CO2 / molar_mass_CO );
        }
        tmass_CO2 *= (*Ap->mantle_mass_ptr); /* total non-dimensional mass */
        p_CO2 = tmass_CO2 / scaling2; /* equivalent surface pressure */
        tmass_CO2 *= scaling; /* total physical mass */
        //tmass_H2O /= 1.4E21; /* non-dimensionalise according to ocean mass of H2O */
        /* below I have taken the Olson and Sharp value for H2, and then scaled by the typical
           molar mass for H2 and H2O, such that by construction, the ocean mass estimates from
           both H2 and H2O are identical */
        tmass_CO2 /= 3.383906780165486e+21;
        p_CO2 *= SC->PRESSURE / 1.0E5; /* to bar */ 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %-15.6g\n","Equivalent present-day mass of ocean water from CO2 (non-dimensional)",(double)tmass_CO2);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"%-30s %-15.6g\n","Equivalent atmospheric pressure of CO2 (bar)",(double)p_CO2);CHKERRQ(ierr);
    }

    if( FLAG_H2 || FLAG_H2O || FLAG_CO || FLAG_CO2 ){
         ierr = PetscPrintf(PETSC_COMM_WORLD,"************************************************\n\n");CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_atmosphere_from_ocean_moles( Ctx *E, Vec sol )
{
    PetscErrorCode       ierr;
    PetscInt             i;
    PetscScalar          abund;
    AtmosphereParameters Ap = E->parameters->atmosphere_parameters;
    Parameters const     P = E->parameters;
    FundamentalConstants const FC = P->fundamental_constants;
    ScalingConstants     const SC = P->scaling_constants;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_atmosphere_from_ocean_moles()\n");CHKERRQ(ierr);

    /* convert from ocean moles to mass and update total abundance */
    for (i=0; i<Ap->n_volatiles; ++i) {
        /* mass of volatile */
        abund = Ap->volatile_parameters[i]->initial_ocean_moles * FC->OCEAN_MOLES;
        abund *= Ap->volatile_parameters[i]->molar_mass;
        abund /= (*Ap->mantle_mass_ptr) * SC->VOLATILE * 4.0 * PETSC_PI;
        Ap->volatile_parameters[i]->initial_total_abundance = abund;
    }

    ierr = set_ic_atmosphere_from_initial_total_abundance( E, sol );CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_atmosphere_from_initial_total_abundance( Ctx *E, Vec sol )
{
    PetscErrorCode       ierr;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_atmosphere_from_initial_total_abundance()\n");CHKERRQ(ierr);

    ierr = solve_for_initial_partial_pressure( E );CHKERRQ(ierr);

    /* below will also update mass reaction terms, which can be
       non-zero */
    ierr = set_solution_from_partial_pressures( E, sol );CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_atmosphere_from_partial_pressure( Ctx *E, Vec sol )
{
    PetscErrorCode       ierr;
    PetscInt             i;
    Atmosphere           *A = &E->atmosphere;
    AtmosphereParameters Ap = E->parameters->atmosphere_parameters;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_atmosphere_from_partial_pressure()\n");CHKERRQ(ierr);

    /* set initial partial pressure to A->volatiles[i].p */
    for (i=0; i<Ap->n_volatiles; ++i) {
        A->volatiles[i].p = Ap->volatile_parameters[i]->initial_atmos_pressure;
    }

    /* mass reaction has not been updated, so should still be zero */
    /* this initial condition is not compatible with reactions,
       but nor does it need to be, since this initial condition is
       used when codes are coupled (e.g. to VULCAN) */
    ierr = set_solution_from_partial_pressures( E, sol );CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode solve_for_initial_partial_pressure( Ctx *E )
{

    /* returns both volatile partial pressures and reaction
       masses, both of which can be non-zero */

    PetscErrorCode ierr;
    SNES           snes;
    Vec            x,r;
    PetscScalar    *xx;
    PetscInt       i;
    Atmosphere                 *A = &E->atmosphere;
    Parameters           const P = E->parameters;
    AtmosphereParameters const Ap = P->atmosphere_parameters;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"solve_for_initial_partial_pressure()\n");CHKERRQ(ierr);

    ierr = SNESCreate( PETSC_COMM_WORLD, &snes );CHKERRQ(ierr);

    /* Use this to address this specific SNES (nonlinear solver) from the command
       line or options file, e.g. -atmosic_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"atmosic_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    /* Vector with one dof per volatile species, and one dof per volatile reaction */
    ierr = VecSetSizes( x, PETSC_DECIDE, Ap->n_volatiles + Ap->n_reactions );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    ierr = SNESSetFunction(snes,r,objective_function_initial_partial_pressure,E);CHKERRQ(ierr);

    /* initialise vector x with initial guess */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    for (i=0; i<Ap->n_volatiles; ++i) {
        xx[i] = Ap->volatile_parameters[i]->initial_atmos_pressure;
    }
    /* initial guesses for reaction masses */
    for (i=Ap->n_volatiles; i<Ap->n_volatiles + Ap->n_reactions; ++i) {
      xx[i] = 0.0; /* assume we close to equilibrium */
    }
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* it's often tricky for the user to know how the mass is partitioned between
       species, and the solver fails if the initial guess of the mass distribution
       is way off.  The newtontr seems to help, and particularly increasing the size
       of delta0 to 10.0.  But I'm sure further optimisations are possible */
    ierr = PetscOptionsSetValue(NULL,"-atmosic_snes_type","newtontr");CHKERRQ(ierr);
    /* Inform the nonlinear solver to generate a finite-difference approximation
       to the Jacobian */
    ierr = PetscOptionsSetValue(NULL,"-atmosic_snes_mf",NULL);CHKERRQ(ierr);
    /* Turn off convergence based on step size */
    ierr = PetscOptionsSetValue(NULL,"-atmosic_snes_stol","0");CHKERRQ(ierr);
    /* Turn off convergenced based on trust region tolerance */
    ierr = PetscOptionsSetValue(NULL,"-atmosic_snes_trtol","0");CHKERRQ(ierr);
    /* https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/SNES/SNESNEWTONTR.html */
    /* initial size of the trust region is delta0*norm2(x) */
    ierr = PetscOptionsSetValue(NULL,"-atmosic_snes_tr_delta0","10.0");CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-atmosic_snes_rtol","1.0e-6");CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-atmosic_snes_atol","1.0e-6");CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-atmosic_ksp_rtol","1.0e-6");CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-atmosic_ksp_rtol","1.0e-6");CHKERRQ(ierr);

    /* below previously suggested by PS, but not currently used */
    //ierr = PetscOptionsSetValue(NULL,"-atmosic_snes_linesearch_damping","0.01");CHKERRQ(ierr);
    //ierr = PetscOptionsSetValue(NULL,"-atmosic_snes_max_it","10000");CHKERRQ(ierr);

    /* For solver analysis/debugging/tuning, activate a custom monitor with a flag */
    {
      PetscBool flg = PETSC_FALSE;

      ierr = PetscOptionsGetBool(NULL,NULL,"-atmosic_snes_verbose_monitor",&flg,NULL);CHKERRQ(ierr);
      if (flg) {
        ierr = SNESMonitorSet(snes,SNESMonitorVerbose,NULL,NULL);CHKERRQ(ierr);
      }
    }

    /* Solve */
    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr); /* Picks up any additional options (note prefix) */
    ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);
    {
      SNESConvergedReason reason;
      ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
      if (reason < 0) SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
          "Nonlinear solver didn't converge: %s\n",SNESConvergedReasons[reason]);
    }

    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    for (i=0; i<Ap->n_volatiles; ++i) {
        if( A->volatiles[i].p < 0.0 ){
            /* Sanity check on solution (since it's non-unique) */
            SETERRQ2(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
                "Unphysical initial volatile partial pressure: volatile %d, x: %g",i,A->volatiles[i].p);
        }
        else{
            A->volatiles[i].p = xx[i];
        }
    }

    for (i=0; i<Ap->n_reactions; ++i) {
      A->mass_reaction[i] = xx[Ap->n_volatiles+i];
    }

    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode objective_function_initial_partial_pressure( SNES snes, Vec x, Vec f, void *ptr)
{
    PetscErrorCode             ierr;
    const PetscScalar          *xx;
    PetscScalar                *ff;
    PetscInt                   i;
    Ctx                        *E = (Ctx*) ptr;
    Atmosphere                 *A = &E->atmosphere;
    Parameters           const P = E->parameters;
    FundamentalConstants const FC = P->fundamental_constants;
    ScalingConstants     const SC = P->scaling_constants;
    AtmosphereParameters const Ap = P->atmosphere_parameters;

    PetscFunctionBeginUser;

    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

    for (i=0; i<Ap->n_volatiles; ++i) {
        A->volatiles[i].p = xx[i];
    }
    for (i=0; i<Ap->n_reactions; ++i) {
        A->mass_reaction[i] = xx[Ap->n_volatiles+i];
    }

    ierr = set_reservoir_volatile_content( A, Ap, FC, SC ); CHKERRQ(ierr);

    /* mass conservation for each volatile */
    for (i=0; i<Ap->n_volatiles; ++i) {
        ff[i] = get_residual_volatile_mass( A, Ap, Ap->volatile_parameters[i], &A->volatiles[i]);
    }

    for (i=0; i<Ap->n_reactions; ++i) {

        PetscScalar Qp,Qr,log10G,G;

        /* fO2 is not accounted for here */
        Qp = get_reaction_quotient_products( Ap->reaction_parameters[i], A );
        Qr = get_reaction_quotient_reactants( Ap->reaction_parameters[i], A );
        /* (Modified) equilibrium constant that accommodates fO2 */
        log10G = get_log10_modified_equilibrium_constant( Ap->reaction_parameters[i], A->tsurf, A );
        G = PetscPowScalar( 10.0, log10G );

        ff[Ap->n_volatiles + i] = Qp/Qr - G;

    }

    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
