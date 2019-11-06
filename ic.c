#include "atmosphere.h"
#include "energy.h"
#include "ic.h"
#include "monitor.h"
#include "parameters.h"
#include "reaction.h"
#include "util.h"

/* interior ic */
static PetscErrorCode set_ic_interior( Ctx *, Vec );
static PetscErrorCode set_ic_interior_default( Ctx *, Vec );
static PetscErrorCode set_ic_interior_entropy( Ctx *, Vec );
static PetscErrorCode set_ic_interior_from_file( Ctx *, Vec );
static PetscErrorCode set_ic_interior_from_solidus( Ctx *, Vec );
static PetscErrorCode set_ic_interior_conform_to_bcs( Ctx *, Vec );
/* atmosphere ic */
static PetscErrorCode set_ic_atmosphere( Ctx *, Vec );
static PetscErrorCode set_ic_atmosphere_default( Ctx *, Vec );
static PetscErrorCode set_ic_atmosphere_from_initial_total_abundance( Ctx *E, Vec sol );
static PetscErrorCode set_ic_atmosphere_from_file( Ctx *, Vec );
static PetscErrorCode set_ic_atmosphere_from_partial_pressure( Ctx *, Vec );
static PetscErrorCode solve_for_initial_melt_abundance( Ctx * );
/* general */
static PetscErrorCode set_ic_from_file( Ctx *, Vec, const char *, const PetscInt *, PetscInt );
static PetscErrorCode objective_function_initial_melt_abundance( SNES, Vec, Vec, void *);
static PetscErrorCode conform_parameters_to_initial_condition( Ctx * );

/* main function to set initial condition */
PetscErrorCode set_initial_condition( Ctx *E, Vec sol)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_initial_condition()\n");CHKERRQ(ierr);

    /* interior initial condition */
    ierr = set_ic_interior( E, sol ); CHKERRQ(ierr);

    /* atmosphere initial condition */
    /* must be after interior initial condition */
    ierr = set_ic_atmosphere( E, sol ); CHKERRQ(ierr);

    /* conform parameters structs to initial condition */
    ierr = conform_parameters_to_initial_condition( E );CHKERRQ(ierr);

    /* this also updates mass reaction terms */
    ierr = set_solution_from_volatile_abundances( E, sol ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* initial condition of interior */

static PetscErrorCode set_ic_interior( Ctx *E, Vec sol)
{
    PetscErrorCode ierr;
    Parameters const *P = &E->parameters;
    PetscInt IC = P->IC_INTERIOR;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior()\n");CHKERRQ(ierr);

    if (IC==1){
        ierr = set_ic_interior_default( E, sol ); CHKERRQ(ierr);
    }
    else if (IC==2){
        ierr = set_ic_interior_from_file( E, sol ); CHKERRQ(ierr);
        /* in parameters.c, this condition also sets the start time
           P->t0 from the interior file */
    }
    else if (IC==3){
        ierr = set_ic_interior_from_solidus( E, sol ); CHKERRQ(ierr);
    }
    else{
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported IC_INTERIOR value %d provided",P->IC_INTERIOR);
    }

    ierr = set_ic_interior_conform_to_bcs( E, sol ); CHKERRQ(ierr);

    /* P->t0 is set in parameters.c */
    /* time is needed for the radiogenic energy input */
    /* below sets interior structure, since possibly required for
       atmosphere initial condition */
    ierr = set_interior_structure_from_solution( E, P->t0, sol ); CHKERRQ(ierr);

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
    Parameters const *P = &E->parameters;
    /* subdomains to use, i.e. dS/dr and S0 */
    PetscInt const arr[2] = {0, 1};

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior_from_file()\n");CHKERRQ(ierr);

    ierr = set_ic_from_file( E, sol, P->ic_interior_filename, arr, 2 ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_atmosphere_from_file( Ctx *E, Vec sol )
{
    PetscErrorCode ierr;
    Parameters const *P = &E->parameters;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    /* subdomains to use, i.e. volatile abundances */
    PetscInt const arr[2] = {2, 3};

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_atmosphere_from_file()\n");CHKERRQ(ierr);

    ierr = set_ic_from_file( E, sol, Ap->ic_atmosphere_filename, arr, 2 ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_interior_entropy( Ctx *E, Vec sol )
{
    /* set initial entropy gradient to constant for all nodes */

    PetscErrorCode   ierr;
    PetscInt         i;
    PetscScalar      S0;
    Parameters const *P = &E->parameters;
    Vec              dSdr_b;
    Vec              *subVecs;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior_entropy()\n");CHKERRQ(ierr);

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    dSdr_b = subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDR_B]];

    ierr = VecSet(dSdr_b, P->ic_dsdr );CHKERRQ(ierr);

    /* these next two lines are simply convenient reminders that the
       first and last values are meaningless because the fluxes here
       are controlled by boundary conditions.  But for debugging and
       clarity it is convenient to explicitly set these values to
       zero */
    ierr = VecSetValue( dSdr_b, 0,             0.0, INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue( dSdr_b, P->numpts_b-1, 0.0, INSERT_VALUES);CHKERRQ(ierr);

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
    /* TODO: can this be improved? */
    // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
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

        /* FIXME: could break if ordering of subdomains changes */
        if (subdomain_num == 0){
            invec = subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDR_B]];
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

static PetscErrorCode set_ic_interior_from_solidus( Ctx *E, Vec sol )
{

    /* entropy tracks the solidus below a cutoff value, and is
       equal to the cutoff value above.  It is simplest to construct
       the ic using S->S_s, and then map the values to the
       solution Vec */

    PetscErrorCode   ierr;
    Parameters const *P = &E->parameters;
    PetscInt         i, numpts_s;
    PetscScalar      S_i;
    Solution         *S = &E->solution;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior_from_solidus()\n");CHKERRQ(ierr);

    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    ierr = VecCopy( S->solidus_s, S->S_s ); CHKERRQ(ierr);

    for(i=1; i<numpts_s-1;++i) {
        ierr = VecGetValues(S->S_s,1,&i,&S_i);CHKERRQ(ierr);
        if(P->ic_adiabat_entropy < S_i){
            ierr = VecSetValues(S->S_s,1,&i,&P->ic_adiabat_entropy,INSERT_VALUES);CHKERRQ(ierr);
        }
    }

    VecAssemblyBegin( S->S_s );
    VecAssemblyEnd( S->S_s );

    /* map S->S_s to the solution Vec (derivatives are taken) */
    ierr = set_solution_from_entropy_at_staggered_nodes( E, sol ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_interior_conform_to_bcs( Ctx *E, Vec sol )
{
    /* ensure that the interior ic is compatible with boundary
       conditions, and that the entropy vecs in the solution
       struct are consistent with the sol Vec (and vice-versa) */

    PetscErrorCode ierr;
    PetscInt ind,numpts_s;
    Solution *S = &E->solution;
    Parameters const *P = &E->parameters;
    PetscInt const ind0=0;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior_conform_to_bcs()\n");CHKERRQ(ierr);

    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ind = numpts_s-1;

    /* conform entropy-related Vecs in solution struct to sol Vec */
    ierr = set_entropy_from_solution( E, sol ); CHKERRQ(ierr);

    /* ensure that initial condition conforms to boundary conditions */
    /* TODO: these are just for isothermal, and also first order, since
       we set the top and bottom staggered node to the surface and
       CMB entropy, respectively */

    if( (P->ic_surface_entropy > 0.0) || (P->ic_core_entropy > 0.0) ){
        if( P->ic_surface_entropy > 0.0 ){
            ierr = VecSetValues( S->S_s,1,&ind0,&P->ic_surface_entropy,INSERT_VALUES );CHKERRQ(ierr);
        }
        if( P->ic_core_entropy > 0.0 ){
            ierr = VecSetValues( S->S_s,1,&ind,&P->ic_core_entropy,INSERT_VALUES );CHKERRQ(ierr);
        }
        ierr = VecAssemblyBegin(S->S_s);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(S->S_s);CHKERRQ(ierr);
    }

    /* conform sol Vec to (potentially) modified S->S_s */
    ierr = set_solution_from_entropy_at_staggered_nodes( E, sol );CHKERRQ(ierr);

    /* conform other entropy-related Vecs in sol struct, since thus
       far only S->S_s is guaranteed to be consistent with sol Vec */
    ierr = set_entropy_from_solution( E, sol );CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* initial condition of atmosphere */

static PetscErrorCode set_ic_atmosphere( Ctx *E, Vec sol )
{

    PetscErrorCode             ierr;
    Parameters const           *P  = &E->parameters;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    Atmosphere                 *A = &E->atmosphere;
    Constants const            *C = &P->constants;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_atmosphere()\n");CHKERRQ(ierr);

    if(Ap->n_volatiles){

        if(Ap->IC_ATMOSPHERE==1){
            ierr = set_ic_atmosphere_default( E, sol );CHKERRQ(ierr);
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

        else{
            SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported IC_ATMOSPHERE value %d provided",Ap->IC_ATMOSPHERE);
        }

        /* all of the above are guaranteed to set Vec sol, for the
           volatile abundances in the melt and the mass reaction
           terms.  Now also ensure that A->volatiles[i].x and
           A->mass_reaction[i] conform to the Vec sol, since these
           are subsequently used to "correct" the mass reaction
           to zero */
        ierr = set_volatile_abundances_from_solution( E, sol );CHKERRQ(ierr);

        /* ensure all atmosphere quantities are consistent with current
           solution */
        ierr = set_atmosphere_volatile_content( A, Ap, C );CHKERRQ(ierr);

        /* again, note that mass reaction terms are not yet
           (necessarily) zero! */

    }

    /* with Ap->n_volatiles=0 there are no entries in
       the solution vector to initialise to zero */

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_atmosphere_default( Ctx *E, Vec sol )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_atmosphere_default()\n");CHKERRQ(ierr);

    ierr = set_ic_atmosphere_from_initial_total_abundance( E, sol ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode conform_parameters_to_initial_condition( Ctx *E )
{

    /* this function has authorisation to update the parameters struct
       to ensure that parameters obey constraints imposed by the
       initial condition */

    PetscErrorCode       ierr;
    PetscInt             i;
    PetscScalar          mass;
    Parameters           *P = &E->parameters;
    Atmosphere           *A = &E->atmosphere;
    AtmosphereParameters *Ap = &P->atmosphere_parameters;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"conform_parameters_to_initial_condition()\n");CHKERRQ(ierr);

    /* prior to this function, A->volatiles[i].x and A->mass_reaction[i]
       are updated */

    /* conform initial_total_abundance */
    for (i=0; i<Ap->n_volatiles; ++i){
        if( Ap->IC_ATMOSPHERE==1 ){
            /* this is equivalent to the operation below for Ap->IC_ATMOSPHERE==3, but a shortcut since
               we do not need to sum all reservoirs to get to the initial total abundance */
            /* below we correct with -= */
            Ap->volatile_parameters[i].initial_total_abundance -= A->volatiles[i].mass_reaction / (*Ap->mantle_mass_ptr);
        }

        /* do not conform if Ap->IC_ATMOSPHERE==2, since we assume the user wants to resume from the exact state as
           defined in the restart file */

        else if( Ap->IC_ATMOSPHERE==3 ){
            mass = A->volatiles[i].mass_liquid + A->volatiles[i].mass_solid + A->volatiles[i].mass_atmos + A->volatiles[i].mass_reaction;
            /* below we set with = */
            Ap->volatile_parameters[i].initial_total_abundance = mass / (*Ap->mantle_mass_ptr);
        }
    }

    /* re-solve to get volatile abundances in the melt (A->volatiles[i].x) */
    ierr = solve_for_initial_melt_abundance( E ); CHKERRQ(ierr);

    /* Ap_>IC_ATMOSPHERE==1 will also set mass reactions close to (but not exactly)
       to zero.  Explicitly zero the entries here */
    /* FIXME: add tolerance check to ensure that mass reactions are
       actually close to zero before overwriting them here */
    for(i=0; i<Ap->n_reactions; ++i){
        A->mass_reaction[i] = 0.0;
    }

    for(i=0; i<Ap->n_volatiles; ++i){
        Ap->volatile_parameters[i].initial_atmos_pressure = A->volatiles[i].p;
    }

    /* the relevant updated values in the structs are mapped back to the sol
       Vec in the next function call following this return */

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_atmosphere_from_initial_total_abundance( Ctx *E, Vec sol )
{
    PetscErrorCode       ierr;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_atmosphere_from_initial_total_abundance()\n");CHKERRQ(ierr);

    ierr = solve_for_initial_melt_abundance( E );CHKERRQ(ierr);

    /* below will also update mass reaction terms, which can be
       non-zero.  We "correct" the mass reaction to be zero in
       conform_parameters_to_initial_condition() */
    ierr = set_solution_from_volatile_abundances( E, sol );CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_atmosphere_from_partial_pressure( Ctx *E, Vec sol )
{
    PetscErrorCode       ierr;
    PetscInt             i;
    Atmosphere           *A = &E->atmosphere;
    AtmosphereParameters *Ap = &E->parameters.atmosphere_parameters;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_atmosphere_from_partial_pressure()\n");CHKERRQ(ierr);

    /* set initial partial pressure to A->volatiles[i].p */
    for (i=0; i<Ap->n_volatiles; ++i) {
        A->volatiles[i].p = Ap->volatile_parameters[i].initial_atmos_pressure;
    }

    /* compute initial volatile abundance in melt */
    ierr = set_volatile_abundances_from_partial_pressure( A, Ap );CHKERRQ(ierr);

    /* mass reaction has not been updated, so should still be zero */
    /* TODO: this initial condition is not compatible with reactions,
       but nor does it need to be, since this initial condition is
       used when codes are coupled (e.g. to VULCAN) */
    ierr = set_solution_from_volatile_abundances( E, sol );CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode solve_for_initial_melt_abundance( Ctx *E )
{

    /* returns both volatile abundances in the melt and reaction
       masses, both of which can be non-zero */

    PetscErrorCode ierr;
    SNES           snes;
    Vec            x,r;
    PetscScalar    *xx;
    PetscInt       i;
    Atmosphere                 *A = &E->atmosphere;
    Parameters           const *P = &E->parameters;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"solve_for_initial_melt_abundance()\n");CHKERRQ(ierr);

    ierr = SNESCreate( PETSC_COMM_WORLD, &snes );CHKERRQ(ierr);

    /* Use this to address this specific SNES (nonlinear solver) from the command
       line or options file, e.g. -atmosic_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"atmosic_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    /* Vector with one dof per volatile species, and one dof per volatile reaction */
    ierr = VecSetSizes( x, PETSC_DECIDE, Ap->n_volatiles + Ap->n_reactions );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    ierr = SNESSetFunction(snes,r,objective_function_initial_melt_abundance,E);CHKERRQ(ierr);

    /* initialise vector x with initial guess */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    for (i=0; i<Ap->n_volatiles; ++i) {
        xx[i] = Ap->volatile_parameters[i].initial_total_abundance;
    }
    /* Initial guesses for reaction masses */
    for (i=Ap->n_volatiles; i<Ap->n_volatiles + Ap->n_reactions; ++i) {
      xx[i] = 0.0; // assume we are at equilibrium
    }
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* Inform the nonlinear solver to generate a finite-difference approximation
       to the Jacobian */
    ierr = PetscOptionsSetValue(NULL,"-atmosic_snes_mf",NULL);CHKERRQ(ierr);

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
        if( A->volatiles[i].x < 0.0 ){
            /* Sanity check on solution (since it's non-unique) */
            SETERRQ2(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
                "Unphysical initial total volatile abundance: volatile %d, x: %g",i,A->volatiles[i].x);
        }
        else{
            A->volatiles[i].x = xx[i];
        }
    }
    /* Save mass offset to reset initial volatile */
    for (i=0; i<Ap->n_reactions; ++i) {
      A->mass_reaction[i] = xx[Ap->n_volatiles + i ];
    }

    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode objective_function_initial_melt_abundance( SNES snes, Vec x, Vec f, void *ptr)
{
    PetscErrorCode             ierr;
    const PetscScalar          *xx;
    PetscScalar                *ff;
    PetscInt                   i;
    Ctx                        *E = (Ctx*) ptr;
    Atmosphere                 *A = &E->atmosphere;
    Parameters           const *P = &E->parameters;
    Constants            const *C = &P->constants;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;

    PetscFunctionBeginUser;

    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

    for (i=0; i<Ap->n_volatiles; ++i) {
        A->volatiles[i].x = xx[i];
    }
    for (i=0; i<Ap->n_reactions; ++i) {
        A->mass_reaction[i] = xx[Ap->n_volatiles+i];
    }

    ierr = set_atmosphere_volatile_content( A, Ap, C ); CHKERRQ(ierr);

    /* mass conservation for each volatile */
    for (i=0; i<Ap->n_volatiles; ++i) {
        ff[i] = get_residual_volatile_mass( A, Ap, &Ap->volatile_parameters[i], &A->volatiles[i]);
    }

    for (i=0; i<Ap->n_reactions; ++i) {

        PetscScalar Qp,Qr,K;

        ReactionParameters const * reaction_parameters_ptr = &Ap->reaction_parameters[i];

        /* return the numerator and denominator separately to retain scalings,
           otherwise the non-linear solver has problems */
        /* since scaling by A->psurf, it may not be necessary any longer
           to split the reaction quotient into two parts, but nevertheless
           why fix something that is not broken? */
        Qp = get_reaction_quotient_products( reaction_parameters_ptr, A );
        Qr = get_reaction_quotient_reactants( reaction_parameters_ptr, A );

        K = get_equilibrium_constant( reaction_parameters_ptr, A->tsurf, C );

        /* the form below was used for testing the simplewater reactions, but seems
           to fail when fO2 is included (which changes the scaling by around 2-3
           orders of magnitude) */
        //ff[Ap->n_volatiles + i] = Qp - K * Qr;
        /* the latest form, below, seems to work for all oxygen fugacity models
           as well as the simplewater tests.  Perhaps a scaling is retained? */
        ff[Ap->n_volatiles + i] = Qp/Qr - K;

    }

#if 0
    /* Objective function, "simple" reactions (2 species, constant epsilon) only */
    for (i=0; i<Ap->n_reactions; ++i) {
      PetscBool is_simple;

      ierr = PetscStrcmp(Ap->reaction_parameters[i]->type,"simple",&is_simple);CHKERRQ(ierr);
      if (is_simple) {
        const PetscInt v0 = Ap->reaction_parameters[i]->volatiles[0];
        const PetscInt v1 = Ap->reaction_parameters[i]->volatiles[1];

        ff[Ap->n_volatiles + i] = Ap->reaction_parameters[i]->epsilon[v0] * A->volatiles[v0].p + Ap->reaction_parameters[i]->epsilon[v1] * A->volatiles[v1].p;
      } else{
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Reaction type %s not recognized",Ap->reaction_parameters[i]->type);
      }
    }
#endif

    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
