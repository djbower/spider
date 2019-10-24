#include "atmosphere.h"
#include "bc.h"
#include "cJSON.h"
#include "energy.h"
#include "ic.h"
#include "matprop.h"
#include "monitor.h"
#include "parameters.h"
#include "reaction.h"
#include "twophase.h"
#include "util.h"

static PetscErrorCode set_ic_default( Ctx *, Vec );
static PetscErrorCode set_ic_entropy( Ctx *, Vec );
static PetscErrorCode set_ic_atmosphere( Ctx *, Vec );
static PetscErrorCode set_ic_interior( Ctx *, Vec );
static PetscErrorCode set_ic_from_file( Ctx *, Vec );
static PetscErrorCode set_ic_from_solidus( Ctx *, Vec );
static PetscErrorCode set_ic_interior_conform_to_bcs( Ctx *, Vec );
static PetscErrorCode set_initial_volatile( Ctx * );
static PetscErrorCode FormFunction1( SNES, Vec, Vec, void *); // TODO DJB: rename this function

PetscErrorCode set_initial_condition( Ctx *E, Vec sol)
{
    PetscErrorCode ierr;
    Parameters const *P = &E->parameters;
    PetscInt IC = P->initial_condition;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_initial_condition()\n");CHKERRQ(ierr);

    ierr = set_ic_interior( E, sol ); CHKERRQ(ierr);

    if( (IC==1 || IC==3) ){
        ierr = set_ic_atmosphere( E, sol ); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);

}

/* initial condition of interior */

static PetscErrorCode set_ic_interior( Ctx *E, Vec sol)
{
    PetscErrorCode ierr;
    Parameters const *P = &E->parameters;
    PetscInt IC = P->initial_condition;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_interior()\n");CHKERRQ(ierr);

    if (IC==1){
        ierr = set_ic_default( E, sol ); CHKERRQ(ierr);
    }
    else if (IC==2){
        /* set everything, including the atmosphere, from a JSON */
        /* also sets P->t0 (start time) */
        ierr = set_ic_from_file( E, sol ); CHKERRQ(ierr);
    }
    else if (IC==3){
        ierr = set_ic_from_solidus( E, sol ); CHKERRQ(ierr);
    }

    ierr = set_ic_interior_conform_to_bcs( E, sol ); CHKERRQ(ierr);

    /* P->t0 (start time) is by default zero, unless set by the user.
       The value is set from the restart JSON for IC==2.  Time
       is needed for the radiogenic energy input */
    ierr = set_interior_structure_from_solution( E, P->t0, sol ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_default( Ctx *E, Vec sol )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_default()\n");CHKERRQ(ierr);

    ierr = set_ic_entropy( E, sol ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_ic_entropy( Ctx *E, Vec sol )
{
    /* set initial entropy gradient to constant for all nodes */

    PetscErrorCode   ierr;
    PetscInt         i;
    PetscScalar      S0;
    Parameters const *P = &E->parameters;
    Vec              dSdr_b;
    Vec              *subVecs;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_entropy()\n");CHKERRQ(ierr);

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

static PetscErrorCode set_ic_from_file( Ctx *E, Vec sol )
{
    /* reads an initial condition from a previously output JSON file
       to enable restarting */

    PetscErrorCode   ierr;
    Parameters *P  = &E->parameters;
    FILE             *fp;
    cJSON            *json, *solution, *subdomain, *values, *data, *item, *time;
    long             length;
    char             *item_str;
    PetscInt         i, subdomain_num;
    PetscScalar      val = 0;
    Vec              invec, *subVecs;
    char             *buffer = 0;
#if (defined PETSC_USE_REAL___FLOAT128)
    char             val_str[30];
#endif

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_from_file()\n");CHKERRQ(ierr);

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    fp = fopen( P->ic_filename, "r" );

    if(fp==NULL) {
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_FILE_OPEN,"Could not open file %s",P->ic_filename);
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

    json = cJSON_Parse( buffer );

    /* time from this restart JSON must be passed to the time stepper
       to continue the integration and keep track of absolute time */
    time = cJSON_GetObjectItem(json,"time");
    P->t0 = time->valuedouble;

    solution = cJSON_GetObjectItem(json,"solution");
    subdomain = cJSON_GetObjectItem(solution,"subdomain data");

    /* loop over subdomains and extract values */
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
        else {
            SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"unexpected number of subdomains");
        }

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

    /* all vec assembly is done here */
    for (i=0; i<E->numFields; ++i) {
      ierr = VecAssemblyBegin(subVecs[i]);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(subVecs[i]);CHKERRQ(ierr);
    }
    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    ierr = PetscFree(subVecs);CHKERRQ(ierr);
    cJSON_Delete( json );
    free( buffer );

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_from_solidus( Ctx *E, Vec sol )
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
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_from_solidus()\n");CHKERRQ(ierr);

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
    ierr = PetscPrintf(PETSC_COMM_WORLD,"ic_conform_to_bcs()\n");CHKERRQ(ierr);

    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ind = numpts_s-1;

    /* conform entropy-related Vecs in solution struct to sol Vec */
    ierr = set_entropy_from_solution( E, sol );

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
    PetscInt                   i,j;
    PetscScalar                x0, factor, massv;
    Vec                        *subVecs;
    Atmosphere                 *A = &E->atmosphere;
    Parameters                 *P  = &E->parameters;
    AtmosphereParameters       *Ap = &P->atmosphere_parameters;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    // FIXME this is horrible - parameters should not ever be changed!! Rather, there should be a dedicated step (with output to the user) which fixes inconsistent ICs
    if(Ap->SOLVE_FOR_VOLATILES){
        ierr = set_initial_volatile( E ); CHKERRQ(ierr);
        /* Resolve for the initial volatile abundances that are necessary to
           satisfy chemical equilibrium.  This operation effectively sets
           the initial chemical reaction mass to zero */
        /* TODO: some overlap here with FormFunction1.  Could have a single
           function that returns the scalings (massv) rather than recomputing
           them */
        for (i=0; i<Ap->n_reactions; ++i){
          const PetscInt v0 = Ap->reaction_parameters[i]->volatiles[0];
          /* by convention, first volatile is a reactant, so stoichiometry (hence factor) will be -ve */
          /* NOTE: introduce scaling by A->psurf to improve scaling for numerical solver (FD Jacobian) */
          /* TODO: swap out A->volatiles[v0].p/A->psurf for the volume mixing ratio? */
          factor = Ap->reaction_parameters[i]->stoichiometry[0] * Ap->volatile_parameters[v0].molar_mass * (A->volatiles[v0].p/A->psurf);
          for (j=0; j<Ap->reaction_parameters[i]->n_volatiles; ++j) {
            const PetscInt v = Ap->reaction_parameters[i]->volatiles[j];
            /* Note: maybe it's possible to just get mass_reaction out of sol, avoiding the state in A */
            /* NOTE: introduce scaling by A->psurf to improve scaling for numerical solver (FD Jacobian) */
            massv = Ap->reaction_parameters[i]->stoichiometry[j] * Ap->volatile_parameters[v].molar_mass * (A->volatiles[v].p/A->psurf);
            massv /= factor;
            Ap->volatile_parameters[v].initial -= massv * A->mass_reaction[i] / (*Ap->mantle_mass_ptr);
          }
        }

        for (i=0; i<Ap->n_volatiles; ++i) {
            x0 = A->volatiles[i].x;
            ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],i,x0,INSERT_VALUES);CHKERRQ(ierr);
        }
    }
    else{
        for (i=0; i<Ap->n_volatiles; ++i) {
            x0 = 0.0;
            ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],i,x0,INSERT_VALUES);CHKERRQ(ierr);
        }
    }

    /* Initialize reaction amounts to zero */
    for (i=0; i<Ap->n_reactions; ++i) {
        ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]],i,0.0,INSERT_VALUES);CHKERRQ(ierr);
    }

    for (i=0; i<E->numFields; ++i) {
      ierr = VecAssemblyBegin(subVecs[i]);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(subVecs[i]);CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    ierr = PetscFree(subVecs);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_initial_volatile( Ctx *E )
{
    PetscErrorCode ierr;
    SNES           snes;
    Vec            x,r;
    PetscScalar    *xx;
    PetscInt       i;
    Atmosphere                 *A = &E->atmosphere;
    Parameters           const *P = &E->parameters;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;

    PetscFunctionBeginUser;

    ierr = SNESCreate( PETSC_COMM_WORLD, &snes );CHKERRQ(ierr);

    /* Use this to address this specific SNES (nonlinear solver) from the command
       line or options file, e.g. -atmosic_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"atmosic_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    /* Vector with one dof per volatile species, and one dof per volatile reaction */
    ierr = VecSetSizes( x, PETSC_DECIDE, Ap->n_volatiles + Ap->n_reactions );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    ierr = SNESSetFunction(snes,r,FormFunction1,E);CHKERRQ(ierr);

    /* initialise vector x with initial guess */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    for (i=0; i<Ap->n_volatiles; ++i) {
        xx[i] = Ap->volatile_parameters[i].initial;
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
                "Unphysical initial volatile concentration: volatile %d, x: %g",i,A->volatiles[i].x);
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

/* Non-linear solver for initial volatile abundance */
static PetscErrorCode FormFunction1( SNES snes, Vec x, Vec f, void *ptr)
{
    PetscErrorCode             ierr;
    const PetscScalar          *xx;
    PetscScalar                *ff;
    const PetscScalar          *mass_r;
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
    mass_r = &xx[Ap->n_volatiles]; /* reaction masses (also to be determined by solver) */

    ierr = set_atmosphere_volatile_content( A, Ap, C ); CHKERRQ(ierr);

    /* Balance equation for volatile abundance */
    for (i=0; i<Ap->n_volatiles; ++i) {
        /* TODO: DJB rename this function */
        ff[i] = get_initial_volatile_abundance( A, Ap, &Ap->volatile_parameters[i], &A->volatiles[i]);
    }

    /* Subtract or add reaction masses */
    for (i=0; i<Ap->n_reactions; ++i) {
      PetscInt j;
      PetscScalar factor; /* normalisation, using the first volatile (reactant) in this chemical reaction */
      const PetscInt v0 = Ap->reaction_parameters[i]->volatiles[0];
      /* by convention, first volatile is a reactant, so stoichiometry (hence factor) will be -ve */
      /* NOTE: introduce scaling by A->psurf to improve scaling for numerical solver (FD Jacobian) */
      /* TODO: swap out A->volatiles[v0].p/A->psurf for the volume mixing ratio? */
      factor = Ap->reaction_parameters[i]->stoichiometry[0] * Ap->volatile_parameters[v0].molar_mass * (A->volatiles[v0].p/A->psurf);
      for (j=0; j<Ap->reaction_parameters[i]->n_volatiles; ++j) {
        const PetscInt v = Ap->reaction_parameters[i]->volatiles[j];
        PetscScalar massv;
        /* NOTE: introduce scaling by A->psurf to improve scaling for numerical solver (FD Jacobian) */
        massv = Ap->reaction_parameters[i]->stoichiometry[j] * Ap->volatile_parameters[v].molar_mass * (A->volatiles[v].p/A->psurf);
        massv /= factor; /* +ve for reactants, -ve for products */
        ff[v] += massv * mass_r[i]; /* convention is mass loss of reactants, mass gain of products */
        /* but regarding above, convention is arbitrary since mass_r[i] will switch sign to balance */
      }
    }

    /* Objective function */
    for (i=0; i<Ap->n_reactions; ++i) {

        PetscScalar Qp,Qr,K;

        ReactionParameters const * reaction_parameters_ptr = &Ap->reaction_parameters[i];

        /* return the numerator and denominator separately to retain scalings,
           otherwise the non-linear solver has problems */
        /* NOTE: since scaling by A->psurf, it may not be necessary any longer
           to split the reaction quotient into two parts, but nevertheless
           why fix something that is not broken? */
        Qp = get_reaction_quotient_products( reaction_parameters_ptr, A );
        Qr = get_reaction_quotient_reactants( reaction_parameters_ptr, A );

        K = get_equilibrium_constant( reaction_parameters_ptr, A->tsurf, C );

        /* this general form retains the pressure scaling and is preferred by the solver:
           passing in the reaction quotient rather than decomposed into a numerator and
           denominator causes problems in computing the FD Jacobian */
        ff[Ap->n_volatiles + i] = Qp - K * Qr;

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
