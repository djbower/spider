#include "ic.h"
#include "bc.h"
#include "util.h"

static PetscErrorCode set_ic_default( Ctx *, Vec );
static PetscErrorCode set_ic_constant_dSdr( Ctx *, Vec );
static PetscErrorCode set_ic_extra( Ctx *, Vec );
static PetscErrorCode set_ic_from_file( Ctx *, Vec );
#if 0
static PetscErrorCode set_ic_from_impact( Ctx *, Vec );
#endif

PetscErrorCode set_initial_condition( Ctx *E, Vec sol)
{
    PetscErrorCode ierr;
    Parameters const *P = &E->parameters;
    PetscInt IC = P->initial_condition;

    PetscFunctionBeginUser;

    /* TODO: these functions do not typically (consistently)
       update entries in the vectors stored in Solution, they
       simply write the data to the solution vector (sol) */

    if( IC==1 ){
        ierr = set_ic_default( E, sol ); CHKERRQ(ierr);
    }
    else if( IC==2 ){
        ierr = set_ic_from_file( E, sol ); CHKERRQ(ierr);
    }

#if 0
    else if( IC==3 ){
        ierr = set_ic_from_impact( E, sol ); CHKERRQ(ierr);
    }
#endif

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_default( Ctx *E, Vec sol )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = set_ic_constant_dSdr( E, sol ); CHKERRQ(ierr);
    ierr = set_ic_extra( E, sol ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_ic_constant_dSdr( Ctx *E, Vec sol )
{
    /* set initial entropy gradient to constant for all nodes */

    PetscErrorCode   ierr;
    Parameters const *P = &E->parameters;
    Vec              dSdr_b;
    Vec              *subVecs;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_constant_dSdr:\n");CHKERRQ(ierr);
#endif

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    dSdr_b = subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDR_B]];

    ierr = VecSet(dSdr_b, P->ic_dsdr ); CHKERRQ( ierr );

    /* these next two lines are simply convenient reminders that the
       first and last values are meaningless because the fluxes here
       are controlled by boundary conditions.  But for debugging and
       clarity it is convenient to explicitly set these values to
       zero */
    ierr = VecSetValue( dSdr_b, 0,             0.0, INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue( dSdr_b, P->numpts_b-1, 0.0, INSERT_VALUES);CHKERRQ(ierr);

    ierr = VecAssemblyBegin( dSdr_b );CHKERRQ(ierr);
    ierr = VecAssemblyEnd( dSdr_b );CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    ierr = PetscFree(subVecs);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


static PetscErrorCode set_ic_extra( Ctx *E, Vec sol )
{
    /* set initial condition for the extra points */

    PetscErrorCode             ierr;
    PetscInt                   i;
    Vec                        *subVecs;
    Parameters           const *P  = &E->parameters;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    VolatileParameters const   *CO2 = &Ap->CO2_volatile_parameters;
    VolatileParameters const   *H2O = &Ap->H2O_volatile_parameters;

    PetscScalar x0, x1;

    PetscFunctionBeginUser;

    /* initial volatile content */
    x0 = 0.0;
    x1 = 0.0;

    /* turn on volatiles for these conditions */
    if(Ap->SOLVE_FOR_VOLATILES || Ap->MODEL==3){
      /* CO2 */
      if( CO2->initial > 0.0 ){
        x0 = get_initial_volatile( E, CO2 );
      }
      /* H2O */
      if( H2O->initial > 0.0 ){
        x1 = get_initial_volatile( E, H2O );
      }
    }

    if( x0 < 0.0 || x1 < 0 ){
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Initial volatile content cannot be negative: %d",x0);
    }

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_CO2]],0,x0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_H2O]],0,x1,INSERT_VALUES);CHKERRQ(ierr);

    /* include initial entropy at staggered node */
    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],0,P->sinit,INSERT_VALUES);CHKERRQ(ierr);

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

    /* reads in the output from sol_[timestep].m to use as the
       initial condition to enable restarting

       Warning: this assumes a particular ordering of the fields in the solution
       vector, so can only assumed to be valid when the file is generated
       by this exact version of SPIDER.
       */

    PetscErrorCode   ierr;
    Parameters const *P  = &E->parameters;
    FILE             *fp;
    PetscInt         i=0, j=0;
    char             string[PETSC_MAX_PATH_LEN];
#if (defined PETSC_USE_REAL___FLOAT128)
    char             xtemp[30];
#endif
    PetscScalar      x;
    const PetscInt   head=3;

    PetscFunctionBeginUser;
    fp = fopen( P->ic_filename, "r" );

    if(fp==NULL) {
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_FILE_OPEN,"Could not open file %s",P->ic_filename);
    }

    while(fgets(string, sizeof(string), fp) != NULL) {
        if( (i>=head) && (i<=P->numpts_b+head+2) ){ /* 3 header lines in sol_[timestep].m */
#if (defined PETSC_USE_REAL___FLOAT128)
            sscanf( string, "%s", xtemp );
            x = strtoflt128(xtemp, NULL);
#else
            sscanf(string, "%lf", &x );
#endif
            j = i-3;
            ierr = VecSetValue(sol,j,x,INSERT_VALUES); CHKERRQ(ierr);
        }
        ++i;
    }

    VecAssemblyBegin( sol );
    VecAssemblyEnd( sol );

    PetscFunctionReturn(0);
}

// PDS: not sure if this is current, so commenting out
#if 0
static PetscErrorCode set_ic_from_impact( Ctx *E, Vec sol )
{
   PetscErrorCode ierr;
   Solution *S = &E->solution;
   PetscScalar S0;
   PetscInt    ind;

   PetscFunctionBeginUser;

   ierr = set_ic_from_file( E, sol ); CHKERRQ(ierr);
   /* for simplicity, offset the tie point by some fixed fraction of
      entropy */
   //ierr = VecSetValue( sol,2,0.2,ADD_VALUES ); CHKERRQ(ierr);

   ind = 2;
   ierr = VecGetValues( sol,1,&ind,&S0); CHKERRQ(ierr);

   /* need to transfer to S->dSdr otherwise set_entropy returns
      constant profile (since dSdr=0 everywhere initially) */
   ierr = FromAug( sol,S->dSdr); CHKERRQ(ierr);

   /* integrate to give S_s and S_b */
   ierr = set_entropy( E, S0 );

   /* TODO: set the entropy to above the liquidus for a certain depth
      range to mimick an impactor that melts a portion of the mantle */

   /* function to go the other way, i.e. compute dSdr from S_s */
   ierr = set_dSdr_b_from_S_s( E );

   PetscFunctionReturn(0);
}
#endif
