#include "ic.h"
#include "bc.h"
#include "util.h"

static PetscErrorCode set_ic_default( Ctx *, Vec );
static PetscErrorCode set_ic_entropy( Ctx *, Vec );
static PetscErrorCode set_ic_atmosphere( Ctx *, Vec );
static PetscErrorCode set_ic_from_file( Ctx *, Vec );

PetscErrorCode set_initial_condition( Ctx *E, Vec sol)
{
    PetscErrorCode ierr;
    PetscMPIInt rank,size;
    PetscInt ind,numpts_s;
    Solution *S = &E->solution;
    Parameters const *P = &E->parameters;
    PetscInt IC = P->initial_condition;
    PetscInt const ind0=0;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ind = numpts_s-1;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

    /* TODO: these functions do not typically (consistently)
       update entries in the vectors stored in Solution, they
       simply write the data to the solution vector (sol) */

    if(IC==1){
        ierr = set_ic_default( E, sol ); CHKERRQ(ierr);
    }
    else if(IC==2){
        ierr = set_ic_from_file( E, sol ); CHKERRQ(ierr);
    }

    /* FIXME: will break in parallel */
    if( (P->ic_surface_entropy > 0.0) || (P->ic_core_entropy > 0.0) ){
        ierr = set_entropy_from_solution( E, sol );
        if( P->ic_surface_entropy > 0.0 ){
            ierr = VecSetValues( S->S_s,1,&ind0,&P->ic_surface_entropy,INSERT_VALUES );CHKERRQ(ierr);
        }
        if( P->ic_core_entropy > 0.0 ){
            ierr = VecSetValues( S->S_s,1,&ind,&P->ic_core_entropy,INSERT_VALUES );CHKERRQ(ierr);
        }
        ierr = VecAssemblyBegin(S->S_s);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(S->S_s);CHKERRQ(ierr);
        ierr = set_solution_from_entropy( E, sol );CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_default( Ctx *E, Vec sol )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = set_ic_entropy( E, sol ); CHKERRQ(ierr);
    ierr = set_ic_atmosphere( E, sol ); CHKERRQ(ierr);

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
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_constant_dSdr:\n");CHKERRQ(ierr);

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


static PetscErrorCode set_ic_atmosphere( Ctx *E, Vec sol )
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
    if(Ap->SOLVE_FOR_VOLATILES || Ap->SURFACE_BC==3){
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

    /* FIXME: currently broken since change to new output format (json) */

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
