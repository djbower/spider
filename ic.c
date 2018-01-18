#include "ic.h"
#include "bc.h"
#include "aug.h"

static PetscErrorCode set_ic_default( Ctx *, Vec );
static PetscErrorCode set_ic_constant_dSdr( Ctx *, Vec );
static PetscErrorCode set_ic_extra( Ctx *, Vec );
static PetscErrorCode set_ic_from_file( Ctx *, Vec );
static PetscErrorCode set_ic_from_impact( Ctx *, Vec );

PetscErrorCode set_initial_condition( Ctx *E, Vec dSdr_b_aug)
{
    PetscErrorCode ierr;
    Parameters const *P = &E->parameters;
    PetscInt IC = P->initial_condition;

    PetscFunctionBeginUser;

    if( IC==1 ){
        ierr = set_ic_default( E, dSdr_b_aug ); CHKERRQ(ierr);
    }
    else if( IC==2 ){
        ierr = set_ic_from_file( E, dSdr_b_aug ); CHKERRQ(ierr);
    }
    else if( IC==3 ){
        ierr = set_ic_from_impact( E, dSdr_b_aug ); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);

}

static PetscErrorCode set_ic_default( Ctx *E, Vec dSdr_b_aug )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = set_ic_constant_dSdr( E, dSdr_b_aug ); CHKERRQ(ierr);
    ierr = set_ic_extra( E, dSdr_b_aug ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_ic_constant_dSdr( Ctx *E, Vec dSdr_b_aug )
{
    /* set initial entropy gradient to constant for all nodes */

    PetscErrorCode   ierr;
    Parameters const *P = &E->parameters;
    Vec dSdr_b;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_constant_dSdr:\n");CHKERRQ(ierr);
#endif

    ierr = DMCreateGlobalVector( E->da_b, &dSdr_b );CHKERRQ(ierr);

    ierr = VecSet( dSdr_b, P->ic_dsdr ); CHKERRQ( ierr );

    /* these next two lines are simply convenient reminders that the
       first and last values are meaningless because the fluxes here
       are controlled by boundary conditions.  But for debugging and
       clarity it is convenient to explicitly set these values to
       zero */
    ierr = VecSetValue( dSdr_b, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue( dSdr_b, P->numpts_b-1, 0.0, INSERT_VALUES); CHKERRQ(ierr);

    VecAssemblyBegin( dSdr_b );
    VecAssemblyEnd( dSdr_b );

    ierr = ToAug( dSdr_b, dSdr_b_aug ); CHKERRQ(ierr);

    ierr = VecDestroy( &dSdr_b ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


static PetscErrorCode set_ic_extra( Ctx *E, Vec dSdr_b_aug )
{
    /* set initial condition for the extra points */

    PetscErrorCode ierr;
    Parameters           const *P  = &E->parameters;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    VolatileParameters const *CO2 = &Ap->CO2_volatile_parameters;
    VolatileParameters const *H2O = &Ap->H2O_volatile_parameters;

    PetscScalar x0, x1;

    PetscFunctionBeginUser;

    /* augmented vector organised as follows:
         0.     CO2 (dissolved content of CO2 in magma ocean)
         1.     H2O (dissolved content of H2O in magma ocean)
         2.     S0 (entropy at uppermost staggered node)
         3-end  dS/dr at basic nodes (see set_ic_dSdr)
    */

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

    ierr = VecSetValue(dSdr_b_aug,0,x0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(dSdr_b_aug,1,x1,INSERT_VALUES);CHKERRQ(ierr);
   
    /* include initial entropy at staggered node */
    ierr = VecSetValue(dSdr_b_aug,2,P->sinit,INSERT_VALUES);CHKERRQ(ierr);

    VecAssemblyBegin( dSdr_b_aug );
    VecAssemblyEnd( dSdr_b_aug );

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_from_file( Ctx *E, Vec dSdr_b_aug )
{

    /* reads in the output from dSdr_b_aug_[timestep].m to use as the
       initial condition to enable restarting */

    PetscErrorCode ierr;
    Parameters const *P  = &E->parameters;
    FILE *fp;
    const PetscInt head=3;
    PetscInt i=0, j=0;
    char string[PETSC_MAX_PATH_LEN];
#if (defined PETSC_USE_REAL___FLOAT128)
    char xtemp[30]
#endif
    PetscScalar x;

    PetscFunctionBeginUser;
    fp = fopen( P->ic_filename, "r" );

    if(fp==NULL) {
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_FILE_OPEN,"Could not open file %s",P->ic_filename);
    }

    while(fgets(string, sizeof(string), fp) != NULL) {
        if( (i>=head) && (i<=P->numpts_b+head+2) ){ /* 3 header lines in dSdr_b_aug_[timestep].m */
#if (defined PETSC_USE_REAL___FLOAT128)
            sscanf( string, "%s", xtemp );
            x = strtoflt182(xtemp, NULL);
#else
            sscanf(string, "%lf", &x );
#endif            
            j = i-3;
            ierr = VecSetValue(dSdr_b_aug,j,x,INSERT_VALUES); CHKERRQ(ierr);
        }
        ++i;
    }

    VecAssemblyBegin( dSdr_b_aug );
    VecAssemblyEnd( dSdr_b_aug );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_ic_from_impact( Ctx *E, Vec dSdr_b_aug )
{
   PetscFunctionBeginUser;

   /* place holder */

   PetscFunctionReturn(0);
}
