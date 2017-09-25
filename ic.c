#include "ic.h"
#include "bc.h"

static PetscErrorCode set_ic_dSdr_constant( Ctx *, Vec );

PetscErrorCode set_ic_aug( Ctx *E, Vec dSdr_b_aug )
{

    PetscErrorCode ierr;
    Mesh                 *M  = &E->mesh;
    Parameters           *P  = &E->parameters;
    AtmosphereParameters *Ap = &P->atmosphere_parameters;

    PetscScalar x0, x1;

    PetscFunctionBeginUser;

    /* augmented vector organised as follows:
         0.     CO2 (dissolved content of CO2 in magma ocean)
         1.     H2O (dissolved content of H2O in magma ocean)
         2.     S0 (entropy at uppermost staggered node)
         3-end  dS/dr at basic nodes (see set_ic_dSdr)
    */

    /* initial mass content of CO2 in the magma ocean */
    x0 = get_initial_volatile( Ap, &Ap->CO2_volatile_parameters, M->mantle_mass );
    ierr = VecSetValue(dSdr_b_aug,0,x0,INSERT_VALUES);CHKERRQ(ierr);
    
    /* initial mass content of H2O in the magma ocean */
    x1 = get_initial_volatile( Ap, &Ap->H2O_volatile_parameters, M->mantle_mass );
    ierr = VecSetValue(dSdr_b_aug,1,x1,INSERT_VALUES);CHKERRQ(ierr);
   
    /* include initial entropy at staggered node */
    ierr = VecSetValue(dSdr_b_aug,2,P->sinit,INSERT_VALUES);CHKERRQ(ierr);

    VecAssemblyBegin( dSdr_b_aug );
    VecAssemblyEnd( dSdr_b_aug );

    PetscFunctionReturn(0);
}


PetscErrorCode set_ic_dSdr( Ctx *E, Vec dSdr_in) 
{
    /* initial entropy gradient for basic nodes */

    PetscFunctionBeginUser;

    set_ic_dSdr_constant( E, dSdr_in );

    /* add more initial condition options here if desired */

    PetscFunctionReturn(0);
}

static PetscErrorCode set_ic_dSdr_constant( Ctx *E, Vec dSdr_in )
{
    /* set initial entropy gradient to constant for basic nodes */

    PetscErrorCode   ierr;
    Parameters const *P = &E->parameters;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_dSdr_constant:\n");CHKERRQ(ierr);
#endif

    ierr = VecSet( dSdr_in, P->ic_dsdr ); CHKERRQ( ierr );

    /* these next two lines are simply convenient reminders that the
       first and last values are meaningless because the fluxes here
       are controlled by boundary conditions.  But for debugging and
       clarity it is convenient to explicitly set these values to
       zero */
    ierr = VecSetValue( dSdr_in, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue( dSdr_in, P->numpts_b-1, 0.0, INSERT_VALUES); CHKERRQ(ierr);

    VecAssemblyBegin( dSdr_in );
    VecAssemblyEnd( dSdr_in );

    PetscFunctionReturn(0);
}
