#include "atmosphere.h"
#include "ic.h"

static PetscErrorCode set_ic_dSdr_constant( Ctx *, Vec );

PetscErrorCode set_ic_aug( Ctx *E, Vec dSdr_b_aug )
{

    PetscErrorCode ierr;
    PetscScalar    xCO2, xH2O;

    PetscFunctionBeginUser;

    /* augmented vector organised as follows:
         0.     CO2 wt. % (dissolved content of CO2 in magma ocean)
         1.     H2O wt. % (dissolved content of H2O in magma ocean)
         2.     S0 (entropy at uppermost staggered node)
         3-end  dS/dr at basic nodes (see set_ic_dSdr)
    */

    /* TODO: this should be updated.  Should be consistent with the
      mass balance of volatiles, so the total initial volatile content
      is actually partitioned betweem the (liquid) magma ocean and the
      atmosphere. */
    /* initial mass content of CO2 in the magma ocean */
    xCO2 = get_initial_xCO2( &E->atmosphere );
    ierr = VecSetValue(dSdr_b_aug,0,xCO2,INSERT_VALUES);CHKERRQ(ierr);
    
    /* initial mass content of H2O in the magma ocean */
    xH2O = get_initial_xH2O( &E->atmosphere );
    ierr = VecSetValue(dSdr_b_aug,1,xH2O,INSERT_VALUES);CHKERRQ(ierr);
   
   /* include initial entropy at staggered node */
   /* TODO: is this over-rideable from the command line? */
    ierr = VecSetValue(dSdr_b_aug,2,SINIT_DEFAULT,INSERT_VALUES);CHKERRQ(ierr);

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

    PetscErrorCode     ierr;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_dSdr_constant:\n");CHKERRQ(ierr);
#endif

    ierr = VecSet( dSdr_in, IC_DSDR ); CHKERRQ( ierr );

    /* these next two lines are simply convenient reminders that the
       first and last values are meaningless because the fluxes here
       are controlled by boundary conditions.  But for debugging and
       clarity it is convenient to explicitly set these values to
       zero */
    ierr = VecSetValue( dSdr_in, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
    // TODO: is next line compatible with command line input? */
    ierr = VecSetValue( dSdr_in, NUMPTS_B_DEFAULT-1, 0.0, INSERT_VALUES); CHKERRQ(ierr);

    VecAssemblyBegin( dSdr_in );
    VecAssemblyEnd( dSdr_in );

    PetscFunctionReturn(0);
}
