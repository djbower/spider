#include "constants.h"

PetscErrorCode FundamentalConstantsSet( FundamentalConstants FC, ScalingConstants const SC )
{
    /* fundamental physical and chemical constants are all here, and should
       always be accessed from this struct to avoid duplication */

    PetscFunctionBeginUser;

    FC->AVOGADRO = 6.02214076E23; /* 1/mol */

    FC->GAS = 8.3144598; /* J/K/mol */
    FC->GAS *= SC->TEMP / SC->ENERGY; /* 1/mol */

    /* BOLTZMANN is 1.380649E-23 J/K */
    /* note that BOLTZMANN is also GAS constant per particle */
    FC->BOLTZMANN = FC->GAS / FC->AVOGADRO;

    FC->GRAVITATIONAL = 6.67408E-11; /* m^3/kg/s^2 */
    FC->GRAVITATIONAL *= SC->DENSITY * PetscPowScalar( SC->TIME, 2.0 );

    FC->STEFAN_BOLTZMANN = 5.670367e-08; /* W/m^2/K^4 */
    FC->STEFAN_BOLTZMANN /= SC->SIGMA;

    PetscFunctionReturn(0);

}

PetscErrorCode ScalingConstantsCreate( ScalingConstants* scaling_constants_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(1,scaling_constants_ptr);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode ScalingConstantsDestroy( ScalingConstants* scaling_constants_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
 
    ierr = PetscFree(*scaling_constants_ptr);CHKERRQ(ierr);
    *scaling_constants_ptr = NULL;
    PetscFunctionReturn(0);
}

PetscErrorCode FundamentalConstantsCreate( FundamentalConstants* fundamental_constants_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(1,fundamental_constants_ptr);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode FundamentalConstantsDestroy( FundamentalConstants* fundamental_constants_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscFree(*fundamental_constants_ptr);CHKERRQ(ierr);
    *fundamental_constants_ptr = NULL;

    PetscFunctionReturn(0);
}
