#include "parameters.h"

/* Note: this could logically be included in parameters.c, but that file was getting crowded */

/* Create a simple reaction, involving 2 volatiles, described by "exchange rate"
   factors gamma0 and gamma1, and an equilibrium ratio factors epsilon0 and epsilon1 */
PetscErrorCode ReactionParametersCreateSimple(ReactionParameters* reaction_parameters_ptr,PetscInt v0,PetscInt v1,PetscReal gamma0,PetscReal gamma1,PetscReal epsilon0, PetscReal epsilon1)
{
  PetscErrorCode     ierr;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->gamma,2,&reaction_parameters->volatiles,2,&reaction_parameters->epsilon);CHKERRQ(ierr);
  reaction_parameters->gamma[0] = gamma0;
  reaction_parameters->gamma[1] = gamma1;
  reaction_parameters->epsilon[0] = epsilon0;
  reaction_parameters->epsilon[1] = epsilon1;
  reaction_parameters->volatiles[0] = v0;
  reaction_parameters->volatiles[1] = v1;
  PetscFunctionReturn(0);
}

PetscErrorCode ReactionParametersDestroy(ReactionParameters* reaction_parameters_ptr)
{
  PetscErrorCode     ierr;

  PetscFunctionBeginUser;
  ierr = PetscFree3((*reaction_parameters_ptr)->gamma,(*reaction_parameters_ptr)->volatiles,(*reaction_parameters_ptr)->epsilon);CHKERRQ(ierr);
  ierr = PetscFree(*reaction_parameters_ptr);CHKERRQ(ierr);/* must be last */
  *reaction_parameters_ptr = NULL;
  PetscFunctionReturn(0);
}
