#include "output.h"

PetscErrorCode add_vector_to_binary_output( Vec vec, PetscViewer viewer)
{
  /* simple wrapper to add a vector to a PetscViewer */

  PetscErrorCode ierr;
  //char vecname[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  // convenient for some output formats, but not Petsc binary
  //ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"phi_s_%lld",step);CHKERRQ(ierr);
  //ierr = PetscObjectSetName((PetscObject)x_aug,vecname);CHKERRQ(ierr);
  ierr = VecView( vec, viewer ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


