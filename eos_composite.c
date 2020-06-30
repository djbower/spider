#include "eos_composite.h"

PetscErrorCode EOSCompositeCreate(EOSComposite* p) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, p);CHKERRQ(ierr);
  // TODO other stuff to init? (to zero?)
  PetscFunctionReturn(0);
}
