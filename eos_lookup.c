#include "eos_lookup.h"

PetscErrorCode EOSLookupPureCreate(EOSLookupPure* p) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, p);CHKERRQ(ierr);
  // TODO other stuff to init? (to zero?)
  PetscFunctionReturn(0);
}
