#include "eos_rtpress.h"

PetscErrorCode EOSRTpressPureCreate(EOSRTpressPure* p) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, p);CHKERRQ(ierr);
  // TODO other stuff to init? (to zero?)
  PetscFunctionReturn(0);
}
