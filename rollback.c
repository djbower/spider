#include "rollback.h"

/* Private header for use ONLY in TSRollBackGenericActivate() */
#include <petsc/private/tsimpl.h>

static PetscErrorCode TSRollBack_Generic(TS ts)
{
  PetscErrorCode ierr;
  Vec            X,Xprev;

  PetscFunctionBegin;
  ierr = PetscObjectQuery((PetscObject)ts,"RollBackGeneric_Xprev",(PetscObject*)(&Xprev));CHKERRQ(ierr);
  if (!Xprev) SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_ARG_WRONGSTATE,"Rollback data not stored before TSRollBack() called (no steps taken?)");
  ierr = TSGetSolution(ts,&X);CHKERRQ(ierr);
  ierr = VecCopy(Xprev,X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PreStep_RollBackGeneric(TS ts)
{
  PetscErrorCode ierr;
  Vec            X,Xprev;

  PetscFunctionBegin;
  ierr = TSGetSolution(ts,&X);CHKERRQ(ierr);
  ierr = PetscObjectQuery((PetscObject)ts,"RollBackGeneric_Xprev",(PetscObject*)(&Xprev));CHKERRQ(ierr);
  if (!Xprev) {
    ierr = VecDuplicate(X,&Xprev);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject)ts,"RollBackGeneric_Xprev",(PetscObject)Xprev);CHKERRQ(ierr);
  }
  ierr = VecCopy(X,Xprev);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TSRollBackGenericActivate(TS ts)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (ts->prestep) SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_ARG_WRONGSTATE,"Generic prestep implementation requires setting the prestep function, which has already been set");
  ierr = TSSetPreStep(ts,PreStep_RollBackGeneric);CHKERRQ(ierr);
  if (!ts->setupcalled) SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_ARG_WRONGSTATE,"Cannot activate generic rollback before calling TSSetUp()");
  if (ts->ops->rollback) SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_ARG_WRONGSTATE,"Refusing to activate generic rollback for a TS implementation that already implements TSRollBack()");
  ts->ops->rollback = TSRollBack_Generic;
  PetscFunctionReturn(0);
}

PetscErrorCode TSRollBackGenericDestroy(TS ts)
{
  PetscErrorCode ierr;
  Vec            Xprev;

  ierr = PetscObjectQuery((PetscObject)ts,"RollBackGeneric_Xprev",(PetscObject*)(&Xprev));CHKERRQ(ierr);
  if (Xprev) {
    ierr = VecDestroy(&Xprev);CHKERRQ(ierr);
    ierr = PetscObjectCompose((PetscObject)ts,"RollBackGeneric_Xprev",NULL);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

