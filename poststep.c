#include "poststep.h"
#include "ctx.h"

PetscErrorCode PostStepDataInitialize(Ctx *pctx, Vec X)
{
  PetscErrorCode ierr;
  PetscInt i = 5; // FIXME just plucking out arbitrary entry..
  PostStepData *data;

  PetscFunctionBeginUser;
  ierr = PetscMalloc(sizeof(PostStepData),&(pctx->postStepData));CHKERRQ(ierr);
  // FIXME should use the usual methods with DMs, dimensionalisable fields etc. to get the correct entry out..
  ierr = PetscPrintf(PetscObjectComm((PetscObject)X),"PLACEHOLDER I'm %s in %s:%d, Populating ctx->postStepData\n",__func__,__FILE__,__LINE__);
  data = (PostStepData*) pctx->postStepData;
  ierr = VecGetValues(X,1,&i,&(data->ref_concentration));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PostStep(TS ts)
{
  PetscErrorCode ierr;
  PetscInt       stepNumber;
  Ctx            *ctx;
  PostStepData   *data;

  PetscFunctionBeginUser;
  ierr = TSGetTimeStepNumber(ts,&stepNumber);CHKERRQ(ierr);
  ierr = TSGetApplicationContext(ts,&ctx);CHKERRQ(ierr);
  if (!ctx->parameters.rollBackActive) SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_SUP,"You must run with -activate_rollback");
  data = (PostStepData*) ctx->postStepData;
  ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, ref_concentration is %g\n",__func__,__FILE__,__LINE__,data->ref_concentration);
  // FIXME djbower add actual check (pull out current data with TSGetSolution(), use usual methods to access fields..)
  if (stepNumber > 4) {
    ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, rolling back solution and setting flags\n",__func__,__FILE__,__LINE__);
    ierr = TSRollBack(ts);CHKERRQ(ierr);
    ierr = TSSetConvergedReason(ts,TS_CONVERGED_USER);CHKERRQ(ierr);
    ctx->stopEarly = PETSC_TRUE;
  }
  PetscFunctionReturn(0);
}
