#include "poststep.h"
#include "ctx.h"

PetscErrorCode PostStepDataInitialize(Ctx *pctx, Vec X)
{
  PetscErrorCode ierr;
  PetscInt i;
  PostStepData *data;
  Parameters const *P = &pctx->parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc(sizeof(PostStepData),&(pctx->postStepData));CHKERRQ(ierr);
  // FIXME should use the usual methods with DMs, dimensionalisable fields etc. to get the correct entry out..
  ierr = PetscPrintf(PetscObjectComm((PetscObject)X),"PLACEHOLDER I'm %s in %s:%d, Populating ctx->postStepData\n",__func__,__FILE__,__LINE__);
  data = (PostStepData*) pctx->postStepData;

  // FIXME: better way to find the indices?
  /* save current volatile concentrations in the melt */
  i = P->numpts_b + 1;
  ierr = VecGetValues(X,1,&i,&(data->x_CO2));CHKERRQ(ierr);
  i = P->numpts_b + 2;
  ierr = VecGetValues(X,1,&i,&(data->x_H2O));CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode PostStep(TS ts)
{
  PetscErrorCode ierr;
  PetscInt       stepNumber;
  Ctx            *ctx;
  PostStepData   *data;
  AtmosphereParameters *Ap;

  PetscFunctionBeginUser;
  ierr = TSGetTimeStepNumber(ts,&stepNumber);CHKERRQ(ierr);
  ierr = TSGetApplicationContext(ts,&ctx);CHKERRQ(ierr);
  Ap = &ctx->parameters.atmosphere_parameters;
  if (!ctx->parameters.rollBackActive) SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_SUP,"You must run with -activate_rollback");
  data = (PostStepData*) ctx->postStepData;
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, x_CO2 is %g\n",__func__,__FILE__,__LINE__,data->x_CO2);
  // FIXME djbower add actual check (pull out current data with TSGetSolution(), use usual methods to access fields..)
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, x_H2O is %g\n",__func__,__FILE__,__LINE__,data->x_H2O);
  // FIXME djbower add actual check (pull out current data with TSGetSolution(), use usual methods to access fields..)

  if ((data->x_CO2 > 1.1 * Ap->CO2_parameters.initial) || (data->x_H2O > 1.1* Ap->H2O_parameters.initial)){
    ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, rolling back solution and setting flags\n",__func__,__FILE__,__LINE__);
    ierr = TSRollBack(ts);CHKERRQ(ierr);
    ierr = TSSetConvergedReason(ts,TS_CONVERGED_USER);CHKERRQ(ierr);
    ctx->stopEarly = PETSC_TRUE;
  }

  //if (stepNumber > 4) {
  //  ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, rolling back solution and setting flags\n",__func__,__FILE__,__LINE__);
  //  ierr = TSRollBack(ts);CHKERRQ(ierr);
   // ierr = TSSetConvergedReason(ts,TS_CONVERGED_USER);CHKERRQ(ierr);
   // ctx->stopEarly = PETSC_TRUE;
  //}

  PetscFunctionReturn(0);
}
