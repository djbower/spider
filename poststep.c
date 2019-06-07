#include "poststep.h"
#include "ctx.h"

PetscErrorCode PostStepDataInitialize(Ctx *E, Vec sol_in)
{
  PetscErrorCode   ierr;
  PostStepData     *data;
  Vec              *subVecs;
  PetscScalar      x0, x1;

  PetscFunctionBeginUser;
  ierr = PetscMalloc(sizeof(PostStepData),&(E->postStepData));CHKERRQ(ierr);
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)X),"PLACEHOLDER I'm %s in %s:%d, Populating ctx->postStepData\n",__func__,__FILE__,__LINE__);
  data = (PostStepData*) E->postStepData;

  ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
  ierr = DMCompositeGetAccessArray(E->dm_sol,sol_in,E->numFields,NULL,subVecs);CHKERRQ(ierr);

  /* concentration of volatiles in melt
     currently these are the only criteria that are used to halt the code */
  {
    const PetscInt ind = SPIDER_VOLATILE_CO2;
    ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],1,&ind,&x0);CHKERRQ(ierr);
  }
  data->x_CO2 = x0;
  {
    const PetscInt ind = SPIDER_VOLATILE_H2O;
    ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],1,&ind,&x1);CHKERRQ(ierr);
  }
  data->x_H2O = x1;

  ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol_in,E->numFields,NULL,subVecs);CHKERRQ(ierr);
  ierr = PetscFree(subVecs);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode PostStep(TS ts)
{
  PetscErrorCode ierr;
  PetscInt       stepNumber;
  Ctx            *E;
  PostStepData   *data;
  Parameters     *P;
  Vec            sol_in;
  PetscScalar    x_CO2, x_H2O, rel_CO2, rel_H2O, max_CO2, max_H2O;
  Vec            *subVecs;
  PetscBool      EVENT = PETSC_FALSE;

  PetscFunctionBeginUser;
  ierr = TSGetTimeStepNumber(ts,&stepNumber);CHKERRQ(ierr);
  ierr = TSGetApplicationContext(ts,&E);CHKERRQ(ierr);
  P = &E->parameters;
  if (!E->parameters.rollBackActive) SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_SUP,"You must run with -activate_rollback");
  data = (PostStepData*) E->postStepData;
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, x_CO2 is %g\n",__func__,__FILE__,__LINE__,data->x_CO2);
  // FIXME djbower add actual check (pull out current data with TSGetSolution(), use usual methods to access fields..)
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, x_H2O is %g\n",__func__,__FILE__,__LINE__,data->x_H2O);
  // FIXME djbower add actual check (pull out current data with TSGetSolution(), use usual methods to access fields..)

  ierr = TSGetSolution(ts,&sol_in);CHKERRQ(ierr);

  ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
  ierr = DMCompositeGetAccessArray(E->dm_sol,sol_in,E->numFields,NULL,subVecs);CHKERRQ(ierr);
  {
    const PetscInt ind = SPIDER_VOLATILE_CO2;
    ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],1,&ind,&x_CO2);CHKERRQ(ierr);
  }
  {
    const PetscInt ind = SPIDER_VOLATILE_H2O;
    ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],1,&ind,&x_H2O);CHKERRQ(ierr);
  }

  ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol_in,E->numFields,NULL,subVecs);CHKERRQ(ierr);

  /* check CO2 */
  max_CO2 = P->atmosphere_parameters.volatile_parameters[SPIDER_VOLATILE_CO2].poststep_change;
  rel_CO2 = PetscAbsReal( (x_CO2-data->x_CO2) / data->x_CO2);
  if ( rel_CO2 > max_CO2 ){
    //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, rolling back solution and setting flags\n",__func__,__FILE__,__LINE__);
    ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"CO2 change exceeded (= %f, max= %f)\n",rel_CO2,max_CO2);
    EVENT = PETSC_TRUE;
  }

  /* check H2O */
  max_H2O = P->atmosphere_parameters.volatile_parameters[SPIDER_VOLATILE_H2O].poststep_change;
  rel_H2O = PetscAbsReal( (x_H2O-data->x_H2O) / data->x_H2O);
  if ( rel_H2O > max_H2O ){
    //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, rolling back solution and setting flags\n",__func__,__FILE__,__LINE__);
    ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"H2O change exceeded (= %f, max= %f)\n",rel_H2O,max_H2O);
    EVENT = PETSC_TRUE;
  }

  if (EVENT){
    ierr = TSRollBack(ts);CHKERRQ(ierr);
    ierr = TSSetConvergedReason(ts,TS_CONVERGED_USER);CHKERRQ(ierr);
    E->stopEarly = PETSC_TRUE;
  }

  PetscFree(subVecs);CHKERRQ(ierr);

  PetscFunctionReturn(0);

}
