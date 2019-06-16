#include "poststep.h"
#include "ctx.h"
#include "util.h"

PetscErrorCode PostStepDataInitialize(Ctx *E, Vec sol_in)
{
  PetscErrorCode   ierr;
  Atmosphere       *A = &E->atmosphere;
  PostStepData     *data;
  PetscInt         i;

  PetscFunctionBeginUser;
  ierr = PetscMalloc(sizeof(PostStepData),&(E->postStepData));CHKERRQ(ierr);
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)X),"PLACEHOLDER I'm %s in %s:%d, Populating ctx->postStepData\n",__func__,__FILE__,__LINE__);
  data = (PostStepData*) E->postStepData;

  /* update the volatile abundances in the atmosphere struct */
  ierr = set_volatile_abundances_from_solution( E, sol_in );CHKERRQ(ierr);

  /* store volatile abundances to the poststep data struct */
  for( i=0; i<SPIDER_MAX_VOLATILE_SPECIES; ++i) {
     data->volatile_abundances[i] = A->volatiles[i].x;
  }

  PetscFunctionReturn(0);
}

PetscErrorCode PostStep(TS ts)
{
  PetscErrorCode ierr;
  PetscInt       stepNumber;
  Ctx            *E;
  PostStepData   *data;
  Atmosphere     *A;
  Parameters     *P;
  Vec            sol_in;
  PetscScalar    maxx, relx;
  PetscInt       i;
  PetscBool      EVENT = PETSC_FALSE;

  PetscFunctionBeginUser;
  ierr = TSGetTimeStepNumber(ts,&stepNumber);CHKERRQ(ierr);
  ierr = TSGetApplicationContext(ts,&E);CHKERRQ(ierr);
  A = &E->atmosphere;
  P = &E->parameters;
  if (!E->parameters.rollBackActive) SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_SUP,"You must run with -activate_rollback");
  data = (PostStepData*) E->postStepData;
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, x_CO2 is %g\n",__func__,__FILE__,__LINE__,data->x_CO2);
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, x_H2O is %g\n",__func__,__FILE__,__LINE__,data->x_H2O);

  ierr = TSGetSolution(ts,&sol_in);CHKERRQ(ierr);

  /* update the volatile abundances in the atmosphere struct */
  ierr = set_volatile_abundances_from_solution( E, sol_in );CHKERRQ(ierr);

  for( i=0; i<SPIDER_MAX_VOLATILE_SPECIES; ++i) {
     Volatile *V = &A->volatiles[i];
     maxx = P->atmosphere_parameters.volatile_parameters[i].poststep_change;
     relx = PetscAbsReal( (V->x-data->volatile_abundances[i]) / data->volatile_abundances[i] );
     if ( relx > maxx ){
         //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, rolling back solution and setting flags\n",__fun    c__,__FILE__,__LINE__);
         ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"volatile %d change exceeded (= %f, max= %f)\n",i,relx,maxx);
         EVENT = PETSC_TRUE;
     }
  }

  if (EVENT){
    ierr = TSRollBack(ts);CHKERRQ(ierr);
    ierr = TSSetConvergedReason(ts,TS_CONVERGED_USER);CHKERRQ(ierr);
    E->stopEarly = PETSC_TRUE;
  }

  PetscFunctionReturn(0);

}
