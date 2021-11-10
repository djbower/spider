#include "poststep.h"
#include "ctx.h"
#include "util.h"

PetscErrorCode PostStepDataInitialize(Ctx *E)
{
  PetscErrorCode             ierr;
  Atmosphere                 *A = &E->atmosphere;
  Parameters const           P = E->parameters;
  AtmosphereParameters const Ap = P->atmosphere_parameters;
  PostStepData               *data;
  PetscInt                   i;

  PetscFunctionBeginUser;
  ierr = PetscMalloc(sizeof(PostStepData),&(E->postStepData));CHKERRQ(ierr);
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)X),"PLACEHOLDER I'm %s in %s:%d, Populating ctx->postStepData\n",__func__,__FILE__,__LINE__);
  data = (PostStepData*) E->postStepData;

  /* PostStepDataInitialize is called after TSCustomMonitor, which means the ctx has
     been updated to be consistent with the solution.  Therefore we can pull data
     from the ctx */

  /* store volatile partial pressures to the poststep data struct */
  for( i=0; i<Ap->n_volatiles; ++i) {
     data->volatile_partialp[i] = A->volatiles[i].p;
  }

  /* store surface temperature to the poststep data struct */
  data->tsurf = A->tsurf;

  PetscFunctionReturn(0);

}

PetscErrorCode PostStep(TS ts)
{
  PetscErrorCode       ierr;
  PetscInt             stepNumber;
  Ctx                  *E;
  PostStepData         *data;
  Atmosphere           *A;
  Parameters           P;
  AtmosphereParameters Ap;
  Vec                  sol_in;
  PetscScalar          maxx, relx;
  PetscInt             i;
  PetscBool            EVENT = PETSC_FALSE;

  PetscFunctionBeginUser;
  ierr = TSGetStepNumber(ts,&stepNumber);CHKERRQ(ierr);
  ierr = TSGetApplicationContext(ts,&E);CHKERRQ(ierr);
  A = &E->atmosphere;
  P = E->parameters;
  Ap = P->atmosphere_parameters;
  if (!E->parameters->rollBackActive) SETERRQ(PetscObjectComm((PetscObject)ts),PETSC_ERR_SUP,"You must run with -activate_rollback");
  data = (PostStepData*) E->postStepData;
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, x_CO2 is %g\n",__func__,__FILE__,__LINE__,data->x_CO2);
  //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, x_H2O is %g\n",__func__,__FILE__,__LINE__,data->x_H2O);

  ierr = TSGetSolution(ts,&sol_in);CHKERRQ(ierr);

  for( i=0; i<Ap->n_volatiles; ++i) {
     Volatile *V = &A->volatiles[i];
     maxx = P->atmosphere_parameters->volatile_parameters[i]->poststep_change;
     if( maxx > 0 ){
         /* remember that the - is for minus! */
         relx = PetscAbsReal( (V->p-data->volatile_partialp[i]) / data->volatile_partialp[i] );
         if ( relx > maxx ){
             //ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"PLACEHOLDER I'm %s in %s:%d, rolling back solution and setting flags\n",__fun    c__,__FILE__,__LINE__);
             ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"volatile %d partial pressure change exceeded (= %f, max= %f)\n",i,relx,maxx);
             EVENT = PETSC_TRUE;
         }
     }
  }

  /* now test tsurf */
  /* next is just a switch to decide whether to test tsurf or not */
  if( Ap->tsurf_poststep_change > 0 ){
      maxx = data->tsurf - Ap->tsurf_poststep_change;
      if ( A->tsurf < maxx ) {
          ierr = PetscPrintf(PetscObjectComm((PetscObject)ts),"tsurf change exceeded (max= %f)\n", maxx);
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
