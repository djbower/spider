#include "ctx.h"
#include "bc.h"
#include "energy.h"
#include "matprop.h"
#include "twophase.h"
#include "util.h"
#include "composition.h"

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec sol_in,Vec rhs,void *ptr)
{
  PetscErrorCode    ierr;
  Ctx                  *E = (Ctx*) ptr;
  Parameters           *P = &E->parameters;
  AtmosphereParameters *Ap = &P->atmosphere_parameters;
  Mesh                 *M = &E->mesh;
  Solution             *S = &E->solution;
  PetscScalar          *arr_dSdt_s, *arr_rhs_b;
  const PetscScalar    *arr_Etot, *arr_lhs_s, *arr_temp_s, *arr_Htot_s, *arr_radius_s, *arr_radius_b;
  PetscMPIInt          rank,size;
  PetscInt             i,ihi_b,ilo_b,w_b,numpts_b;
  DM                   da_s = E->da_s, da_b=E->da_b;
  Vec                  rhs_b;
  PetscScalar          S0, dS0dt;
  PetscScalar          x0, x1;
  Vec                  *subVecs;

  const PetscInt       ind0 = 0;

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  /* for looping over basic nodes */
  ierr = DMDAGetCorners(da_s,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
  ihi_b = ilo_b + w_b;

  ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
  ierr = DMCompositeGetAccessArray(E->dm_sol,sol_in,E->numFields,NULL,subVecs);CHKERRQ(ierr);

  /* Transfer from the input vector to "S->dSdr", which is the same, minus the
     extra points */
  ierr = VecCopy(subVecs[0],S->dSdr);CHKERRQ(ierr);

  /* extract other necessary quantities from sol array */
  /* volatiles have already been initialised to zero in the
     initial condition, so reading them here is fine even if
     we are not explicitly using volatiles in a model run */

  /* Get first staggered node value (stored as S0) */
  ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],1,&ind0,&S0);CHKERRQ(ierr);

  /* CO2 content of magma ocean (liquid phase) */
  ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_CO2]],1,&ind0,&x0);CHKERRQ(ierr);

  /* H2O content of magma ocean (liquid phase) */
  ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_H2O]],1,&ind0,&x1);CHKERRQ(ierr);

  ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol_in,E->numFields,NULL,subVecs);CHKERRQ(ierr);

  /* Create rhs vector of "normal" size (no extra point)
     this is initialised with zeros (potential for bug?) */
  ierr = VecDuplicate(S->dSdr,&rhs_b);CHKERRQ(ierr);

  ierr = set_entropy_from_solution( E, sol_in );CHKERRQ(ierr);

  ierr = set_gphi_smooth( E );CHKERRQ(ierr);

  ierr = set_melt_fraction_staggered( E ); CHKERRQ(ierr);

  // TODO: move later, end of time step?
  ierr = set_rheological_front( E ); CHKERRQ(ierr);

  ierr = set_capacitance_staggered( E );CHKERRQ(ierr);

  ierr = set_matprop_basic( E );CHKERRQ(ierr);

  ierr = set_Etot( E );CHKERRQ(ierr);

  ierr = set_Htot( E, t );CHKERRQ(ierr);

  /* will populate A->p?, A->dp?dx, and A->m? with zeros if
     x0 and x1 and/or are zero */
  ierr = set_atmosphere_volatile_content( E, x0, x1 );CHKERRQ(ierr);

  /* boundary conditions must be after all arrays are set */
  ierr = set_surface_flux( E );CHKERRQ(ierr);

  ierr = set_core_mantle_flux( E );CHKERRQ(ierr);

  /* loop over basic nodes except last node */
  ierr = DMDAVecGetArray(da_b,rhs_b,&arr_rhs_b);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(E->da_b,S->Etot,INSERT_VALUES,E->work_local_b);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(E->da_b,S->Etot,INSERT_VALUES,E->work_local_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_b,E->work_local_b,&arr_Etot);CHKERRQ(ierr);
  ierr = DMDAVecGetArray    (da_s,S->dSdt_s,&arr_dSdt_s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->Htot_s,&arr_Htot_s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->lhs_s,&arr_lhs_s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->temp_s,&arr_temp_s); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);

  /* first staggered node */
  arr_dSdt_s[0] = ( arr_Etot[1] - arr_Etot[0] ) / arr_lhs_s[0];
  arr_dSdt_s[0] += arr_Htot_s[0] / arr_temp_s[0];
  dS0dt = arr_dSdt_s[0];

  for(i=ilo_b+1; i<ihi_b; ++i){
    /* dSdt at staggered nodes */
    /* need this quantity for coupling to atmosphere evollution */
    arr_dSdt_s[i] = ( arr_Etot[i+1] - arr_Etot[i] ) / arr_lhs_s[i];
    arr_dSdt_s[i] += arr_Htot_s[i] / arr_temp_s[i];
    /* d/dt(dS/dr) at internal basic nodes */
    arr_rhs_b[i] = arr_dSdt_s[i] - arr_dSdt_s[i-1];
    arr_rhs_b[i] /= arr_radius_s[i] - arr_radius_s[i-1]; // note dr is negative
  }

  ierr = DMDAVecRestoreArrayRead(da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da_b,rhs_b,&arr_rhs_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_b,E->work_local_b,&arr_Etot);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da_s,S->dSdt_s,&arr_dSdt_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->Htot_s,&arr_Htot_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->lhs_s,&arr_lhs_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);


  /* Transfer back  */
  ierr = DMCompositeGetAccessArray(E->dm_sol,rhs,E->numFields,NULL,subVecs);CHKERRQ(ierr);
  ierr = VecCopy(rhs_b,subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDR_B]]);CHKERRQ(ierr);

  /* these update A->Mliq, A->Msol, and A->dMliqdt */
  ierr = set_Mliq( E );CHKERRQ(ierr);
  ierr = set_Msol( E );CHKERRQ(ierr);
  ierr = set_dMliqdt( E );CHKERRQ(ierr); /* must be after dS/dt computation */

  /* S0 */
  /* TODO: I think this breaks for parallel */
  ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],0,dS0dt,INSERT_VALUES);CHKERRQ(ierr);

  /* time-dependence of additional quantities */
  if (Ap->SOLVE_FOR_VOLATILES || Ap->SURFACE_BC==3){
    PetscScalar dx0dt, dx1dt;
    dx0dt = get_dx0dt( E, x0 );
    dx1dt = get_dx1dt( E, x1 );
    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_CO2]],0,dx0dt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_H2O]],0,dx1dt,INSERT_VALUES);CHKERRQ(ierr);
  }
  else{
    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_CO2]],0,0.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_H2O]],0,0.0,INSERT_VALUES);CHKERRQ(ierr);
  }
  for (i=1; i<E->numFields; ++i) {
    ierr = VecAssemblyBegin(subVecs[i]);CHKERRQ(ierr);
  }
  for (i=1; i<E->numFields; ++i) {
    ierr = VecAssemblyEnd(subVecs[i]);CHKERRQ(ierr);
  }

  ierr = DMCompositeRestoreAccessArray(E->dm_sol,rhs,E->numFields,NULL,subVecs);CHKERRQ(ierr);
  ierr = PetscFree(subVecs);CHKERRQ(ierr);

  ierr = VecDestroy(&rhs_b);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
