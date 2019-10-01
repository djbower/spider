#include "ctx.h"
#include "bc.h"
#include "energy.h"
#include "matprop.h"
#include "rheologicalfront.h"
#include "twophase.h"
#include "util.h"
#include "atmosphere.h"
// FIXME
//#include "composition.h"

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec sol_in,Vec rhs,void *ptr)
{
  PetscErrorCode    ierr;
  Ctx                  *E = (Ctx*) ptr;
  Parameters           *P = &E->parameters;
  Constants            *C = &P->constants;
  AtmosphereParameters *Ap = &P->atmosphere_parameters;
  Atmosphere           *A = &E->atmosphere;
  Mesh                 *M = &E->mesh;
  Solution             *S = &E->solution;
  PetscScalar          *arr_dSdt_s, *arr_rhs_b;
  const PetscScalar    *arr_Etot, *arr_lhs_s, *arr_temp_s, *arr_Htot_s, *arr_radius_s, *arr_radius_b;
  PetscMPIInt          rank,size;
  PetscInt             i,v,ihi_b,ilo_b,w_b,numpts_b;
  DM                   da_s = E->da_s, da_b=E->da_b;
  Vec                  rhs_b, *subVecs;

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  /* for looping over basic nodes */
  ierr = DMDAGetCorners(da_s,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
  ihi_b = ilo_b + w_b;

  /* allocate memory for RHS vector */
  ierr = VecCreate( PETSC_COMM_WORLD, &rhs_b );
  ierr = VecSetSizes( rhs_b, PETSC_DECIDE, numpts_b );CHKERRQ(ierr);
  ierr = VecSetFromOptions( rhs_b );CHKERRQ(ierr);
  ierr = VecSetUp( rhs_b );CHKERRQ(ierr);

  /* set solution in the relevant structs */
  ierr = set_entropy_from_solution( E, sol_in );CHKERRQ(ierr);
  ierr = set_volatile_abundances_from_solution( E, sol_in );CHKERRQ(ierr);

  /* set material properties and energy fluxes and sources */
  ierr = set_gphi_smooth( E );CHKERRQ(ierr);
  ierr = set_melt_fraction_staggered( E ); CHKERRQ(ierr);
  ierr = set_capacitance_staggered( E );CHKERRQ(ierr);
  ierr = set_matprop_basic( E );CHKERRQ(ierr);
  ierr = set_Etot( E );CHKERRQ(ierr);
  ierr = set_Htot( E, t );CHKERRQ(ierr);

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

  /* these update A->Mliq, A->Msol, and A->dMliqdt */
  ierr = set_Mliq( E );CHKERRQ(ierr);
  ierr = set_Msol( E );CHKERRQ(ierr);
  ierr = set_dMliqdt( E );CHKERRQ(ierr); /* must be after dS/dt computation */

  ierr = set_rheological_front( E ); CHKERRQ(ierr);

  /* transfer d/dt to solution */
  /* dS/dr at basic nodes */
  ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
  ierr = DMCompositeGetAccessArray(E->dm_sol,rhs,E->numFields,NULL,subVecs);CHKERRQ(ierr);
  ierr = VecCopy(rhs_b,subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDR_B]]);CHKERRQ(ierr);

  /* S0, TODO: I think this breaks for parallel */
  ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],0,arr_dSdt_s[0],INSERT_VALUES);CHKERRQ(ierr);

  /* volatiles and reactions */
  /* TODO: only necessary for chemical reactions, but since this is a cheap operation, compute regardless */
  ierr = set_oxygen_fugacity( A, Ap, C );CHKERRQ(ierr);

  if (Ap->SOLVE_FOR_VOLATILES) {
    ierr = solve_dxdts( E );
    for (v=0; v<Ap->n_volatiles; ++v) {
      ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],v,A->volatiles[v].dxdt,INSERT_VALUES);CHKERRQ(ierr);
    }
    for (v=0; v<Ap->n_reactions; ++v) {
      ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]],v,A->reactions[v].dmrdt,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  else{
    for (v=0; v<Ap->n_volatiles; ++v) {
      ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],v,0.0,INSERT_VALUES);CHKERRQ(ierr);
    }
    for (v=0; v<Ap->n_reactions; ++v) {
      ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]],v,0.0,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  for (i=1; i<E->numFields; ++i) {
    ierr = VecAssemblyBegin(subVecs[i]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(subVecs[i]);CHKERRQ(ierr);
  }

  ierr = DMCompositeRestoreAccessArray(E->dm_sol,rhs,E->numFields,NULL,subVecs);CHKERRQ(ierr);

  ierr = PetscFree(subVecs);CHKERRQ(ierr);
  ierr = VecDestroy(&rhs_b);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
