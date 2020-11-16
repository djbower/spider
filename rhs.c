#include "atmosphere.h"
#include "bc.h"
#include "ctx.h"
#include "energy.h"
#include "matprop.h"
#include "rheologicalfront.h"
#include "twophase.h"
#include "util.h"

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec sol_in,Vec rhs,void *ptr)
{
  PetscErrorCode       ierr;
  Ctx                  *E = (Ctx*) ptr;
  Parameters           P = E->parameters;
  AtmosphereParameters Ap = P->atmosphere_parameters;
  Atmosphere           *A = &E->atmosphere;
  Mesh                 *M = &E->mesh;
  Solution             *S = &E->solution;
  PetscScalar          *arr_dSdt_s, *arr_rhs_b;
  const PetscScalar    *arr_Etot, *arr_capacitance_s, *arr_temp_s, *arr_cp_s, *arr_Htot_s, *arr_xi_s, *arr_xi_b;
  PetscMPIInt          rank,size;
  DM                   da_s = E->da_s, da_b=E->da_b;
  PetscInt             i,v,ihi_s,ilo_s,w_s,numpts_b;
  Vec                  rhs_b, *subVecs;

  PetscFunctionBeginUser;
  ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  /* for looping over basic nodes */
  ierr = DMDAGetCorners(E->da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
  ihi_s = ilo_s + w_s;

  /* allocate memory for RHS vector */
  ierr = VecCreate( PETSC_COMM_WORLD, &rhs_b );
  ierr = VecSetSizes( rhs_b, PETSC_DECIDE, numpts_b );CHKERRQ(ierr);
  ierr = VecSetFromOptions( rhs_b );CHKERRQ(ierr);
  ierr = VecSetUp( rhs_b );CHKERRQ(ierr);

  /* sets everything possible (and consistently) from entropy */
  ierr = set_interior_structure_from_solution( E, t, sol_in );CHKERRQ(ierr);

  /* below also sets reaction masses in Atmosphere struct (required
     for time-stepping) */
  ierr = set_partial_pressures_from_solution( E, sol_in );CHKERRQ(ierr);

  /* boundary conditions must be after all arrays are set */
  ierr = set_surface_flux( E );CHKERRQ(ierr);

  /* legacy function for core mantle boundary sets Jtot and Etot
     at the CMB to adhere to boundary condition */
  if(1){
    ierr = set_core_mantle_flux_legacy( E );CHKERRQ(ierr);
  }

  /* now we compute the time-dependent quantities */

  /* loop over basic nodes except last node */
  ierr = DMDAVecGetArray(da_b,rhs_b,&arr_rhs_b);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(E->da_b,S->Etot,INSERT_VALUES,E->work_local_b);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(E->da_b,S->Etot,INSERT_VALUES,E->work_local_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_b,E->work_local_b,&arr_Etot);CHKERRQ(ierr);
  ierr = DMDAVecGetArray    (da_s,S->dSdt_s,&arr_dSdt_s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->Htot_s,&arr_Htot_s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->capacitance_s,&arr_capacitance_s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->temp_s,&arr_temp_s); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->cp_s,&arr_cp_s); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,M->xi_s,&arr_xi_s);CHKERRQ(ierr);

  /* first staggered node */
  arr_dSdt_s[0] = ( arr_Etot[1] - arr_Etot[0] ) / arr_capacitance_s[0];
  arr_dSdt_s[0] += arr_Htot_s[0] / arr_temp_s[0];

  for(i=ilo_s+1; i<ihi_s; ++i){
    /* dSdt at staggered nodes */
    /* need this quantity for coupling to atmosphere evolution */
    arr_dSdt_s[i] = ( arr_Etot[i+1] - arr_Etot[i] ) / arr_capacitance_s[i];
    arr_dSdt_s[i] += arr_Htot_s[i] / arr_temp_s[i];
    /* d/dt(dS/dxi) at internal basic nodes */
    arr_rhs_b[i] = arr_dSdt_s[i] - arr_dSdt_s[i-1];
    arr_rhs_b[i] /= arr_xi_s[i] - arr_xi_s[i-1]; // note dxi is negative
  }

  /* dTsurf/dr */
  /* A->dtsurfdt already contains contribution of dTsurf/dT */
  /* By chain rule, just need dT/dt */
  A->dtsurfdt *= arr_dSdt_s[0] * arr_temp_s[0] / arr_cp_s[0];
  /* TODO, add effect of gradient to above 0.5*d/dt (dS/dxi) */
  /* but this should be a minor effect */

  ierr = DMDAVecRestoreArrayRead(da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,M->xi_s,&arr_xi_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da_b,rhs_b,&arr_rhs_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_b,E->work_local_b,&arr_Etot);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da_s,S->dSdt_s,&arr_dSdt_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->Htot_s,&arr_Htot_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->capacitance_s,&arr_capacitance_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->cp_s,&arr_cp_s);CHKERRQ(ierr);

  /* must be here since must be after dS/dt computation */

  /* TODO: only relevant for 2 phases.  Perhaps need to tidy up similar
     functionality, like output of rheological front, etc. */
  if( P->n_phases == 2 ){
      ierr = set_dMliqdt( E );CHKERRQ(ierr);
  }

  /* transfer d/dt to solution */
  /* dS/dr at basic nodes */
  ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
  ierr = DMCompositeGetAccessArray(E->dm_sol,rhs,E->numFields,NULL,subVecs);CHKERRQ(ierr);
  ierr = VecCopy(rhs_b,subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDXI_B]]);CHKERRQ(ierr);

  /* S0, TODO: I think this breaks for parallel */
  ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],0,arr_dSdt_s[0],INSERT_VALUES);CHKERRQ(ierr);

  /* volatiles and reactions */

  if (Ap->n_volatiles){
    ierr = solve_dpdts( E );
    for (v=0; v<Ap->n_volatiles; ++v) {
      ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],v,A->volatiles[v].dpdt,INSERT_VALUES);CHKERRQ(ierr);
    }
    for (v=0; v<Ap->n_reactions; ++v) {
      ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]],v,A->reactions[v].dmrdt,INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  /* TODO: check, but with Ap->n_volatiles=0 there are no entries in
     the rhs vector to initialise to zero */

  for (i=1; i<E->numFields; ++i) {
    ierr = VecAssemblyBegin(subVecs[i]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(subVecs[i]);CHKERRQ(ierr);
  }

  ierr = DMCompositeRestoreAccessArray(E->dm_sol,rhs,E->numFields,NULL,subVecs);CHKERRQ(ierr);

  ierr = PetscFree(subVecs);CHKERRQ(ierr);
  ierr = VecDestroy(&rhs_b);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
