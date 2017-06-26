#include "ctx.h"
#include "aug.h"
#include "bc.h"
#include "energy.h"
#include "matprop.h"
#include "twophase.h"

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec dSdr_b_aug_in,Vec rhs_b_aug,void *ptr)
{
  PetscErrorCode ierr;
  Ctx               *E = (Ctx*) ptr;
  Mesh              *M = &E->mesh;
  Solution          *S = &E->solution;
  PetscScalar       *arr_rhs_b, *arr_radius_s, *arr_radius_b;
  PetscScalar       *arr_S_b, *arr_S_s, *arr_dSdr_b;
  const PetscScalar *arr_Etot,*arr_lhs_s;
  PetscMPIInt       rank,size;
  PetscInt          i,ihi_b,ilo_b,w_b,numpts_b;
  DM                da_s = E->da_s, da_b=E->da_b;
  Vec               rhs_b;
  PetscInt          ind;
  PetscScalar       S0,dS0dt;

  PetscFunctionBeginUser;
#if (defined VERBOSE)
  ierr = PetscPrintf(PETSC_COMM_WORLD,"rhs:\n");CHKERRQ(ierr);
#endif
  ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  /* for looping over basic nodes */
  ierr = DMDAGetCorners(da_s,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
  ihi_b = ilo_b + w_b;

  /* Transfer from the input vector to "S->dSdr", which is the same, minus the 
     extra point */
  /* TODO: check notation.  I think memory address is right for first line,
     and pointer to vec for second? */
  ierr = FromAug(dSdr_b_aug_in,S->dSdr);CHKERRQ(ierr);

  /* Create rhs vector of "normal" size (no extra point) */
  ierr = VecDuplicate(S->dSdr,&rhs_b);CHKERRQ(ierr);

  /* Get first staggered node value (stored as S0) */
  ind = 0;
  ierr = VecGetValues(dSdr_b_aug_in,1,&ind,&S0);CHKERRQ(ierr);

  /* next section involves integrating to get the S profile, and this
     is coded entirely for a serial run */
  /* TODO: is this going to slow down the code if we run this check everytime? */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
  if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");

  /* integrate to get S profile */
  ierr = DMDAVecGetArrayRead(da_b,S->dSdr,&arr_dSdr_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da_b,S->S,&arr_S_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);

  /* S (absolute) at staggered and basic nodes */
  arr_S_s[0] = 0.0;
  /* note plus one for start of loop */
  for(i=ilo_b+1; i<ihi_b; ++i){
    /* S (absolute) at staggered nodes */
    arr_S_s[i] = arr_dSdr_b[i] * (arr_radius_s[i] - arr_radius_s[i-1] );
    arr_S_s[i] += arr_S_s[i-1]; // dS relative to first staggered value
    arr_S_b[i] = arr_dSdr_b[i] * 0.5 * (arr_radius_b[i] - arr_radius_b[i-1] );
    arr_S_b[i] += arr_S_s[i-1];
    arr_S_s[i-1] += S0; // add large constant at end to try and retain precision
    arr_S_b[i-1] += S0; // add large constant at end to try and retain precision
  }

  /* perhaps not required because first and last points are determined by
     boundary conditions, but by setting to the value of the first staggered
     node we might avoid falling off the end of the lookup and hence avoid
     some potential bugs.  This is also the same as how the python code
     currently operates */
  arr_S_b[ihi_b] = arr_dSdr_b[ihi_b] * 0.5 * (arr_radius_b[ihi_b] - arr_radius_b[ihi_b-1] );
  arr_S_b[ihi_b] += arr_S_s[ihi_b-1];
  arr_S_s[ihi_b-1] += S0; // add large constant at end to try and retain precision
  arr_S_b[ihi_b-1] += S0; // add large constant at end to try and retain precision
  arr_S_b[ihi_b] += S0; // add large constant at end to try and retain precision

  ierr = DMDAVecRestoreArrayRead(da_b,S->dSdr,&arr_dSdr_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da_b,S->S,&arr_S_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);

  set_gphi_smooth( E );

  /* loop over staggered nodes and populate E struct */
  set_capacitance_staggered( E );

  /* loop over basic (internal) nodes and populate E struct */
  set_matprop_basic( E );

  //set_matprop_and_flux( E );

  set_Etot( E );

  set_surface_flux( E );

  set_core_mantle_flux( E );

  /* loop over basic nodes except last node */
  ierr = DMDAVecGetArray(da_b,rhs_b,&arr_rhs_b);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(E->da_b,S->Etot,INSERT_VALUES,E->work_local_b);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(E->da_b,S->Etot,INSERT_VALUES,E->work_local_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_b,E->work_local_b,&arr_Etot);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->lhs_s,&arr_lhs_s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);

  /* note plus one for start of loop */
  for(i=ilo_b+1; i<ihi_b; ++i){
    arr_rhs_b[i] = arr_Etot[i+1] * ( 1.0 / arr_lhs_s[i] );
    arr_rhs_b[i] += arr_Etot[i] * ( -1.0 / arr_lhs_s[i] - 1.0 / arr_lhs_s[i-1] );
    arr_rhs_b[i] += arr_Etot[i-1] * ( 1.0 / arr_lhs_s[i-1] );
    arr_rhs_b[i] /= arr_radius_s[i] - arr_radius_s[i-1];
  }

  /* TODO: this will break in parallel*/
  dS0dt = ( arr_Etot[1] - arr_Etot[0] ) / arr_lhs_s[0];

  ierr = DMDAVecRestoreArrayRead(da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da_b,rhs_b,&arr_rhs_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_b,E->work_local_b,&arr_Etot);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->lhs_s,&arr_lhs_s);CHKERRQ(ierr);

  /* Transfer back  */
  ierr = ToAug(rhs_b,rhs_b_aug);CHKERRQ(ierr);

  /* Set zero in first position */
  /* TODO: I think this breaks for parallel */
  ierr = VecSetValue(rhs_b_aug,0,dS0dt,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(rhs_b_aug);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs_b_aug);CHKERRQ(ierr);

  ierr = VecDestroy(&rhs_b);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
