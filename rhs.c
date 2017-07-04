#include "ctx.h"
#include "aug.h"
#include "bc.h"
#include "energy.h"
#include "matprop.h"
#include "twophase.h"
#include "util.h"
#include "atmosphere.h"

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec dSdr_b_aug_in,Vec rhs_b_aug,void *ptr)
{
  PetscErrorCode    ierr;
  Ctx               *E = (Ctx*) ptr;
  Mesh              *M = &E->mesh;
  Solution          *S = &E->solution;
  PetscScalar       *arr_rhs_b, *arr_radius_s, *arr_radius_b;
  const PetscScalar *arr_Etot, *arr_lhs_s, *arr_temp_s, *arr_Htot_s;
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
  /* this is initialised with zeros (potential for bug?) */
  ierr = VecDuplicate(S->dSdr,&rhs_b);CHKERRQ(ierr);

  /* Get first staggered node value (stored as S0) */
  ind = 0;
  ierr = VecGetValues(dSdr_b_aug_in,1,&ind,&S0);CHKERRQ(ierr);

  set_entropy( E, S0 ); CHKERRQ(ierr);

  set_gphi_smooth( E );

  set_capacitance_staggered( E );

  set_matprop_basic( E );

  set_Etot( E );

  /* note pass in current time in years here */
  set_Htot( E, t );

  /* note pass in current time in years here */
  set_surface_flux( E, t );

  set_core_mantle_flux( E );

  /* loop over basic nodes except last node */
  ierr = DMDAVecGetArray(da_b,rhs_b,&arr_rhs_b);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(E->da_b,S->Etot,INSERT_VALUES,E->work_local_b);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(E->da_b,S->Etot,INSERT_VALUES,E->work_local_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_b,E->work_local_b,&arr_Etot);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->Htot_s,&arr_Htot_s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->lhs_s,&arr_lhs_s);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->temp_s,&arr_temp_s); CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);

  /* note plus one for start of loop */
  for(i=ilo_b+1; i<ihi_b; ++i){
    /* fluxes */
    arr_rhs_b[i] = arr_Etot[i+1] * ( 1.0 / arr_lhs_s[i] );
    arr_rhs_b[i] += arr_Etot[i] * ( -1.0 / arr_lhs_s[i] - 1.0 / arr_lhs_s[i-1] );
    arr_rhs_b[i] += arr_Etot[i-1] * ( 1.0 / arr_lhs_s[i-1] );
    /* internal heat generation */
    arr_rhs_b[i] += arr_Htot_s[i] / arr_temp_s[i];
    arr_rhs_b[i] -= arr_Htot_s[i] / arr_temp_s[i-1];
    /* d/dr term */
    arr_rhs_b[i] /= arr_radius_s[i] - arr_radius_s[i-1]; // note this is negative
  }

  /* TODO: this will break in parallel*/
  /* fluxes */
  dS0dt = ( arr_Etot[1] - arr_Etot[0] ) / arr_lhs_s[0];
  /* internal heat generation */
  dS0dt += arr_Htot_s[0] / arr_temp_s[0];

  ierr = DMDAVecRestoreArrayRead(da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da_b,rhs_b,&arr_rhs_b);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_b,E->work_local_b,&arr_Etot);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->Htot_s,&arr_Htot_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->lhs_s,&arr_lhs_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);

  /* Transfer back  */
  ierr = ToAug(rhs_b,rhs_b_aug);CHKERRQ(ierr);

  /* Set zero in first position */
  /* TODO: I think this breaks for parallel */
  ierr = VecSetValue(rhs_b_aug,0,dS0dt,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(rhs_b_aug);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs_b_aug);CHKERRQ(ierr);

  /* convert time from per second to per year */
  ierr = VecScale( rhs_b_aug, SECSINYR );

  ierr = VecDestroy(&rhs_b);CHKERRQ(ierr);

  /* DJB testing atmosphere here */
  atmosphere_test( E );

  /* end of testing atmosphere */

  PetscFunctionReturn(0);
}
