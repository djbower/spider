#include "ctx.h"
#include "aug.h"
#include "bc.h"
#include "energy.h"
#include "matprop.h"
#include "twophase.h"
#include "util.h"

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec dSdr_b_aug_in,Vec rhs_b_aug,void *ptr)
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
  PetscInt             ind;
  PetscScalar          S0, dS0dt;
  PetscScalar          x0, x1;

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
     extra points */
  ierr = FromAug(dSdr_b_aug_in,S->dSdr);CHKERRQ(ierr);

  /* Create rhs vector of "normal" size (no extra point) */
  /* this is initialised with zeros (potential for bug?) */
  ierr = VecDuplicate(S->dSdr,&rhs_b);CHKERRQ(ierr);

  /* extract other necessary quantities from augmented array */
  /* FIXME: clean up */
  if( Ap->MODEL == 3 ){ //MO_ATMOSPHERE_TYPE_VOLATILES){
    /* C02 content of magma ocean (solid and liquid phase) */
    ind = 0;
    ierr = VecGetValues(dSdr_b_aug_in,1,&ind,&x0);CHKERRQ(ierr);

    /* H20 content of magma ocean (solid and liquid phase) */
    ind = 1;
    ierr = VecGetValues(dSdr_b_aug_in,1,&ind,&x1);CHKERRQ(ierr);
  }

  /* Get first staggered node value (stored as S0) */
  ind = 2;
  ierr = VecGetValues(dSdr_b_aug_in,1,&ind,&S0);CHKERRQ(ierr);

  ierr = set_entropy( E, S0 );CHKERRQ(ierr);

  ierr = set_gphi_smooth( E );CHKERRQ(ierr);

  ierr = set_capacitance_staggered( E );CHKERRQ(ierr);

  ierr = set_matprop_basic( E );CHKERRQ(ierr);

  /* updates A->Mliq, A->Msol, A->dMliqdt */
  ierr = set_Mliq( E );CHKERRQ(ierr);
  ierr = set_Msol( E );CHKERRQ(ierr);
  ierr = set_dMliqdt( E );CHKERRQ(ierr);

  /* FIXME: add switch depending on atmosphere model and
     passive atmosphere tracking */
  ierr = set_atmosphere_volatile_content( E, x0, x1 );

  ierr = set_Etot( E );CHKERRQ(ierr);

  /* FIXME: make sure time is correct in years depending on
     non dimensionalisation scheme */
  /* note pass in current time in years here */
  ierr = set_Htot( E, t );CHKERRQ(ierr);

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
  ierr = ToAug(rhs_b,rhs_b_aug);CHKERRQ(ierr);

  /* time-dependence of additional quantities at the top of the augmented array */

  /* now that we have dS/dt, compute change in volatile concentrations */
  if (Ap->MODEL == 3){//MO_ATMOSPHERE_TYPE_VOLATILES){
    PetscScalar dx0dt, dx1dt;
    dx0dt = get_dx0dt( E, x0 );
    dx1dt = get_dx1dt( E, x1 );
    ierr = VecSetValue(rhs_b_aug,0,dx0dt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(rhs_b_aug,1,dx1dt,INSERT_VALUES);CHKERRQ(ierr);
  }
  else{
    ierr = VecSetValue(rhs_b_aug,0,0.0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(rhs_b_aug,1,0.0,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* S0 */
  /* TODO: I think this breaks for parallel */
  ierr = VecSetValue(rhs_b_aug,2,dS0dt,INSERT_VALUES);CHKERRQ(ierr);

  /* VecAssembly */
  ierr = VecAssemblyBegin(rhs_b_aug);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs_b_aug);CHKERRQ(ierr);

  /* FIXME: this must blow up noise so get rid of this */
  /* convert time from per second to per year */
  //ierr = VecScale( rhs_b_aug, SECSINYR );

  ierr = VecDestroy(&rhs_b);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
