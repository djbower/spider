#include "petsc.h"
#include "ctx.h"
#include "aug.h"

#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec S_s_aug_in,Vec rhs_s_aug,void *ptr)
{
  PetscErrorCode ierr;
  Ctx               *E = (Ctx*) ptr;
  Mesh              *M = &E->mesh;
  Solution          *S = &E->solution;
  PetscScalar       *arr_rhs_s;
  const PetscScalar *arr_Etot,*arr_lhs_s;
  PetscMPIInt       rank,size;
  PetscInt          i,ihi,ilo,w_s,numpts_b;
  DM                da_s = E->da_s, da_b=E->da_b;
  Vec               S_s_in,rhs_s;

  PetscFunctionBeginUser;
#if (defined VERBOSE)
  ierr = PetscPrintf(PETSC_COMM_WORLD,"rhs:\n");CHKERRQ(ierr);
#endif
  ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);

  /* Transfer from the input vector to "S_s_in", which is the same, minus the 
     extra point */
  ierr = CreateUnAug(S_s_aug_in,&S_s_in);CHKERRQ(ierr);
  ierr = FromAug(S_s_aug_in,S_s_in);CHKERRQ(ierr);

  /* Create rhs vector of "normal" size (no extra point) */
  ierr = VecDuplicate(S_s_in,&rhs_s);CHKERRQ(ierr);

  /* S_s_in is the solution array.  It's easiest to store this in
     the E struct for future access, and this step is done at the
     top of set_capacitance */

  /* DJB FOR CONSISTENCY WITH PYTHON CODE, AND TRANSPARENCY, SHOULD
     HERE TAKE DXDR AND INTEGRATE ALONG PROFILE TO GET S ETC.  THEN
     SAVE THESE PROFILES IN STRUCT SO SUBSEQUENT FUNCTIONS CAN READ
     THEM */

  /* loop over staggered nodes and populate E struct */
  set_capacitance( E, S_s_in );

  /* loop over basic (internal) nodes and populate E struct */
  set_matprop_and_flux( E, S_s_in );

  /* surface radiative boundary condition
     parameterised ultra-thin thermal boundary layer
     constants given by fit (python script) */

  // NOTE: here we assume that the first rank has the first point
  if (!rank) {
    PetscScalar temp_s_0,val;
    PetscInt ind = 0;
    ierr = VecGetValues(S->temp_s,1,&ind,&temp_s_0);CHKERRQ(ierr);
    val = radiative_flux_with_dT( temp_s_0 );

    // DJB set to isothermal for testing
    //ierr = VecGetValues(S->Etot,1,&ind,&temp_s_0);CHKERRQ(ierr);
    //val = temp_s_0;

    /* Note - below is true because non-dim geom is exactly equal
       to one, so do not need to multiply by area of surface */
    ierr = VecSetValue(S->Etot,0,val,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(S->Jtot,0,val,INSERT_VALUES);CHKERRQ(ierr);
  }

  // NOTE: here, we somewhat dangerously assume that the last proc has the last point 
  if (rank == size-1) {
    PetscScalar val, val2;
    PetscInt    ind, ind2;
    //PetscScalar temp_s_0;

    ind  = numpts_b-2; // penultimate basic node index
    ind2 = numpts_b-1; // last basic node index

    /* energy flux */
    ierr = VecGetValues(S->Jtot,1,&ind,&val);CHKERRQ(ierr);
    val *= E->BC_BOT_FAC;

    // DJB set to isothermal for testing
    //ierr = VecGetValues(S->Jtot,1,&ind,&temp_s_0);CHKERRQ(ierr);
    //val = temp_s_0;

    ierr = VecSetValue(S->Jtot,ind2,val,INSERT_VALUES);CHKERRQ(ierr);

    /* energy flow */
    ierr = VecGetValues(M->area_b,1,&ind2,&val2);CHKERRQ(ierr);
    val2 *= val;
    ierr = VecSetValue(S->Etot,ind2,val2,INSERT_VALUES);CHKERRQ(ierr);

  }
  ierr = VecAssemblyBegin(S->Etot);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(S->Etot);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);

  /* loop over staggered nodes except last node */
  ierr = DMDAGetCorners(da_s,&ilo,0,0,&w_s,0,0);CHKERRQ(ierr);
  ihi = ilo + w_s;
  ierr = DMDAVecGetArray(da_s,rhs_s,&arr_rhs_s);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(E->da_b,S->Etot,INSERT_VALUES,E->work_local_b);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(E->da_b,S->Etot,INSERT_VALUES,E->work_local_b);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_b,E->work_local_b,&arr_Etot);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayRead(da_s,S->lhs_s,&arr_lhs_s);CHKERRQ(ierr);
  for(i=ilo; i<ihi; ++i){
    arr_rhs_s[i] = arr_Etot[i+1] - arr_Etot[i];
    arr_rhs_s[i] /= arr_lhs_s[i];
  }
  ierr = DMDAVecRestoreArray(da_s,rhs_s,&arr_rhs_s);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_b,E->work_local_b,&arr_Etot);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayRead(da_s,S->lhs_s,&arr_lhs_s);CHKERRQ(ierr);

  /* Transfer back  */
  ierr = ToAug(rhs_s,rhs_s_aug);CHKERRQ(ierr);

  /* Set zero in first position */
  ierr = VecSetValue(rhs_s_aug,0,0.0,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(rhs_s_aug);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs_s_aug);CHKERRQ(ierr);

  ierr = VecDestroy(&S_s_in);CHKERRQ(ierr);
  ierr = VecDestroy(&rhs_s);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
