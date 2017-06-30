/* Helper functions to convert back and forth between standard and
   "augmented" vectors */

#include "aug.h"

PetscErrorCode CreateAug(Vec in, Vec *out_aug){
  PetscErrorCode ierr;
  PetscInt       global_size,local_size;
  PetscMPIInt    rank;
  VecType        vectype;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,out_aug);CHKERRQ(ierr);
  ierr = VecGetSize(in,&global_size);CHKERRQ(ierr);
  ierr = VecGetType(in,&vectype);CHKERRQ(ierr);
  ierr = VecGetLocalSize(in,&local_size);CHKERRQ(ierr);
  ierr = VecSetType(*out_aug,vectype);CHKERRQ(ierr);
  if (!rank) {
    ierr = VecSetSizes(*out_aug,local_size+1,global_size+1);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Created aug vec %d %d\n",local_size+1,global_size+1);CHKERRQ(ierr);
  }else{
    ierr = VecSetSizes(*out_aug,local_size,global_size+1);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode CreateUnAug(Vec in_aug, Vec *out)
{
  PetscErrorCode ierr;
  PetscInt       global_size_aug,local_size_aug;
  PetscMPIInt    rank;
  VecType        vectype;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD,out);CHKERRQ(ierr);
  ierr = VecGetSize(in_aug,&global_size_aug);CHKERRQ(ierr);
  ierr = VecGetLocalSize(in_aug,&local_size_aug);CHKERRQ(ierr);
  ierr = VecGetType(in_aug,&vectype);CHKERRQ(ierr);
  ierr = VecSetType(*out,vectype);CHKERRQ(ierr);
  if (!rank) {
    ierr = VecSetSizes(*out,local_size_aug-1,global_size_aug-1);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"Created unaug vec %d %d\n",local_size_aug-1,global_size_aug-1);CHKERRQ(ierr);
  } else {
    ierr = VecSetSizes(*out,local_size_aug,global_size_aug-1);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ToAug(Vec in,Vec out_aug)
{
  PetscErrorCode     ierr;
  PetscInt           global_size_aug,local_size_aug;
  PetscMPIInt        rank;
  PetscScalar        *out_aug_arr;
  const PetscScalar  *in_arr;
  PetscInt           i;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = VecGetSize(out_aug,&global_size_aug);CHKERRQ(ierr);
  ierr = VecGetLocalSize(out_aug,&local_size_aug);CHKERRQ(ierr);
  ierr = VecGetArray(out_aug,&out_aug_arr);CHKERRQ(ierr);
  ierr = VecGetArrayRead(in,&in_arr);CHKERRQ(ierr);
  if(!rank){
    for(i=1;i<local_size_aug;++i) out_aug_arr[i]=in_arr[i-1];
  }else{
    for(i=0;i<local_size_aug;++i) out_aug_arr[i]=in_arr[i];
  }
  ierr = VecRestoreArray(out_aug,&out_aug_arr);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(in,&in_arr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode FromAug(Vec in_aug,Vec out)
{
  PetscErrorCode     ierr;
  PetscInt           global_size_aug,local_size_aug;
  PetscMPIInt        rank;
  const PetscScalar  *in_aug_arr;
  PetscScalar        *out_arr;
  PetscInt           i;

  PetscFunctionBeginUser;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = VecGetSize(in_aug,&global_size_aug);CHKERRQ(ierr);
  ierr = VecGetLocalSize(in_aug,&local_size_aug);CHKERRQ(ierr);
  ierr = VecGetArray(out,&out_arr);CHKERRQ(ierr);
  ierr = VecGetArrayRead(in_aug,&in_aug_arr);CHKERRQ(ierr);
  if(!rank){
    for(i=0;i<local_size_aug-1;++i) out_arr[i]=in_aug_arr[i+1];
  }else{
    for(i=0;i<local_size_aug;++i) out_arr[i]=in_aug_arr[i];
  }
  ierr = VecRestoreArray(out,&out_arr);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(in_aug,&in_aug_arr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
