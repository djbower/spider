#include "output.h"

PetscErrorCode add_vector_to_binary_output( Vec vec, PetscViewer viewer)
{
  /* simple wrapper to add a vector to a PetscViewer */

  PetscErrorCode ierr;
  //char vecname[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  // convenient for some output formats, but not Petsc binary
  //ierr = PetscSNPrintf(vecname,PETSC_MAX_PATH_LEN,"phi_s_%lld",step);CHKERRQ(ierr);
  //ierr = PetscObjectSetName((PetscObject)x_aug,vecname);CHKERRQ(ierr);
  ierr = VecView( vec, viewer ); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode atmosphere_struct_to_vec( Atmosphere *A, Vec vec )
{

    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecSetValue(vec,0,A->M0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,1,A->Mliq,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,2,A->Msol,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,3,A->dMliqdt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,4,A->tau,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,5,A->emissivity,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,6,A->x0init,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,7,A->x0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,8,A->dx0dt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,9,A->p0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,10,A->dp0dx,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,11,A->m0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,12,A->tau0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,13,A->x1init,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,14,A->x1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,15,A->dx1dt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,16,A->p1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,17,A->dp1dx,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,18,A->m1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,19,A->tau1,INSERT_VALUES);CHKERRQ(ierr);

    ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
