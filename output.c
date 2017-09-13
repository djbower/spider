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

    ierr = VecSetValue(vec,0,A->VOLSCALE,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,1,A->H2O_INITIAL,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,2,A->CO2_INITIAL,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,3,A->EMISSIVITY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,4,A->SIGMA,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,5,A->TEQM,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,6,A->P0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,7,A->H2O_KDIST,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,8,A->H2O_KABS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,9,A->H2O_HENRY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,10,A->H2O_HENRY_POW,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,11,A->CO2_KDIST,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,12,A->CO2_KABS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,13,A->CO2_HENRY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,14,A->CO2_HENRY_POW,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,15,A->RADIUS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,16,A->GRAVITY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,17,A->M0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,18,A->Mliq,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,19,A->Msol,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,20,A->dMliqdt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,21,A->tau,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,22,A->emissivity,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,23,A->x0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,24,A->dx0dt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,25,A->p0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,26,A->dp0dx,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,27,A->m0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,28,A->tau0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,29,A->x1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,30,A->dx1dt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,31,A->p1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,32,A->dp1dx,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,33,A->m1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,34,A->tau1,INSERT_VALUES);CHKERRQ(ierr);

    ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
