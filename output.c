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

PetscErrorCode atmosphere_structs_to_vec( Atmosphere const *A, AtmosphereParameters const * Ap, Vec vec )
{

    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecSetValue(vec,0,Ap->MODEL,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,1,Ap->HYBRID,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,2,Ap->EMISSIVITY0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,3,Ap->SIGMA,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,4,Ap->TEQM,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,5,Ap->CONSTBC,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,6,Ap->VOLSCALE,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,7,Ap->P0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,8,Ap->H2O_INITIAL,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,9,Ap->H2O_KDIST,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,10,Ap->H2O_KABS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,11,Ap->H2O_HENRY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,12,Ap->H2O_HENRY_POW,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,13,Ap->CO2_INITIAL,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,14,Ap->CO2_KDIST,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,15,Ap->CO2_KABS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,16,Ap->CO2_HENRY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,17,Ap->CO2_HENRY_POW,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,18,Ap->RADIUS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,19,Ap->GRAVITY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,20,A->M0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,21,A->Mliq,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,22,A->Msol,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,23,A->dMliqdt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,24,A->tau,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,25,A->x0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,26,A->dx0dt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,27,A->p0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,28,A->dp0dx,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,29,A->m0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,30,A->tau0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,31,A->x1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,32,A->dx1dt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,33,A->p1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,34,A->dp1dx,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,35,A->m1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,36,A->tau1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,37,A->emissivity,INSERT_VALUES);CHKERRQ(ierr);

    ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
