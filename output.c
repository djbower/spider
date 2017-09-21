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
    ierr = VecSetValue(vec,2,Ap->emissivity0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,3,Ap->sigma,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,4,Ap->teqm,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,5,Ap->PARAM_UTBL,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,6,Ap->param_utbl_const,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,7,Ap->volscale,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,8,Ap->P0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,9,Ap->H2O_volatile_parameters.initial,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,10,Ap->H2O_volatile_parameters.kdist,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,11,Ap->H2O_volatile_parameters.kabs,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,12,Ap->H2O_volatile_parameters.henry,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,13,Ap->H2O_volatile_parameters.henry_pow,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,14,Ap->CO2_volatile_parameters.initial,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,15,Ap->CO2_volatile_parameters.kdist,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,16,Ap->CO2_volatile_parameters.kabs,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,17,Ap->CO2_volatile_parameters.henry,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,18,Ap->CO2_volatile_parameters.henry_pow,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,19,Ap->RADIUS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,20,Ap->GRAVITY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,21,A->M0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,22,A->Mliq,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,23,A->Msol,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,24,A->dMliqdt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,25,A->tau,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,26,A->x0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,27,A->dx0dt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,28,A->p0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,29,A->dp0dx,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,30,A->m0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,31,A->tau0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,32,A->x1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,33,A->dx1dt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,34,A->p1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,35,A->dp1dx,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,36,A->m1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,37,A->tau1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,38,A->emissivity,INSERT_VALUES);CHKERRQ(ierr);

    ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
