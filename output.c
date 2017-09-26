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

PetscErrorCode atmosphere_structs_to_vec( Ctx *E, Vec vec )
{

    PetscErrorCode ierr;

    Atmosphere const *A = &E->atmosphere;
    Parameters const *P = &E->parameters;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    VolatileParameters const *CO2 = &Ap->CO2_volatile_parameters;
    VolatileParameters const *H2O = &Ap->H2O_volatile_parameters;
    Mesh const *M = &E->mesh;

    PetscFunctionBeginUser;

    ierr = VecSetValue(vec,0,Ap->MODEL,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,1,Ap->HYBRID,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,2,Ap->emissivity0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,3,Ap->sigma,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,4,Ap->teqm,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,5,Ap->PARAM_UTBL,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,6,Ap->param_utbl_const,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,7,Ap->SOLVE_FOR_VOLATILES,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,8,Ap->P0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,9,CO2->initial,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,10,CO2->kdist,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,11,CO2->kabs,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,12,CO2->henry,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,13,CO2->henry_pow,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,14,H2O->initial,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,15,H2O->kdist,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,16,H2O->kabs,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,17,H2O->henry,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,18,H2O->henry_pow,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,19,M->mantle_mass,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,20,A->Mliq,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,21,A->Msol,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,22,A->dMliqdt,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,23,A->tsurf,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,24,A->tau,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,25,A->p0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,26,A->dp0dx,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,27,A->m0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,28,A->tau0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,29,A->p1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,30,A->dp1dx,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,31,A->m1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,32,A->tau1,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,33,A->emissivity,INSERT_VALUES);CHKERRQ(ierr);

    ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode constants_struct_to_vec( Constants const *C, Vec vec )
{

    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecSetValue(vec,0,C->RADIUS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,1,C->TEMP,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,2,C->ENTROPY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,3,C->DENSITY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,4,C->AREA,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,5,C->VOLUME,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,6,C->MASS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,7,C->TIME,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,8,C->TIMEYRS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,9,C->SENERGY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,10,C->ENERGY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,11,C->PRESSURE,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,12,C->POWER,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,13,C->FLUX,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,14,C->DPDR,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,15,C->GRAVITY,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,16,C->KAPPA,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,17,C->DTDP,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,18,C->DSDR,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,19,C->DTDR,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,20,C->GSUPER,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,21,C->VISC,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,22,C->LOG10VISC,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,23,C->COND,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,24,C->SIGMA,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,25,C->LHS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,26,C->RHS,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(vec,27,C->VOLSCALE,INSERT_VALUES);CHKERRQ(ierr);

    ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
