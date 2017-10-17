#include "output.h"

static PetscErrorCode add_vector_to_viewer( Vec, PetscViewer );

PetscErrorCode SetScalingsForOutput(Ctx *E)
{
    /* scalings for all output quantities */

    PetscFunctionBeginUser;

    Constants const *C = &E->parameters.constants;

    PetscScalar *meshScalings_b = E->mesh.meshScalings_b;
    PetscScalar *meshScalings_s = E->mesh.meshScalings_s;
    PetscScalar *solutionScalings_b = E->solution.solutionScalings_b;
    PetscScalar *solutionScalings_s = E->solution.solutionScalings_s;

    /* must ensure these entries correspond positionally to the order
       of the vectors as given in set_mesh */

    /* basic vecs */
    meshScalings_b[0] = C->AREA * 4.0 * PETSC_PI;
    meshScalings_b[1] = C->DPDR;
    meshScalings_b[2] = C->PRESSURE;
    meshScalings_b[3] = C->RADIUS;
    meshScalings_b[4] = C->RADIUS;

    /* staggered vecs */
    meshScalings_s[0] = C->PRESSURE;
    meshScalings_s[1] = C->RADIUS;
    meshScalings_s[2] = C->VOLUME * 4.0 * PETSC_PI;
    meshScalings_s[3] = C->DPDR;
    meshScalings_s[4] = C->AREA * 4.0 * PETSC_PI;
    meshScalings_s[5] = C->DENSITY;
    meshScalings_s[6] = C->MASS * 4.0 * PETSC_PI;

    /* must ensure these entries correspond positionally to the order
       of the vectors as given in SetupCtx in ctx.c */

    /* basic vecs */
    solutionScalings_b[0] = 1.0 / C->TEMP;
    solutionScalings_b[1] = 1.0 / C->TEMP;
    solutionScalings_b[2] = C->COND;
    solutionScalings_b[3] = C->ENTROPY;
    solutionScalings_b[4] = C->ENTROPY;
    solutionScalings_b[5] = C->DSDR;
    solutionScalings_b[6] = C->DTDR;
    solutionScalings_b[7] = C->DSDR;
    solutionScalings_b[8] = C->DSDR;
    solutionScalings_b[9] = C->DSDR;
    solutionScalings_b[10] = C->DTDR;
    solutionScalings_b[11] = C->DTDR;
    solutionScalings_b[12] = C->POWER * 4.0 * PETSC_PI; // total energy flow over spherical surface
    solutionScalings_b[13] = C->ENTROPY;
    solutionScalings_b[14] = C->ENTROPY;
    solutionScalings_b[15] = C->TEMP;
    solutionScalings_b[16] = C->DENSITY;
    solutionScalings_b[17] = C->TEMP;
    solutionScalings_b[18] = 1.0; // weight is non-dimensional
    solutionScalings_b[19] = 1.0; // weight is non-dimensional
    solutionScalings_b[20] = 1.0; // (generalised) melt fraction is non-dimensional
    solutionScalings_b[21] = C->GSUPER;
    solutionScalings_b[22] = C->FLUX;
    solutionScalings_b[23] = C->FLUX;
    solutionScalings_b[24] = C->FLUX;
    solutionScalings_b[25] = C->FLUX;
    solutionScalings_b[26] = C->FLUX;
    solutionScalings_b[27] = C->KAPPA;
    solutionScalings_b[28] = C->ENTROPY;
    solutionScalings_b[29] = C->DENSITY;
    solutionScalings_b[30] = C->TEMP;
    solutionScalings_b[31] = C->KAPPA;
    solutionScalings_b[32] = 1.0; // melt fraction is non-dimensional
    solutionScalings_b[33] = C->DENSITY;
    solutionScalings_b[34] = C->ENTROPY;
    solutionScalings_b[35] = C->ENTROPY;
    solutionScalings_b[36] = C->DENSITY;
    solutionScalings_b[37] = C->TEMP;
    solutionScalings_b[38] = C->TEMP;
    solutionScalings_b[39] = C->VISC;

    /* staggered vecs */
    solutionScalings_s[0] = C->ENTROPY;
    solutionScalings_s[1] = C->ENTROPY;
    solutionScalings_s[2] = C->RHS; // note: C->RHS is 1.0 (see parameters.c)
    solutionScalings_s[3] = C->ENTROPY;
    solutionScalings_s[4] = C->ENTROPY;
    solutionScalings_s[5] = C->TEMP;
    solutionScalings_s[6] = C->TEMP;
    solutionScalings_s[7] = 1.0; // non-dimensional weight
    solutionScalings_s[8] = 1.0; // non-dimenisonal weight
    solutionScalings_s[9] = 1.0; // generalised melt fraction is non-dimensional
    solutionScalings_s[10] = C->SENERGY / C->TIME; // W/kg
    solutionScalings_s[11] = C->SENERGY / C->TIME; // W/kg
    solutionScalings_s[12] = C->SENERGY / C->TIME; // W/kg
    solutionScalings_s[13] = C->LHS * 4.0 * PETSC_PI;
    solutionScalings_s[14] = C->DENSITY;
    solutionScalings_s[15] = C->ENTROPY;
    solutionScalings_s[16] = C->TEMP;
    solutionScalings_s[17] = 1.0; // melt fraction is non-dimensional
    solutionScalings_s[18] = C->DENSITY;
    solutionScalings_s[19] = C->ENTROPY;
    solutionScalings_s[20] = C->ENTROPY;
    solutionScalings_s[21] = C->DENSITY;
    solutionScalings_s[22] = C->TEMP;
    solutionScalings_s[23] = C->TEMP;

    PetscFunctionReturn(0);

}

PetscErrorCode scale_vectors_and_output( Vec *vec, PetscScalar *arr_scale, PetscInt NUM, PetscViewer viewer)
{

  PetscErrorCode ierr;
  Vec vec_scaled;
  PetscScalar scale;
  PetscInt i;

  PetscFunctionBeginUser;

  for (i=0;i<NUM;++i){
    ierr = VecDuplicate( vec[i], &vec_scaled); CHKERRQ(ierr);
    ierr = VecCopy( vec[i], vec_scaled ); CHKERRQ(ierr);
    scale = arr_scale[i];
    ierr = VecScale( vec_scaled, scale );
    ierr = add_vector_to_viewer( vec_scaled, viewer ); CHKERRQ(ierr);
    ierr = VecDestroy(&vec_scaled); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);

}

static PetscErrorCode add_vector_to_viewer( Vec vec, PetscViewer viewer)
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
    ierr = VecSetValue(vec,27,C->VOLATILE,INSERT_VALUES);CHKERRQ(ierr);

    ierr = VecAssemblyBegin(vec);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(vec);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
