#include "output.h"

static PetscErrorCode add_vector_to_viewer( Vec, PetscViewer );

// TODO PDS : this will be taken care of by ScalableField eventually
PetscErrorCode scale_vectors_and_output( DimensionalisableField *f, PetscInt NUM, PetscViewer viewer)
{
  PetscErrorCode ierr;
  Vec            vec,vec_scaled;
  PetscScalar    scale;
  PetscInt       i,numFields;

  PetscFunctionBeginUser;
  for (i=0; i<NUM; ++i) {
    ierr = DimensionalisableFieldGetGlobalVec(f[i],&vec);CHKERRQ(ierr);
    ierr = VecDuplicate( vec, &vec_scaled); CHKERRQ(ierr);
    ierr = VecCopy( vec, vec_scaled ); CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetScaling(f[i],&numFields,&scale);CHKERRQ(ierr);
    if (numFields != 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Not supported - only expecting to output single fields with this function");
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

PetscErrorCode atmosphere_structs_to_vec( Vec x_aug, Ctx *E, Vec vec )
{

    PetscErrorCode ierr;

    Atmosphere           const *A = &E->atmosphere;
    Parameters           const *P = &E->parameters;
    Constants            const *C  = &P->constants;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    VolatileParameters   const *CO2 = &Ap->CO2_volatile_parameters;
    VolatileParameters   const *H2O = &Ap->H2O_volatile_parameters;
    Mesh const *M = &E->mesh;
    PetscInt ind;
    PetscScalar Msol,Mliq;
    PetscScalar sol0,liq0,atm0,tot0,sol1,liq1,atm1,tot1;
    PetscScalar x0,x1;

    PetscScalar FAC, MASS;

    PetscFunctionBeginUser;

    /* CO2 content of magma ocean (liquid phase) */
    ind = 0;
    ierr = VecGetValues(x_aug,1,&ind,&x0);CHKERRQ(ierr);

    /* H2O content of magma ocean (liquid phase) */
    ind = 1;
    ierr = VecGetValues(x_aug,1,&ind,&x1);CHKERRQ(ierr);

    /* scalings */
    MASS = 4.0 * PETSC_PI * C->MASS; // includes 4*PI for spherical geometry
    FAC = C->VOLATILE / 1.0E6;

    Msol = A->Msol * MASS;
    Mliq = A->Mliq * MASS;

    // CO2
    sol0 = FAC * x0 * CO2->kdist * Msol; // solid
    liq0 = FAC * x0 * Mliq; // liquid
    atm0 = A->m0 * MASS; // atmosphere
    tot0 = FAC * CO2->initial * M->mantle_mass * MASS; // total

    // H2O
    sol1 = FAC * x1 * H2O->kdist * Msol; // solid
    liq1 = FAC * x1 * Mliq; // liquid
    atm1 = A->m1 * MASS; // atmosphere
    tot1 = FAC * H2O->initial * M->mantle_mass * MASS; // total

    // total liquid mass of mantle, kg
    ierr = VecSetValue(vec,0,Mliq,INSERT_VALUES);CHKERRQ(ierr);
    // total solid mass of mantle, kg
    ierr = VecSetValue(vec,1,Msol,INSERT_VALUES);CHKERRQ(ierr);
    // total mass of mantle, kg (for sanity check)
    ierr = VecSetValue(vec,2,M->mantle_mass*MASS,INSERT_VALUES);CHKERRQ(ierr);
    // surface temperature, K
    ierr = VecSetValue(vec,3,A->tsurf*C->TEMP,INSERT_VALUES);CHKERRQ(ierr);
    // optical depth, non-dimensional
    ierr = VecSetValue(vec,4,A->tau,INSERT_VALUES);CHKERRQ(ierr);
    // (effective) emissivity, non-dimensional
    ierr = VecSetValue(vec,5,A->emissivity,INSERT_VALUES);CHKERRQ(ierr);
    // CO2 related
      // volatile mass in liquid mantle, kg
    ierr = VecSetValue(vec,6,liq0,INSERT_VALUES);CHKERRQ(ierr);
      // volatile mass in solid mantle, kg
    ierr = VecSetValue(vec,7,sol0,INSERT_VALUES);CHKERRQ(ierr);
      // volatile mass in atmosphere, kg
    ierr = VecSetValue(vec,8,atm0,INSERT_VALUES);CHKERRQ(ierr);
      // volatile mass in all reservoirs (for sanity check)
    ierr = VecSetValue(vec,9,tot0,INSERT_VALUES);CHKERRQ(ierr);
      // volatile partial pressure, Pa
    ierr = VecSetValue(vec,10,A->p0*C->PRESSURE,INSERT_VALUES);CHKERRQ(ierr);
      // volatile derivative, Pa/ppm
    ierr = VecSetValue(vec,11,A->dp0dx*C->PRESSURE/C->VOLATILE,INSERT_VALUES);CHKERRQ(ierr);
      // volatile optical depth, non-dimensional
    ierr = VecSetValue(vec,12,A->tau0,INSERT_VALUES);CHKERRQ(ierr);
    // H2O related
      // volatile mass in liquid mantle, kg
    ierr = VecSetValue(vec,13,liq1,INSERT_VALUES);CHKERRQ(ierr);
      // volatile mass in solid mantle, kg
    ierr = VecSetValue(vec,14,sol1,INSERT_VALUES);CHKERRQ(ierr);
      // volatile mass in atmosphere, kg
    ierr = VecSetValue(vec,15,atm1,INSERT_VALUES);CHKERRQ(ierr);
      // volatile mass in all reservoirs (for sanity check)
    ierr = VecSetValue(vec,16,tot1,INSERT_VALUES);CHKERRQ(ierr);
      // volatile partial pressure, Pa
    ierr = VecSetValue(vec,17,A->p1*C->PRESSURE,INSERT_VALUES);CHKERRQ(ierr);
      // volatile derivative, Pa/ppm
    ierr = VecSetValue(vec,18,A->dp1dx*C->PRESSURE/C->VOLATILE,INSERT_VALUES);CHKERRQ(ierr);
      // volatile optical depth, non-dimensional
    ierr = VecSetValue(vec,19,A->tau1,INSERT_VALUES);CHKERRQ(ierr);

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
