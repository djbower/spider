#ifndef EOS_H_
#define EOS_H_

#include "constants.h"
#include "interp.h"

/* A (temporary) struct that is used to hold  EOS properties at a given V,T or P, S.  */
typedef struct EOSEvalData_ {
  PetscScalar P; /* pressure */
  PetscScalar S; /* entropy */
  PetscScalar V; /* volume */
  PetscScalar T; /* temperature */
  PetscScalar Cp; /* heat capacity at constant pressure */
  PetscScalar Cv; /* heat capacity at constant volume */
  PetscScalar alpha; /* thermal expansion */
  PetscScalar rho; /* density */
  PetscScalar dTdPs; /* adiabatic temperature gradient */
  PetscScalar cond;
  PetscScalar log10visc;
  PetscScalar phase_fraction; // by definition unity for single species, but can be 0<x<1 for a composite EOS (with two phases) */
  PetscScalar porosity;
  PetscScalar fusion; // only relevant for EOSComposite (with two phases)
} EOSEvalData;

typedef const char* EOSType;

typedef struct data_EOS_ {
  EOSType type;  /* Implementation type */

  PetscScalar cond; /* thermal conductivity, W/m/K */
  PetscScalar log10visc; /* log base 10 of viscosity */
  PetscScalar activation_energy;
  PetscScalar activation_volume;
  PetscScalar activation_volume_pressure_scale;
  PetscScalar visc_ref_temp;
  PetscScalar visc_ref_pressure;
  PetscScalar visc_comp;
  PetscScalar visc_ref_comp;

  /* phase boundary which is evaluated using this EOS */
  PetscBool PHASE_BOUNDARY; /* is a phase boundary for this EOS defined? */
  char phase_boundary_filename[PETSC_MAX_PATH_LEN]; // filename of phase boundary
  /* in generality, the phase boundary should be a function pointer as well.  Currently,
     it is always a 1D lookup, but could be any function */
  Interp1d phase_boundary; /* pressure-entropy space, J/kg/K */

  /* Pointer to implementation-specific data */
  void *impl_data;

  /* Flag to signify that set up has been completed and the object can be used */
  PetscBool is_setup;

  /* Pointers to implementation-specific functions */
  // Note: no "create" pointer here, since we have a factory method (EOSCreate())
  PetscErrorCode (*eval)(struct data_EOS_*, PetscScalar, PetscScalar, EOSEvalData*);
  PetscErrorCode (*destroy)(struct data_EOS_*);
  PetscErrorCode (*setupfromoptions)(struct data_EOS_*, const char*, const FundamentalConstants, const ScalingConstants);
} data_EOS;
typedef data_EOS *EOS;

PetscErrorCode EOSCheckType(EOS,EOSType,PetscBool*);
PetscErrorCode EOSCreate(EOS*, EOSType);
PetscErrorCode EOSDestroy(EOS*);
PetscErrorCode EOSEval(EOS, PetscScalar, PetscScalar, EOSEvalData*);
PetscErrorCode EOSSetUpFromOptions(EOS, const char*, const FundamentalConstants, const ScalingConstants);
PetscErrorCode EOSGetPhaseBoundary(EOS,PetscScalar, PetscScalar*, PetscScalar*);
PetscErrorCode EOSGetType(EOS,EOSType*);
PetscErrorCode EOSEvalSetViscosity(EOS,EOSEvalData*);

/* Note that this is the only place in this header that anything
   related to specific types is mentioned. These must correspond to all available
   values for EOSType, giving a string and a Creation function. */
#define SPIDER_EOS_LOOKUP "lookup"
#define SPIDER_EOS_RTPRESS "rtpress"
#define SPIDER_EOS_COMPOSITE "composite"
#define SPIDER_EOS_ADAMSWILLIAMSON "adamswilliamson"

PetscErrorCode EOSCreate_Lookup(EOS);
// never fully tested or debugged.  Now in dev/
//PetscErrorCode EOSCreate_RTpress(EOS);
PetscErrorCode EOSCreate_Composite(EOS);
PetscErrorCode EOSCreate_AdamsWilliamson(EOS);

#endif
