#ifndef EOS_H_
#define EOS_H_

#include "parameters.h" // TODO remove once new EOS class in place
#include "constants.h"
#include "interp.h"

/* A (temporary) struct that is used to hold  EOS properties at a given V,T or P, S.  */
typedef struct EosEval_ {
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
  PetscScalar fusion; // only relevant for EosComposite (with two phases)
} EosEval; // TODO capitalize to EOSEval

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
  PetscErrorCode (*eval)(struct data_EOS_*, PetscScalar, PetscScalar, EosEval*);
  PetscErrorCode (*destroy)(struct data_EOS_*);
  PetscErrorCode (*setupfromoptions)(struct data_EOS_*, const char*, const FundamentalConstants, const ScalingConstants);
} data_EOS;
typedef data_EOS *EOS;

PetscErrorCode EOSCreate(EOS*, EOSType);
PetscErrorCode EOSDestroy(EOS*);
PetscErrorCode EOSEval(EOS, PetscScalar, PetscScalar, EosEval*);
PetscErrorCode EOSSetUpFromOptions(EOS, const char*, const FundamentalConstants, const ScalingConstants);
PetscErrorCode EOSGetPhaseBoundary(EOS,PetscScalar, PetscScalar*, PetscScalar*);
PetscErrorCode EOSEvalSetViscosity(EOS,EosEval*);

/* Note that this is the only place in this header that anything
   related to specific types is mentioned. These must correspond to all available
   values for EOSType, giving a string and a Creation function. */
#define SPIDER_EOS_LOOKUP "lookup"
#define SPIDER_EOS_RTPRESS "rtpress"
#define SPIDER_EOS_COMPOSITE "composite"

PetscErrorCode EOSCreate_Lookup(EOS);
PetscErrorCode EOSCreate_RTpress(EOS);
PetscErrorCode EOSCreate_Composite(EOS);

// TODO below here, remove once new classes are tested ---------------------------------
typedef struct RTpressEval_ {
  PetscScalar P; /* pressure */
  PetscScalar S; /* entropy */
  RTpressParameters rtp;
} RTpressEval;

PetscErrorCode EosParametersCreate( EosParameters * );
PetscErrorCode EosParametersDestroy( EosParameters * );
PetscErrorCode EosParametersSetFromOptions( EosParameters, const FundamentalConstants, const ScalingConstants, PetscInt );
PetscErrorCode SetEosEval( const EosParameters, PetscScalar, PetscScalar, EosEval * );

PetscErrorCode SetPhaseBoundary( const EosParameters, PetscScalar, PetscScalar *, PetscScalar * );
PetscErrorCode SetEosEvalViscosity( const EosParameters, EosEval * );

#endif
