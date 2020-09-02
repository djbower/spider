#ifndef EOS_H_
#define EOS_H_

#include "ctx.h"
#include "parameters.h"

/* A (temporary) struct that is used to set the eos properties at a
   given V,T or P, S.  There should be as  */
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
  PetscScalar phase_fraction; // by definition unity for single species, but can be 0<x<1 for the EosComposite (with two phases) */
  PetscScalar fusion; // only relevant for EosComposite (with two phases)
} EosEval;

typedef struct RTpressEval_ {
  PetscScalar P; /* pressure */
  PetscScalar S; /* entropy */
  RTpressParameters rtp;
} RTpressEval;

PetscErrorCode EosParametersCreate( EosParameters * );
PetscErrorCode EosParametersDestroy( EosParameters * );
PetscErrorCode EosParametersSetFromOptions( EosParameters, const FundamentalConstants, const ScalingConstants );
PetscErrorCode SetPhaseBoundary( const EosParameters, PetscScalar, PetscScalar *, PetscScalar * );

PetscErrorCode SetEosEval( const EosParameters, PetscScalar, PetscScalar, EosEval * );
PetscErrorCode SetEosCompositeEval( const EosComposite, PetscScalar, PetscScalar, EosEval * );
PetscErrorCode SetTwoPhasePhaseFractionNoTruncation( const EosComposite eos_composite, PetscScalar, PetscScalar, PetscScalar * );

PetscErrorCode EosCompositeCreateTwoPhase( EosComposite *, const EosParameters[], PetscInt );
PetscErrorCode EosCompositeDestroy( EosComposite * );

PetscErrorCode JSON_add_phase_boundary( const Ctx *, const EosParameters, const char *, cJSON * );

/////////////////// WIP TODO

typedef const char* EOSType;

typedef struct data_EOS_ {
  EOSType type;  /* Implementation type */

  char prefix[128];  /* Maximum prefix length */
  // TODO do we even need this? It's only used once

  PetscScalar cond; /* thermal conductivity, W/m/K */
  PetscScalar log10visc; /* log base 10 of viscosity */
  PetscScalar activation_energy;
  PetscScalar activation_volume;
  PetscScalar activation_volume_pressure_scale;
  PetscScalar visc_ref_temp;
  PetscScalar visc_ref_pressure;
  PetscScalar visc_comp;
  PetscScalar visc_ref_comp;

  /* Pointer to implementation-specific data */
  void * impl_data;

  /* Pointers to implementation-specific functions */
  // Note: no "create" pointer here, since we have a factory method (EOSCreate())
  PetscErrorCode (*eval)(const struct data_EOS_, PetscScalar, PetscScalar, EosEval*);
  PetscErrorCode (*destroy)(struct data_EOS_*);
  PetscErrorCode (*setupfromoptions)(struct data_EOS_*, const char*);
  // TODO setup from options function

} data_EOS;
typedef data_EOS *EOS;

PetscErrorCode EOSCreate(EOS*, EOSType);
PetscErrorCode EOSDestroy(EOS*);
PetscErrorCode EOSEval(const EOS, PetscScalar, PetscScalar, EosEval*);
PetscErrorCode EOSSetUpFromOptions(EOS, const char*);

#endif
