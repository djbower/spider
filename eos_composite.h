#if !defined(EOS_COMPOSITE_H_)
#define EOS_COMPOSITE_H_

#include "eos.h"

typedef struct {
    EOS *eos_parameters;
    PetscInt n_eos;
    /* it's only for composite structures that we need to know
       if/how to blend together the material parameters across
       phase boundaries.  Most are blended using a simple linear
       weighting and a tanh functon, but the viscosity smoothing
       is a bit more complicated.  To reiterate, these quantities
       are not relevant for single-phase EOS implementations */
    /* these are all smoothing-related */
    PetscScalar matprop_smooth_width; // numerical reasons only
    PetscScalar phi_critical; // physical transition between melt and solid for viscosity
    PetscScalar phi_width; // physical transition between melt and solid for viscosity

    char phase_boundary_filename[PETSC_MAX_PATH_LEN]; // filename of phase boundary
    /* in generality, the phase boundary should be a function pointer as well.  Currently,
       it is always a 1D lookup, but could be any function */
    Interp1d phase_boundary; /* pressure-entropy space, J/kg/K */
} data_EOSComposite;

// TODO existing EosParameters-based functions, to refactor into proper class hierarchy

PetscErrorCode SetPhaseBoundary( const EosParameters, PetscScalar, PetscScalar *, PetscScalar * );

PetscErrorCode SetEosCompositeEval( const EosComposite, PetscScalar, PetscScalar, EosEval * );
PetscErrorCode SetTwoPhasePhaseFractionNoTruncation( const EosComposite eos_composite, PetscScalar, PetscScalar, PetscScalar * );

PetscErrorCode EosCompositeCreateTwoPhase( EosComposite *, const EosParameters[], PetscInt );
PetscErrorCode EosCompositeDestroy( EosComposite * );

#endif
