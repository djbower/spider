#if !defined(EOS_COMPOSITE_H_)
#define EOS_COMPOSITE_H_

#include "eos.h"

typedef struct {
    EOS *eos;
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

} data_EOSComposite;

PetscErrorCode EOSCompositeSetSubEOS(EOS,EOS*,PetscInt);

// TODO  --- remove below once EOS class is finished -----


PetscErrorCode SetEosCompositeEval( const EosComposite, PetscScalar, PetscScalar, EosEval * );
PetscErrorCode SetTwoPhasePhaseFractionNoTruncation( const EosComposite eos_composite, PetscScalar, PetscScalar, PetscScalar * );

PetscErrorCode EosCompositeCreateTwoPhase( EosComposite *, const EosParameters[], PetscInt );
PetscErrorCode EosCompositeDestroy( EosComposite * );

#endif
