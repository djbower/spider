#if !defined(EOS_COMPOSITE_H_)
#define EOS_COMPOSITE_H_

#include "eos.h"

typedef struct
{
   EOS *eos;
   PetscInt n_eos;
   /* For composite structures that we need to know if/how to blend together the material
      parameters across phase boundaries.  Most are blended using a simple linear weighting and a
      tanh function, but the viscosity smoothing is a bit more complicated. */
   PetscScalar matprop_smooth_width; // numerical reasons only
   PetscScalar phi_critical;         // physical transition between melt and solid for viscosity
   PetscScalar phi_width;            // physical transition between melt and solid for viscosity
   PetscScalar entropy_of_fusion;    // entropy of fusion

   /* For brittle logic which specifies which slots to look in for special purposes */
   PetscInt melt_slot, solid_slot, liquidus_slot, solidus_slot;
} data_EOSComposite;

PetscErrorCode EOSCompositeGetMatpropSmoothWidth(EOS, PetscScalar *);
PetscErrorCode EOSCompositeGetSubEOS(EOS, EOS **, PetscInt *);
PetscErrorCode EOSCompositeGetTwoPhasePhaseFractionNoTruncation(EOS, PetscScalar, PetscScalar, PetscScalar *);
PetscErrorCode EOSCompositeSetSubEOS(EOS, EOS *, PetscInt);

#endif
