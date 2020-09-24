#if !defined(EOS_ADAMSWILLIAMSON_H_)
#define EOS_ADAMSWILLIAMSON_H_

#include "eos.h"

typedef struct {
    PetscScalar gravity;
    PetscScalar rhos;
    PetscScalar beta;
} data_EOSAdamsWilliamson;

#endif
