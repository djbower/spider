#ifndef REACTION_H_
#define REACTION_H_

#include "parameters.h"

PetscScalar get_equilibrium_constant( ReactionParameters *, PetscScalar, const Constants * );

#endif
