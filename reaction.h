#ifndef REACTION_H_
#define REACTION_H_

#include "atmosphere.h"
#include "parameters.h"

PetscScalar get_equilibrium_constant( ReactionParameters *, PetscScalar, const Constants * );
PetscScalar get_reaction_quotient( const ReactionParameters *, const Atmosphere * );

#endif
