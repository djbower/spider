#ifndef REACTION_H_
#define REACTION_H_

#include "atmosphere.h"
#include "parameters.h"

PetscScalar get_equilibrium_constant( const ReactionParameters *, PetscScalar, const Constants * );
PetscScalar get_equilibrium_constant_derivative( const ReactionParameters *, PetscScalar, const Constants * );
PetscScalar get_reaction_quotient_products( const ReactionParameters *, const Atmosphere * );
PetscScalar get_reaction_quotient_reactants( const ReactionParameters *, const Atmosphere * );

#endif
