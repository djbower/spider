#ifndef REACTION_H_
#define REACTION_H_

#include "atmosphere.h"
#include "parameters.h"

PetscScalar get_log10_modified_equilibrium_constant( const ReactionParameters, PetscScalar, const Atmosphere *, const ScalingConstants );
PetscScalar get_dlog10GdT( const ReactionParameters, PetscScalar, const Atmosphere * );
PetscScalar get_reaction_quotient_products( const ReactionParameters, const Atmosphere * );
PetscScalar get_reaction_quotient_reactants( const ReactionParameters, const Atmosphere * );
PetscScalar get_reaction_quotient_products_time_derivative( const ReactionParameters, const Atmosphere *, const AtmosphereParameters );
PetscScalar get_reaction_quotient_reactants_time_derivative( const ReactionParameters, const Atmosphere *, const AtmosphereParameters );

#endif
