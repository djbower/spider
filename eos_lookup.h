#if !defined(EOS_LOOKUP_H_)
#define EOS_LOOKUP_H_

#include "eos.h"
#include "interp.h"
#include "parameters.h" // TODO bad (for Lookup, which should not be there)

#define SPIDER_EOS_LOOKUP_PURE "lookup"
typedef struct {
    /* lookup data filenames */
    char        rho_filename[PETSC_MAX_PATH_LEN];
    char        dTdPs_filename[PETSC_MAX_PATH_LEN];
    char        cp_filename[PETSC_MAX_PATH_LEN];
    char        temp_filename[PETSC_MAX_PATH_LEN];
    char        alpha_filename[PETSC_MAX_PATH_LEN];
    /* lookup objects to evaluate */
    /* memory is not necessarily allocated to these lookup
       objects if they are not required */
    Interp2d rho; /* density, kg/m^3 */
    Interp2d dTdPs; /* adiabatic temperature gradient, K/Pa */
    Interp2d cp; /* heat capacity, J/kg/K */
    Interp2d temp; /* temperature, K */
    Interp2d alpha; /* thermal expansion, 1/K */
} data_EOSLookupPure;
typedef data_EOSLookupPure* EOSLookupPure;

PetscErrorCode EOSLookupPureCreate(EOSLookupPure*);

// TODO these may ultimately not need to be in this header, just used in eos_lookup.c, if at all
PetscErrorCode LookupCreate( Lookup * );
PetscErrorCode LookupDestroy( Lookup * );
PetscErrorCode LookupFilenameSet( const char *, const char *, char *, PetscBool * );
PetscErrorCode SetEosEvalFromLookup( const Lookup, PetscScalar, PetscScalar, EosEval * );

#endif
