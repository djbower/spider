#if !defined(EOS_LOOKUP_H_)
#define EOS_LOOKUP_H_

#include "eos.h"
#include "interp.h"

typedef struct
{
    /* lookup data filenames */
    char rho_filename[PETSC_MAX_PATH_LEN];
    char dTdPs_filename[PETSC_MAX_PATH_LEN];
    char cp_filename[PETSC_MAX_PATH_LEN];
    char temp_filename[PETSC_MAX_PATH_LEN];
    char alpha_filename[PETSC_MAX_PATH_LEN];
    /* lookup objects to evaluate */
    /* memory is not necessarily allocated to these lookup
       objects if they are not required */
    Interp2d rho;   /* density, kg/m^3 */
    Interp2d dTdPs; /* adiabatic temperature gradient, K/Pa */
    Interp2d cp;    /* heat capacity, J/kg/K */
    Interp2d temp;  /* temperature, K */
    Interp2d alpha; /* thermal expansion, 1/K */
} data_EOSLookup;

#endif
