#ifndef LOOKUP_H_
#define LOOKUP_H_

#include <petsc.h>

typedef struct _Interp1d {
    PetscInt    NX;
    PetscScalar *xa;
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar *ya;
    PetscScalar ymin;
    PetscScalar ymax;
} Interp1d;

typedef struct _Interp2d {
    PetscInt    NX;
    PetscScalar *xa;
    PetscScalar xmin;
    PetscScalar xmax;
    PetscScalar dx;
    PetscInt    NY;
    PetscScalar *ya;
    PetscScalar ymin;
    PetscScalar ymax;
    PetscScalar dy;
    PetscScalar **za;
} Interp2d;

/* lookup for a single phase */
typedef struct _Lookup {
    Interp2d rho; /* density, kg/m^3 */
    Interp2d dTdPs; /* adiabatic temperature gradient, K/Pa */
    Interp2d cp; /* heat capacity, J/kg/K */
    Interp2d temp; /* temperature, K */
    Interp2d alpha; /* thermal expansion, 1/K */
    PetscScalar cond; /* thermal conductivity, W/m/K */
    PetscScalar log10visc; /* log base 10 of viscosity */
    /* for single component system these are the same for both solid
      and melt phases, but formulating it this way might provide a
      a way forward for a multicomponent system */
    Interp1d liquidus; /* liquidus, J/kg/K */
    Interp1d solidus; /* solidus, J/kg/K */
} Lookup;

PetscErrorCode set_interp1d( const char *, Interp1d *, PetscScalar, PetscScalar );
PetscScalar get_val1d( Interp1d const *, PetscScalar );
void free_interp1d( Interp1d * );

PetscErrorCode set_interp2d( const char *, Interp2d *, PetscScalar, PetscScalar, PetscScalar );
PetscScalar get_val2d( Interp2d const *, PetscScalar, PetscScalar );
void free_interp2d( Interp2d * );

#endif
