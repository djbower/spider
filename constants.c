#include "constants.h"

PetscErrorCode set_constants( Ctx *E )
{
    Constants *C = &E->constants;
    PetscScalar SQRTST;

    PetscFunctionBeginUser;

    // 35 constants to set
    SQRTST = PetscSqrtScalar( ENTROPY0 * TEMPERATURE0 );

    C->RADIUS   = RADIUS0;
    C->TEMP     = TEMPERATURE0;
    C->ENTROPY  = ENTROPY0;
    C->DENSITY  = DENSITY0;
    C->AREA     = PetscSqr( C->RADIUS );
    C->AREAG    = C->AREA * 4.0 * PETSC_PI;
    C->VOLUME   = C->AREA * C->RADIUS;
    C->VOLUMEG  = C->VOLUME * 4.0 * PETSC_PI;
    C->MASS     = C->DENSITY * C->VOLUME;
    C->MASSG    = C->MASS * 4.0 * PETSC_PI;
    C->TIME     = C->RADIUS / SQRTST;
    C->TIMEYRS  = C->TIME / (60.0*60.0*24.0*365.25);
    C->SENERGY  = C->ENTROPY * C->TEMP;
    C->ENERGY   = C->SENERGY * C->MASS;
    C->ENERGYG  = C->ENERGY * 4.0 * PETSC_PI;
    C->PRESSURE = C->ENTROPY * C->TEMP * C->DENSITY;
    C->POWER    = C->ENERGY / C->TIME;
    C->POWERG   = C->POWER * 4.0 * PETSC_PI;
    C->FLUX     = C->POWER / C->AREA;
    C->DPDR     = C->PRESSURE / C->RADIUS;
    C->ALPHA    = 1.0 / C->TEMP;
    C->GRAVITY  = (C->ENTROPY * C->TEMP) / C->RADIUS;
    C->KAPPA    = C->RADIUS * SQRTST;
    C->DTDP     = 1.0 / (C->DENSITY * C->ENTROPY);
    C->DSDR     = C->ENTROPY / C->RADIUS;
    C->DTDR     = C->TEMP / C->RADIUS;
    C->GSUPER   = C->GRAVITY * C->DTDR;
    C->ETA      = C->DENSITY * C->KAPPA;
    C->LOG10ETA = PetscLog10Real( C->ETA );
    C->NU       = C->KAPPA;
    C->COND     = C->ENTROPY * C->DENSITY * C->KAPPA;
    C->SIGMA    = C->FLUX * 1.0 / PetscPowScalar( C->TEMP, 4.0 );
    C->LHS      = C->DENSITY * C->VOLUME * C->TEMP;
    C->LHSG     = C->LHS * 4.0 * PETSC_PI;
    C->RHS      = C->ENTROPY / C->TIME;

    PetscFunctionReturn(0);
}
