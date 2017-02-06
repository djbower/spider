#ifndef MONITOR_H_
#define MONITOR_H_

#include "ctx.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscInt step, PetscReal time, PetscReal time0, PetscReal timeprev, Vec x, void * ptr, double walltime0, double* walltimeprev);

#endif
