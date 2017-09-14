#ifndef MONITOR_H_
#define MONITOR_H_

#include "ctx.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscReal dtmacro, PetscInt step, PetscReal time, Vec x, void * ptr, double walltime0, double* walltimeprev);

#endif
