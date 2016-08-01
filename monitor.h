#ifndef MONITOR_H_
#define MONITOR_H_

#include "ctx.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscInt step, PetscReal time, Vec x, void * ptr);

#endif
