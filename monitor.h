#ifndef MONITOR_H_
#define MONITOR_H_

#include "ctx.h"

PetscErrorCode TSCustomMonitor(TS ts, PetscReal dtmacro, PetscInt step, PetscReal time, Vec x, void * ptr, double walltime0, double* walltimeprev);

typedef struct TSMonitorWalltimedCtx {
  double walltime0,walltimeprev;
} TSMonitorWalltimedCtx;

PetscErrorCode TSMonitorWalltimed(TS ts,PetscInt steps,PetscReal time,Vec x,void *mctx);

#endif
