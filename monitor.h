#ifndef MONITOR_H_
#define MONITOR_H_

#include "ctx.h"

typedef struct MonitorCtx
{
  double walltime0, walltimeprev;
  PetscBool outputDirectoryExistenceConfirmed;
} MonitorCtx;

PetscErrorCode SNESMonitorVerbose(SNES snes, PetscInt its, PetscReal norm, void *mctx);
PetscErrorCode TSCustomMonitor(TS ts, PetscInt step, PetscReal time, Vec x, void *ptr, MonitorCtx *mctx);
PetscErrorCode TSMonitorWalltimed(TS ts, PetscInt steps, PetscReal time, Vec x, void *mctx);

#endif
