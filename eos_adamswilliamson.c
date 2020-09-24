#include "eos_adamswilliamson.h"
#include "util.h"

/* Prototypes for local functions used in EOS interface functions */
static PetscErrorCode EOSAdamsWilliamson_GetRho( const data_EOSAdamsWilliamson*,PetscScalar,PetscScalar,PetscScalar*);

/* EOS interface functions */
static PetscErrorCode EOSEval_AdamsWilliamson(EOS eos, PetscScalar P, PetscScalar S, EOSEvalData *eval)
{
  PetscErrorCode  ierr;
  data_EOSAdamsWilliamson *adams = (data_EOSAdamsWilliamson*) eos->impl_data;

  PetscFunctionBegin;
  eval->P = P;
  eval->S = S;
  eval->T = 0.0;
  ierr = EOSAdamsWilliamson_GetRho( adams, P, S, &eval->rho );CHKERRQ(ierr);
  eval->Cp = 0.0;
  eval->dTdPs = 0.0;
  eval->alpha = 0.0;
  eval->Cv = 0.0;
  eval->V = 0.0;
  eval->log10visc = 0.0;
  eval->cond = 0.0;
  eval->phase_fraction = 1.0; // by definition, since only one phase
  PetscFunctionReturn(0);
}

static PetscErrorCode EOSDestroy_AdamsWilliamson(EOS eos)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(eos->impl_data);CHKERRQ(ierr);
  eos->impl_data = NULL;
  PetscFunctionReturn(0);
}


PetscErrorCode EOSSetUpFromOptions_AdamsWilliamson(EOS eos, const char *prefix, const FundamentalConstants FC, const ScalingConstants SC)
{
  PetscErrorCode  ierr;
  char            buf[PETSC_MAX_PATH_LEN]; /* max size */
  data_EOSAdamsWilliamson *data = (data_EOSAdamsWilliamson*) eos->impl_data;

  PetscFunctionBegin;
  (void) FC; // unused

  /* gravity (m/s^2), must be negative */
  data->gravity = -10.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-gravity",&data->gravity,NULL);CHKERRQ(ierr);
  data->gravity /= SC->GRAVITY;

  /* surface density (kg/m^3) */
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_rhos");CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf,&data->rhos,4078.95095544,NULL);CHKERRQ(ierr);
  data->rhos /= SC->DENSITY;

  /* parameter (1/m) */
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_beta");CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf,&data->beta,1.1115348931000002e-07,NULL);CHKERRQ(ierr);
  data->beta *= SC->RADIUS;

  PetscFunctionReturn(0);
}

/* Creation Function */
PetscErrorCode EOSCreate_AdamsWilliamson(EOS eos) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, (data_EOSAdamsWilliamson**) (&eos->impl_data));CHKERRQ(ierr);
  eos->eval = EOSEval_AdamsWilliamson;
  eos->destroy = EOSDestroy_AdamsWilliamson;
  eos->setupfromoptions = EOSSetUpFromOptions_AdamsWilliamson;
  PetscFunctionReturn(0);
}

/* Helper functions */

static PetscErrorCode EOSAdamsWilliamson_GetRho( const data_EOSAdamsWilliamson *adams,PetscScalar P,PetscScalar S, PetscScalar *rho_ptr)
{
    PetscScalar rho;
    
    PetscFunctionBeginUser;

    rho = adams->rhos - P * adams->beta / adams->gravity;

    *rho_ptr = rho;

    PetscFunctionReturn(0);
}
