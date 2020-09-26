#include "eos_adamswilliamson.h"
#include "util.h"

/* Prototypes for local functions used in EOS interface functions */
static PetscErrorCode EOSAdamsWilliamson_GetRho(const data_EOSAdamsWilliamson*,PetscScalar,PetscScalar,PetscScalar*);
static PetscErrorCode EOSAdamsWilliamson_GetRadiusFromPressure( const data_EOSAdamsWilliamson*, PetscScalar, PetscScalar * );
static PetscErrorCode EOSAdamsWilliamson_GetPressureFromRadius( const data_EOSAdamsWilliamson*, PetscScalar, PetscScalar * );
static PetscErrorCode EOSAdamsWilliamson_GetMassWithinRadius( const data_EOSAdamsWilliamson*, PetscScalar, PetscScalar *);
static PetscErrorCode EOSAdamsWilliamson_GetMassWithinShell( const data_EOSAdamsWilliamson*, PetscScalar, PetscScalar, PetscScalar *);

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

  /* radius of planet (m) */
  ierr = PetscOptionsGetPositiveScalar("-radius",&data->radius,6371000.0,NULL);CHKERRQ(ierr); // m
  data->radius /= SC->RADIUS;

  /* core radius relative to physical radius i.e. radius
     therefore, scaled (code) core radius is P->coresize * P->radius
     and actual physical core radius is P->coresize * P->radius * SC->RADIUS */
  ierr = PetscOptionsGetPositiveScalar("-coresize",&data->radius_core,0.55,NULL);CHKERRQ(ierr); // Earth core radius
  /* already non-dimensonal, but more convenient to have core radius directly */
  data->radius_core *= data->radius;

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
   /* Adams-Williamson density is a simple function of depth (radius)
      Sketch derivation:
          dP/dr = dP/drho * drho/dr = -rho g
          dP/drho \sim (dP/drho)_s (adiabatic)
          drho/dr = -rho g / Si
          then integrate to give the form rho(r) = k * exp(-(g*r)/c)
          (g is positive)
          apply the limit that rho = rhos at r=R
          gives:
              rho(z) = rhos * exp( beta * z )
          where z = R-r

    this is arguably the simplest relation to get rho directly from r, but other
    EOSs can be envisaged */

    PetscScalar rho;
    
    PetscFunctionBeginUser;
    (void) S;

    /* using pressure, expression is simpler than sketch derivation above */
    rho = adams->rhos - P * adams->beta / adams->gravity;

    *rho_ptr = rho;

    PetscFunctionReturn(0);
}

/* Mass coordinate mapping and static structure */

/* Currently, these functions are accessed directly, not through an interface, and hence they
   take EOS as an argument (which must be a AdamsWilliamson EOS obviously) */

static PetscErrorCode EOSAdamsWilliamson_GetRadiusFromPressure( const data_EOSAdamsWilliamson *adams, PetscScalar P, PetscScalar *R_ptr )
{
    PetscScalar R;

    PetscFunctionBeginUser;

    R = adams->radius - (1.0/adams->beta)*PetscLogScalar(1.0-(P*adams->beta)/(adams->rhos*adams->gravity));

    *R_ptr = R;

    PetscFunctionReturn(0);
}

static PetscErrorCode EOSAdamsWilliamson_GetPressureFromRadius( const data_EOSAdamsWilliamson *adams, PetscScalar R, PetscScalar *P_ptr )
{
    PetscScalar P;

    PetscFunctionBeginUser;

    P = -adams->rhos * adams->gravity / adams->beta;
    P *= PetscExpScalar( adams->beta*(adams->radius-R) ) - 1.0;

    *P_ptr = P;

    PetscFunctionReturn(0);
}

static PetscErrorCode EOSAdamsWilliamson_GetMassWithinRadius( const data_EOSAdamsWilliamson *adams, PetscScalar R, PetscScalar *mass_ptr)
{
  /* return integral from 0 to r of r^2 * rho dr */

  PetscErrorCode  ierr;
  PetscScalar mass, P, rho; /* note mass without 4*pi scaling, as convention in SPIDER */
  PetscScalar const beta = adams->beta;
  PetscScalar S = 0.0; // S not used in this function;

  PetscFunctionBeginUser;

  mass = -2.0/PetscPowScalar(beta,3) - PetscPowScalar(R,2)/beta -2*R/PetscPowScalar(beta,2);
  ierr = EOSAdamsWilliamson_GetPressureFromRadius( adams, R, &P );CHKERRQ(ierr);
  ierr = EOSAdamsWilliamson_GetRho( adams, P, S, &rho );CHKERRQ(ierr);
  mass *= rho;

  *mass_ptr = mass;

  PetscFunctionReturn(0);
}

static PetscErrorCode EOSAdamsWilliamson_GetMassWithinShell( const data_EOSAdamsWilliamson *adams, PetscScalar Rout, PetscScalar Rin, PetscScalar *mass_ptr)
{
  PetscErrorCode ierr;
  PetscScalar massout, massin; /* note mass without 4*pi scaling, as convention in SPIDER */

  PetscFunctionBeginUser;

  ierr = EOSAdamsWilliamson_GetMassWithinRadius( adams, Rout, &massout );CHKERRQ(ierr);
  ierr = EOSAdamsWilliamson_GetMassWithinRadius( adams, Rin, &massin );CHKERRQ(ierr);

  *mass_ptr = massout - massin;

  PetscFunctionReturn(0);
}

PetscErrorCode EOSAdamsWilliamson_GetMassCoordinateAverageRho( EOS eos, PetscScalar *rho_ptr )
{
  PetscErrorCode ierr;
  PetscScalar rho, mass;
  data_EOSAdamsWilliamson *adams = (data_EOSAdamsWilliamson*) eos->impl_data;

  PetscFunctionBeginUser;

  ierr = EOSAdamsWilliamson_GetMassWithinShell( adams, adams->radius, adams->radius_core, &mass ); CHKERRQ(ierr);

  /* a choice is made here to define rho based on a mass coordinate
     0 < xi < radius.  Recall that radius is non-dimensional, so it
     can always be chosen (if desired) such that 0 < xi < 1.0 */
  rho = mass * 3.0 / PetscPowScalar( adams->radius, 3.0 );

  *rho_ptr = rho;

  PetscFunctionReturn(0);
}
