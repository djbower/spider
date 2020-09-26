#include "eos_adamswilliamson.h"
#include "util.h"

/* Prototypes for local functions used in EOS interface functions */
static PetscErrorCode EOSAdamsWilliamson_GetRho(const data_EOSAdamsWilliamson*,PetscScalar,PetscScalar,PetscScalar*);
static PetscErrorCode EOSAdamsWilliamson_GetRadiusFromPressure( const data_EOSAdamsWilliamson*, PetscScalar, PetscScalar * );
static PetscErrorCode EOSAdamsWilliamson_GetPressureFromRadius( const data_EOSAdamsWilliamson*, PetscScalar, PetscScalar * );
static PetscErrorCode EOSAdamsWilliamson_GetMassWithinRadius( const data_EOSAdamsWilliamson*, PetscScalar, PetscScalar *);
static PetscErrorCode EOSAdamsWilliamson_GetMassWithinShell( const data_EOSAdamsWilliamson*, PetscScalar, PetscScalar, PetscScalar *);
static PetscErrorCode EOSAdamsWilliamson_GetMassCoordinateAverageRho( const data_EOSAdamsWilliamson*, PetscScalar * );
static PetscErrorCode EOSAdamsWilliamson_GetMassElement( const data_EOSAdamsWilliamson*, PetscScalar, PetscScalar * );

/* EOS interface functions */
static PetscErrorCode EOSEval_AdamsWilliamson(EOS eos, PetscScalar P, PetscScalar S, EOSEvalData *eval)
{
  PetscErrorCode  ierr;
  data_EOSAdamsWilliamson *adams = (data_EOSAdamsWilliamson*) eos->impl_data;
  (void) S; //unused

  PetscFunctionBegin;
  eval->P = P;
  ierr = EOSAdamsWilliamson_GetRho( adams, P, S, &eval->rho );CHKERRQ(ierr);
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
  ierr = PetscOptionsGetPositiveScalar(buf,&data->density_surface,4078.95095544,NULL);CHKERRQ(ierr);
  data->density_surface /= SC->DENSITY;

  /* parameter (1/m) */
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_beta");CHKERRQ(ierr);
  ierr = PetscOptionsGetPositiveScalar(buf,&data->beta,1.1115348931000002e-07,NULL);CHKERRQ(ierr);
  data->beta *= SC->RADIUS;

  /* we need the average density for the mass coordinate to radius mapping */
  ierr = EOSAdamsWilliamson_GetMassCoordinateAverageRho( data, &data->density_average );CHKERRQ(ierr);

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
    rho = adams->density_surface - P * adams->beta / adams->gravity;

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

    R = adams->radius - (1.0/adams->beta)*PetscLogScalar(1.0-(P*adams->beta)/(adams->density_surface*adams->gravity));

    *R_ptr = R;

    PetscFunctionReturn(0);
}

static PetscErrorCode EOSAdamsWilliamson_GetPressureFromRadius( const data_EOSAdamsWilliamson *adams, PetscScalar R, PetscScalar *P_ptr )
{
    PetscScalar P;

    PetscFunctionBeginUser;

    P = -adams->density_surface * adams->gravity / adams->beta;
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

static PetscErrorCode EOSAdamsWilliamson_GetMassElement( const data_EOSAdamsWilliamson *adams, PetscScalar R, PetscScalar *mass_ptr )
{
  PetscErrorCode ierr;
  PetscScalar mass, P, S=0.0; // S unused

  PetscFunctionBeginUser;

  ierr = EOSAdamsWilliamson_GetPressureFromRadius( adams, R, &P );CHKERRQ(ierr);
  ierr = EOSAdamsWilliamson_GetRho( adams, P, S, &mass );CHKERRQ(ierr);
  mass *= PetscPowScalar( R, 2.0 );

  *mass_ptr = mass;

  PetscFunctionReturn(0);
}

static PetscErrorCode EOSAdamsWilliamson_GetMassCoordinateAverageRho( const data_EOSAdamsWilliamson *adams, PetscScalar *rho_ptr )
{
  PetscErrorCode ierr;
  PetscScalar rho, mass;

  PetscFunctionBeginUser;

  /* radius_core < r < radius: mass within mantle */
  ierr = EOSAdamsWilliamson_GetMassWithinShell( adams, adams->radius, adams->radius_core, &mass ); CHKERRQ(ierr);

  /* use the same range applied to construct the mass coordinate mesh
     radius < xi < radius_core
     hence rho is the actual average density of the mantle */
  rho = mass * 3.0 / ( PetscPowScalar( adams->radius, 3.0 ) - PetscPowScalar( adams->radius_core, 3.0) );

  *rho_ptr = rho;

  PetscFunctionReturn(0);
}

PetscErrorCode EOSAdamsWilliamson_ObjectiveFunctionRadius( SNES snes, Vec x, Vec f, void *ptr )
{
  /* uses definition of mass coordinates */

  PetscErrorCode    ierr;
  const PetscScalar *xx, *xi_b, *xi_s;
  PetscScalar       *ff, xi; 
  Ctx               *E = (Ctx*) ptr;
  Parameters  const P = E->parameters;
  Mesh        const *M = &E->mesh;
  PetscInt          i,numpts_b,numpts_s;
  EOS               eos = P->eos_mesh;
  data_EOSAdamsWilliamson *adams = (data_EOSAdamsWilliamson*) eos->impl_data;

  PetscFunctionBeginUser;

  /* TODO: replace with composite DM */
  /* below is simple loop over all nodes (basic and staggered), but will clearly break for parallel */
  ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

  ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr); /* initial guess of r */
  /* below instead should use composite DM */
  ierr = VecGetArrayRead(M->xi_b,&xi_b);CHKERRQ(ierr); /* target mass coordinate */
  ierr = VecGetArrayRead(M->xi_s,&xi_s);CHKERRQ(ierr);
  ierr = VecGetArray(f,&ff);CHKERRQ(ierr); /* residual function */

  for(i=0; i< numpts_b+numpts_s; ++i){
    /* mass within radial shell of mantle */
    ierr = EOSAdamsWilliamson_GetMassWithinShell( adams, xx[i], adams->radius_core, &ff[i] ); CHKERRQ(ierr);

    /* get mass coordinate */
    if (i<numpts_b){
      xi = xi_b[i]; // basic nodes
    }   
    else{
      xi = xi_s[i-numpts_b]; // staggered nodes
    }   

    /* difference from mass coordinate mass */
    ff[i] -= (adams->density_average / 3.0) * (PetscPowScalar(xi,3.0)-PetscPowScalar(adams->radius_core,3.0));

  }

  ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(M->xi_b,&xi_b);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(M->xi_s,&xi_s);CHKERRQ(ierr);
  ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

PetscErrorCode EOSAdamsWilliamson_JacobianRadius( SNES snes, Vec x, Mat jac, Mat B, void *ptr)
{
  PetscErrorCode    ierr;
  const PetscScalar *xx;
  Ctx               *E = (Ctx*) ptr;
  Parameters  const P = E->parameters;
  PetscInt          i,numpts_b,numpts_s;
  Vec               diag;
  PetscScalar       *arr_diag;
  EOS               eos = P->eos_mesh;
  data_EOSAdamsWilliamson *adams = (data_EOSAdamsWilliamson*) eos->impl_data;

  PetscFunctionBeginUser;

  ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

  /* work vector for diagonal of Jacobian */
  ierr = VecCreate( PETSC_COMM_WORLD, &diag );CHKERRQ(ierr);
  ierr = VecSetSizes( diag, PETSC_DECIDE, numpts_b+numpts_s );CHKERRQ(ierr);
  ierr = VecSetFromOptions(diag);CHKERRQ(ierr);

  ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
  ierr = VecGetArray(diag,&arr_diag);CHKERRQ(ierr);

  /* set diagonal */
  for(i=0; i < numpts_b+numpts_s; ++i){
    ierr = EOSAdamsWilliamson_GetMassElement( adams, xx[i], &arr_diag[i] );CHKERRQ(ierr); 
  }

  ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArray(diag,&arr_diag);CHKERRQ(ierr);

  ierr = MatDiagonalSet( B, diag, INSERT_VALUES );CHKERRQ(ierr);

  /* Assemble matrix */
  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
  if (jac != B) {
    MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);
  }

  ierr = VecDestroy(&diag);

  PetscFunctionReturn(0);

}

/* below could be generalised to accommodate any EOS */
PetscErrorCode EOSAdamsWilliamson_MassCoordinateSpatialDerivative( const data_EOSAdamsWilliamson *adams, PetscScalar R, PetscScalar xi, PetscScalar *dxidr_ptr )
{
  PetscErrorCode ierr;
  PetscScalar dxidr, P, S=0.0; // S unused

  PetscFunctionBeginUser;

  ierr = EOSAdamsWilliamson_GetPressureFromRadius( adams, R, &P );CHKERRQ(ierr);
  ierr = EOSAdamsWilliamson_GetRho( adams, P, S, &dxidr );CHKERRQ(ierr);
  dxidr /= adams->density_average;
  dxidr *= (R/xi) * (R/xi);

  *dxidr_ptr = dxidr;

  PetscFunctionReturn(0);

}
