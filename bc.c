#include "bc.h"
#include "util.h"

static PetscScalar hybrid( Ctx *, PetscScalar );
static PetscScalar tsurf_param( PetscScalar, AtmosphereParameters const * );
static PetscScalar grey_body( Atmosphere const *, AtmosphereParameters const * );
static PetscScalar steam_atmosphere_zahnle_1988( Atmosphere const *, Constants const *C );
static PetscScalar get_emissivity_abe_matsui( Parameters const *, Atmosphere * );
static PetscScalar get_emissivity_from_flux( Atmosphere const *, AtmosphereParameters const *, PetscScalar );
static PetscScalar get_atmosphere_mass( Parameters const *, PetscScalar );
static PetscScalar get_optical_depth( Parameters const *, PetscScalar, VolatileParameters const * );
static PetscScalar get_dxdt( Ctx const *, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_partial_pressure_volatile( PetscScalar, VolatileParameters const *, Constants const * );
static PetscScalar get_partial_pressure_derivative_volatile( PetscScalar, VolatileParameters const *, Constants const * );

/* to solve for the initial volatile content of the magma ocean (liquid)
   we use Newton's method */
static PetscScalar f( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar f_prim( PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar newton( PetscScalar, PetscScalar, PetscScalar );

PetscErrorCode set_surface_flux( Ctx *E )
{
    PetscErrorCode       ierr;
    PetscMPIInt          rank;
    PetscScalar          temp0, Qout;
    PetscInt             ind;
    Atmosphere           *A  = &E->atmosphere;
    Parameters           const *P  = &E->parameters;
    Constants            const *C  = &P->constants;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    Solution             *S  = &E->solution;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_surface_flux:\n");CHKERRQ(ierr);
#endif
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

    if (!rank){
      /* get surface values to pass to flux boundary condition
         function that is chosen by user */
      ind = 0;

      /* temperature (potential temperature if coarse mesh is used)
         this is taken from the top staggered node */
      ierr = VecGetValues(S->temp,1,&ind,&temp0); CHKERRQ(ierr);

      /* correct for ultra-thin thermal boundary layer at the surface */
      if( Ap->PARAM_UTBL ){
        A->tsurf = tsurf_param( temp0, Ap); // parameterised boundary layer
      }
      else{
        A->tsurf = temp0; // surface temperature is potential temperature
      }

      /* determine surface flux */
      switch( Ap->MODEL ){
        case 1:
          // grey-body with constant emissivity
          A->emissivity = Ap->emissivity0;
          Qout = grey_body( A, Ap );
          break;
        case 2:
          /* trying to pass the Constants struct resulted in a circular
             dependency, which was easiest to address by just passing
             in the two required constants instead */
          Qout = steam_atmosphere_zahnle_1988( A, C );
          break;
        case 3:
          // atmosphere evolution
          A->emissivity = get_emissivity_abe_matsui( P, A );
          Qout = grey_body( A, Ap );
          break;
        default:
          SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported model value %d provided",Ap->MODEL);
      }

      /* for legacy purposes */
      if( Ap->HYBRID ){
          Qout = hybrid( E, Qout );
      }

      /* some atmosphere models do not explicitly set an emissivity,
         so we should back-compute it here to always ensure that the
         emissivity is consistent with our chosen atmosphere scenario
         For cases where the emissivity is set, this should yield the
         same answer of course, and adds little computational overhead
         in those cases */
      A->emissivity = get_emissivity_from_flux( A, Ap, Qout );

      ierr = VecSetValue(S->Jtot,0,Qout,INSERT_VALUES);CHKERRQ(ierr);
      /* TODO: probably OK, but remember this has a scaling dictated
         by the dS/dr non-dimensionalisation scheme, and not the
         atmosphere scheme */
      Qout *= PetscSqr( P->radius );
      ierr = VecSetValue(S->Etot,0,Qout,INSERT_VALUES);CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscScalar hybrid( Ctx *E, PetscScalar Qin )
{
    PetscErrorCode ierr;
    PetscScalar    Qout;
    PetscScalar    G0, R0, R1, R2, E0, E1, E2, Q2, fwt, phi0;
    PetscInt       ind;
    Mesh           *M = &E->mesh;
    Parameters     *P = &E->parameters;
    Solution       *S = &E->solution;

    /* for legacy purposes, enable the ability for the magma ocean
       to cool at a rate dictated by the upper mantle cooling rate,
       This was originally a hack to prevent a viscous lid from
       forming at the top */

    /* for weight of different fluxes */
    ind = 0;
    ierr = VecGetValues(S->phi,1,&ind,&phi0); CHKERRQ(ierr);
    /* SWIDTH or PHI_WIDTH most appropriate choice here? */
    fwt = tanh_weight( phi0, P->phi_critical, P->swidth );

    // energy flux from energy gradient
    ierr = VecGetValues(M->area_b,1,&ind,&G0); CHKERRQ(ierr);
    ierr = VecGetValues(M->radius_b,1,&ind,&R0); CHKERRQ(ierr);
    ind = 1;
    ierr = VecGetValues(M->radius_b,1,&ind,&R1); CHKERRQ(ierr);
    ierr = VecGetValues(S->Etot,1,&ind,&E1); CHKERRQ(ierr);
    ind = 2;
    ierr = VecGetValues(M->radius_b,1,&ind,&R2); CHKERRQ(ierr);
    ierr = VecGetValues(S->Etot,1,&ind,&E2); CHKERRQ(ierr);
    E0 = E1 - (E2-E1)*(R2-R1)/(R1-R0); // energy at surface
    Q2 = E0 / G0;
    // weight Q1 and Q2 to give total flux
    Qout = Qin * fwt + Q2 * (1.0 - fwt);

    return Qout;
}

PetscErrorCode set_core_mantle_flux( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          ix;             // index of last basic node
    PetscInt          ix2;            // index of penultimate basic node
    PetscInt          numpts_b;
    PetscScalar       fac,vol,vol_core,area1,area2,rho_cmb,cp_cmb;
    PetscScalar       val;
    PetscMPIInt       rank, size;

    Mesh              *M = &E->mesh;
    Parameters        *P = &E->parameters;
    Solution          *S = &E->solution;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_core_mantle_flux:\n");CHKERRQ(ierr);
#endif
    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ix  = numpts_b-1;
    ix2 = numpts_b-2;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

    /* Assume that the last rank contains the last two points */
    if (rank == size-1 ){
        ierr = VecGetValues( M->area_b,1,&ix,&area1);CHKERRQ(ierr); // without 4*pi
        ierr = VecGetValues( M->area_b,1,&ix2,&area2);CHKERRQ(ierr); // without 4*pi
        ierr = VecGetValues( M->volume_s,1,&ix2,&vol);CHKERRQ(ierr); // without 4*pi
        ierr = VecGetValues( S->rho_s,1,&ix2,&rho_cmb);CHKERRQ(ierr);
        ierr = VecGetValues( S->cp_s,1,&ix2,&cp_cmb);CHKERRQ(ierr);
        vol_core = 1.0/3.0 * PetscPowScalar(P->coresize,3.0) * PetscPowScalar(P->radius,3.0); // without 4*pi
        fac = vol / vol_core; // excluding 4*pi is OK since here we take the ratio of two volumes
        fac *= rho_cmb / P->rho_core;
        fac *= cp_cmb / P->cp_core;
        fac /= P->tfac_core_avg;

        fac = 1.0 / (1.0 + fac);
        fac *= area2 / area1; // excluding 4*pi is OK since here we take the ratio of two areas

        /* energy flux (Jtot) */
        ierr = VecGetValues(S->Jtot,1,&ix2,&val);CHKERRQ(ierr);
        val *= fac;
        ierr = VecSetValue(S->Jtot,ix,val,INSERT_VALUES);CHKERRQ(ierr);

        /* energy flow (Etot) */
        val *= area1; // without 4*pi is correct
        ierr = VecSetValue(S->Etot,ix,val,INSERT_VALUES);CHKERRQ(ierr);
    }

    /* this is the only place where VecAssembly calls operate
       on Jtot and Etot */
    ierr = VecAssemblyBegin(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

////////////////////////
/* atmosphere related */
////////////////////////

PetscErrorCode set_atmosphere_volatile_content( Ctx *E , PetscScalar x0, PetscScalar x1 )
{
    Atmosphere           *A  = &E->atmosphere;
    Parameters           const *P  = &E->parameters;
    Constants            const *C  = &P->constants;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    VolatileParameters   const *CO2 = &Ap->CO2_volatile_parameters;
    VolatileParameters   const *H2O = &Ap->H2O_volatile_parameters;

    PetscFunctionBeginUser;

    /* if x0 and/or x1 are zero, the quantities below will also all
       be set to zero */

    /* CO2 */
    A->p0 = get_partial_pressure_volatile( x0, CO2, C );
    A->dp0dx = get_partial_pressure_derivative_volatile( x0, CO2, C );
    A->m0 = get_atmosphere_mass( P, A->p0 );

    /* H2O */
    A->p1 = get_partial_pressure_volatile( x1, H2O, C );
    A->dp1dx = get_partial_pressure_derivative_volatile( x1, H2O, C );
    A->m1 = get_atmosphere_mass( P, A->p1 );

    PetscFunctionReturn(0);
}

static PetscScalar get_partial_pressure_volatile( PetscScalar x, VolatileParameters const *V, Constants const *C)
{
    /* partial pressure of volatile */
    PetscScalar p;

    p = 1.0 / PetscPowScalar( V->henry, V->henry_pow );
    p *= PetscPowScalar( x, V->henry_pow );
    p /= C->PRESSURE; // non-dimensionalise

    return p;
}

static PetscScalar get_partial_pressure_derivative_volatile( PetscScalar x, VolatileParameters const *V, Constants const *C )
{
    /* partial pressure derivative of volatile */
    PetscScalar dpdx;

    dpdx = 1.0 / PetscPowScalar( V->henry, V->henry_pow );
    dpdx *= V->henry_pow * PetscPowScalar( x, V->henry_pow-1.0 );
    dpdx /= C->PRESSURE; // non-dimensionalise

    return dpdx;
}

static PetscScalar get_atmosphere_mass( Parameters const *P, PetscScalar p )
{
    /* mass of atmosphere */ 
    PetscScalar mass_atm;

    mass_atm = PetscSqr(P->radius) * p / -P->gravity;

    return mass_atm;
}


PetscScalar tsurf_param( PetscScalar temp, AtmosphereParameters const *Ap )
{
    PetscScalar Ts, c, fac, num, den;
    c = Ap->param_utbl_const;

    fac = 3.0*PetscPowScalar(c,3.0)*(27.0*PetscPowScalar(temp,2.0)*c+4.0);
    fac = PetscPowScalar( fac, 1.0/2.0 );
    fac += 9.0*temp*PetscPowScalar(c,2.0);
    // numerator
    num = PetscPowScalar(2.0,1.0/3)*PetscPowScalar(fac,2.0/3)-2.0*PetscPowScalar(3.0,1.0/3)*c;
    // denominator
    den = PetscPowScalar(6.0,2.0/3)*c*PetscPowScalar(fac,1.0/3);
    // surface temperature
    Ts = num / den;

    return Ts; 
}

PetscScalar grey_body( Atmosphere const *A, AtmosphereParameters const *Ap )
{
    PetscScalar Fsurf;

    Fsurf = PetscPowScalar(A->tsurf,4.0)-PetscPowScalar(Ap->teqm,4.0);
    Fsurf *= Ap->sigma * A->emissivity; /* Note emissivity may vary */

    return Fsurf;
}

PetscScalar steam_atmosphere_zahnle_1988( Atmosphere const *A, Constants const *C )
{
    PetscScalar       Tsurf, Fsurf;

    /* fit to Zahnle et al. (1988) from Solomatov and Stevenson (1993)
       Eqn. 40.  Expressed dimensionally so must convert here using
       TEMP and FLUX scalings */

    Tsurf = A->tsurf * C->TEMP;
    Fsurf = 1.5E2 + 1.02E-5 * PetscExpScalar(0.011*Tsurf);
    Fsurf /= C->FLUX; // non-dimensionalise

    return Fsurf;

}

PetscScalar get_emissivity_from_flux( Atmosphere const *A, AtmosphereParameters const *Ap, PetscScalar flux )
{
    PetscScalar emissivity;

    emissivity = flux / Ap->sigma;
    emissivity /= PetscPowScalar(A->tsurf,4.0)-PetscPowScalar(Ap->teqm,4.0);

    return emissivity;

}

PetscScalar get_emissivity_abe_matsui( Parameters const *P, Atmosphere *A )
{
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;

    PetscScalar emissivity;

    /* CO2 */
    A->tau0 = get_optical_depth( P, A->m0, &Ap->CO2_volatile_parameters );

    /* H2O */
    A->tau1 = get_optical_depth( P, A->m1, &Ap->H2O_volatile_parameters );

    /* total */
    A->tau = A->tau0 + A->tau1;
    emissivity = 2.0 / (A->tau + 2.0);

    return emissivity;

}

static PetscScalar get_optical_depth( Parameters const *P, PetscScalar mass_atm, VolatileParameters const *V )
{
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    PetscScalar tau;

    // note negative gravity
    tau = 0.5 * mass_atm / PetscSqr( P->radius );
    tau *= PetscSqrtScalar( 3.0 * V->kabs * -P->gravity / Ap->P0 );

    return tau; // dimensionless (by definition)
}

/////////////////////////////////////
/* general functions for volatiles */
/////////////////////////////////////
static PetscScalar get_dxdt( Ctx const *E, PetscScalar x, PetscScalar dpdx, PetscScalar kdist )
{
    Parameters           const *P = &E->parameters;
    Atmosphere           const *A = &E->atmosphere;
    Constants            const *C = &P->constants;
    Mesh                 const *M = &E->mesh;

    PetscScalar dxdt;
    PetscScalar num, den;

    num = x * (kdist-1.0) * A->dMliqdt;
    den = kdist * M->mantle_mass + (1.0-kdist) * A->Mliq;
    den += C->VOLSCALE * PetscSqr(P->radius) * dpdx / -P->gravity;

    dxdt = num / den;

    return dxdt;
}

PetscScalar get_dx0dt( Ctx const *E, PetscScalar x0 )
{
    /* update for dissolved CO2 content in the magma ocean */
    Atmosphere const *A = &E->atmosphere;
    AtmosphereParameters const *Ap = &E->parameters.atmosphere_parameters;
    VolatileParameters const *CO2 = &Ap->CO2_volatile_parameters;

    PetscScalar dx0dt;

    dx0dt = get_dxdt( E, x0, A->dp0dx, CO2->kdist );

    return dx0dt;
}

PetscScalar get_dx1dt( Ctx const *E, PetscScalar x1 )
{
    /* update for dissolved H2O content in the magma ocean */
    Atmosphere const *A = &E->atmosphere;
    AtmosphereParameters const *Ap = &E->parameters.atmosphere_parameters;
    VolatileParameters const *H2O = &Ap->H2O_volatile_parameters;

    PetscScalar dx1dt;

    dx1dt = get_dxdt( E, x1, A->dp1dx, H2O->kdist );

    return dx1dt;
}

PetscScalar get_initial_volatile( Ctx const *E, VolatileParameters const *V )
{
    /* initial volatile in the aqueous phase */
    Parameters           const *P  = &E->parameters;
    Constants            const *C  = &P->constants;
    Mesh                 const *M  = &E->mesh;

    PetscScalar fac, x;

    /* remember that V->henry contains another factor of VOLSCALE */

    fac = PetscSqr( P->radius );
    fac /= -P->gravity * M->mantle_mass;
    fac /= C->PRESSURE;

    fac *= C->VOLSCALE;
    fac /= PetscPowScalar(V->henry,V->henry_pow);

    x = newton( fac, V->henry_pow, V->initial );

    return x;
}


/////////////////////
/* Newton's method */
/////////////////////

/* for determining the initial mass fraction of volatiles in the
   melt.  The initial condition can be expressed as:

       x + A * x ** B = C */

static PetscScalar f( PetscScalar x, PetscScalar A, PetscScalar B, PetscScalar C )
{
    PetscScalar result;

    result = x + A * PetscPowScalar( x, B ) - C;

    return result;

}

static PetscScalar f_prim( PetscScalar x, PetscScalar A, PetscScalar B, PetscScalar C )
{
    PetscScalar result;

    result = 1.0 + A*B*PetscPowScalar( x, B-1.0 );

    return result;

}

static PetscScalar newton( PetscScalar A, PetscScalar B, PetscScalar xinit )
{
    PetscInt i=0;
    PetscScalar x;
    x = xinit;
    while(i < 50){
        x = x - f( x, A, B, xinit ) / f_prim( x, A, B, xinit );
        i++;
    }
    return x;

}
