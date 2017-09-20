#include "bc.h"
#include "util.h"

static PetscScalar grey_body( PetscScalar, Atmosphere *, AtmosphereParameters * );
static PetscScalar zahnle( PetscScalar,  AtmosphereParameters * );
static PetscScalar tsurf_param( PetscScalar, AtmosphereParameters * );
static PetscScalar hybrid( Ctx *, PetscScalar );

PetscErrorCode set_surface_flux( Ctx *E )
{
    PetscErrorCode       ierr;
    PetscMPIInt          rank;
    // initialise Qout    to avoid compiler warning
    PetscScalar          temp0, Tsurf, Qout=0.0;
    PetscInt             ind;
    Atmosphere           *A  = &E->atmosphere;
    Parameters           *P  = &E->parameters;
    AtmosphereParameters *Ap = &P->atmosphere_parameters;
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
        Tsurf = tsurf_param( temp0, Ap); // parameterised boundary layer
      }
      else{
        Tsurf = temp0; // surface temperature is potential temperature
      }

      /* determine flux */
      switch( Ap->MODEL ){
        case 1:
          // grey-body
          A->emissivity = Ap->emissivity0;
          Qout = grey_body( Tsurf, A, Ap );
          break;
        case 2:
          // zahnle
          Qout = zahnle( Tsurf, Ap );
          break;
        case 3:
          // atmosphere evolution
          set_emissivity_abe_matsui( A, Ap ); // updates A->emissivity
          Qout = grey_body( Tsurf, A, Ap );
          break;
      }

      /* TODO: for legacy purposes, perhaps to remove at some point */
      if( Ap->HYBRID ){
          Qout = hybrid( E, Qout );
      }

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

///////////////////////
/* atmosphere models */
///////////////////////

static PetscScalar grey_body( PetscScalar Tsurf, Atmosphere *A, AtmosphereParameters *Ap )
{
    PetscScalar Fsurf;

    Fsurf = PetscPowScalar(Tsurf,4.0)-PetscPowScalar(Ap->teqm,4.0);
    Fsurf *= Ap->sigma * A->emissivity; /* Note emissivity may vary */

    return Fsurf;
}

static PetscScalar zahnle( PetscScalar Tsurf, AtmosphereParameters *Ap )
{
    PetscScalar       Fsurf;

    /* fit to Zahnle et al. (1988) from Solomatov and Stevenson (1993)
       Eqn. 40 */

    /* FIXME: will break for non-dimensional
       see commit 780b1dd to reverse this */

    Fsurf = 1.5E2 + 1.02E-5 * PetscExpScalar(0.011*Tsurf);

    return Fsurf;

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

    /* it is of interest to know what the effective emissivity is,
       particularly for the hybrid method */
    /* TODO: not correct for parameterised UTBL */
    /* TODO: this overwrites the previous version of emissivity calculated
       from the grey-body model.  This could mess things up depending on the order
       of function calls */
    //A->emissivity = Q1 / ( PetscPowScalar(temp0,4.0) - PetscPowScalar(TEQM,4.0) );
    //A->emissivity /= SIGMA;

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
        ierr = VecGetValues( M->area_b,1,&ix,&area1);CHKERRQ(ierr);
        ierr = VecGetValues( M->area_b,1,&ix2,&area2);CHKERRQ(ierr);
        ierr = VecGetValues( M->volume_s,1,&ix2,&vol);CHKERRQ(ierr);
        ierr = VecGetValues( S->rho_s,1,&ix2,&rho_cmb);CHKERRQ(ierr);
        ierr = VecGetValues( S->cp_s,1,&ix2,&cp_cmb);CHKERRQ(ierr);
        fac = 4.0 * M_PI * vol;
        vol_core = 1.0/3.0 * PetscPowScalar(P->coresize,3.0) * PetscPowScalar(P->radius,3.0);
        fac = vol / vol_core;
        fac *= rho_cmb / P->rho_core;
        fac *= cp_cmb / P->cp_core;
        fac /= P->tfac_core_avg;

        fac = 1.0 / (1.0 + fac);
        fac *= area2 / area1;

        /* energy flux (Jtot) */
        ierr = VecGetValues(S->Jtot,1,&ix2,&val);CHKERRQ(ierr);
        val *= fac;
        ierr = VecSetValue(S->Jtot,ix,val,INSERT_VALUES);CHKERRQ(ierr);

        /* energy flow (Etot) */
        val *= area1;
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

static PetscScalar tsurf_param( PetscScalar temp, AtmosphereParameters *Ap )
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
