#include "bc.h"
#include "util.h"

static PetscScalar hybrid( Ctx *, PetscScalar, PetscReal );
#if defined GREYBODY
static PetscScalar greybody_with_dT( PetscScalar, PetscReal );
static PetscScalar greybody( PetscScalar, PetscReal );
static PetscScalar tsurf_param( PetscScalar );
#endif
#ifdef ZAHNLE
static PetscScalar zahnle( PetscScalar, PetscReal );
#endif
#ifdef HAMANO
static PetscScalar hamano( PetscScalar, PetscReal );
#endif

PetscErrorCode set_surface_flux( Ctx *E, PetscReal tyrs )
{
    PetscErrorCode    ierr;
    PetscMPIInt       rank;
    PetscScalar       temp0, Qout;
    PetscInt          ind;
    Solution          *S = &E->solution;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_surface_flux:\n");CHKERRQ(ierr);
#endif
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

    if (!rank){
      /* get surface values to pass to flux boundary condition
         function that is chosen by user */
      ind = 0;

      /* surface temperature */
      ierr = VecGetValues(S->temp,1,&ind,&temp0); CHKERRQ(ierr);

      Qout = hybrid( E, temp0, tyrs );

      ierr = VecSetValue(S->Jtot,0,Qout,INSERT_VALUES);CHKERRQ(ierr);
      Qout *= PetscSqr( RADIUS );
      ierr = VecSetValue(S->Etot,0,Qout,INSERT_VALUES);CHKERRQ(ierr);

    }

    ierr = VecAssemblyBegin(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscScalar hybrid( Ctx *E, PetscScalar temp0, PetscReal tyrs )
{
    PetscScalar    Q1;
#ifdef HYBRID
    PetscErrorCode ierr;
    PetscScalar    G0, R0, R1, R2, E0, E1, E2, Q2, fwt, phi0;
    PetscInt       ind;
    Mesh           *M = &E->mesh;
    Solution       *S = &E->solution;
#endif

#ifdef HAMANO
    Q1 = hamano( temp0, tyrs );
#endif
#ifdef ZAHNLE
    Q1 = zahnle( temp0, tyrs );
#endif
#ifdef GREYBODY
    Q1 = greybody_with_dT( temp0, tyrs );
#endif
#ifdef HYBRID
    /* for weight of different fluxes */
    ind = 0;
    ierr = VecGetValues(S->phi,1,&ind,&phi0); CHKERRQ(ierr);
    /* SWIDTH or PHI_WIDTH most appropriate choice here? */
    fwt = tanh_weight( phi0, PHI_CRITICAL, SWIDTH );

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
    Q2 = E0 / G0; // G0 should be 1.0 by definition
    Q1 = Q1 * fwt + Q2 * (1.0 - fwt);
#endif

    return Q1;

}

#ifdef ZAHNLE
static PetscScalar zahnle( PetscScalar Tsurf, PetscReal tyrs )
{
    PetscScalar       Fsurf;

    /* fit to Zahnle et al. (1988) from Solomatov and Stevenson (1993)
       Eqn. 40 */

    Fsurf = 1.5E2 + 1.02E-5 * PetscExpScalar(0.011*Tsurf);

    return Fsurf;

}
#endif

#ifdef HAMANO
static PetscScalar hamano( PetscScalar Tsurf, PetscReal tyrs )
{
    PetscScalar       Fsurf;

    /* add code here */
    /* build your own function of surface temperature (Tsurf) and
       time in years (tyrs) */
    Fsurf = 0.0;

    return Fsurf;

}
#endif

PetscErrorCode set_core_mantle_flux( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          ix;             // index of last basic node
    PetscInt          ix2;            // index of penultimate basic node
    PetscInt          numpts_b;
    PetscScalar       fac,vol,vol_core,area1,area2,rho_cmb,cp_cmb;
    PetscScalar       val;
    PetscMPIInt       rank, size;

    Mesh              *M; 
    Solution          *S;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_core_mantle_flux:\n");CHKERRQ(ierr);
#endif
    M = &E->mesh;
    S = &E->solution;
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
        vol_core = 1.0/3.0 * PetscPowScalar(CORESIZE,3.0) * PetscPowScalar(RADIUS,3.0);
        fac = vol / vol_core;
        fac *= rho_cmb / RHO_CORE;
        fac *= cp_cmb / CP_CORE;
        fac /= TFAC_CORE_AVG;

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

#if defined(GREYBODY)
static PetscScalar tsurf_param( PetscScalar temp )
{
    PetscScalar Ts, c, fac, num, den;
    c = CONSTBC;

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

static PetscScalar greybody( PetscScalar Tsurf, PetscReal tyrs )
{
    PetscScalar Fsurf;

    Fsurf = PetscPowScalar(Tsurf,4.0)-PetscPowScalar(TEQM,4.0);
    Fsurf *= SIGMA * EMISSIVITY;

    return Fsurf;
}

static PetscScalar greybody_with_dT( PetscScalar Tsurf, PetscReal tyrs )
{
    PetscScalar Ts, Fsurf;

    /* no parameterised ultra-thin thermal boundary layer */
    if( CONSTBC == 0.0 ){
        Ts = Tsurf;
    }
    /* parameterised ultra-thin thermal boundary layer */
    else{
        Ts = tsurf_param( Tsurf );
    }
    Fsurf = greybody( Ts, tyrs );

    return Fsurf;
}
#endif
