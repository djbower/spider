#include "bc.h"
#include "util.h"

static PetscScalar radiative_flux_with_dT( PetscScalar );
static PetscScalar radiative_flux( PetscScalar );
static PetscScalar tsurf_param( PetscScalar );

PetscErrorCode set_surface_flux( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscMPIInt       rank;
    PetscScalar       temp_s_0, phi_s_0, Q1, Q2, Qout, fwt;
    PetscInt          ind;
    PetscScalar       G0, R0, R1, R2, E0, E1, E2;
    Mesh              *M; 
    Solution          *S;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_surface_flux:\n");CHKERRQ(ierr);
#endif
    M = &E->mesh;
    S = &E->solution;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

    if (!rank){
      ind = 0;
      /* could possibly formulate below using surface values, as now we
         extrapolate to get values at the surface */
      ierr = VecGetValues(S->phi_s,1,&ind,&phi_s_0);CHKERRQ(ierr);
      /* SWIDTH or PHI_WIDTH most appropriate choice here? */
      fwt = tanh_weight( phi_s_0, PHI_CRITICAL, SWIDTH );

      // radiative surface
      ierr = VecGetValues(S->temp_s,1,&ind,&temp_s_0);CHKERRQ(ierr);
      Q1 = radiative_flux_with_dT( temp_s_0 );

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

      Qout = Q1 * fwt + Q2 * (1.0 - fwt);

      /* Note - below is true because non-dim geom is exactly equal
         to one, so do not need to multiply by area of surface */
      ierr = VecSetValue(S->Etot,0,Qout,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(S->Jtot,0,Qout,INSERT_VALUES);CHKERRQ(ierr);

    }

    ierr = VecAssemblyBegin(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);

    PetscFunctionReturn(0);

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
        /* previous value used when constant */
        //rho_cmb = 1.1922543609124883;
        //cp_cmb = 0.40093215388392944;
        // gives a fac of:
        // 1.00069301288
        vol_core = 1.0/3.0 * PetscPowScalar(RADIN,3.0);
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

    ierr = VecAssemblyBegin(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

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

static PetscScalar radiative_flux( PetscScalar temp )
{
    PetscScalar flux;

    flux = PetscPowScalar(temp,4.0)-PetscPowScalar(TEQM,4.0);
    flux *= SIGMA * EMISSIVITY;

    return flux;
}


static PetscScalar radiative_flux_with_dT( PetscScalar temp )
{
    PetscScalar Ts, flux;

    if( CONSTBC == 0.0 ){
        Ts = temp;
    }
    else{
        Ts = tsurf_param( temp );
    }
    flux = radiative_flux( Ts );

    return flux;
}
