#include "bc.h"

static PetscScalar radiative_flux_with_dT( PetscScalar );
static PetscScalar radiative_flux( PetscScalar );
static PetscScalar utbl_temp_drop( PetscScalar );

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
    PetscScalar       fac,vol,area1,area2,rho_cmb,cp_cmb;
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
        fac *= rho_cmb * cp_cmb;
        fac /= CP_CORE * MCORE * TFAC_CORE_AVG;
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

static PetscScalar utbl_temp_drop( PetscScalar temp )
{
    PetscScalar dT;

    dT = CONSTBC * PetscPowScalar( temp, EXPBC );

    return dT;
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
    PetscScalar dT, temp0, flux;

    dT = utbl_temp_drop( temp );
    temp0 = temp - dT;
    flux = radiative_flux( temp0 );

    return flux;
}
