#include "bc.h"

static PetscScalar radiative_flux( PetscScalar );
static PetscScalar utbl_temp_drop( PetscScalar );

PetscErrorCode get_core_cooling( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          ix;             // index of last basic node
    PetscInt          ix2;            // index of penultimate basic node
    PetscInt          numpts_b;
    PetscScalar       fac,vol,area1,area2,rho_cmb,cp_cmb;
    PetscMPIInt       rank, size;

    Mesh              *M; 
    Solution          *S;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"get_core_cooling:\n");CHKERRQ(ierr);
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
        fac *= rho_cmb * cp_cmb;
        fac /= CP_CORE * MCORE * TFAC_CORE_AVG;
        fac = 1.0 / (1.0 + fac);
        fac *= area2 / area1;
        E->BC_BOT_FAC = fac;
    }

    PetscFunctionReturn(0);

}

static PetscScalar utbl_temp_drop( PetscScalar temp )
{
    PetscScalar dT;

    dT = CONSTBC * PetscPowScalar( temp, EXPBC );

    return dT;
}

PetscScalar radiative_flux_with_dT( PetscScalar temp )
{
    PetscScalar dT, temp0, flux;

    dT = utbl_temp_drop( temp );
    temp0 = temp - dT;
    flux = radiative_flux( temp0 );

    return flux;
}

static PetscScalar radiative_flux( PetscScalar temp )
{
    PetscScalar flux;

    flux = PetscPowScalar(temp,4.0)-PetscPowScalar(TEQM,4.0);
    flux *= SIGMA * EMISSIVITY;

    return flux;
}
