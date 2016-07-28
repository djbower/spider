#include "bc.h"

static PetscScalar radiative_flux( PetscScalar );
static PetscScalar utbl_temp_drop( PetscScalar );

PetscErrorCode set_core_cooling( Ctx *E )
{
    PetscErrorCode    ierr;
    const PetscInt    ix0 = 0;        // index of first node
    PetscInt          ix;             // index of last basic node
    PetscInt          ix2;            // index of penultimate basic node
    PetscInt          numpts;
    PetscScalar       fac,radi,rado,vol,area1,area2;
    PetscMPIInt       rank, size;

    Mesh              *M; 

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_core_cooling:\n");CHKERRQ(ierr);
#endif
    M = &E->mesh;
    ierr = DMDAGetInfo(E->da_b,NULL,&numpts,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ix  = numpts-1;
    ix2 = numpts-2;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

    /* Assume that rank 0 contains the first point, and grab the 
       radius at the first point to broadcast */
    if (!rank) {
        ierr = VecGetValues( M->radius_b,1,&ix0,&rado);CHKERRQ(ierr);
    } else {
        rado = 0.0; // just for warning suppression
    }   
    ierr = MPI_Bcast(&rado,1,MPIU_SCALAR,0,PETSC_COMM_WORLD);CHKERRQ(ierr);

    /* Assume that the last rank contains the last two points */
    if (rank == size-1 ){
        ierr = VecGetValues( M->area_b,1,&ix,&area1);CHKERRQ(ierr);
        ierr = VecGetValues( M->area_b,1,&ix2,&area2);CHKERRQ(ierr);
        ierr = VecGetValues( M->radius_b,1,&ix,&radi);CHKERRQ(ierr);
        ierr = VecGetValues( M->volume_s,1,&ix2,&vol);CHKERRQ(ierr);
        fac = 4.0*M_PI*PetscPowScalar(rado, 3.0) * vol;
        fac *= RHO_CMB * CP_CMB;
        fac /= CP_CORE * MCORE * TFAC_CORE_AVG;
        fac += 1.0;
        fac = area2 / ( area1 * fac );
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
