#include "ic.h"

static PetscErrorCode make_super_adiabatic( Ctx *E ) 
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    PetscScalar       *arr_S_s,pres_b_last;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Vec               pres_b,pres_s;
    Mesh              *M;  
    Solution          *S;  
    PetscMPIInt       rank,size;
    DM                da_s=E->da_s, da_b=E->da_b;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    printf("make_super_adiabatic:\n");
#endif
    M = &E->mesh;
    S = &E->solution;
    pres_b = M->pressure_b;
    pres_s = M->pressure_s;
    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s; 

    /* Scatter the last value to all procs */
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (rank == size-1) { /* Assume that the last processor contains the last value */
      const PetscInt ix = NUMPTS-1; // Dangerous if PetscInt is not int!
      ierr = VecGetValues(pres_b,1,&ix,&pres_b_last);
#if (defined DEBUGOUTPUT)
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%d]   make_super_adiabatic: scattering value %f\n",rank,pres_b_last);CHKERRQ(ierr);
#endif
    }    
    MPI_Bcast(&pres_b_last,1,MPIU_SCALAR,size-1,PETSC_COMM_WORLD);

#if (defined DEBUGOUTPUT)
  ierr = PetscPrintf(PETSC_COMM_SELF,"[%d]   make_super_adiabatic: value of last point is %f\n",rank,pres_b_last);CHKERRQ(ierr);
#endif
    ierr = DMDAVecGetArray(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        arr_S_s[i] *= 1.0 + 0.01 * SUPERFAC * arr_pres_s[i]/pres_b_last;
    }
    ierr = DMDAVecRestoreArray(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_initial_condition( Ctx *E ) 
{
    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    printf("set_initial_condition:\n");
#endif
    S = &E->solution;

    ierr = VecSet(S->S_s,SINIT);CHKERRQ(ierr);

    make_super_adiabatic( E ); 

    PetscFunctionReturn(0);
}
