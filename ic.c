#include "ic.h"

static PetscErrorCode make_super_adiabatic( Ctx *, Vec );

PetscErrorCode set_initial_condition(Ctx *E, Vec S_in) 
{
# if (defined VERBOSE)
    PetscErrorCode ierr;
#endif

    PetscFunctionBeginUser;
#if (defined VERBOSE)
  ierr = PetscPrintf(PETSC_COMM_WORLD,"set_initial_condition : basing i.c. on a value of %f\n",S_init);CHKERRQ(ierr);
#endif
    /* Note that we don't initialize the vector to anything here, as we used to,
       as this is a pertubative approach and the values are set in the following 
       function */

    make_super_adiabatic( E, S_in ); 

    PetscFunctionReturn(0);
}

static PetscErrorCode make_super_adiabatic( Ctx *E, Vec S_in ) 
{
    /* linear increase of entropy with pressure to make the initial
       entropy profile slightly superadiabatic */

    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s,numpts_b;
    PetscScalar       *arr_S_s,pres_b_last;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Vec               pres_b,pres_s;
    Mesh              *M;  
    PetscMPIInt       rank,size;
    DM                da_s=E->da_s, da_b=E->da_b;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"make_super_adiabatic:\n");CHKERRQ(ierr);
#endif
    M = &E->mesh;
    pres_b = M->pressure_b;
    pres_s = M->pressure_s;
    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s; 
    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* Scatter the last value to all procs */
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (rank == size-1) { /* Assume that the last processor contains the last value */
      const PetscInt ix = numpts_b-1; // Dangerous if PetscInt is not int!
      ierr = VecGetValues(pres_b,1,&ix,&pres_b_last);
#if (defined DEBUGOUTPUT)
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%d]   make_super_adiabatic: scattering value %f\n",rank,pres_b_last);CHKERRQ(ierr);
#endif
    }    
    MPI_Bcast(&pres_b_last,1,MPIU_SCALAR,size-1,PETSC_COMM_WORLD);

#if (defined DEBUGOUTPUT)
  ierr = PetscPrintf(PETSC_COMM_SELF,"[%d]   make_super_adiabatic: value of last point is %f\n",rank,pres_b_last);CHKERRQ(ierr);
#endif
    ierr = DMDAVecGetArray(da_s,S_in,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        arr_S_s[i] = PERTURB * arr_pres_s[i]/pres_b_last;
    }
    ierr = DMDAVecRestoreArray(da_s,S_in,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
