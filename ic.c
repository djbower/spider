#include "ic.h"

static PetscErrorCode make_ic_from_adiabat( Ctx *, Vec );
static PetscErrorCode make_ic_from_melt_fraction( Ctx *, Vec );
static PetscErrorCode make_super_adiabatic( Ctx *, Vec );

PetscErrorCode set_initial_condition(Ctx *E, Vec S_in) 
{

    PetscFunctionBeginUser;

    /* ic from adiabat */
    //make_ic_from_adiabat( E, S_in );

    /* ic from melt fraction */
    make_ic_from_melt_fraction( E, S_in );

    PetscFunctionReturn(0);
}

static PetscErrorCode make_ic_from_adiabat( Ctx *E, Vec S_in )
{

    PetscErrorCode    ierr;
    PetscScalar       dS;

    PetscFunctionBeginUser;

    /* perturbation relative to reference entropy, which is defined
       at a non-dimensional entropy of 1.0 */
    dS = E->S_init - 1.0;
    ierr = VecSet( S_in, dS ); CHKERRQ( ierr );

    /* now add small gradient to initial perturbed value */
    make_super_adiabatic( E, S_in );

    PetscFunctionReturn(0);
}

static PetscErrorCode make_ic_from_melt_fraction( Ctx *E, Vec S_in )
{

    PetscErrorCode    ierr;
    PetscScalar       meltf, val, maxval, step;
    PetscInt          i, ilo;
    Vec               liq_s, dfus_s;
    Solution          *S;
    PetscScalar       *arr;

    PetscFunctionBeginUser;

    S = &E->solution;
    dfus_s = S->fusion_s;
    liq_s = S->liquidus_s;

    /* melt fraction contour to follow */
    meltf = 0.99;
    val = 1.0 - meltf;

    ierr = VecCopy( liq_s, S_in ); CHKERRQ(ierr);
    ierr = VecAXPY( S_in, -val, dfus_s ); CHKERRQ(ierr);
    ierr = VecShift( S_in, -1.0 ); CHKERRQ(ierr); // (1.0 is reference)

    /* find overturn point and set everything to right to overturn value */
    ierr = VecMax( S_in, &ilo, &maxval );
    ierr = VecGetArray( S_in, &arr );
    for(i=ilo; i<NUMPTS_S_DEFAULT; ++i){
        arr[i] = maxval;
    }

    /* entropy drop.  Everything to left cannot drop by more than this
       value relative to the overturn value */
    step = 0.0075;
    for(i=0; i<ilo; ++i){
        if(arr[i] < maxval-step){
            arr[i] = maxval-step;
        }
    }

    ierr = VecRestoreArray( S_in, &arr );

    /* now add small gradient to initial perturbed value */
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
      const PetscInt ix = numpts_b-1; // FIXME: Dangerous if PetscInt is not int!
      ierr = VecGetValues(pres_b,1,&ix,&pres_b_last);
#if (defined DEBUGOUTPUT)
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%d]   make_super_adiabatic: scattering value %f\n",rank,(double)pres_b_last);CHKERRQ(ierr);
#endif
    }    
    MPI_Bcast(&pres_b_last,1,MPIU_SCALAR,size-1,PETSC_COMM_WORLD);

#if (defined DEBUGOUTPUT)
  ierr = PetscPrintf(PETSC_COMM_SELF,"[%d]   make_super_adiabatic: value of last point is %f\n",rank,(double)pres_b_last);CHKERRQ(ierr);
#endif
    ierr = DMDAVecGetArray(da_s,S_in,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        arr_S_s[i] += PERTURB * arr_pres_s[i]/pres_b_last;
    }
    ierr = DMDAVecRestoreArray(da_s,S_in,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
