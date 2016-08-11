#include "util.h"

PetscScalar average( PetscScalar a, PetscScalar b ) 
{
    /* arithmetic average of two numbers */

    PetscScalar out;

    out = 0.5 * (a+b);

    return out;
}

PetscScalar combine_matprop( PetscScalar weight, PetscScalar mat1, PetscScalar mat2 )
{
    /* linear weighting of two quantities */

    PetscScalar out;

    out = weight * mat1 + (1.0-weight) * mat2;

    return out;

}

/* this is still needed for fusion curve derivative */
PetscErrorCode d_dr( Ctx *E, Vec in_s, Vec out_b )
{

    /* use staggered nodes to compute spatial derivative at basic
       (internal) nodes.  This is 2nd order accurate for a mesh with
       constant spacing */

    PetscErrorCode     ierr;
    PetscInt           i,ilo_b,ihi_b,w_b,ilo,ihi,numpts_b;
    DM                 da_s=E->da_s,da_b=E->da_b;
    PetscScalar        dr,*arr_out_b;
    const PetscScalar *arr_in_s;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"d_dr\n");CHKERRQ(ierr);
#endif
    dr = E->mesh.dx_s;

    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ilo = ilo_b==0        ? 1          : ilo_b;
    ihi = ihi_b==numpts_b ? numpts_b-1 : ihi_b;

    // TODO: here and elsewhere, we are a little glib about assuming things abou the way the DA's partition things. We should introduce checks for any function which involves both DAs at once, that the expected ranges apply.

    /* Scatter to local vectors, since we may require potentially off-processor ghost values */
    ierr = DMGlobalToLocalBegin(da_s,in_s,INSERT_VALUES,E->work_local_s);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_s,in_s,INSERT_VALUES,E->work_local_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,E->work_local_s,&arr_in_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,out_b,&arr_out_b);CHKERRQ(ierr);
    for(i=ilo; i<ihi; i++){
        arr_out_b[i] = 1.0/dr * ( arr_in_s[i]-arr_in_s[i-1] );
    }
    ierr = DMDAVecRestoreArrayRead(da_s,E->work_local_s,&arr_in_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,out_b,&arr_out_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

//PetscErrorCode set_d_dr2( Ctx *E )
//{
    /* create a matrix containing the 2nd order accurate stencil for
       computing the 1st order derivative of a quantity at the 
       staggered nodes

       d/dr will be given by MatMult( A, x, y )

       where:
           x: (Vec) input: quantity at staggered nodes (size numpts_s)
           y: (Vec) output: d/dr at staggered nodes (size numpts_s)

       currently this matrix is global */

/*    PetscErrorCode ierr;
    PetscInt i, numpts_s, col[3], rstart, rend, ilo_s, ihi_s, w_s;
    PetscScalar dr, value[3];
    Mat A;

    PetscFunctionBeginUser;

    dr = E->mesh.dx_s;

    ierr = DMDAGetCorners(E->da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

#if (defined DEBUGOUTPUT)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_d_dr2 : creating a matrix of size %D\n",numpts_s);CHKERRQ(ierr);
    {
       PetscMPIInt rank;
       ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_WORLD,"[%D] set_d_dr2 : local rows %d<=i<%d\n",rank,ilo_s,ihi_s);CHKERRQ(ierr);
    }
#endif

    ierr = MatCreate( PETSC_COMM_WORLD, &E->d_dr2 );CHKERRQ(ierr);
    A = E->d_dr2;
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,numpts_s,numpts_s);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);*/

    /* Set first and last row values, if appropriate */
    /*if (!ilo_s) {
      i = 0;
      col[0]=0; col[1]=1; col[2]=2;
      value[0]=-3.0/2.0; value[1]=2.0; value[2]=-0.5;
      ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }
    if (ihi_s == numpts_s) {
      i = numpts_s-1;
      col[0]=numpts_s-3; col[1]=numpts_s-2, col[2]=numpts_s-1;
      value[0]=0.5; value[1]=-2.0; value[2]=3.0/2.0;
      ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }*/

    /* set values corresponding to the mesh interior */
    /*rstart = (ilo_s == 0      )  ? 1          : ilo_s;
    rend   = (ihi_s == numpts_s) ? numpts_s-1 : ihi_s;
    value[0]=-0.5; value[1]=0.0; value[2]=0.5;
    for (i=rstart; i<rend; i++) {
        col[0] = i-1; col[1] = i; col[2] = i+1;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }*/

    /* Assemble the matrix */
    /*ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);*/

    /* Must occur after assembly */
    /*ierr = MatScale( A, 1.0/dr);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}*/
