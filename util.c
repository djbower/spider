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
    PetscInt           i,ilo_b,ihi_b,w_b,ilo,ihi;
    DM                 da_s=E->da_s,da_b=E->da_b;
    PetscScalar        dr,*arr_out_b;
    const PetscScalar *arr_in_s;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    {
      PetscMPIInt rank;
      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] d_dr\n");CHKERRQ(ierr);
    }
#endif
    dr = E->mesh.dx_s;

    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ilo = ilo_b==0      ? 1        : ilo_b;
    ihi = ihi_b==NUMPTS ? NUMPTS-1 : ihi_b;

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

PetscErrorCode set_d_dr2( Ctx *E )
{
    /* create a matrix containing the 2nd order accurate stencil for
       computing the 1st order derivative of a quantity at the 
       staggered nodes

       d/dr will be given by MatMult( A, x, y )

       where:
           x: (Vec) input: quantity at staggered nodes (size NUMPTSS)
           y: (Vec) output: d/dr at staggered nodes (size NUMPTSS)

       currently this matrix is global */

    PetscErrorCode ierr;
    PetscInt i, n=NUMPTSS, col[3], rstart, rend;
    PetscScalar dr, value[3];
    Mat A;

    PetscFunctionBeginUser;

    dr = E->mesh.dx_s;

    ierr = MatCreate( PETSC_COMM_WORLD, &A );CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,NUMPTSS,NUMPTSS);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);

    //if (!rstart) {
    rstart = 1; i = 0;
    col[0]=0; col[1]=1; col[2]=2;
    value[0]=-3.0/2.0; value[1]=2.0; value[2]=-0.5;
    ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    //}
    //if (rend == n) {
    rend = n-1; i =n-1;
    col[0]=n-3; col[1]=n-2, col[2]=n-1;
    value[0]=0.5; value[1]=-2.0; value[2]=3.0/2.0;
    ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    //}

    /* set values corresponding to the mesh interior */
    value[0]=-0.5; value[1]=0.0; value[2]=0.5;
    for (i=rstart; i<rend; i++) {
        col[0] = i-1; col[1] = i; col[2] = i+1;
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }

    /* Assemble the matrix */
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    /* Must occur after assembly */
    ierr = MatScale( A, 1.0/dr);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
