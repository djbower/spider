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

PetscErrorCode set_d_dr( Ctx *E )
{
    /* interpolate quadratic functions to get val and dval/dr at
       basic nodes.  For overlapping regions can get an estimate
       from the left and an estimate from the right, which are
       combined by the arithmetic average

       Y     = Ax**2 + Bx + C
       dY/dr = 2Ax + B

    */

    PetscErrorCode ierr;
    PetscInt i, col[3], rstart, rend;
    PetscInt numpts_s, ilo_s, ihi_s, w_s;
    PetscInt numpts_b, ilo_b, ihi_b, w_b;
    PetscScalar value[3];
    PetscScalar *arr_radius_s, *arr_radius_b, dh, h1, h2;
    Mat A, B, C, S1a, S1b;
    Mesh *M;
    Vec dx, dxsq;

    PetscFunctionBeginUser;

    M = &E->mesh;

    // staggered
    ierr = DMDAGetCorners(E->da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);

    // basic
    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);

    /* the first and last row are all zeros (implicitly, I think),
       to enable us to return a vector with length numpts_b */

    // FIXME TODO this initialises the 'final' matrix
    //ierr = MatCreate( PETSC_COMM_WORLD, &E->d_dr );CHKERRQ(ierr);

    // initialise A coefficient matrix
    ierr = MatCreate( PETSC_COMM_WORLD, &A );CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,numpts_b,numpts_s);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);

    // initialise B coefficient matrix
    ierr = MatCreate( PETSC_COMM_WORLD, &B );CHKERRQ(ierr);
    ierr = MatSetSizes(B,PETSC_DECIDE,PETSC_DECIDE,numpts_b,numpts_s);CHKERRQ(ierr);
    ierr = MatSetFromOptions(B);CHKERRQ(ierr);
    ierr = MatSetUp(B);CHKERRQ(ierr);

    // initialise C coefficient matrix
    ierr = MatCreate( PETSC_COMM_WORLD, &C );CHKERRQ(ierr);
    ierr = MatSetSizes(C,PETSC_DECIDE,PETSC_DECIDE,numpts_b,numpts_s);CHKERRQ(ierr);
    ierr = MatSetFromOptions(C);CHKERRQ(ierr);
    ierr = MatSetUp(C);CHKERRQ(ierr);

    // now build coefficient matrices

    /* set values corresponding to the mesh interior */
    rstart = (ilo_s == 0      )  ? 1          : ilo_s;
    rend   = (ihi_s == numpts_s) ? numpts_s-1 : ihi_s;
    for (i=rstart; i<rend; i++) {
        col[0]=i-1; col[1]=i; col[2]=i+1;
        h1 = arr_radius_s[i] - arr_radius_s[i-1];
        h2 = arr_radius_s[i+1] - arr_radius_s[i];
        // A
        value[0] = 1.0 / (h1*(h1+h2));
        value[1] = -1.0 / (h1*h2);
        value[2] = 1.0 / (h2*(h1+h2));
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
        // B
        value[0] = (-2.0*h1-h2)/(h1*(h1+h2));
        value[1] = (h1+h2)/(h1*h2);
        value[2] = -h1/(h2*(h1+h2));
        ierr = MatSetValues(B,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
        // C
        value[0] = 1.0;
        value[1] = 0.0;
        value[2] = 0.0;
        ierr = MatSetValues(C,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }

    /* set last row values, if appropriate */
    /* use backward difference for last point, to retain ability to dot
       product with half the basic node spacing when sequentially stepping
       forward */
    if (ihi_s == numpts_s) {
        i = numpts_s-1;
        col[0]=numpts_s-3; col[1]=numpts_s-2, col[2]=numpts_s-1;
        h1 = arr_radius_s[rend-1] - arr_radius_s[rend-2];
        h2 = arr_radius_s[rend] - arr_radius_s[rend-1];
        // A
        value[0] = 1.0 / (h1*(h1+h2)); // same A as above
        value[1] = -1.0 / (h1*h2); // same as A above
        value[2] = 1.0 / (h2*(h1+h2)); // same as A above
        ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
        // B
        value[0] = -h2 / (h1*(h1+h2));
        value[1] = (h2-h1)/(h1*h2);
        value[2] = h1 / (h2*(h1+h2));
        ierr = MatSetValues(B,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
        // C
        value[0] = 0.0;
        value[1] = 1.0;
        value[2] = 0.0;
        ierr = MatSetValues(C,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }

    /* this vector gives the distance of the basic node from the
       nearest staggered node to its left */

    ierr = VecCreate( PETSC_COMM_WORLD, &dx );CHKERRQ(ierr);
    ierr = VecSetSizes( dx, PETSC_DECIDE, numpts_b );CHKERRQ(ierr);
    ierr = VecSetFromOptions( dx );CHKERRQ(ierr);
    ierr = VecSetUp( dx );CHKERRQ(ierr);

    /* set values corresponding to the vector interior */
    /* remember that top and bottom entries are always zero */
    rstart = (ilo_b == 0      )  ? 1          : ilo_b;
    rend   = (ihi_b == numpts_b) ? numpts_b-1 : ihi_b;
    for (i=rstart; i<rend; i++) {
        dh = 0.5 * (arr_radius_b[i] - arr_radius_b[i-1]);
        VecSetValues( dx, 1, &i, &dh, INSERT_VALUES );CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(dx);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(dx);CHKERRQ(ierr);

    // new vector with square of dx
    ierr = VecDuplicate(dx,&dxsq);CHKERRQ(ierr);
    ierr = VecCopy(dx,dxsq);CHKERRQ(ierr);
    ierr = VecPow(dxsq,2.0);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayRead(E->da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(E->da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);

    /* assemble the matrices */
    // A
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    // B
    ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    // C
    ierr = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    /* for computing first estimate of S */
    // S1a: ax**2 term
    ierr = MatDuplicate( A, MAT_COPY_VALUES, &S1a );CHKERRQ(ierr);
    ierr = MatDiagonalScale( S1a, dxsq, NULL); CHKERRQ(ierr);
    // S1b: bx term
    ierr = MatDuplicate( B, MAT_COPY_VALUES, &S1b );CHKERRQ(ierr);
    ierr = MatDiagonalScale( S1b, dx, NULL); CHKERRQ(ierr);
    // add (final matrix is S1a)
    ierr = MatAXPY( S1a, 1.0, S1b, SAME_NONZERO_PATTERN ); CHKERRQ(ierr); 
    ierr = MatAXPY( S1a, 1.0, C, SAME_NONZERO_PATTERN ); CHKERRQ(ierr);





    /* clean up temporary matrices and vectors */
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = MatDestroy(&B); CHKERRQ(ierr);
    ierr = MatDestroy(&C); CHKERRQ(ierr);
    ierr = MatDestroy(&S1a); CHKERRQ(ierr);
    ierr = MatDestroy(&S1b); CHKERRQ(ierr);
    ierr = VecDestroy(&dx); CHKERRQ(ierr);
    ierr = VecDestroy(&dxsq); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/*PetscErrorCode d_dr_regular( Ctx *E, Vec in_s, Vec out_b )
{*/

    /* use staggered nodes to compute spatial derivative at basic
       (internal) nodes.  This is 2nd order accurate for a mesh with
       constant spacing */

   /* PetscErrorCode     ierr;
    PetscInt           i,ilo_b,ihi_b,w_b,ilo,ihi,numpts_b;
    DM                 da_s=E->da_s,da_b=E->da_b;
    PetscScalar        dr,*arr_out_b;
    const PetscScalar *arr_in_s;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"d_dr\n");CHKERRQ(ierr);
#endif*/
    // TODO FIXME
    //dr = E->mesh.dx_s;

    /*ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ilo = ilo_b==0        ? 1          : ilo_b;
    ihi = ihi_b==numpts_b ? numpts_b-1 : ihi_b;*/

    // TODO: here and elsewhere, we are a little glib about assuming things abou the way the DA's partition things. We should introduce checks for any function which involves both DAs at once, that the expected ranges apply.

    /* Scatter to local vectors, since we may require potentially off-processor ghost values */
    /*ierr = DMGlobalToLocalBegin(da_s,in_s,INSERT_VALUES,E->work_local_s);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_s,in_s,INSERT_VALUES,E->work_local_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,E->work_local_s,&arr_in_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,out_b,&arr_out_b);CHKERRQ(ierr);
    for(i=ilo; i<ihi; i++){
        arr_out_b[i] = 1.0/dr * ( arr_in_s[i]-arr_in_s[i-1] );
    }
    ierr = DMDAVecRestoreArrayRead(da_s,E->work_local_s,&arr_in_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,out_b,&arr_out_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}*/

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
