#include "util.h"

static PetscErrorCode set_d_dr_linear( Ctx * );
//static PetscErrorCode set_d_dr_quadratic( Ctx * );

PetscErrorCode set_entropy_from_solution( Ctx *E, Vec sol )
{
    /* Set entropy at the basic (S_b) and staggered (S_s) nodes
       as well as the entropy gradient at the basic nodes using 
       the solution fields SPIDER_SOLUTION_FIELD_S0 and
       SPIDER_SOLUTION_FIELD_DSDR_B */

    PetscErrorCode ierr;
    Mesh           *M = &E->mesh;
    Solution       *S = &E->solution;
    PetscScalar    S0;
    PetscScalar    *arr_S_b, *arr_S_s, *arr_dSdr_b, *arr_radius_s, *arr_radius_b;
    PetscInt       i, ihi_b, ilo_b, w_b;
    PetscMPIInt    size;
    DM             da_s = E->da_s, da_b=E->da_b;
    Vec            *subVecs;

    const PetscInt ind0 = 0;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    /* set S->dS/dr at basic nodes */
    ierr = VecCopy( subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDR_B]], S->dSdr );CHKERRQ(ierr);
    /* get first staggered node value (store as S0) */
    ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],1,&ind0,&S0);CHKERRQ(ierr);

    /* for looping over basic nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    ierr = DMDAVecGetArray(da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dSdr,&arr_dSdr_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);

    /* S (absolute) at all staggered and basic internal nodes */
    arr_S_s[0] = 0.0;
    /* note plus one for start of loop */
    for(i=ilo_b+1; i<ihi_b-1; ++i){
      /* S (absolute) at staggered nodes */
      arr_S_s[i] = arr_dSdr_b[i] * (arr_radius_s[i] - arr_radius_s[i-1] );
      arr_S_s[i] += arr_S_s[i-1]; // dS relative to first staggered value
      arr_S_b[i] = arr_dSdr_b[i] * 0.5 * (arr_radius_b[i] - arr_radius_b[i-1] );
      arr_S_b[i] += arr_S_s[i-1];
      arr_S_s[i-1] += S0; // add at end to try and retain precision
      arr_S_b[i-1] += S0; // add at end to try and retain precision
    }
    /* loop above terminates before we have added the constant offset
       to the last staggered node value.  So do that here */
    arr_S_s[ihi_b-2] += S0; // add at end to try and retain precision
    arr_S_b[ihi_b-2] += S0; // add at end to try and retain precision

    /* now deal with top and bottom surfaces (outermost basic nodes)
       use gradients to give estimates of entropy at the top
       and bottom surfaces by extrapolation.  Remember that S0 has already
       been included by the loop above over the basic internal nodes */
    arr_S_b[0] = -arr_dSdr_b[1] * 0.5 * (arr_radius_b[1] - arr_radius_b[0]);
    arr_S_b[0] += arr_S_s[0];
    arr_S_b[ihi_b-1] = arr_dSdr_b[ihi_b-2] * 0.5 * (arr_radius_b[ihi_b-1]-arr_radius_b[ihi_b-2]);
    arr_S_b[ihi_b-1] += arr_S_s[ihi_b-2];

    ierr = DMDAVecRestoreArray(da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dSdr,&arr_dSdr_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    PetscFree(subVecs);

    PetscFunctionReturn(0);
}

PetscErrorCode set_partial_pressures_from_solution( Ctx *E, Vec sol )
{
    PetscErrorCode             ierr;
    Atmosphere                 *A = &E->atmosphere;
    Parameters                 *P = &E->parameters;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    PetscInt                    i;
    PetscMPIInt                 size;
    Vec                        *subVecs;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    /* partial pressures */
    for( i=0; i<Ap->n_volatiles; ++i) {
        ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],1,&i,&A->volatiles[i].p);CHKERRQ(ierr);
    }

    /* mass reaction terms */
    for( i=0; i<Ap->n_reactions; ++i) {
        ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]],1,&i,&A->mass_reaction[i]);CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    PetscFree(subVecs);

    PetscFunctionReturn(0);

}

PetscErrorCode set_solution_from_partial_pressures( Ctx *E, Vec sol )
{
    PetscErrorCode             ierr;
    Parameters                 *P = &E->parameters;
    Atmosphere                 *A = &E->atmosphere;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    PetscInt                   i;
    Vec                        *subVecs;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    /* partial pressures */
    for( i=0; i<Ap->n_volatiles; ++i) {
        ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],i,A->volatiles[i].p,INSERT_VALUES);CHKERRQ(ierr);
    }

    /* mass reaction terms */
    for( i=0; i<Ap->n_reactions; ++i) {
        ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]],i,A->mass_reaction[i],INSERT_VALUES);CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]]);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]]);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    PetscFree(subVecs);

    PetscFunctionReturn(0);

}

PetscErrorCode set_solution_from_entropy_at_staggered_nodes( Ctx *E, Vec sol )
{
    /* Set the solution fields SPIDER_SOLUTION_FIELD_S0 and
       SPIDER_SOLUTION_FIELD_DSDR_B using the entropy at
       staggered nodes (S_s) */

    PetscErrorCode ierr;
    Mesh           *M = &E->mesh;
    Solution       *S = &E->solution;
    PetscScalar    S0;
    PetscScalar    *arr_dSdr_b, *arr_S_s, *arr_radius_s;
    PetscInt       i, ihi_b, ilo_b, w_b;
    PetscMPIInt    size;
    DM             da_s = E->da_s, da_b=E->da_b;
    Vec            dSdr_b;
    Vec            *subVecs;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    dSdr_b = subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDR_B]];

    /* for looping over basic nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    ierr = DMDAVecGetArray(da_b,dSdr_b,&arr_dSdr_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);

    arr_dSdr_b[0] = 0.0; // first point constrained by boundary conditions
    arr_dSdr_b[ihi_b-1] = 0.0; // last point constrained by boundary conditions
    for(i=ilo_b+1; i<ihi_b-1; ++i){
        arr_dSdr_b[i] = arr_S_s[i] - arr_S_s[i-1];
        arr_dSdr_b[i] /= arr_radius_s[i] - arr_radius_s[i-1];
    }

    /* set entropy at top staggered node */
    S0 = arr_S_s[0];
    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],0,S0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]]);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArray(da_b,dSdr_b,&arr_dSdr_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);
    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    ierr = PetscFree(subVecs);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscScalar combine_matprop( PetscScalar weight, PetscScalar mat1, PetscScalar mat2 )
{
    /* linear weighting of two quantities */

    PetscScalar out;

    out = weight * mat1 + (1.0-weight) * mat2;

    return out;
}

PetscScalar tanh_weight( PetscScalar qty, PetscScalar threshold, PetscScalar width )
{
    /* tanh weight for viscosity profile and smoothing */

    PetscScalar    fwt, z;

    z = ( qty - threshold ) / width;
    fwt = 0.5 * ( 1.0 + PetscTanhScalar( z ) );
    return fwt;
}

PetscErrorCode set_d_dr( Ctx *E )
{

    PetscFunctionBeginUser;

    /* FIXME: can manually select derivative order here */

    set_d_dr_linear( E );
    //set_d_dr_quadratic( E );

    PetscFunctionReturn(0);
}


static PetscErrorCode set_d_dr_linear( Ctx *E )
{

    /* piece-wise linear interpolation between neighbouring
       staggered nodes to get the basic node value and the
       derivative at the basic node */

    PetscErrorCode ierr;
    PetscInt i, rstart, rend, col[2];
    PetscInt numpts_s, ilo_s, ihi_s, w_s;
    PetscInt numpts_b, ilo_b, w_b; // ihi_b, w_b;
    PetscScalar value[2], dh, dh2;
    PetscScalar *arr_radius_s, *arr_radius_b;
    Mat A, B;
    Mesh *M;

    PetscFunctionBeginUser;

    M = &E->mesh;

    // staggered
    ierr = DMDAGetCorners(E->da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);

    // basic
    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    //ihi_b = ilo_b + w_b;
    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);

    /* the first and last row are all zeros (implicitly, I think),
       to enable us to return a vector with length numpts_b */

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

    // add coefficients to matrices

    /* set values corresponding to the mesh interior */
    rstart = (ilo_s == 0      )  ? 1          : ilo_s;
    rend   = (ihi_s == numpts_s) ? numpts_s-1 : ihi_s;
    for (i=rstart; i<=rend; i++) {
        col[0]=i-1; col[1]=i;
        dh = arr_radius_s[i] - arr_radius_s[i-1];
        // A
        value[0] = -1.0 / dh;
        value[1] = 1.0 / dh;
        ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES); CHKERRQ(ierr);
        // B
        dh2 = arr_radius_b[i] - arr_radius_s[i-1];
        value[0] = 1.0 - dh2/dh;
        value[1] = dh2/dh;
        ierr = MatSetValues(B,1,&i,2,col,value,INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    /* store to context */
    ierr = MatDuplicate( B, MAT_COPY_VALUES, &E->qty_at_b ); CHKERRQ(ierr);
    ierr = MatDuplicate( A, MAT_COPY_VALUES, &E->ddr_at_b ); CHKERRQ(ierr);

    /* clean up */
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = MatDestroy(&B); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* TODO: can this function be converted to PetscMalloc calls?  Also, the pointer
   should be read in as an argument, rather than returned, to be consistent with
   other functions */
PetscScalar** Make2DPetscScalarArray( PetscInt arraySizeX, PetscInt arraySizeY)
{

    PetscScalar** theArray;
    PetscInt i;

    theArray = (PetscScalar**) malloc(arraySizeX*sizeof(PetscScalar*));
    for( i=0; i<arraySizeX; i++){
        theArray[i] = (PetscScalar *) malloc(arraySizeY*sizeof(PetscScalar));
    }

    return theArray;
}

PetscErrorCode make_vec_mask( DM dm, PetscInt index, Vec * mask_ptr )
{
    PetscErrorCode ierr;
    PetscInt       numpts,i,ilo,ihi,w;
    const PetscScalar one=1.0, zero=0.0;

    ierr = DMDAGetInfo(dm,NULL,&numpts,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* create mask vector */
    ierr = VecCreate( PETSC_COMM_WORLD, mask_ptr );CHKERRQ(ierr);
    ierr = VecSetSizes( *mask_ptr, PETSC_DECIDE, numpts );CHKERRQ(ierr);
    ierr = VecSetFromOptions( *mask_ptr );CHKERRQ(ierr);
    ierr = VecSetUp( *mask_ptr );CHKERRQ(ierr);

    ierr = VecSet( *mask_ptr, 0.0 );CHKERRQ(ierr);

    ierr = DMDAGetCorners(dm,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;

    for(i=ilo; i<ihi; ++i){
        if( i < index ){
            ierr = VecSetValues( *mask_ptr, 1, &i, &one, INSERT_VALUES );CHKERRQ(ierr);
        }
        else{
            ierr = VecSetValues( *mask_ptr, 1, &i, &zero, INSERT_VALUES );CHKERRQ(ierr);
            break;
        }
    }

    ierr = VecAssemblyBegin(*mask_ptr);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*mask_ptr);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

PetscErrorCode average_by_mass_staggered( Ctx *E, Vec in_vec, Vec * in_mask_ptr, PetscScalar *out )
{
    PetscErrorCode    ierr;
    Mesh              *M = &E->mesh;
    Vec               qty_s, mass_s;
    PetscScalar       mass;
    PetscInt          numpts_s;
    Vec               in_mask = *in_mask_ptr;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* create work vectors */
    ierr = VecCreate( PETSC_COMM_WORLD, &qty_s ); CHKERRQ(ierr);
    ierr = VecSetSizes( qty_s, PETSC_DECIDE, numpts_s ); CHKERRQ(ierr);
    ierr = VecSetFromOptions( qty_s ); CHKERRQ(ierr);
    ierr = VecSetUp( qty_s ); CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &mass_s ); CHKERRQ(ierr);
    ierr = VecSetSizes( mass_s, PETSC_DECIDE, numpts_s ); CHKERRQ(ierr);
    ierr = VecSetFromOptions( mass_s ); CHKERRQ(ierr);
    ierr = VecSetUp( mass_s ); CHKERRQ(ierr);

    /* compute mass-weighted average */
    ierr = VecPointwiseMult(qty_s,in_vec,in_mask); CHKERRQ(ierr);
    ierr = VecPointwiseMult(mass_s,M->mass_s,in_mask); CHKERRQ(ierr);
    ierr = VecPointwiseMult(qty_s,qty_s,mass_s); CHKERRQ(ierr);
    ierr = VecSum(qty_s,out); CHKERRQ(ierr);
    ierr = VecSum(mass_s,&mass); CHKERRQ(ierr);

    *out = *out / mass;

    /* this block would prevent 'nans' from being returned, but 'nans'
       are helpful because it indicates the rheological front is not
       yet moving */
#if 0
    if (mass > 0.0){
        *out = *out / mass;
    }
    else{
        *out = 0.0;
    }
#endif

    ierr = VecDestroy( &qty_s ); CHKERRQ(ierr);
    ierr = VecDestroy( &mass_s ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

PetscErrorCode invert_vec_mask( Vec * in_vec_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecShift( *in_vec_ptr, -1 ); CHKERRQ(ierr);
    ierr = VecScale( *in_vec_ptr, -1 ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

#if 0
static PetscErrorCode set_d_dr_quadratic( Ctx *E )
{
    /* interpolate quadratic functions to get val and dval/dr at
       basic nodes.  For overlapping regions can get an estimate
       from the left and an estimate from the right, which are
       combined by the arithmetic average

       Y     = Ax**2 + Bx + C
       dY/dr = 2Ax + B

    */

    PetscErrorCode ierr;
    PetscInt i, i2, col[3], rstart, rend;
    PetscInt numpts_s, ilo_s, ihi_s, w_s;
    PetscInt numpts_b, ilo_b, ihi_b, w_b;
    PetscScalar value[3], cc;
    PetscScalar *arr_radius_s, *arr_radius_b, dh, h1, h2;
    Mat A1, B1, C1, S1a, S1b, dS1, A2, B2, C2, S2a, S2b, dS2;
    Mesh *M;
    Vec dx1, dx1sq, dx2, dx2sq, count;
    //Vec test;

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

    // initialise A coefficient matrices
    ierr = MatCreate( PETSC_COMM_WORLD, &A1 );CHKERRQ(ierr);
    ierr = MatSetSizes(A1,PETSC_DECIDE,PETSC_DECIDE,numpts_b,numpts_s);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A1);CHKERRQ(ierr);
    ierr = MatSetUp(A1);CHKERRQ(ierr);

    ierr = MatCreate( PETSC_COMM_WORLD, &A2 );CHKERRQ(ierr);
    ierr = MatSetSizes(A2,PETSC_DECIDE,PETSC_DECIDE,numpts_b,numpts_s);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A2);CHKERRQ(ierr);
    ierr = MatSetUp(A2);CHKERRQ(ierr);

    // initialise B coefficient matrices
    ierr = MatCreate( PETSC_COMM_WORLD, &B1 );CHKERRQ(ierr);
    ierr = MatSetSizes(B1,PETSC_DECIDE,PETSC_DECIDE,numpts_b,numpts_s);CHKERRQ(ierr);
    ierr = MatSetFromOptions(B1);CHKERRQ(ierr);
    ierr = MatSetUp(B1);CHKERRQ(ierr);

    ierr = MatCreate( PETSC_COMM_WORLD, &B2 );CHKERRQ(ierr);
    ierr = MatSetSizes(B2,PETSC_DECIDE,PETSC_DECIDE,numpts_b,numpts_s);CHKERRQ(ierr);
    ierr = MatSetFromOptions(B2);CHKERRQ(ierr);
    ierr = MatSetUp(B2);CHKERRQ(ierr);

    // initialise C coefficient matrices
    ierr = MatCreate( PETSC_COMM_WORLD, &C1 );CHKERRQ(ierr);
    ierr = MatSetSizes(C1,PETSC_DECIDE,PETSC_DECIDE,numpts_b,numpts_s);CHKERRQ(ierr);
    ierr = MatSetFromOptions(C1);CHKERRQ(ierr);
    ierr = MatSetUp(C1);CHKERRQ(ierr);

    ierr = MatCreate( PETSC_COMM_WORLD, &C2 );CHKERRQ(ierr);
    ierr = MatSetSizes(C2,PETSC_DECIDE,PETSC_DECIDE,numpts_b,numpts_s);CHKERRQ(ierr);
    ierr = MatSetFromOptions(C2);CHKERRQ(ierr);
    ierr = MatSetUp(C2);CHKERRQ(ierr);

    // add coefficients to matrices

    /* set values corresponding to the mesh interior */
    rstart = (ilo_s == 0      )  ? 1          : ilo_s;
    rend   = (ihi_s == numpts_s) ? numpts_s-1 : ihi_s;
    for (i=rstart; i<rend; i++) {
        col[0]=i-1; col[1]=i; col[2]=i+1;
        i2 = i+1;
        h1 = arr_radius_s[i] - arr_radius_s[i-1];
        h2 = arr_radius_s[i+1] - arr_radius_s[i];
        // A
        value[0] = 1.0 / (h1*(h1+h2));
        value[1] = -1.0 / (h1*h2);
        value[2] = 1.0 / (h2*(h1+h2));
        ierr = MatSetValues(A1,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr); // A1
        if (i<rend-1){
        ierr = MatSetValues(A2,1,&i2,3,col,value,INSERT_VALUES);CHKERRQ(ierr); // A2
        }
        // B
        value[0] = (-2.0*h1-h2)/(h1*(h1+h2));
        value[1] = (h1+h2)/(h1*h2);
        value[2] = -h1/(h2*(h1+h2));
        ierr = MatSetValues(B1,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr); // B1
        if (i<rend-1){
        ierr = MatSetValues(B2,1,&i2,3,col,value,INSERT_VALUES);CHKERRQ(ierr); // B2
        }
        // C
        value[0] = 1.0;
        value[1] = 0.0;
        value[2] = 0.0;
        ierr = MatSetValues(C1,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr); // C1
        if (i<rend-1){
        ierr = MatSetValues(C2,1,&i2,3,col,value,INSERT_VALUES);CHKERRQ(ierr); // C2
        }
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
        ierr = MatSetValues(A1,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
        // B
        value[0] = -h2 / (h1*(h1+h2));
        value[1] = (h2-h1)/(h1*h2);
        value[2] = h1 / (h2*(h1+h2));
        ierr = MatSetValues(B1,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
        // C
        value[0] = 0.0;
        value[1] = 1.0;
        value[2] = 0.0;
        ierr = MatSetValues(C1,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    }

    /* this vector gives the distance of the basic node from the
       nearest staggered node to its left */

    ierr = VecCreate( PETSC_COMM_WORLD, &dx1 );CHKERRQ(ierr);
    ierr = VecSetSizes( dx1, PETSC_DECIDE, numpts_b );CHKERRQ(ierr);
    ierr = VecSetFromOptions( dx1 );CHKERRQ(ierr);
    ierr = VecSetUp( dx1 );CHKERRQ(ierr);

    /* this vector gives distance of the basic node from the
       second staggerd node to its left */

    ierr = VecCreate( PETSC_COMM_WORLD, &dx2 );CHKERRQ(ierr);
    ierr = VecSetSizes( dx2, PETSC_DECIDE, numpts_b );CHKERRQ(ierr);
    ierr = VecSetFromOptions( dx2 );CHKERRQ(ierr);
    ierr = VecSetUp( dx2 );CHKERRQ(ierr);

    /* counter to keep track of number of estimates of val and dval/dr */
    ierr = VecCreate( PETSC_COMM_WORLD, &count );CHKERRQ(ierr);
    ierr = VecSetSizes( count, PETSC_DECIDE, numpts_b );CHKERRQ(ierr);
    ierr = VecSetFromOptions( count );CHKERRQ(ierr);
    ierr = VecSetUp( count );CHKERRQ(ierr);

    /* set values corresponding to the vector interior */
    /* remember that top and bottom entries are always zero */
    rstart = (ilo_b == 0      )  ? 1          : ilo_b;
    rend   = (ihi_b == numpts_b) ? numpts_b-1 : ihi_b;
    for (i=rstart; i<rend; i++) {
        dh = 0.5 * (arr_radius_b[i] - arr_radius_b[i-1]);
        VecSetValues( dx1, 1, &i, &dh, INSERT_VALUES );CHKERRQ(ierr);
        cc = 1.0;
        VecSetValues( count, 1, &i, &cc, INSERT_VALUES ); CHKERRQ(ierr);
    }

    for (i=rstart+1; i<rend-1; i++) {
        dh = 0.5 * (arr_radius_b[i-1] - arr_radius_b[i-2]);
        dh += arr_radius_b[i] - arr_radius_b[i-1];
        VecSetValues( dx2, 1, &i, &dh, INSERT_VALUES );CHKERRQ(ierr);
        cc = 1.0;
        VecSetValues( count, 1, &i, &cc, ADD_VALUES ); CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(dx1);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(dx1);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(dx2);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(dx2);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(count);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(count);CHKERRQ(ierr);

    // square of dx1
    ierr = VecDuplicate(dx1,&dx1sq);CHKERRQ(ierr);
    ierr = VecCopy(dx1,dx1sq);CHKERRQ(ierr);
    ierr = VecPow(dx1sq,2.0);CHKERRQ(ierr);

    // square of dx2
    ierr = VecDuplicate(dx2,&dx2sq);CHKERRQ(ierr);
    ierr = VecCopy(dx2,dx2sq);CHKERRQ(ierr);
    ierr = VecPow(dx2sq,2.0);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayRead(E->da_s,M->radius_s,&arr_radius_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(E->da_b,M->radius_b,&arr_radius_b);CHKERRQ(ierr);

    /* assemble the matrices */
    // A
    ierr = MatAssemblyBegin(A1, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A1, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(A2, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A2, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    // B
    ierr = MatAssemblyBegin(B1, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B1, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(B2, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B2, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    // C
    ierr = MatAssemblyBegin(C1, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C1, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyBegin(C2, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C2, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    /* for computing first estimate of S */
    // S1a: ax**2 term
    ierr = MatDuplicate( A1, MAT_COPY_VALUES, &S1a );CHKERRQ(ierr);
    ierr = MatDiagonalScale( S1a, dx1sq, NULL);CHKERRQ(ierr);
    // S1b: bx term
    ierr = MatDuplicate( B1, MAT_COPY_VALUES, &S1b );CHKERRQ(ierr);
    ierr = MatDiagonalScale( S1b, dx1, NULL);CHKERRQ(ierr);
    // add (final matrix is S1a)
    ierr = MatAXPY( S1a, 1.0, S1b, SAME_NONZERO_PATTERN );CHKERRQ(ierr); 
    ierr = MatAXPY( S1a, 1.0, C1, SAME_NONZERO_PATTERN );CHKERRQ(ierr);

    /* for computing first estimate of dS/dr */
    ierr = MatDuplicate( A1, MAT_COPY_VALUES, &dS1 );CHKERRQ(ierr);
    ierr = MatDiagonalScale( dS1, dx1, NULL);CHKERRQ(ierr);
    ierr = MatScale( dS1, 2.0);CHKERRQ(ierr);
    ierr = MatAXPY( dS1, 1.0, B1, SAME_NONZERO_PATTERN );

    /* for computing second estimate of S */
    // S2a: ax**2 term
    ierr = MatDuplicate( A2, MAT_COPY_VALUES, &S2a );CHKERRQ(ierr);
    ierr = MatDiagonalScale( S2a, dx2sq, NULL); CHKERRQ(ierr);
    // S2b: bx term
    ierr = MatDuplicate( B2, MAT_COPY_VALUES, &S2b );CHKERRQ(ierr);
    ierr = MatDiagonalScale( S2b, dx2, NULL); CHKERRQ(ierr);
    // add (final matrix is S2a)
    ierr = MatAXPY( S2a, 1.0, S2b, SAME_NONZERO_PATTERN ); CHKERRQ(ierr); 
    ierr = MatAXPY( S2a, 1.0, C2, SAME_NONZERO_PATTERN ); CHKERRQ(ierr);

    /* for computing second estimate of dS/dr */
    ierr = MatDuplicate( A2, MAT_COPY_VALUES, &dS2 );CHKERRQ(ierr);
    ierr = MatDiagonalScale( dS2, dx2, NULL);CHKERRQ(ierr);
    ierr = MatScale( dS2, 2.0);CHKERRQ(ierr);
    ierr = MatAXPY( dS2, 1.0, B2, SAME_NONZERO_PATTERN );

    /* combine estimates by arithmetic average */
    ierr = VecReciprocal( count );CHKERRQ(ierr);
    ierr = MatAXPY( S1a, 1.0, S2a, DIFFERENT_NONZERO_PATTERN ); CHKERRQ(ierr);
    ierr = MatDiagonalScale( S1a, count, NULL); CHKERRQ(ierr); 
    ierr = MatAXPY( dS1, 1.0, dS2, DIFFERENT_NONZERO_PATTERN ); CHKERRQ(ierr);
    ierr = MatDiagonalScale( dS1, count, NULL); CHKERRQ(ierr);

    /* store to context */
    ierr = MatDuplicate( S1a, MAT_COPY_VALUES, &E->qty_at_b ); CHKERRQ(ierr);
    ierr = MatDuplicate( dS1, MAT_COPY_VALUES, &E->ddr_at_b ); CHKERRQ(ierr);

    /* test by taking d/dr of radius at staggered nodes */
    //ierr = VecDuplicate(M->radius_b,&test);CHKERRQ(ierr);
    //ierr = MatMult( E->ddr_at_b, M->radius_s, test );CHKERRQ(ierr);

    /* clean up temporary matrices and vectors */
    ierr = MatDestroy(&A1); CHKERRQ(ierr);
    ierr = MatDestroy(&A2); CHKERRQ(ierr);
    ierr = MatDestroy(&B1); CHKERRQ(ierr);
    ierr = MatDestroy(&B2); CHKERRQ(ierr);
    ierr = MatDestroy(&C1); CHKERRQ(ierr);
    ierr = MatDestroy(&C2); CHKERRQ(ierr);
    ierr = MatDestroy(&S1a); CHKERRQ(ierr);
    ierr = MatDestroy(&S1b); CHKERRQ(ierr);
    ierr = MatDestroy(&dS1); CHKERRQ(ierr);
    ierr = MatDestroy(&S2a); CHKERRQ(ierr);
    ierr = MatDestroy(&S2b); CHKERRQ(ierr);
    ierr = MatDestroy(&dS2); CHKERRQ(ierr);
    ierr = VecDestroy(&dx1); CHKERRQ(ierr);
    ierr = VecDestroy(&dx1sq); CHKERRQ(ierr);
    ierr = VecDestroy(&dx2); CHKERRQ(ierr);
    ierr = VecDestroy(&dx2sq); CHKERRQ(ierr);
    ierr = VecDestroy(&count); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
#endif
