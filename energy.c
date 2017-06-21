#include "energy.h"

static PetscErrorCode set_Jconv( Ctx * );
static PetscErrorCode set_Jcond( Ctx * );
static PetscErrorCode set_Jgrav( Ctx * );

/* total heat flux */
PetscErrorCode set_Jtot( Ctx *E )
{
    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;

    S = &E->solution;

    ierr = set_Jconv( E );
    ierr = set_Jcond( E );
    ierr = set_Jgrav( E );

    /* total heat flux by summing terms */
    ierr = VecAYPX( S->Jtot, 1.0, S->Jconv );CHKERRQ(ierr);
    ierr = VecAYPX( S->Jtot, 1.0, S->Jcond );CHKERRQ(ierr);
    ierr = VecAYPX( S->Jtot, 1.0, S->Jgrav );CHKERRQ(ierr);

    PetscFunctionReturn(0);

}


/* convective heat flux */
static PetscErrorCode set_Jconv( Ctx *E )
{
    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;

    S = &E->solution;

    /* convective heat flux */
    //   arr_Jconv[i] = -arr_dSdr[i] * arr_kappah[i] * arr_rho[i] * arr_temp[i];

    ierr = VecPointwiseMult(S->Jconv,S->dSdr, S->kappah);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jconv,S->Jconv, S->rho);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jconv, S->Jconv, S->temp);CHKERRQ(ierr);
    ierr = VecScale(S->Jconv, -1.0);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

/* conductive heat flux */
static PetscErrorCode set_Jcond( Ctx *E )
{
    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;

    S = &E->solution;

    /* conductive heat flux */
    //   arr_Jcond[i] = arr_temp[i] / arr_cp[i] * arr_dSdr[i] + arr_dTdrs[i];
    //   arr_Jcond[i] *= -arr_cond[i];

    ierr = VecPointwiseDivide(S->Jcond, S->temp, S->cp);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jcond, S->Jcond, S->dSdr);CHKERRQ(ierr);
    ierr = VecAYPX(S->Jcond, 1.0, S->dTdrs);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jcond, S->Jcond, S->cond);CHKERRQ(ierr);
    ierr = VecScale(S->Jcond, -1.0);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

/* gravitational separation heat flux */
static PetscErrorCode set_Jgrav( Ctx *E )
{
    PetscErrorCode ierr;
    Solution *S;
    Vec cond1, cond2, rho, rhol, rhos, F;
    PetscInt i,ilo_b,ihi_b,w_b,ilo,ihi,numpts_b;
    DM da_b=E->da_b;
    const PetscScalar *arr_cond1, *arr_cond2, *arr_liquidus_rho, *arr_phi, *arr_solidus_rho;
    PetscScalar icond1, icond2, irhos, irhol, iphi;
    PetscScalar *arr_F;

    PetscFunctionBeginUser;

    S = &E->solution;
    numpts_b = NUMPTS_B_DEFAULT;
    rho = S->rho;
    rhol = S->liquidus_rho;
    rhos = S->solidus_rho;

    ierr = VecCreate( PETSC_COMM_WORLD, &F );CHKERRQ(ierr);
    ierr = VecSetSizes( F, PETSC_DECIDE, numpts_b );CHKERRQ(ierr);
    ierr = VecSetFromOptions( F );CHKERRQ(ierr);
    ierr = VecSetUp( F );CHKERRQ(ierr);

    /* these are actually time-independent, so could be precomputed
       outside of the time loop */
    ierr = VecCreate( PETSC_COMM_WORLD, &cond1 );CHKERRQ(ierr);
    ierr = VecSetSizes( cond1, PETSC_DECIDE, numpts_b );CHKERRQ(ierr);
    ierr = VecSetFromOptions( cond1 );CHKERRQ(ierr);
    ierr = VecSetUp( cond1 );CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &cond2 );CHKERRQ(ierr);
    ierr = VecSetSizes( cond2, PETSC_DECIDE, numpts_b );CHKERRQ(ierr);
    ierr = VecSetFromOptions( cond2 );CHKERRQ(ierr);
    ierr = VecSetUp( cond2 );CHKERRQ(ierr);

    //PetscScalar cond1 = rhol / (11.993*rhos + rhol);
    VecWAXPY(cond1, 11.993, rhos, rhol);
    VecPointwiseDivide( cond1, rhol, cond1 );
    //PetscScalar cond2 = rhol / (0.29624*rhos + rhol);
    VecWAXPY(cond2, 0.29624, rhos, rhol);
    VecPointwiseDivide( cond2, rhol, cond2 );

    /* loop over all basic internal nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ilo = ilo_b == 0        ? 1            : ilo_b;
    ihi = ihi_b == numpts_b ? numpts_b - 1 : ihi_b;

    ierr = DMDAVecGetArrayRead(da_b, cond1, &arr_cond1);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b, cond2, &arr_cond2);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b, F, &arr_F);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->phi,&arr_phi);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);

    /* (I think) unavoidably have to loop over array to build F, since
       Petsc Vecs do not support logic operations? */
    for(i=ilo; i<ihi; ++i){
        icond1 = arr_cond1[i];
        icond2 = arr_cond2[i];
        iphi = arr_phi[i];
        irhol = arr_liquidus_rho[i];
        irhos = arr_solidus_rho[i];

        if(iphi < icond1){
            arr_F[i] = 0.001*PetscPowScalar(irhos,2)*PetscPowScalar(iphi,3);
            arr_F[i] /= PetscPowScalar(irhol,2)*(1.0-iphi);
        } else if(iphi > icond2){
            arr_F[i] = 2.0/9.0 * iphi * (1.0-iphi);
        } else{
            arr_F[i] = 5.0/7.0*PetscPowScalar(irhos,4.5)*PetscPowScalar(iphi,5.5)*(1.0-iphi);
            arr_F[i] /= PetscPowScalar( irhol+(irhos-irhol)*iphi, 4.5 );
        }
    }

    ierr = DMDAVecRestoreArrayRead(da_b, cond1, &arr_cond1);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, cond2, &arr_cond2);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b, F, &arr_F);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, S->liquidus_rho, &arr_liquidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, S->phi, &arr_phi);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, S->solidus_rho, &arr_solidus_rho);CHKERRQ(ierr);

    // arr_Jgrav[i] = (rhol-rhos) * rho;
    ierr = VecWAXPY( S->Jgrav, -1.0, rhos, rhol );
    ierr = VecPointwiseMult( S->Jgrav, S->Jgrav, rho );
    // arr_Jgrav[i] *= pref * PetscPowScalar(GRAIN,2) * GRAVITY * F;
    ierr = VecPointwiseMult( S->Jgrav, S->Jgrav, S->fusion );CHKERRQ(ierr);
    ierr = VecPointwiseMult( S->Jgrav, S->Jgrav, S->temp );CHKERRQ(ierr);
    ierr = VecScale( S->Jgrav, PetscPowScalar(GRAIN,2) );CHKERRQ(ierr);
    ierr = VecScale( S->Jgrav, GRAVITY );CHKERRQ(ierr);
    ierr = VecPointwiseMult( S->Jgrav, S->Jgrav, F );CHKERRQ(ierr);
    // arr_Jgrav[i] /= PetscPowScalar(10.0, LOG10VISC_MEL);
    ierr = VecScale( S->Jgrav, 1.0/PetscPowScalar(10.0, LOG10VISC_MEL));CHKERRQ(ierr);

    ierr = VecDestroy(&cond1);CHKERRQ(ierr);
    ierr = VecDestroy(&cond2);CHKERRQ(ierr);
    ierr = VecDestroy(&F);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
