#include "energy.h"

static PetscErrorCode set_Jtot( Ctx * );
#ifdef CONDUCTION
static PetscErrorCode append_Jcond( Ctx * );
#endif
#ifdef CONVECTION
static PetscErrorCode append_Jconv( Ctx * );
#endif
#ifdef MIXING
static PetscErrorCode append_Jmix( Ctx * );
#endif
#ifdef SEPARATION
static PetscErrorCode append_Jgrav( Ctx * );
#endif
#ifdef HRADIO
static PetscErrorCode append_Hradio( Ctx *, PetscReal t );
#endif

///////////////////////////
/* internal heat sources */
///////////////////////////
/* total internal heat generation */
PetscErrorCode set_Htot( Ctx *E, PetscReal t )
{
    PetscErrorCode ierr;
    //Mesh           *M = &E->mesh;
    Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    /* Htot = int_V rho H dV */

    /* initialise to zero */
    ierr = VecSet( S->Htot_s, 0.0 ); CHKERRQ(ierr);

    /* total internal heat generation by summing terms */
#ifdef HRADIO
    ierr = append_Hradio( E, t ); CHKERRQ(ierr);
#endif

    PetscFunctionReturn(0);
}

#ifdef HRADIO
/* radiogenic heat generation */
static PetscErrorCode append_Hradio( Ctx *E, PetscReal t )
{
    PetscErrorCode ierr;
    //Mesh           *M = &E->mesh;
    Solution       *S = &E->solution;
    PetscReal      tyrs = t * TIME0YEARS; // time in years
    PetscScalar    Alheat;
    PetscScalar    thalf = 7.17E5; // years (from Tim's spreadsheet)
    /* initial heating rate of Al26 */
    /* way higher than it should be, but useful for checking
       the code */
    PetscScalar    Alinit = 1.526E-7; // W/kg (from Tim's spreadsheet)
    PetscScalar    decayconst = log(2)/thalf; // 1/years

    PetscFunctionBeginUser;

    /* a simple calculation for Earth (r=6371000, core size=0.55)
       tells you that to balance a heat flux of 10^6 W requires
       internal heating around 1.41E-4 W/kg, which non-dimensionalised
       is about 2.14E-8 */

    /* do stuff here */

    /* this applies constant heating throughout the domain */
    //ierr = VecSet( S->Hradio_s, 1.0E-8 ); CHKERRQ(ierr);

    // no heating example */
    //ierr = VecSet( S->Hradio_s, 0.0 ); CHKERRQ(ierr);

    /* decay Al26 radiogenic heat production with time */
    Alheat = Alinit*exp(-decayconst*tyrs);
    Alheat /= 6584.128807671617; // non-dimensionalise as below
    ierr = VecSet( S->Hradio_s, Alheat ); CHKERRQ(ierr);

    /* the dimensional scaling for heat sources is:
           6584.128807671617 W / kg
       so once you compute your desired dimensional heating,
       divide by the above number to get the quantity in
       non-dimensional units (required for the code) */

    ierr = VecAXPY( S->Htot_s, 1.0, S->Hradio_s ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
#endif

///////////////////
/* energy fluxes */
///////////////////

/* total energy flow (flux*area) at basic nodes */
PetscErrorCode set_Etot( Ctx *E )
{
    PetscErrorCode ierr;
    Mesh *M;
    Solution *S;

    PetscFunctionBeginUser;

    M = &E->mesh;
    S = &E->solution;

    ierr = set_Jtot(E); CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Etot,S->Jtot,M->area_b); CHKERRQ(ierr);

    // I don't think this has to be done here?
    //ierr = VecAssemblyBegin(S->Etot);CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(S->Etot);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

/* total heat flux at basic nodes */
static PetscErrorCode set_Jtot( Ctx *E )
{
    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;

    S = &E->solution;

    /* must initialise to zero */
    ierr = VecSet( S->Jtot, 0.0 ); CHKERRQ(ierr);

    /* total heat flux by summing terms */
#ifdef CONDUCTION
    ierr = append_Jcond( E );
#endif
#ifdef CONVECTION
    ierr = append_Jconv( E );
#endif
#ifdef MIXING
    ierr = append_Jmix( E );
#endif
#ifdef SEPARATION
    ierr = append_Jgrav( E );
#endif

    // I don't think this has to be done here?
    //ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
    //ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

#ifdef CONVECTION
/* convective heat flux at basic nodes */
static PetscErrorCode append_Jconv( Ctx *E )
{
    PetscErrorCode ierr;
    Solution       *S;

    PetscFunctionBeginUser;

    S = &E->solution;

    /* convective heat flux */
    //   arr_Jconv[i] = -arr_dSdr[i] * arr_kappah[i] * arr_rho[i] * arr_temp[i];

    ierr = VecPointwiseMult(S->Jconv,S->dSdr, S->kappah);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jconv,S->Jconv, S->rho);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jconv, S->Jconv, S->temp);CHKERRQ(ierr);
    ierr = VecScale(S->Jconv, -1.0);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jconv ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
#endif

#ifdef MIXING
/* mixing heat flux (latent heat transport) at basic nodes */
static PetscErrorCode append_Jmix( Ctx *E )
{
    PetscErrorCode ierr;
    DM             da_b=E->da_b;
    Solution       *S;
    PetscInt       i, ilo, ihi, w;
    PetscScalar    *arr_gphi, *arr_fwtl, *arr_fwts, *arr_Jmix;

    PetscFunctionBeginUser;

    S = &E->solution;

    /* convective mixing */
    // these first two lines give F_Jmix (in the python script)
    // arr_Jmix[i] = arr_dSdr[i] - arr_phi[i] * arr_dSliqdr[i];
    // arr_Jmix[i] += (arr_phi[i]-1.0) * arr_dSsoldr[i];
    ierr = VecWAXPY(S->Jmix,-1.0,S->dSliqdr,S->dSsoldr);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jmix,S->Jmix,S->phi);CHKERRQ(ierr);
    ierr = VecAYPX(S->Jmix,1.0,S->dSdr);CHKERRQ(ierr);
    ierr = VecAXPY(S->Jmix,-1.0,S->dSsoldr);CHKERRQ(ierr);
    // arr_Jmix[i] *= -arr_kappah[i] * arr_rho[i] * arr_temp[i];
    ierr = VecPointwiseMult(S->Jmix,S->Jmix,S->kappah);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jmix,S->Jmix,S->rho);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jmix,S->Jmix,S->temp);CHKERRQ(ierr);
    ierr = VecScale(S->Jmix,-1.0);CHKERRQ(ierr);

    ierr = DMDAGetCorners(da_b,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;

    ierr = DMDAVecGetArrayRead(da_b,S->fwtl,&arr_fwtl);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fwts,&arr_fwts);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->gphi,&arr_gphi);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->Jmix,&arr_Jmix);CHKERRQ(ierr);

    /* to smoothly blend in convective mixing across the liquidus
       and solidus */
    for(i=ilo; i<ihi; ++i){
        if(arr_gphi[i] > 0.5){
            arr_Jmix[i] *= 1.0 - arr_fwtl[i];
        }
        else if (arr_gphi[i] <= 0.5){
            arr_Jmix[i] *= arr_fwts[i];
        }   
    }

    ierr = DMDAVecRestoreArrayRead(da_b,S->fwtl,&arr_fwtl);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fwts,&arr_fwts);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->gphi,&arr_gphi);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->Jmix,&arr_Jmix);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jmix ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
#endif

#ifdef CONDUCTION
/* conductive heat flux at basic nodes */
static PetscErrorCode append_Jcond( Ctx *E )
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

    /* at the moment, thermophysical properties are only computed at
       the basic internal nodes, so temp/cp = 0/0 = NaN at the top
       and bottom surface */

    ierr = VecAXPY( S->Jtot, 1.0, S->Jcond ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
#endif

#ifdef SEPARATION
/* gravitational separation heat flux at basic nodes */
static PetscErrorCode append_Jgrav( Ctx *E )
{
    PetscErrorCode ierr;
    Solution *S;
    Vec cond1, cond2, rho, rhol, rhos, F;
    PetscInt i,ilo_b,ihi_b,w_b,numpts_b;
    DM da_b=E->da_b;
    const PetscScalar *arr_cond1, *arr_cond2, *arr_liquidus_rho, *arr_phi, *arr_solidus_rho;
    PetscScalar icond1, icond2, irhos, irhol, iphi;
    PetscScalar *arr_F;

    PetscFunctionBeginUser;

    S = &E->solution;
    rho = S->rho;
    rhol = S->liquidus_rho;
    rhos = S->solidus_rho;

    ierr = DMDAGetInfo(da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

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
    ihi_b = ilo_b + w_b;

    ierr = DMDAVecGetArrayRead(da_b, cond1, &arr_cond1);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b, cond2, &arr_cond2);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b, F, &arr_F);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->phi,&arr_phi);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);

    /* (I think) unavoidably have to loop over array to build F, since
       Petsc Vecs do not support logic operations? */
    for(i=ilo_b; i<ihi_b; ++i){
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

    ierr = VecAXPY( S->Jtot, 1.0, S->Jgrav ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
#endif