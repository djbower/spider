#include "atmosphere.h"
#include "energy.h"
#include "matprop.h"
#include "twophase.h"
#include "util.h"

static PetscErrorCode set_Jtot( Ctx * );
static PetscErrorCode append_Jcond( Ctx * );
static PetscErrorCode append_Jconv( Ctx * );
static PetscErrorCode append_Jmix( Ctx * );
static PetscErrorCode append_Jgrav( Ctx * );
static PetscErrorCode append_Hradio( Ctx *, PetscReal );
static PetscErrorCode append_Htidal( Ctx *, PetscReal );
static PetscScalar get_radiogenic_heat_production( RadionuclideParameters const, PetscReal );
static PetscScalar get_tsurf_using_parameterised_boundary_layer( PetscScalar, const AtmosphereParameters );
static PetscScalar get_dtsurf_using_parameterised_boundary_layer( PetscScalar, const AtmosphereParameters );

///////////////////////////
/* internal heat sources */
///////////////////////////
/* total internal heat generation */
PetscErrorCode set_Htot( Ctx *E, PetscReal time )
{
    PetscErrorCode ierr;
    Parameters     const P = E->parameters;
    Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    /* Htot = int_V rho H dV */

    /* initialise to zero */
    ierr = VecSet( S->Hradio_s, 0.0 );CHKERRQ(ierr);
    ierr = VecSet( S->Htot_s, 0.0 );CHKERRQ(ierr);
    ierr = VecSet( S->Hal26_s, 0.0 );CHKERRQ(ierr);
    ierr = VecSet( S->Hk40_s, 0.0 );CHKERRQ(ierr);
    ierr = VecSet( S->Hfe60_s, 0.0 );CHKERRQ(ierr);
    ierr = VecSet( S->Hth232_s, 0.0 );CHKERRQ(ierr);
    ierr = VecSet( S->Hu235_s, 0.0 );CHKERRQ(ierr);
    ierr = VecSet( S->Hu238_s, 0.0 );CHKERRQ(ierr);

    /* total internal heat generation by summing terms */
    if (P->HRADIO){
      ierr = append_Hradio( E, time ); CHKERRQ(ierr);
    }
    if (P->HTIDAL){
      ierr = append_Htidal( E, time ); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

/* radiogenic heat generation */
static PetscErrorCode append_Hradio( Ctx *E, PetscReal time )
{

    PetscErrorCode ierr;
    PetscInt i;
    PetscScalar H;
    Solution       *S = &E->solution;
    Parameters const P = E->parameters;
    RadionuclideParameters Rp;

// FIXME: REMOVE
#if 0
    RadiogenicIsotopeParameters const *al26 = &P->al26_parameters;
    RadiogenicIsotopeParameters const *k40 = &P->k40_parameters;
    RadiogenicIsotopeParameters const *fe60 = &P->fe60_parameters;
    RadiogenicIsotopeParameters const *th232 = &P->th232_parameters;
    RadiogenicIsotopeParameters const *u235 = &P->u235_parameters;
    RadiogenicIsotopeParameters const *u238 = &P->u238_parameters;
#endif

    PetscFunctionBeginUser;

    for (i=0;i<P->n_radionuclides; ++i){

        Rp = P->radionuclide_parameters[i];
        H = get_radiogenic_heat_production( Rp, time );
        ierr = VecShift( S->Hradio_s, H ); CHKERRQ(ierr);
    }

//FIXME: REMOVE
#if 0
    // al26
    H = get_radiogenic_heat_production( al26, time );
    ierr = VecSet(S->Hal26_s,H);CHKERRQ(ierr);
    ierr = VecAXPY( S->Hradio_s, 1.0, S->Hal26_s ); CHKERRQ(ierr);
    // k40
    H = get_radiogenic_heat_production( k40, time );
    ierr = VecSet(S->Hk40_s,H);CHKERRQ(ierr);
    ierr = VecAXPY( S->Hradio_s, 1.0, S->Hk40_s ); CHKERRQ(ierr);
    // fe60
    H = get_radiogenic_heat_production( fe60, time );
    ierr = VecSet(S->Hfe60_s,H);CHKERRQ(ierr);
    ierr = VecAXPY( S->Hradio_s, 1.0, S->Hfe60_s ); CHKERRQ(ierr);
    // th232
    H = get_radiogenic_heat_production( th232, time );
    ierr = VecSet(S->Hth232_s,H);CHKERRQ(ierr);
    ierr = VecAXPY( S->Hradio_s, 1.0, S->Hth232_s ); CHKERRQ(ierr);
    // u235
    H = get_radiogenic_heat_production( u235, time );
    ierr = VecSet(S->Hu235_s,H);CHKERRQ(ierr);
    ierr = VecAXPY( S->Hradio_s, 1.0, S->Hu235_s ); CHKERRQ(ierr);
    // u238
    H = get_radiogenic_heat_production( u238, time );
    ierr = VecSet(S->Hu238_s,H);CHKERRQ(ierr);
    ierr = VecAXPY( S->Hradio_s, 1.0, S->Hu238_s ); CHKERRQ(ierr);
#endif

    // append total of radiogenic heating to total heating vector
    ierr = VecAXPY( S->Htot_s, 1.0, S->Hradio_s ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* tidal heat generation */
static PetscErrorCode append_Htidal( Ctx *E, PetscReal tyrs )
{
    /* TODO: template code snippets below, but tidal heating is
       currently not implemented */

    //PetscErrorCode ierr;
    //Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    /* do stuff here */
    //ierr = VecSet( S->Htidal_s, 0.0 ); CHKERRQ(ierr);

    // final command is always an append call to the Htot_s array
    //ierr = VecAXPY( S->Htot_s, 1.0, S->Htidal_s ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscScalar get_radiogenic_heat_production( RadionuclideParameters const Iso, PetscReal time )
{
    PetscScalar H;

    H = (Iso->t0-time) * PetscLogScalar(2.0);
    H /= Iso->half_life;
    H = PetscExpScalar(H);
    H *= Iso->heat_production * Iso->abundance * Iso->concentration * 1.0E-6; // since ppm

    return H;

}

///////////////////
/* energy fluxes */
///////////////////

/* total energy flow (flux*area) at basic nodes */
PetscErrorCode set_Etot( Ctx *E )
{
    PetscErrorCode ierr;
    Mesh           *M = &E->mesh;
    Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    ierr = set_Jtot(E); CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Etot,S->Jtot,M->area_b); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

/* total heat flux at basic nodes */
static PetscErrorCode set_Jtot( Ctx *E )
{
    PetscErrorCode ierr;
    Parameters     const P = E->parameters;
    Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    /* initialise to zero */
    ierr = VecSet( S->Jtot, 0.0 ); CHKERRQ(ierr);

    /* total heat flux by summing terms */
    if (P->CONDUCTION){
      ierr = append_Jcond( E );
    }
    if (P->CONVECTION){
      ierr = append_Jconv( E );
    }
    if (P->MIXING){
      ierr = append_Jmix( E );
    }
    if (P->SEPARATION){
      ierr = append_Jgrav( E );
    }

    PetscFunctionReturn(0);

}

/* convective heat flux at basic nodes */
static PetscErrorCode append_Jconv( Ctx *E )
{
    PetscErrorCode ierr;
    Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    /* convective heat flux */
    //   arr_Jconv[i] = -arr_dSdr[i] * arr_kappah[i] * arr_rho[i] * arr_temp[i];

    ierr = VecPointwiseMult(S->Jconv,S->dSdr, S->kappah);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jconv,S->Jconv, S->rho);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jconv, S->Jconv, S->temp);CHKERRQ(ierr);
    ierr = VecScale(S->Jconv, -1.0);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jconv ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

/* mixing heat flux (latent heat transport) at basic nodes */
static PetscErrorCode append_Jmix( Ctx *E )
{
    PetscErrorCode ierr;
    DM             da_b = E->da_b;
    Solution       *S = &E->solution;
    PetscInt       i, ilo, ihi, w;
    PetscScalar    *arr_gphi, *arr_fwtl, *arr_fwts, *arr_Jmix;

    PetscFunctionBeginUser;

    /* convective mixing */
    // arr_Jmix[i] = arr_dSdr[i] - arr_phi[i] * arr_dSliqdr[i];
    // arr_Jmix[i] += (arr_phi[i]-1.0) * arr_dSsoldr[i];
    ierr = VecWAXPY(S->Jmix,-1.0,S->dSliqdr,S->dSsoldr);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jmix,S->Jmix,S->phi);CHKERRQ(ierr);
    ierr = VecAYPX(S->Jmix,1.0,S->dSdr);CHKERRQ(ierr);
    ierr = VecAXPY(S->Jmix,-1.0,S->dSsoldr);CHKERRQ(ierr);
    // arr_Jmix[i] *= -arr_kappac[i] * arr_rho[i] * arr_temp[i];
    ierr = VecPointwiseMult(S->Jmix,S->Jmix,S->kappac);CHKERRQ(ierr);
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
            // FIXME: smoothing width is hard-coded
            arr_Jmix[i] *= 1.0 - tanh_weight( arr_gphi[i], 1.0, 1.0E-2 );
        }
        else if (arr_gphi[i] <= 0.5){
            // FIXME: smoothing width is hard-coded
            arr_Jmix[i] *= tanh_weight( arr_gphi[i], 0.0, 1.0E-2 );
        }   
    }

    ierr = DMDAVecRestoreArrayRead(da_b,S->fwtl,&arr_fwtl);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fwts,&arr_fwts);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->gphi,&arr_gphi);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->Jmix,&arr_Jmix);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jmix ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

/* conductive heat flux at basic nodes */
static PetscErrorCode append_Jcond( Ctx *E )
{
    PetscErrorCode ierr;
    Solution *S = &E->solution;

    PetscFunctionBeginUser;

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

/* gravitational separation heat flux at basic nodes */
static PetscErrorCode append_Jgrav( Ctx *E )
{

    /* FIXME: this needs updating for composition! */

    PetscErrorCode ierr;
    Solution *S = &E->solution;
    Vec cond1, cond2, F;
    Vec rho = S->rho;
    Vec rhol = S->liquidus_rho;
    Vec rhos = S->solidus_rho;
    Parameters P = E->parameters;

//rho, rhol, rhos, F;
    PetscInt i,ilo_b,ihi_b,w_b,numpts_b;
    DM da_b=E->da_b;
    const PetscScalar *arr_cond1, *arr_cond2, *arr_liquidus_rho, *arr_phi, *arr_solidus_rho;
    PetscScalar icond1, icond2, irhos, irhol, iphi;
    PetscScalar *arr_F;

    PetscFunctionBeginUser;

    //S = &E->solution;
    //rho = S->rho;
    //rhol = S->liquidus_rho;
    //rhos = S->solidus_rho;

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
    ierr = VecScale( S->Jgrav, PetscPowScalar(P->grain,2) );CHKERRQ(ierr);
    ierr = VecScale( S->Jgrav, P->gravity );CHKERRQ(ierr);
    ierr = VecPointwiseMult( S->Jgrav, S->Jgrav, F );CHKERRQ(ierr);
    // arr_Jgrav[i] /= PetscPowScalar(10.0, LOG10VISC_MEL);
    ierr = VecScale( S->Jgrav, 1.0/PetscPowScalar(10.0, P->eos1_parameters->log10visc));CHKERRQ(ierr);

    ierr = VecDestroy(&cond1);CHKERRQ(ierr);
    ierr = VecDestroy(&cond2);CHKERRQ(ierr);
    ierr = VecDestroy(&F);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jgrav ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

PetscErrorCode set_interior_structure_from_solution( Ctx *E, PetscReal t, Vec sol_in )
{

    /* set all possible quantities for a given entropy structure (i.e. top staggered 
       value S0 and dS/dr at all basic nodes excluding the top and bottom which are 
       controlled by boundary conditions).  This one function ensure that everything
       is set self-consistently. */

    PetscErrorCode       ierr;
    PetscMPIInt          rank;
    PetscScalar          temp0;
    PetscInt             const ind0 = 0;
    Atmosphere           *A  = &E->atmosphere;
    Parameters           const P  = E->parameters;
    Solution             const *S  = &E->solution;
    ScalingConstants     const SC  = P->scaling_constants;
    AtmosphereParameters const Ap = P->atmosphere_parameters;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

    /* set solution in the relevant structs */
    ierr = set_entropy_from_solution( E, sol_in );CHKERRQ(ierr);

    /* set material properties and energy fluxes and sources */
    ierr = set_gphi_smooth( E );CHKERRQ(ierr);
    ierr = set_melt_fraction_staggered( E ); CHKERRQ(ierr);
    ierr = set_capacitance_staggered( E );CHKERRQ(ierr);
    ierr = set_matprop_basic( E );CHKERRQ(ierr);
    ierr = set_Etot( E );CHKERRQ(ierr);
    ierr = set_Htot( E, t );CHKERRQ(ierr);

    ierr = set_Mliq( E );CHKERRQ(ierr);
    ierr = set_Msol( E );CHKERRQ(ierr);
    /* NOTE: we cannot set_dMliqdt, since it must be after dS/dt computation */

    ierr = set_rheological_front( E ); CHKERRQ(ierr);

    /* set surface temperature */
    if (!rank) {
        /* temperature (potential temperature if coarse mesh is used) */
        ierr = VecGetValues(S->temp,1,&ind0,&temp0); CHKERRQ(ierr);

        /* correct for ultra-thin thermal boundary layer at the surface */
        if( Ap->PARAM_UTBL ){
            A->tsurf = get_tsurf_using_parameterised_boundary_layer( temp0, Ap); // parameterised boundary layer
            A->dtsurfdt = get_dtsurf_using_parameterised_boundary_layer( temp0, Ap); // dTsurf/dT
        }
        else{
            A->tsurf = temp0; // surface temperature is potential temperature
            A->dtsurfdt = 1.0; // dTsurf/dT
        }
    }

    /* must be after A->tsurf is set for fO2 calculation */
    if( Ap->OXYGEN_FUGACITY ){
        ierr = set_oxygen_fugacity( A, Ap, SC );CHKERRQ(ierr);
    }
    else{
        /* TODO: maybe initialise these variables elsewhere? */
        A->log10fO2 = 0;
        A->dlog10fO2dT = 0;
    }

    PetscFunctionReturn(0);

}

static PetscScalar get_tsurf_using_parameterised_boundary_layer( PetscScalar temp, const AtmosphereParameters Ap )
{
    PetscScalar Ts, c, fac, num, den;
    c = Ap->param_utbl_const;

    fac = 3.0*PetscPowScalar(c,3.0)*(27.0*PetscPowScalar(temp,2.0)*c+4.0);
    fac = PetscPowScalar( fac, 1.0/2.0 );
    fac += 9.0*temp*PetscPowScalar(c,2.0);
    // numerator
    num = PetscPowScalar(2.0,1.0/3)*PetscPowScalar(fac,2.0/3)-2.0*PetscPowScalar(3.0,1.0/3)*c;
    // denominator
    den = PetscPowScalar(6.0,2.0/3)*c*PetscPowScalar(fac,1.0/3);
    // surface temperature
    Ts = num / den;

    return Ts; 
}

static PetscScalar get_dtsurf_using_parameterised_boundary_layer( PetscScalar temp, const AtmosphereParameters Ap )
{
    PetscScalar dTsdT, c, fac1, fac2, fac3, num1, den1, num2, den2, part1, part2;
    c = Ap->param_utbl_const;

    fac1 = 27*PetscSqrtScalar(3) * PetscPowScalar(c,3) * temp;
    fac2 = PetscSqrtScalar( PetscPowScalar(c,3) * (4.0 + 27.0 * PetscSqr(temp) ) );
    fac3 = PetscSqrtScalar( PetscPowScalar(c,3) * (4.0 + 27.0 * PetscSqr(temp) * c ) );

    num1 = PetscPowScalar(2.0/3.0,2.0/3.0) * (9 * PetscSqr(c) + fac1/fac2);
    den1 = PetscPowScalar( 9.0*PetscSqr(c)*temp + PetscSqrtScalar(3) * fac2, 1.0/3.0 );
    den1 *= PetscPowScalar( 9.0*PetscSqr(c)*temp + PetscSqrtScalar(3) * fac3, 1.0/3.0);
    den1 *= 3.0 * c;

    part1 = num1 / den1;

    num2 = 9.0*PetscSqr(c);
    num2 += (fac1*c) / fac3;
    num2 *= (-2.0*PetscPowScalar(3,1.0/3.0)*c + PetscPowScalar(2.0,1.0/3.0) * PetscPowScalar(9.0*PetscSqr(c)*temp + PetscSqrtScalar(3)*fac2,2.0/3.0));
    den2 = 3.0 * PetscPowScalar(6.0,2.0/3.0)*c * PetscPowScalar(9.0*PetscSqr(c)*temp + PetscSqrtScalar(3)*fac3,4.0/3.0);

    part2 = num2 / den2;

    dTsdT = part1 - part2;

    return dTsdT;

}
