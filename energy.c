#include "atmosphere.h"
#include "energy.h"
#include "eos.h"
#include "eos_composite.h"
#include "matprop.h"
#include "twophase.h"
#include "util.h"
#include "monitor.h"

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
static PetscErrorCode objective_function_surfacebc( SNES, Vec , Vec , void *ptr);
static PetscErrorCode solve_surface_entropy( Ctx * );

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

    /* total internal heat generation by summing terms */
    if (P->n_radionuclides>0){
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

    PetscFunctionBeginUser;

    for (i=0;i<P->n_radionuclides; ++i){

        Rp = P->radionuclide_parameters[i];
        H = get_radiogenic_heat_production( Rp, time );
        ierr = VecShift( S->Hradio_s, H ); CHKERRQ(ierr);
    }

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
    (void) E; // unused for now
    (void) tyrs; // unused for now

    /* do stuff here */
    //ierr = VecSet( S->Htidal_s, 0.0 ); CHKERRQ(ierr);

    // final command is always an append call to the Htot_s array
    //ierr = VecAXPY( S->Htot_s, 1.0, S->Htidal_s ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscScalar get_radiogenic_heat_production( RadionuclideParameters const Rp, PetscReal time )
{
    PetscScalar H;

    H = (Rp->t0-time) * PetscLogScalar(2.0);
    H /= Rp->half_life;
    H = PetscExpScalar(H);
    H *= Rp->heat_production * Rp->abundance * Rp->concentration;

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
    Mesh     const *M = &E->mesh;

    PetscFunctionBeginUser;

    /* convective heat flux */
    //   arr_Jconv[i] = -arr_dSdr[i] * arr_kappah[i] * arr_rho[i] * arr_temp[i];

    /* line below converts mass coordinate to physical dS/dr */
    ierr = VecPointwiseMult(S->Jconv, S->dSdxi, M->dxidr_b );CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jconv,S->Jconv, S->kappah);CHKERRQ(ierr);
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
    Mesh           *M = &E->mesh;
    Solution       *S = &E->solution;
    Parameters const P = E->parameters;
    PetscInt       i, ilo, ihi, w, should_be_two;
    PetscScalar    *arr_Jmix;
    PetscScalar    dSliqdP, dSsoldP, smth, gphi;
    const PetscScalar *arr_phi, *arr_dSdxi, *arr_kappac, *arr_rho, *arr_temp, *arr_pres, *arr_S, *arr_dPdr, *arr_dxidr;
    EOS *sub_eos;

    /* Jmix requires two phases */
    ierr = EOSCompositeGetSubEOS(P->eos, &sub_eos, &should_be_two);CHKERRQ(ierr);
    if (should_be_two!=2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Expecting two sub-EOSs");
    const EOS Ep0 = sub_eos[0];
    const EOS Ep1 = sub_eos[1];

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_b,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;

    ierr = DMDAVecGetArrayRead(da_b,S->phi,&arr_phi);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dSdxi,&arr_dSdxi);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->kappac,&arr_kappac);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->rho,&arr_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->temp,&arr_temp);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->S,&arr_S);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->pressure_b,&arr_pres);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->dPdr_b,&arr_dPdr);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->dxidr_b,&arr_dxidr);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->Jmix,&arr_Jmix);CHKERRQ(ierr);

    for(i=ilo; i<ihi; ++i){
        ierr = EOSGetPhaseBoundary( Ep0, arr_pres[i], NULL, &dSliqdP );CHKERRQ(ierr);
        ierr = EOSGetPhaseBoundary( Ep1, arr_pres[i], NULL, &dSsoldP );CHKERRQ(ierr);
        arr_Jmix[i] = arr_dSdxi[i] * arr_dxidr[i] - arr_phi[i] * dSliqdP * arr_dPdr[i];
        arr_Jmix[i] += (arr_phi[i]-1.0) * dSsoldP * arr_dPdr[i];
        arr_Jmix[i] *= -arr_kappac[i] * arr_rho[i] * arr_temp[i];

        /* (optional) smoothing across phase boundaries for two phase composite */
        ierr = EOSCompositeGetTwoPhasePhaseFractionNoTruncation(P->eos, arr_pres[i], arr_S[i], &gphi);CHKERRQ(ierr);
        {
          PetscScalar matprop_smooth_width;

          ierr = EOSCompositeGetMatpropSmoothWidth(P->eos, &matprop_smooth_width);CHKERRQ(ierr);
          smth = get_smoothing(matprop_smooth_width, gphi );
        }
        arr_Jmix[i] *= smth;

    }

    ierr = DMDAVecRestoreArrayRead(da_b,S->phi,&arr_phi);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dSdxi,&arr_dSdxi);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->kappac,&arr_kappac);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->rho,&arr_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->temp,&arr_temp);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->S,&arr_S);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->pressure_b,&arr_pres);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->dPdr_b,&arr_dPdr);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->dxidr_b,&arr_dxidr);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->Jmix,&arr_Jmix);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jmix ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

/* conductive heat flux at basic nodes */
static PetscErrorCode append_Jcond( Ctx *E )
{
    PetscErrorCode ierr;
    Solution *S = &E->solution;
    Mesh const *M = &E->mesh;
    Vec adiabat;

    PetscFunctionBeginUser;

    /* create work vector */
    ierr = VecDuplicate(M->dxidr_b,&adiabat);CHKERRQ(ierr);
    ierr = VecCopy(M->dxidr_b,adiabat);CHKERRQ(ierr);
    ierr = VecPointwiseMult(adiabat,adiabat,S->dTdxis);CHKERRQ(ierr);

    /* conductive heat flux */
    //   arr_Jcond[i] = arr_temp[i] / arr_cp[i] * arr_dSdr[i] + arr_dTdrs[i];
    //   arr_Jcond[i] *= -arr_cond[i];

    ierr = VecPointwiseDivide(S->Jcond, S->temp, S->cp);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jcond, S->Jcond, S->dSdxi);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jcond, S->Jcond, M->dxidr_b);CHKERRQ(ierr);
    ierr = VecAYPX(S->Jcond, 1.0, adiabat);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Jcond, S->Jcond, S->cond);CHKERRQ(ierr);
    ierr = VecScale(S->Jcond, -1.0);CHKERRQ(ierr);

    /* at the moment, thermophysical properties are only computed at
       the basic internal nodes, so temp/cp = 0/0 = NaN at the top
       and bottom surface */

    ierr = VecAXPY( S->Jtot, 1.0, S->Jcond ); CHKERRQ(ierr);

    ierr = VecDestroy(&adiabat);

    PetscFunctionReturn(0);

}

/* gravitational separation heat flux at basic nodes */
static PetscErrorCode append_Jgrav( Ctx *E )
{

    PetscErrorCode ierr;
    Mesh *M = &E->mesh;
    Solution *S = &E->solution;
    Parameters P = E->parameters;
    PetscScalar *arr_Jgrav, F, cond1, cond2, phi, rhol, rhos, Sliq, Ssol;
    PetscInt i,ilo_b,ihi_b,w_b,numpts_b, should_be_two;
    DM da_b = E->da_b;
    const PetscScalar *arr_phi, *arr_pres, *arr_rho, *arr_temp;
    EOS *sub_eos;

    /* Jgrav requires two phases */
    ierr = EOSCompositeGetSubEOS(P->eos, &sub_eos, &should_be_two);CHKERRQ(ierr);
    if (should_be_two!=2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Expecting two sub-EOSs");
    const EOS Ep0 = sub_eos[0];
    const EOS Ep1 = sub_eos[1];

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* loop over all basic internal nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    ierr = DMDAVecGetArrayRead(da_b,S->phi,&arr_phi);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->pressure_b,&arr_pres);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->Jgrav,&arr_Jgrav);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->rho,&arr_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->temp,&arr_temp);CHKERRQ(ierr);

    for(i=ilo_b; i<ihi_b; ++i){
        EOSEvalData eval_liq, eval_sol;

        ierr = EOSGetPhaseBoundary( Ep0, arr_pres[i], &Sliq, NULL );CHKERRQ(ierr);
        ierr = EOSGetPhaseBoundary( Ep1, arr_pres[i], &Ssol, NULL );CHKERRQ(ierr);

        ierr = EOSEval(Ep0, arr_pres[i], Sliq, &eval_liq);CHKERRQ(ierr);
        rhol = eval_liq.rho;
        ierr = EOSEval(Ep1, arr_pres[i], Ssol, &eval_sol);CHKERRQ(ierr);
        rhos = eval_sol.rho;

        cond1 = rhol / (11.993 * rhos + rhol);
        cond2 = rhol / (0.29624 * rhos + rhol);

        phi = arr_phi[i];

        /* TODO: check that F=0 for phi=0 and phi=1 */
        if(phi < cond1){
            F = 0.001*PetscPowScalar(rhos,2)*PetscPowScalar(phi,3);
            F /= PetscPowScalar(rhol,2)*(1.0-phi);
        } else if(phi > cond2){
            F = 2.0/9.0 * phi * (1.0-phi);
        } else{
            F = 5.0/7.0*PetscPowScalar(rhos,4.5)*PetscPowScalar(phi,5.5)*(1.0-phi);
            F /= PetscPowScalar( rhol+(rhos-rhol)*phi, 4.5 );
        }

        /* changing the order of these operations, or even consolidating the lines
           actually changes the behaviour of the timestepper, and makes direct
           comparison with current test data more tricky */
        // arr_Jgrav[i] = (rhol-rhos) * rho;
        arr_Jgrav[i] = rhol - rhos;
        arr_Jgrav[i] *= arr_rho[i];
        // arr_Jgrav[i] *= pref * PetscPowScalar(GRAIN,2) * GRAVITY * F;
        arr_Jgrav[i] *= Sliq - Ssol; // entropy of fusion
        arr_Jgrav[i] *= arr_temp[i];
        arr_Jgrav[i] *= PetscPowScalar(P->grain,2);
        arr_Jgrav[i] *= P->gravity;
        arr_Jgrav[i] *= F;
        arr_Jgrav[i] /= PetscPowScalar(10.0, P->eos_phases[0]->log10visc);

    }

    ierr = DMDAVecRestoreArrayRead(da_b, S->phi, &arr_phi);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, M->pressure_b, &arr_pres);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, S->Jgrav, &arr_Jgrav);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, S->rho, &arr_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, S->temp, &arr_temp);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jgrav ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

PetscErrorCode solve_surface_entropy( Ctx *E )
{
    PetscErrorCode             ierr;
    SNES                       snes;
    Vec                        x,r;
    PetscScalar                *xx, *arr_xi_b, *arr_S_b, *arr_S_s, *arr_dSdxi_b;
    DM                         da_b=E->da_b,da_s=E->da_s;
    Mesh                       *M = &E->mesh;
    Solution                   *S = &E->solution;

    PetscFunctionBeginUser;

    ierr = SNESCreate( PETSC_COMM_WORLD, &snes );CHKERRQ(ierr);

    /* Use this to address this specific SNES (nonlinear solver) from the command
       line or options file, e.g. -surfacebc_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"surfacebc_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    ierr = VecSetSizes( x, PETSC_DECIDE, 1 );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    ierr = SNESSetFunction(snes,r,objective_function_surfacebc,E);CHKERRQ(ierr);

    /* maybe guess entropy at staggered node? */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    xx[0] = -10;
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* Inform the nonlinear solver to generate a finite-difference approximation
       to the Jacobian */
    ierr = PetscOptionsSetValue(NULL,"-surfacebc_snes_mf",NULL);CHKERRQ(ierr);

    /* Turn off convergence based on step size */
    ierr = PetscOptionsSetValue(NULL,"-surfacebc_snes_stol","0");CHKERRQ(ierr);

    /* Turn off convergenced based on trust region tolerance */
    ierr = PetscOptionsSetValue(NULL,"-surfacebc_snes_trtol","0");CHKERRQ(ierr);

    ierr = PetscOptionsSetValue(NULL,"-surfacebc_snes_rtol","1.0e-9");CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-surfacebc_snes_atol","1.0e-9");CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-surfacebc_ksp_rtol","1.0e-9");CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-surfacebc_ksp_atol","1.0e-9");CHKERRQ(ierr);

    /* For solver analysis/debugging/tuning, activate a custom monitor with a flag */
    {
      PetscBool flg = PETSC_FALSE;

      ierr = PetscOptionsGetBool(NULL,NULL,"-surfacebc_snes_verbose_monitor",&flg,NULL);CHKERRQ(ierr);
      if (flg) {
        ierr = SNESMonitorSet(snes,SNESMonitorVerbose,NULL,NULL);CHKERRQ(ierr);
      }
    }

    /* Solve */
    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr); /* Picks up any additional options (note prefix) */
    ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);
    {
      SNESConvergedReason reason;
      ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
      if (reason < 0) SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
          "Nonlinear solver didn't converge: %s\n",SNESConvergedReasons[reason]);
    }

    /* double check solution */
    //objective_function_surfacebc( NULL, x, r, Ctx );

    ierr = DMDAVecGetArray(da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    arr_dSdxi_b[0] = xx[0];
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* over-rides value calculated in set_entropy_from_solution */
    arr_S_b[0] = -arr_dSdxi_b[0] * 0.5 * (arr_xi_b[1] - arr_xi_b[0]);
    arr_S_b[0] += arr_S_s[0];
 
    ierr = DMDAVecRestoreArray(da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode objective_function_surfacebc( SNES snes, Vec x, Vec f, void *ptr)
{   
    PetscErrorCode             ierr;
    const PetscScalar          *xx;
    PetscScalar                *ff; 
    PetscScalar                Ss0, dxidr0, Sb0, res, dSdxi0, radius0;
    const PetscScalar          *arr_xi_b, *arr_dPdr_b;
    Ctx                        *E = (Ctx*) ptr;
    //Atmosphere                 *A = &E->atmosphere;
    Parameters           const P = E->parameters;
    FundamentalConstants const FC = P->fundamental_constants;
    Mesh                 const *M = &E->mesh;
    Solution                   *S = &E->solution;
    AtmosphereParameters const Ap = P->atmosphere_parameters;
    DM                         da_b = E->da_b;
    EOSEvalData                eos_eval;
    
    const PetscInt ind0 = 0;
    
    PetscFunctionBeginUser;
    
    ierr = DMDAVecGetArrayRead(da_b,M->xi_b,&arr_xi_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->dPdr_b,&arr_dPdr_b); CHKERRQ(ierr);
    
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

    /* gradient we are solving for, dSdxi for the top basic node must adhere to
       the energy balance of radiation out and heat from the interior in */ 
    dSdxi0 = xx[ind0];
    
    /* get first staggered node value (store as Ss0) */
    ierr = VecGetValues(S->S_s,1,&ind0,&Ss0);CHKERRQ(ierr);
    /* surface mapping from mass coordinate to radius */
    ierr = VecGetValues(M->dxidr_b,1,&ind0,&dxidr0);CHKERRQ(ierr);
    ierr = VecGetValues(M->radius_b,1,&ind0,&radius0);CHKERRQ(ierr);

    /* based on surface gradient (which we are solving for), compute surface
       entropy (uppermost basic node) using our reconstruction */
    Sb0 = -dSdxi0 * 0.5 * (arr_xi_b[1] - arr_xi_b[0]) + Ss0;
    
    /* need material properties at this entropy and surface pressure (0 GPa).  Since these
       are lookup quantities it precludes defining a Jacobian */
    ierr = EOSEval( P->eos, 0.0, Sb0, &eos_eval );CHKERRQ(ierr);
    
    /* TODO: test only, but should be emissivity consistent with atmosphere */
    const PetscScalar emissivity = 1.0;
    
    /* radiative flux */
    res = emissivity * FC->STEFAN_BOLTZMANN * ( PetscPowScalar( eos_eval.T, 4.0 ) - PetscPowScalar( Ap->teqm, 4.0 ) );
    
    /* conductive flux (negative by definition, so positive for residual) */
    res += eos_eval.cond * (eos_eval.T / eos_eval.Cp * dSdxi0 * dxidr0 + arr_dPdr_b[ind0] * eos_eval.dTdPs);
  
    /* convective flux (negative by definition, so positive for residual) */
    PetscScalar kappah;
    ierr = GetEddyDiffusivity( eos_eval, P, radius0, dSdxi0, dxidr0, &kappah, NULL, NULL );CHKERRQ(ierr);
    res += dSdxi0 * dxidr0 * kappah * eos_eval.rho * eos_eval.T;

    /* set residual of fluxes */
    ff[ind0] = res;
    
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);
    
    ierr = DMDAVecRestoreArrayRead(da_b,M->dPdr_b,&arr_dPdr_b); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode set_interior_structure_from_solution( Ctx *E, PetscReal t, Vec sol_in )
{

    /* set all possible quantities for a given entropy structure (i.e. top staggered 
       value S0 and dS/dxi at all basic nodes excluding the top and bottom which are 
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

#if 1
    /* TODO: testing, solve for surface entropy based on boundary condition */
    ierr = solve_surface_entropy( E );CHKERRQ(ierr);
#endif

    /* set material properties and energy fluxes and sources */
    ierr = set_phase_fraction_staggered( E ); CHKERRQ(ierr);
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
