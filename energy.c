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
/* permeability laws */
static PetscScalar GetPermeabilityBlakeKozenyCarman( PetscScalar grainsize, PetscScalar porosity, PetscScalar constant );
static PetscScalar GetPermeabilityRumpfGupte( PetscScalar grainsize, PetscScalar porosity, PetscScalar constant );
static PetscScalar GetPermeabilityRudge( PetscScalar grainsize, PetscScalar porosity, PetscScalar constant );

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
    PetscInt       ind_cmb, numpts_b;
    PetscMPIInt    rank, size;
    Parameters     const P = E->parameters;
    Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ind_cmb  = numpts_b-1; // index of last basic node (i.e., cmb)

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

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

    /* assume that the last rank contains the last two points */
    if (rank == size-1){

        /* conform basal mantle flux to core mantle boundary condition */
        switch( P->CORE_BC ){
            case 1:
                // core cooling
                /* do nothing since core cools by mantle heat flux
                   as determined above */
                break;
            case 2:
                // set heat flux to desired value
                ierr = VecSetValue( S->Jtot, ind_cmb, P->core_bc_value, INSERT_VALUES);CHKERRQ(ierr);
                ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
                ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);
                break;
            case 3:
                // isotherm core, i.e. no cooling
                ierr = VecSetValue( S->Jtot, ind_cmb, 0.0, INSERT_VALUES);CHKERRQ(ierr);
                ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
                ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);
                break;
            default:
                SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported CORE_BC value %d provided",P->CORE_BC);
        }
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
    PetscScalar *arr_Jgrav, F, cond1, cond2, Sliq, Ssol, dv;
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
        PetscScalar porosity;

        ierr = EOSGetPhaseBoundary( Ep0, arr_pres[i], &Sliq, NULL );CHKERRQ(ierr);
        ierr = EOSGetPhaseBoundary( Ep1, arr_pres[i], &Ssol, NULL );CHKERRQ(ierr);

        ierr = EOSEval(Ep0, arr_pres[i], Sliq, &eval_liq);CHKERRQ(ierr);
        ierr = EOSEval(Ep1, arr_pres[i], Ssol, &eval_sol);CHKERRQ(ierr);

        porosity = (eval_sol.rho - arr_rho[i]) / ( eval_sol.rho - eval_liq.rho );

        switch( P->SEPARATION ){
            case 1:
                /* Abe formulation */
                /* these switches depend on the functions below, and are constructed to
                   ensure F is a smooth function of porosity.  They are given in Abe in
                   terms of the volume fraction of solid, hence the 1.0 minus in the if
                   statement.  They depend on the choice of constants to the flow laws */
                /* See Eq. 44 in Abe (1995).  But below we use permeability directly to
                   stay connected to the physics */
                cond1 = 0.7714620383592684;
                cond2 = 0.0769618;

                /* solid_volume < cond1 (Abe) is porosity > 1-cond1 (here) */
                if(porosity > cond1){
                    /* Stokes settling factor with grainsize squared */
                    F = (2.0/9.0) * PetscPowScalar( P->grain, 2.0);
                }
                else if(porosity < cond2){
                    /* permeability includes grainsize squared */
                    F = GetPermeabilityBlakeKozenyCarman( P->grain, porosity, 1.0E-3 );
                    F /= porosity;
                }
                else{
                    /* permeability includes grainsize squared */
                    F = GetPermeabilityRumpfGupte( P->grain, porosity, 5.0/7.0 );
                    F /= porosity;
                }
                break;
            case 2:
                /* numerically it seems problematic to have no separation and then
                   try and turn it on for low melt fraction.  Hence this seems
                   to work less well than the Abe smooth approach above */
                /* apply separation below 40% melt volume fraction */
                //if( porosity < 0.4 ){
                    /* large permeability from Rudge (2018) */
                F = GetPermeabilityRudge( P->grain, porosity, 1.0/75 );
                F /= porosity;
                //}
                //else{
                //    F = 0.0;
                //}
                {
                    PetscScalar matprop_smooth_width;
                    ierr = EOSCompositeGetMatpropSmoothWidth(P->eos, &matprop_smooth_width);CHKERRQ(ierr);
                    F *= 1.0 - tanh_weight( porosity, 0.3, 1.0E-2 );
                }
                break;
            default:
                SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported SEPARATION value %d provided",P->SEPARATION);
                break;
        }

        /* relative velocity: velocity_melt - velocity_solid */
        /* works for Stokes and permeability since F is different above */
        /* e.g. Abe (1995) Eq. 39 and 40 */
        dv = ( eval_liq.rho - eval_sol.rho ) * P->gravity * F;
        dv /= PetscPowScalar(10.0, P->eos_phases[0]->log10visc);

        /* mass flux, e.g. Abe (1995) Eq. 8 */
        /* here, clear that Jgrav is zero when phi is single phase (phi=0
           or phi=1) */
        arr_Jgrav[i] = arr_rho[i] * arr_phi[i] * ( 1.0-arr_phi[i] ) * dv;

        /* energy flux */
        arr_Jgrav[i] *= Sliq - Ssol; // entropy of fusion
        arr_Jgrav[i] *= arr_temp[i];

    }

    ierr = DMDAVecRestoreArrayRead(da_b, S->phi, &arr_phi);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, M->pressure_b, &arr_pres);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, S->Jgrav, &arr_Jgrav);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, S->rho, &arr_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, S->temp, &arr_temp);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jgrav ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

PetscErrorCode set_interior_structure_from_solution( Ctx *E, PetscReal t, Vec sol_in )
{

    /* set all possible quantities for a given entropy structure (i.e. top staggered 
       value S0 and dS/dxi at all basic nodes).  This one function ensure that 
       everything is set self-consistently. */

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

static PetscScalar GetPermeabilityBlakeKozenyCarman( PetscScalar grainsize, PetscScalar porosity, PetscScalar constant )
{
    /* Abe (1995) Eq. 42a.  Also see McKenzie (1984) */

    PetscScalar kp;

    kp = PetscPowScalar( grainsize, 2.0 ) * PetscPowScalar( porosity, 3.0 ) / PetscPowScalar( 1.0-porosity, 2.0 );
    kp *= constant;

    return kp;

}

static PetscScalar GetPermeabilityRumpfGupte( PetscScalar grainsize, PetscScalar porosity, PetscScalar constant )
{
    /* Abe (1995) Eq. 42b.  Also see McKenzie (1984) */

    PetscScalar kp;

    kp = PetscPowScalar( grainsize, 2.0 ) * PetscPowScalar( porosity, 5.5 );
    kp *= constant;

    return kp;

}

static PetscScalar GetPermeabilityRudge( PetscScalar grainsize, PetscScalar porosity, PetscScalar constant )
{
    /* Rudge (2018) Eq. 7.4, large permeability for 0.1 < porosity < 0.3 */

    PetscScalar kp;

    kp = PetscPowScalar( grainsize, 2.0 ) * PetscPowScalar( porosity, 3.0 );
    kp *= constant;

    return kp;

}
