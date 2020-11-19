#include "atmosphere.h"
#include "bc.h"
#include "energy.h"
#include "eos.h"
#include "eos_composite.h"
#include "matprop.h"
#include "twophase.h"
#include "util.h"
#include "monitor.h"

/* heat fluxes */
static PetscScalar GetConductiveHeatFlux( Ctx *, PetscInt *);
static PetscScalar GetConvectiveHeatFlux( Ctx *, PetscInt *);
static PetscScalar GetMixingHeatFlux( Ctx *, PetscInt *);
static PetscScalar GetGravitationalHeatFlux( Ctx *, PetscInt * );

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

    /* also ensures the core mantle boundary flux adheres to the
       imposed boundary condition */

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

    ierr = SetCoreMantleFluxBC( E );CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

/* convective heat flux at basic nodes */
static PetscErrorCode append_Jconv( Ctx *E )
{
    PetscErrorCode ierr;
    PetscScalar    Jconv;
    PetscInt       i,ilo_b,ihi_b,w_b;
    Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    /* loop over all basic internal nodes */
    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    for(i=ilo_b; i<ihi_b; ++i){
        Jconv = GetConvectiveHeatFlux( E, &i );
        ierr = VecSetValue(S->Jconv,i,Jconv,INSERT_VALUES);CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(S->Jconv);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jconv);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jconv ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscScalar GetConvectiveHeatFlux( Ctx *E, PetscInt * ind_ptr)
{
    PetscErrorCode ierr;
    PetscScalar    dSdxi,dxidr,temp,rho,kappah,Jconv;
    Solution const *S = &E->solution;
    Mesh const     *M = &E->mesh;

    ierr = VecGetValues(S->dSdxi,1,ind_ptr,&dSdxi);CHKERRQ(ierr);
    ierr = VecGetValues(M->dxidr_b,1,ind_ptr,&dxidr);CHKERRQ(ierr);
    ierr = VecGetValues(S->temp,1,ind_ptr,&temp);CHKERRQ(ierr);
    ierr = VecGetValues(S->rho,1,ind_ptr,&rho);CHKERRQ(ierr);
    ierr = VecGetValues(S->kappah,1,ind_ptr,&kappah);CHKERRQ(ierr);

    Jconv = -dSdxi * dxidr * kappah * rho * temp;

    return Jconv;
}


/* mixing heat flux (latent heat transport) at basic nodes */
static PetscErrorCode append_Jmix( Ctx *E )
{
    PetscErrorCode ierr;
    PetscScalar    Jmix;
    Solution       *S = &E->solution;
    PetscInt       i,ilo_b,ihi_b,w_b;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    for(i=ilo_b; i<ihi_b; ++i){
        Jmix = GetMixingHeatFlux( E, &i );
        ierr = VecSetValue(S->Jmix,i,Jmix,INSERT_VALUES);CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(S->Jmix);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jmix);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jmix ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscScalar GetMixingHeatFlux( Ctx *E, PetscInt * ind_ptr )
{
    PetscErrorCode  ierr;
    PetscScalar     dSdxi,dxidr,temp,rho,kappac,phi,pres,Sval,Jmix,dPdr,dSliqdP,dSsoldP,gphi,smth;
    PetscInt        should_be_two;
    Mesh const      *M = &E->mesh;
    Solution const  *S = &E->solution;
    Parameters const P = E->parameters;
    EOS *sub_eos;

    /* Jmix requires two phases */
    ierr = EOSCompositeGetSubEOS(P->eos, &sub_eos, &should_be_two);CHKERRQ(ierr);
    if (should_be_two!=2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Expecting two sub-EOSs");
    const EOS Ep0 = sub_eos[0];
    const EOS Ep1 = sub_eos[1];

    ierr = VecGetValues(S->dSdxi,1,ind_ptr,&dSdxi);CHKERRQ(ierr);
    ierr = VecGetValues(M->dxidr_b,1,ind_ptr,&dxidr);CHKERRQ(ierr);
    ierr = VecGetValues(M->pressure_b,1,ind_ptr,&pres);CHKERRQ(ierr);
    ierr = VecGetValues(S->temp,1,ind_ptr,&temp);CHKERRQ(ierr);
    ierr = VecGetValues(S->rho,1,ind_ptr,&rho);CHKERRQ(ierr);
    ierr = VecGetValues(S->kappac,1,ind_ptr,&kappac);CHKERRQ(ierr);
    ierr = VecGetValues(S->phi,1,ind_ptr,&phi);CHKERRQ(ierr);
    ierr = VecGetValues(M->dPdr_b,1,ind_ptr,&dPdr);CHKERRQ(ierr);
    ierr = VecGetValues(S->S,1,ind_ptr,&Sval);CHKERRQ(ierr);

    ierr = EOSGetPhaseBoundary( Ep0, pres, NULL, &dSliqdP );CHKERRQ(ierr);
    ierr = EOSGetPhaseBoundary( Ep1, pres, NULL, &dSsoldP );CHKERRQ(ierr);

    Jmix = dSdxi * dxidr - phi * dSliqdP * dPdr;
    Jmix += (phi-1.0) * dSsoldP * dPdr;
    Jmix *= -kappac * rho * temp;

    /* (optional) smoothing across phase boundaries for two phase composite */
    ierr = EOSCompositeGetTwoPhasePhaseFractionNoTruncation(P->eos, pres, Sval, &gphi);CHKERRQ(ierr);
    { 
      PetscScalar matprop_smooth_width;
      
      ierr = EOSCompositeGetMatpropSmoothWidth(P->eos, &matprop_smooth_width);CHKERRQ(ierr);
      smth = get_smoothing(matprop_smooth_width, gphi );
    }

    Jmix *= smth;

    return Jmix;
}

/* conductive heat flux at basic nodes */
static PetscErrorCode append_Jcond( Ctx *E )
{
    PetscErrorCode ierr;
    PetscScalar    Jcond;
    PetscInt       i,ilo_b,ihi_b,w_b;
    Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    /* loop over all basic internal nodes */
    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    for(i=ilo_b; i<ihi_b; ++i){
        Jcond = GetConductiveHeatFlux( E, &i );
        ierr = VecSetValue(S->Jcond,i,Jcond,INSERT_VALUES);CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(S->Jcond);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jcond);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jcond ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscScalar GetConductiveHeatFlux( Ctx *E, PetscInt * ind_ptr)
{
    PetscErrorCode ierr;
    PetscScalar    dSdxi,dxidr,temp,cp,dTdxis,cond,Jcond;
    Solution const *S = &E->solution;
    Mesh const     *M = &E->mesh;

    ierr = VecGetValues(S->dSdxi,1,ind_ptr,&dSdxi);CHKERRQ(ierr);
    ierr = VecGetValues(M->dxidr_b,1,ind_ptr,&dxidr);CHKERRQ(ierr);
    ierr = VecGetValues(S->temp,1,ind_ptr,&temp);CHKERRQ(ierr);
    ierr = VecGetValues(S->cp,1,ind_ptr,&cp);CHKERRQ(ierr);
    ierr = VecGetValues(S->dTdxis,1,ind_ptr,&dTdxis);CHKERRQ(ierr);
    ierr = VecGetValues(S->cond,1,ind_ptr,&cond);CHKERRQ(ierr);

    Jcond = ( temp / cp ) * dSdxi * dxidr + dTdxis * dxidr;
    Jcond *= -cond;

    return Jcond;
}

/* gravitational separation heat flux at basic nodes */
static PetscErrorCode append_Jgrav( Ctx *E )
{

    PetscErrorCode ierr;
    PetscScalar    Jgrav;
    PetscInt       i,ilo_b,ihi_b,w_b;
    Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    /* loop over all basic internal nodes */
    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    for(i=ilo_b; i<ihi_b; ++i){
        Jgrav = GetGravitationalHeatFlux( E, &i );
        ierr = VecSetValue(S->Jgrav,i,Jgrav,INSERT_VALUES);CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(S->Jgrav);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jgrav);CHKERRQ(ierr);

    ierr = VecAXPY( S->Jtot, 1.0, S->Jgrav ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscScalar GetGravitationalHeatFlux( Ctx *E, PetscInt * ind_ptr )
{
    PetscErrorCode  ierr;
    PetscScalar     porosity,cond1,cond2,F,dv,pres,rho,Sliq,Ssol,phi,Jgrav,temp;
    PetscInt        should_be_two;
    Solution const  *S = &E->solution;
    Mesh const      *M = &E->mesh;
    Parameters const P = E->parameters;
    EOS              *sub_eos;
    EOSEvalData     eval_liq, eval_sol;

    /* Jgrav requires two phases */
    ierr = EOSCompositeGetSubEOS(P->eos, &sub_eos, &should_be_two);CHKERRQ(ierr);
    if (should_be_two!=2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Expecting two sub-EOSs");
    const EOS Ep0 = sub_eos[0];
    const EOS Ep1 = sub_eos[1];

    ierr = VecGetValues(M->pressure_b,1,ind_ptr,&pres);CHKERRQ(ierr);
    ierr = VecGetValues(S->rho,1,ind_ptr,&rho);CHKERRQ(ierr);
    ierr = VecGetValues(S->phi,1,ind_ptr,&phi);CHKERRQ(ierr);
    ierr = VecGetValues(S->temp,1,ind_ptr,&temp);CHKERRQ(ierr);

    ierr = EOSGetPhaseBoundary( Ep0, pres, &Sliq, NULL );CHKERRQ(ierr);
    ierr = EOSEval(Ep0, pres, Sliq, &eval_liq);CHKERRQ(ierr);
    ierr = EOSGetPhaseBoundary( Ep1, pres, &Ssol, NULL );CHKERRQ(ierr);
    ierr = EOSEval(Ep1, pres, Ssol, &eval_sol);CHKERRQ(ierr);

    porosity = (eval_sol.rho - rho) / ( eval_sol.rho - eval_liq.rho );

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
    Jgrav = rho * phi * ( 1.0 - phi ) * dv;

    /* energy flux */
    Jgrav *= temp * ( Sliq - Ssol ); // enthalpy

    return Jgrav;
}

PetscErrorCode set_interior_structure_from_solution( Ctx *E, PetscReal t, Vec sol_in )
{

    /* set all possible quantities for a given entropy structure (i.e. top staggered 
       value S0 and dS/dxi at all basic nodes).  This one function ensures that 
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
