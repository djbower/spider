#include "bc.h"
#include "energy.h"
#include "eos.h"
#include "eos_composite.h"
#include "matprop.h"
#include "monitor.h"
#include "twophase.h"
#include "util.h"

/* heat fluxes */
static PetscErrorCode set_Jtot( Ctx * );
static PetscErrorCode append_Jcond( Ctx * );
static PetscErrorCode append_Jconv( Ctx * );
static PetscErrorCode append_Jmix( Ctx * );
static PetscErrorCode append_Jgrav( Ctx * );
static PetscErrorCode append_Hradio( Ctx *, PetscReal );
static PetscErrorCode append_Htidal( Ctx *, PetscReal );
/* next enforces a boundary condition, but to avoid circular
   dependencies it lives here and not in bc.c */
static PetscErrorCode solve_for_surface_radiation_balance( Ctx *, PetscReal );
static PetscErrorCode objective_function_surface_radiation_balance( SNES, Vec, Vec, void *);
static PetscErrorCode set_current_state( Ctx *, PetscReal );
static PetscScalar get_radiogenic_heat_production( RadionuclideParameters const, PetscReal );
/* permeability laws control efficiency of gravitational separation */
static PetscScalar GetPermeabilityBlakeKozenyCarman( PetscScalar grainsize, PetscScalar porosity, PetscScalar constant );
static PetscScalar GetPermeabilityRumpfGupte( PetscScalar grainsize, PetscScalar porosity, PetscScalar constant );
static PetscScalar GetPermeabilityRudge( PetscScalar grainsize, PetscScalar porosity, PetscScalar constant );

PetscErrorCode set_Htot( Ctx *E, PetscReal time )
{
    /* total internal heat generation */

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

static PetscErrorCode append_Hradio( Ctx *E, PetscReal time )
{
    /* radiogenic heat generation */

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

static PetscErrorCode append_Htidal( Ctx *E, PetscReal tyrs )
{
    /* tidal heat generation */

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

PetscErrorCode set_Etot( Ctx *E )
{
    /* total energy flow (flux*area) at basic nodes */

    PetscErrorCode ierr;
    Mesh           *M = &E->mesh;
    Solution       *S = &E->solution;

    PetscFunctionBeginUser;

    ierr = set_Jtot(E); CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->Etot,S->Jtot,M->area_b); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_Jtot( Ctx *E )
{
    /* total heat flux at basic nodes */

    /* also ensures the core mantle boundary flux adheres to the
       imposed boundary condition */

    PetscErrorCode ierr;
    Parameters     const P = E->parameters;
    Solution       *S = &E->solution;
    AtmosphereParameters const Ap = P->atmosphere_parameters;

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

    /* if we are using the simple (less accurate) surface boundary
       condition, then we need to impose the atmospheric flux at the
       top of the interior */
    /* TODO: need an exception for the isothermal bcs? */
    if( !Ap->SURFACE_BC_ACC ){
        ierr = set_surface_flux_from_atmosphere( E );CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

static PetscErrorCode append_Jconv( Ctx *E )
{
    /* convective heat flux at basic nodes */

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

PetscScalar GetConvectiveHeatFlux( Ctx *E, PetscInt * ind_ptr)
{
    PetscErrorCode ierr;
    PetscScalar    dSdxi,dxidr,temp,rho,kappah,Jconv;
    Solution const *S = &E->solution;
    Mesh const     *M = &E->mesh;
    PetscInt       ind_cmb,ind_abv_cmb,ilo_b,ihi_b,w_b;
    PetscInt const ind1=1;

    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b; /* total number of basic nodes */
    ind_cmb = ihi_b-1;
    ind_abv_cmb = ind_cmb-1;

    /* for surface and cmb, use kappah at the next basic node in to avoid
       numerical issues.  kappah is a nasty non-linear function of the entropy
       gradient and other material properties */

    /* surface, use kappah at node below */
    if(!*ind_ptr){
        ierr = VecGetValues(S->kappah,1,&ind1,&kappah);CHKERRQ(ierr);
    }
    /* cmb, use kappah at node above */
    else if(*ind_ptr == ind_cmb){
        ierr = VecGetValues(S->kappah,1,&ind_abv_cmb,&kappah);CHKERRQ(ierr);
    }
    /* otherwise use kappah at this node */
    else{
        ierr = VecGetValues(S->kappah,1,ind_ptr,&kappah);CHKERRQ(ierr);
    }

    ierr = VecGetValues(S->rho,1,ind_ptr,&rho);CHKERRQ(ierr);
    ierr = VecGetValues(S->dSdxi,1,ind_ptr,&dSdxi);CHKERRQ(ierr);
    ierr = VecGetValues(M->dxidr_b,1,ind_ptr,&dxidr);CHKERRQ(ierr);
    ierr = VecGetValues(S->temp,1,ind_ptr,&temp);CHKERRQ(ierr);

    Jconv = -dSdxi * dxidr * kappah * rho * temp;

    return Jconv;
}

static PetscErrorCode append_Jmix( Ctx *E )
{
    /* mixing heat flux (latent heat transport) at basic nodes */

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

PetscScalar GetMixingHeatFlux( Ctx *E, PetscInt * ind_ptr )
{
    PetscErrorCode  ierr;
    PetscScalar     dSdxi,dxidr,temp,rho,kappac,phi,pres,Sval,Jmix,dPdr,dSliqdP,dSsoldP,gphi,smth;
    PetscInt        should_be_two;
    Mesh const      *M = &E->mesh;
    Solution const  *S = &E->solution;
    Parameters const P = E->parameters;
    EOS *sub_eos;
    PetscInt        ilo_b,ihi_b,w_b;

    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b; // this is one more than last index of basic array

    /* no energy transport by mixing at boundaries since no mass transfer */
    if( (!*ind_ptr) || (*ind_ptr == ihi_b-1) ){
        return 0.0; 
    }   

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

static PetscErrorCode append_Jcond( Ctx *E )
{
    /* conductive heat flux at basic nodes */

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

PetscScalar GetConductiveHeatFlux( Ctx *E, PetscInt * ind_ptr)
{
    PetscErrorCode ierr;
    PetscScalar    dSdxi,dxidr,temp,cp,dTdxis,cond,Jcond;
    Solution const *S = &E->solution;
    Mesh const     *M = &E->mesh;
    PetscInt       ilo_b,ihi_b,w_b;

    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b; // this is one more than last index of basic array

    ierr = VecGetValues(S->dSdxi,1,ind_ptr,&dSdxi);CHKERRQ(ierr);
    ierr = VecGetValues(M->dxidr_b,1,ind_ptr,&dxidr);CHKERRQ(ierr);
    ierr = VecGetValues(S->temp,1,ind_ptr,&temp);CHKERRQ(ierr);
    ierr = VecGetValues(S->cp,1,ind_ptr,&cp);CHKERRQ(ierr);
    ierr = VecGetValues(S->dTdxis,1,ind_ptr,&dTdxis);CHKERRQ(ierr);
    /* conductivity is constant, so this is OK for the boundary nodes
       which also take kappah from the nodes inset */
    ierr = VecGetValues(S->cond,1,ind_ptr,&cond);CHKERRQ(ierr);

    Jcond = ( temp / cp ) * dSdxi * dxidr + dTdxis * dxidr;
    Jcond *= -cond;

    return Jcond;
}

static PetscErrorCode append_Jgrav( Ctx *E )
{
    /* gravitational separation heat flux at basic nodes */

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

PetscScalar GetGravitationalHeatFlux( Ctx *E, PetscInt * ind_ptr )
{
    PetscErrorCode  ierr;
    PetscScalar     porosity,cond1,cond2,F,dv,pres,rho,Sliq,Ssol,phi,Jgrav,temp;
    PetscInt        should_be_two;
    Solution const  *S = &E->solution;
    Mesh const      *M = &E->mesh;
    Parameters const P = E->parameters;
    EOS              *sub_eos;
    EOSEvalData     eval_liq, eval_sol;
    PetscInt        ilo_b,ihi_b,w_b;

    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b; // this is one more than last index of basic array

    /* no gravitational separation at boundaries (since no mass transfer possible) */
    if( (!*ind_ptr) || (*ind_ptr == ihi_b-1) ){
        return 0.0;
    }

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

static PetscErrorCode solve_for_surface_radiation_balance( Ctx *E, PetscReal t ) 
{
    /* to formally balance radiation with the interior heat flux at the surface,
       we must solve a coupled system for the surface entropy gradient and
       surface entropy */

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
       line or options file, e.g. -surfrad_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"surfrad_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    ierr = VecSetSizes( x, PETSC_DECIDE, 1 );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    /* set current time to Ctx, since the Ctx is passed to the
       objective function with t updated */
    E->t = t;

    ierr = SNESSetFunction(snes,r,objective_function_surface_radiation_balance,E);CHKERRQ(ierr);

    /* initial guess of surface entropy gradient is gradient at basic node below surface */
    ierr = DMDAVecGetArrayRead(da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    xx[0] = arr_dSdxi_b[1];
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);

    ierr = PetscOptionsSetValue(NULL,"-surfrad_snes_mf",NULL);CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-surfrad_snes_stol","0");CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-surfrad_snes_rtol","1.0E-9");CHKERRQ(ierr);
    /* atol will give accurate result to within 0.001 W/m^2.  Could
       likely relax this further, at least for the magma ocean stage */
    ierr = PetscOptionsSetValue(NULL,"-surfrad_snes_atol","1.0E-3");CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-surfrad_ksp_rtol","1.0E-9");CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL,"-surfrad_ksp_atol","1.0E-9");CHKERRQ(ierr);

    /* For solver analysis/debugging/tuning, activate a custom monitor with a flag */
    {   
      PetscBool flg = PETSC_FALSE;

      ierr = PetscOptionsGetBool(NULL,NULL,"-surfrad_snes_verbose_monitor",&flg,NULL);CHKERRQ(ierr);
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

    ierr = DMDAVecGetArray(da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    /* set entropy gradient at surface */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    arr_dSdxi_b[0] = xx[0];
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* set entropy at surface */
    arr_S_b[0] = -arr_dSdxi_b[0] * 0.5 * (arr_xi_b[1] - arr_xi_b[0]);
    arr_S_b[0] += arr_S_s[0];

    ierr = DMDAVecRestoreArray(da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode objective_function_surface_radiation_balance( SNES snes, Vec x, Vec f, void *ptr)
{
    PetscErrorCode             ierr;
    PetscMPIInt                rank;
    const PetscScalar          *xx;
    PetscScalar                *ff;
    PetscScalar                Ss0, Sb0, dSdxi0, Jtot0, res;
    const PetscScalar          *arr_xi_b;
    Ctx                        *E = (Ctx*) ptr;
    Parameters            const P  = E->parameters;
    Mesh                 const *M = &E->mesh;
    Atmosphere           const *A = &E->atmosphere;
    ScalingConstants      const SC = P->scaling_constants;
    Solution                   *S = &E->solution;
    DM                         da_b = E->da_b;
    PetscInt             const ind0 = 0;
    PetscInt             const ind1 = 1;
    PetscReal                  t;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

    /* set current time from Ctx, since the Ctx is passed to the
       objective function */
    t = E->t;

    ierr = DMDAVecGetArrayRead(da_b,M->xi_b,&arr_xi_b); CHKERRQ(ierr);
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

    /* conform entropy Vecs in struct to our current guess of the 
       entropy gradient at the surface */
    dSdxi0 = xx[ind0];
    ierr = VecSetValue( S->dSdxi, ind0, dSdxi0, INSERT_VALUES );CHKERRQ(ierr);
    ierr = VecAssemblyBegin( S->dSdxi );CHKERRQ(ierr);
    ierr = VecAssemblyEnd( S->dSdxi );CHKERRQ(ierr);

    /* compute surface entropy using the reconstruction */
    ierr = VecGetValues(S->S_s,1,&ind0,&Ss0);CHKERRQ(ierr);
    Sb0 = -dSdxi0 * 0.5 * (arr_xi_b[ind1] - arr_xi_b[ind0]) + Ss0;
    ierr = VecSetValue( S->S, ind0, Sb0, INSERT_VALUES );CHKERRQ(ierr);
    ierr = VecAssemblyBegin( S->S );CHKERRQ(ierr);
    ierr = VecAssemblyEnd( S->S );CHKERRQ(ierr);

    ierr = set_current_state( E, t );CHKERRQ(ierr);

    /* compute residual */
    ierr = VecGetValues(S->Jtot,1,&ind0,&Jtot0);CHKERRQ(ierr);
    res = A->Fatm - Jtot0;

    /* scale residual to physical flux (W/m^2), so the solver
       tolerances effectively enforce a minimum flux difference
       in W/m^2 */
    res *= SC->FLUX;

    /* set residual of fluxes */
    ff[ind0] = res;

    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_current_state( Ctx *E, PetscReal t )
{
    /* for a given entropy profile (all basic and staggered nodes),
       set the current state of the system for the interior and
       atmosphere */

    PetscErrorCode        ierr;
    PetscMPIInt           rank;
    Atmosphere            *A  = &E->atmosphere;
    Parameters            const P  = E->parameters;
    FundamentalConstants  const FC = P->fundamental_constants;
    ScalingConstants      const SC  = P->scaling_constants;
    AtmosphereParameters  const Ap = P->atmosphere_parameters;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

    /* material properties at basic and staggered nodes */
    ierr = set_capacitance_staggered( E );CHKERRQ(ierr);
    ierr = set_matprop_basic( E );CHKERRQ(ierr);

    /* update surface temperature, fO2 */
    ierr = set_interior_atmosphere_interface_from_surface_entropy( E );CHKERRQ(ierr);

    /* update all volatile reservoirs */
    ierr = set_reservoir_volatile_content( A, Ap, FC, SC ); CHKERRQ(ierr);

    /* compute A->emissivity and A->Fatm */
    ierr = set_atmosphere_emissivity_and_flux( A, Ap, FC, SC );CHKERRQ(ierr);

    ierr = set_Etot( E );CHKERRQ(ierr);
    ierr = set_Htot( E, t );CHKERRQ(ierr);

    ierr = set_rheological_front( E ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_interior_atmosphere_interface_from_surface_entropy( Ctx *E )
{
    PetscErrorCode       ierr;
    Parameters           const P = E->parameters;
    Mesh                 const *M  = &E->mesh;
    AtmosphereParameters const Ap = P->atmosphere_parameters;
    ScalingConstants     const SC = P->scaling_constants;
    Atmosphere           *A = &E->atmosphere;
    Solution             const *S  = &E->solution;
    PetscScalar          T0, Sb0, Pres0;
    EOSEvalData          eos_eval;
    PetscInt             const ind0 = 0;

    PetscFunctionBeginUser;

    ierr = set_phase_fraction_staggered( E ); CHKERRQ(ierr);
    /* to set atmosphere. we need to know A->Mliq and A->Msol */
    ierr = set_Mliq( E );CHKERRQ(ierr);
    ierr = set_Msol( E );CHKERRQ(ierr);

    ierr = VecGetValues( S->S, 1, &ind0, &Sb0 );
    ierr = VecGetValues( M->pressure_b, 1, &ind0, &Pres0 );
    ierr = EOSEval( P->eos, Pres0, Sb0, &eos_eval );
    T0 = eos_eval.T;

    /* correct for ultra-thin thermal boundary layer at the surface */
    if( Ap->PARAM_UTBL ){
        A->tsurf = get_tsurf_using_parameterised_boundary_layer( T0, Ap); // parameterised boundary layer
        A->dtsurfdt = get_dtsurf_using_parameterised_boundary_layer( T0, Ap); // dTsurf/dT
    }
    else{
        A->tsurf = T0;
        A->dtsurfdt = 1.0; // dTsurf/dT
    }

    /* must be after A->tsurf is set for fO2 calculation */
    if( Ap->OXYGEN_FUGACITY ){
        ierr = set_oxygen_fugacity( A, Ap, SC );CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode set_current_state_from_solution( Ctx *E, PetscReal t, Vec sol_in )
{
    /* set all possible quantities for a given sol Vec.  This one
       function ensures that everything is set self-consistently. */

    PetscErrorCode       ierr;
    Parameters           const P  = E->parameters;
    AtmosphereParameters const Ap = P->atmosphere_parameters;

    PetscFunctionBeginUser;

    /* set solution in the relevant structs */
    /* entropy */
    ierr = set_entropy_from_solution( E, sol_in );CHKERRQ(ierr);
    /* atmosphere */
    ierr = set_partial_pressures_from_solution( E, sol_in );CHKERRQ(ierr);

    /* set_entropy_from_solution has already extrapolated to get an
       estimate of the surface entropy and entropy gradient. */

    /* we ensure consistency with the surface boundary condition here */

    /* for an isothermal surface, we can simply impose the surface
       entropy and entropy gradient */
    if( Ap->SURFACE_BC == 5 ){
        ierr = set_surface_entropy_constant( E );CHKERRQ(ierr);
        ierr = set_current_state( E, t);CHKERRQ(ierr);
    }

    if( Ap->SURFACE_BC_ACC ){
        /* must solve for entropy and entropy gradient at surface
           to adhere to the radiation-interior flux balance */
        ierr = solve_for_surface_radiation_balance( E, t );CHKERRQ(ierr);
    }
    else{
        ierr = set_current_state( E, t);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
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
