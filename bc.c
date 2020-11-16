#include "atmosphere.h"
#include "bc.h"
#include "monitor.h"
#include "util.h"

static PetscScalar get_viscous_mantle_cooling_rate( const Ctx *, PetscScalar );
static PetscScalar get_isothermal_surface( const Ctx * );
static PetscScalar isothermal_or_cooling_cmb( const Ctx *, PetscScalar );
static PetscScalar get_core_cooling_factor( const Ctx * );

PetscErrorCode set_surface_flux( Ctx *E )
{

    PetscErrorCode       ierr;
    PetscMPIInt          rank;
    PetscScalar          Qout,area0;
    PetscInt             const ind0=0;
    Atmosphere           *A  = &E->atmosphere;
    Mesh                 const *M  = &E->mesh;
    Parameters           const P  = E->parameters;
    FundamentalConstants const FC = P->fundamental_constants;
    ScalingConstants     const SC  = P->scaling_constants;
    AtmosphereParameters const Ap = P->atmosphere_parameters;
    Solution             *S  = &E->solution;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

    if (!rank){

      ierr = VecGetValues(M->area_b,1,&ind0,&area0);CHKERRQ(ierr);

      /* must be after A->tsurf is set for fO2 calculation */
      /* therefore set_surface_flux always called after set_interior_structure_from_solution */
      ierr = set_reservoir_volatile_content( A, Ap, FC, SC ); CHKERRQ(ierr);

      /* determine surface flux */
      /* in all cases, compute flux and emissivity consistently */
      switch( Ap->SURFACE_BC ){
        case 1:
          // grey-body with constant emissivity
          A->emissivity = Ap->emissivity0;
          Qout = get_grey_body_flux( A, Ap, FC );
          break;
        case 2:
          // Zahnle steam atmosphere
          Qout = get_steam_atmosphere_zahnle_1988_flux( A, SC );
          A->emissivity = get_emissivity_from_flux( A, Ap, FC, Qout );
          break;
        case 3:
          // two stream approximation
          A->emissivity = get_emissivity_abe_matsui( A, Ap );
          Qout = get_grey_body_flux( A, Ap, FC );
          break;
        case 4:
          // heat flux
          Qout = Ap->surface_bc_value;
          A->emissivity = get_emissivity_from_flux( A, Ap, FC, Qout );
          break;
        case 5:
          // isothermal (constant entropy)
          // TODO: is this consistent with A->tsurf?
          Qout = get_isothermal_surface( E );
          A->emissivity = get_emissivity_from_flux( A, Ap, FC, Qout );
          break;
        default:
          SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported SURFACE_BC value %d provided",Ap->SURFACE_BC);
          break;
      }

      /* this is always the radiative flux of the atmosphere, by definition, but not
         necessarily the cooling rate of the mantle, which might be influenced by
         a near-surface lid which restricts the cooling rate */
      A->Fatm = Qout;

      /* smoothly transition the atmospheric radiation limit to the 
         viscous mantle cooling rate below the rheological transition */
      if( Ap->VISCOUS_MANTLE_COOLING_RATE ){
          Qout = get_viscous_mantle_cooling_rate( E, Qout );
      }

      /* always honour the emissivity, so ensure consistency by
         adjusting the surface temperature.  This is only relevant
         if VISCOUS_MANTLE_COOLING_RATE is set and should not change
         the surface temperature for the other options, since the
         surface temperature is already consistent with the flux
         and emissivity */
      ierr = set_surface_temperature_from_flux( A, Ap, FC ); CHKERRQ(ierr);

      /* FIXME: if we want atmospheric escape with VISCOUS_MANTLE_COOLING_RATE,
         we would need to now update the escape parameters since the surface
         temperature would have changed */

      // energy flux (Jtot)
      ierr = VecSetValue(S->Jtot,0,Qout,INSERT_VALUES);CHKERRQ(ierr);
      // energy flow (Etot)
      Qout *= area0;
      ierr = VecSetValue(S->Etot,0,Qout,INSERT_VALUES);CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Etot);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscScalar get_viscous_mantle_cooling_rate( const Ctx *E, PetscScalar Qin )
{
    PetscErrorCode ierr;
    PetscScalar    Qout;
    PetscScalar    G0, R0, R1, R2, E0, E1, E2, Q2, fwt, phi0;
    PetscInt       ind;
    Mesh           const *M = &E->mesh;
    Parameters     const P = E->parameters;
    Solution       const *S = &E->solution;

    /* below the rheological transition, the mantle cooling rate is
       not dictated by the atmosphere, since
           interior flux (100 mW) << stellar flux (1000 W)
       instead, the cooling rate of the mantle is restricted by the
       ability of the near surface lid to conduct heat to the surface,
       where it is then ``instantaneously'' radiated away. */

    /* this function helps to account for this, by only allowing the mantle
       to cool at a rate dictated by the near-surface lid.  But in doing so,
       the surface temperature is clearly wrong, because it assumes that the
       surface temperature is only "cooled" by the mantle, when in reality
       it is the atmosphere that determines the surface temperature.  And
       this should tend very quickly to the equilibrium temperature */

    /* So the surface temperature output by the model during this stage of
       evolution is not representative of the true surface temperature */

    /* for weight of different fluxes */
    ind = 0;
    ierr = VecGetValues(S->phi,1,&ind,&phi0); CHKERRQ(ierr);
    /* TODO: SWIDTH or PHI_WIDTH most appropriate choice here? */
    /* FIXME width of 1.0E-2 is hard-coded here */
    fwt = tanh_weight( phi0, P->phi_critical, 1.0E-2 );

    // energy flux from energy gradient
    ierr = VecGetValues(M->area_b,1,&ind,&G0); CHKERRQ(ierr);
    ierr = VecGetValues(M->radius_b,1,&ind,&R0); CHKERRQ(ierr);
    ind = 1;
    ierr = VecGetValues(M->radius_b,1,&ind,&R1); CHKERRQ(ierr);
    ierr = VecGetValues(S->Etot,1,&ind,&E1); CHKERRQ(ierr);
    ind = 2;
    ierr = VecGetValues(M->radius_b,1,&ind,&R2); CHKERRQ(ierr);
    ierr = VecGetValues(S->Etot,1,&ind,&E2); CHKERRQ(ierr);
    /* TODO: does not account for spherical geometry, but should
       not make a noticable difference */
    E0 = E1 - (E2-E1)*(R2-R1)/(R1-R0); // energy at surface
    Q2 = E0 / G0;
    // weight Q1 and Q2 to give total flux
    Qout = Qin * fwt + Q2 * (1.0 - fwt);

    return Qout;
}

PetscErrorCode set_core_mantle_flux_legacy( Ctx *E )
{
    /* computes the necessary core flux such that the lowermost staggered node approximates
       the desired behaviour of the core mantle boundary condition */
    /* updates Jtot and Etot */
    /* TODO: replace with improved bcs at the CMB */

    PetscErrorCode    ierr;
    PetscInt          ix,numpts_b;
    PetscScalar       fac,area1,Qin;
    PetscMPIInt       rank,size;

    Mesh        const *M = &E->mesh;
    Parameters  const P = E->parameters;
    Solution          *S = &E->solution;

    PetscFunctionBeginUser;
    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ix  = numpts_b-1; // index of last basic node

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

    /* assume that the last rank contains the last two points */
    if (rank == size-1){

      ierr = VecGetValues(M->area_b,1,&ix,&area1);CHKERRQ(ierr); // area of cmb

      switch( P->CORE_BC ){
        case 1:
          // core cooling
          fac = get_core_cooling_factor( E );
          Qin = isothermal_or_cooling_cmb( E, fac );
          break;
        case 2:
          // heat flux
          Qin = P->core_bc_value; // CMB heat flux
          break;
        case 3:
          // constant entropy (equivalent to isothermal)
          fac = 1.0;
          Qin = isothermal_or_cooling_cmb( E, fac );
          break;
        default:
          SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported CORE_BC value %d provided",P->CORE_BC);
      }

      /* energy flux (Jtot) */
      ierr = VecSetValue(S->Jtot,ix,Qin,INSERT_VALUES);CHKERRQ(ierr);
      /* energy flow (Etot) */
      Qin *= area1;
      ierr = VecSetValue(S->Etot,ix,Qin,INSERT_VALUES);CHKERRQ(ierr);

      ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);
      ierr = VecAssemblyBegin(S->Etot);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(S->Etot);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

static PetscScalar get_core_cooling_factor( const Ctx *E )
{
    PetscErrorCode ierr;
    PetscInt ix2,numpts_b;
    PetscScalar fac,vol,vol_core,rho_cmb,cp_cmb;
    Mesh const *M = &E->mesh;
    Parameters const P = E->parameters;
    Solution const *S = &E->solution;

    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ix2 = numpts_b-2; // penultimate basic node (also last staggered node)

    // staggered (last staggered node)
    ierr = VecGetValues( M->volume_s,1,&ix2,&vol);CHKERRQ(ierr);
    ierr = VecGetValues( S->rho_s,1,&ix2,&rho_cmb);CHKERRQ(ierr);
    ierr = VecGetValues( S->cp_s,1,&ix2,&cp_cmb);CHKERRQ(ierr);

    vol_core = 1.0/3.0 * PetscPowScalar(P->coresize,3.0) * PetscPowScalar(P->radius,3.0);
    fac = vol / vol_core;
    fac *= rho_cmb / P->rho_core;
    fac *= cp_cmb / P->cp_core;
    fac /= P->tfac_core_avg;
    fac = 1.0 / (1.0 + fac);

    return fac;
}

static PetscScalar get_isothermal_surface( const Ctx *E )
{
    PetscErrorCode ierr;
    PetscScalar Qout;
    PetscScalar area0,area1,jtot1; // basic
    PetscScalar htot0,vol0,rho0; // staggered
    PetscInt const ind0=0, ind1=1;
    Mesh const *M  = &E->mesh;
    Solution const *S  = &E->solution;

    // basic
    ierr = VecGetValues(M->area_b,1,&ind0,&area0);CHKERRQ(ierr);
    ierr = VecGetValues(M->area_b,1,&ind1,&area1);CHKERRQ(ierr);
    ierr = VecGetValues(S->Jtot,1,&ind1,&jtot1);CHKERRQ(ierr);
    // staggered
    ierr = VecGetValues(S->Htot_s,1,&ind0,&htot0);CHKERRQ(ierr);
    ierr = VecGetValues(M->volume_s,1,&ind0,&vol0);CHKERRQ(ierr);
    ierr = VecGetValues(S->rho_s,1,&ind0,&rho0);CHKERRQ(ierr);

    Qout = jtot1*area1 + vol0*rho0*htot0;
    Qout /= area0;

    return Qout;
}

static PetscScalar isothermal_or_cooling_cmb( const Ctx *E, PetscScalar cooling_factor )
{
    PetscErrorCode ierr;
    PetscInt ix,ix2,numpts_b;
    PetscScalar area1,area2,jtot,vol,rho_cmb,htot_cmb;
    PetscScalar Qin;
    Mesh const *M = &E->mesh;
    Solution const *S = &E->solution;

    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ix  = numpts_b-1; // last basic node
    ix2 = numpts_b-2; // penultimate basic node (also last staggered node)

    // basic
    ierr = VecGetValues(M->area_b,1,&ix,&area1);CHKERRQ(ierr);
    ierr = VecGetValues(M->area_b,1,&ix2,&area2);CHKERRQ(ierr);
    ierr = VecGetValues(S->Jtot,1,&ix2,&jtot);CHKERRQ(ierr); // flux at penultimate basic node

    // staggered
    ierr = VecGetValues( M->volume_s,1,&ix2,&vol);CHKERRQ(ierr);
    ierr = VecGetValues( S->rho_s,1,&ix2,&rho_cmb);CHKERRQ(ierr);
    ierr = VecGetValues( S->Htot_s,1,&ix2,&htot_cmb);CHKERRQ(ierr);

    Qin = jtot*area2 - vol*rho_cmb*htot_cmb;
    Qin /= area1;
    Qin *= cooling_factor;

    return Qin;
}

PetscErrorCode solve_dpdts( Ctx *E )
{
    PetscErrorCode             ierr;
    SNES                       snes;
    Vec                        x,r;
    PetscScalar                *xx;
    PetscInt                   i;
    Atmosphere                 *A  = &E->atmosphere;
    Parameters const           P  = E->parameters;
    AtmosphereParameters const Ap = P->atmosphere_parameters;

    PetscFunctionBeginUser;

    ierr = SNESCreate( PETSC_COMM_WORLD, &snes );CHKERRQ(ierr);

    /* Use this to address this specific SNES (nonlinear solver) from the command
       line or options file, e.g. -atmosic_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"atmosts_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    /* one dof per volatile, and one dof per reaction */
    ierr = VecSetSizes( x, PETSC_DECIDE, Ap->n_volatiles + Ap->n_reactions );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    ierr = SNESSetFunction(snes,r,objective_function_volatile_evolution,E);CHKERRQ(ierr);

    /* initialise vector x with initial guess */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    for (i=0; i<Ap->n_volatiles; ++i) {
        /* initial guess of dp/dt, which we solve for.  TODO: currently must be non-zero
           otherwise the time stepper results in NaN for Paolo Sossi solubility */
        /* FIXME: for "standard" power law this value works well as zero (not 0.1) */
        /* FIXME: but for Sossi solubility, this cannot be zero otherwise it gives
           NaNs in the time-stepper */
        xx[i] = 0.0;
    }
    for (i=0; i<Ap->n_reactions; ++i) {
        xx[Ap->n_volatiles + i] = 0.0;
    }
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* Inform the nonlinear solver to generate a finite-difference approximation
       to the Jacobian */
    ierr = PetscOptionsSetValue(NULL,"-atmosts_snes_mf",NULL);CHKERRQ(ierr);

    /* For solver analysis/debugging/tuning, activate a custom monitor with a flag */
    {   
      PetscBool flg = PETSC_FALSE;

      ierr = PetscOptionsGetBool(NULL,NULL,"-atmosts_snes_verbose_monitor",&flg,NULL);CHKERRQ(ierr);
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

    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    for (i=0; i<Ap->n_volatiles; ++i) {
        A->volatiles[i].dpdt = xx[i];
    }
    for (i=0; i<Ap->n_reactions; ++i) {
        A->reactions[i].dmrdt = xx[Ap->n_volatiles + i];
    }
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
