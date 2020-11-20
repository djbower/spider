#include "atmosphere.h"
#include "bc.h"
#include "monitor.h"
#include "util.h"

static PetscScalar get_viscous_mantle_cooling_rate( const Ctx *, PetscScalar );
static PetscScalar get_isothermal_surface( const Ctx * );

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

#if 0
PetscErrorCode SetCoreMantleFluxBC( Ctx *E )
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

    /* conform heat flux at the core mantle boundary to the
       imposed boundary condition */

    /* assume that the last rank contains the last basic node */
    if (rank == size-1){

        /* conform basal mantle flux to core mantle boundary condition */
        switch( P->CORE_BC ){
            case 1:
                // core cooling
                /* do nothing since core cools by mantle heat flux
                   as determined above */
                break;
            case 2:
                // constant heat flux
                ierr = VecSetValue( S->Jtot, ind_cmb, P->core_bc_value, INSERT_VALUES);CHKERRQ(ierr);
                break;
            case 3:
                // isothermal core, i.e. CMB entropy/temperature fixed
                ierr = VecSetValue( S->Jtot, ind_cmb, 0.0, INSERT_VALUES);CHKERRQ(ierr);
                break;
            default:
                SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported CORE_BC value %d provided",P->CORE_BC);
        }

        ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
        ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);

    }

    PetscFunctionReturn(0);

}
#endif

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
