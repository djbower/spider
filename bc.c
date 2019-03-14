#include "bc.h"
#include "atmosphere.h"
#include "util.h"

static PetscScalar get_viscous_mantle_cooling_rate( Ctx *, PetscScalar );
static PetscScalar tsurf_param( PetscScalar, const AtmosphereParameters * );
static PetscScalar get_isothermal_surface( const Ctx * );
static PetscScalar isothermal_or_cooling_cmb( const Ctx *, PetscScalar );
static PetscScalar get_core_cooling_factor( const Ctx * );

PetscErrorCode set_surface_flux( Ctx *E )
{
    PetscErrorCode       ierr;
    PetscMPIInt          rank;
    PetscScalar          temp0,Qout,area0;
    PetscInt             const ind0=0;
    Atmosphere           *A  = &E->atmosphere;
    Mesh                 const *M  = &E->mesh;
    Parameters           const *P  = &E->parameters;
    Constants            const *C  = &P->constants;
    AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    Solution             *S  = &E->solution;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

    if (!rank){

      ierr = VecGetValues(M->area_b,1,&ind0,&area0);CHKERRQ(ierr);

      /* temperature (potential temperature if coarse mesh is used) */
      ierr = VecGetValues(S->temp,1,&ind0,&temp0); CHKERRQ(ierr);

      /* correct for ultra-thin thermal boundary layer at the surface */
      if( Ap->PARAM_UTBL ){
        A->tsurf = tsurf_param( temp0, Ap); // parameterised boundary layer
      }
      else{
        A->tsurf = temp0; // surface temperature is potential temperature
      }

      /* must be after A->tsurf is set */
      ierr = set_atmosphere_volatile_content( Ap, A ); CHKERRQ(ierr);

      /* determine surface flux */
      /* in all cases, compute flux and emissivity consistently */
      switch( Ap->SURFACE_BC ){
        case 1:
          // grey-body with constant emissivity
          A->emissivity = Ap->emissivity0;
          Qout = get_grey_body_flux( A, Ap );
          break;
        case 2:
          // Zahnle steam atmosphere
          Qout = get_steam_atmosphere_zahnle_1988_flux( A, C );
          A->emissivity = get_emissivity_from_flux( A, Ap, Qout );
          break;
        case 3:
          // two stream approximation
          A->emissivity = get_emissivity_abe_matsui( Ap, A );
          Qout = get_grey_body_flux( A, Ap );
          break;
        case 4:
          // heat flux
          Qout = Ap->surface_bc_value;
          A->emissivity = get_emissivity_from_flux( A, Ap, Qout );
          break;
        case 5:
          // isothermal (constant entropy)
          // TODO: is this consistent with A->tsurf?
          Qout = get_isothermal_surface( E );
          A->emissivity = get_emissivity_from_flux( A, Ap, Qout );
          break;
        default:
          SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported SURFACE_BC value %d provided",Ap->SURFACE_BC);
          break;
      }

      /* smoothly transition the cooling rate to the viscous mantle
         cooling rate below the rheological transition */
      if( Ap->VISCOUS_MANTLE_COOLING_RATE ){
          Qout = get_viscous_mantle_cooling_rate( E, Qout );
      }

      /* to ensure conservation of energy at the interface of the
         interior and atmosphere, the fluxes must be equal */
      A->Fatm = Qout;

      /* always honour the emissivity, so ensure consistency by
         adjusting the surface temperature.  This is only relevant
         if VISCOUS_MANTLE_COOLING_RATE is set and should not change
         the surface temperature for the other options, since the
         surface temperature is already consistent with the flux
         and emissivity */
      ierr = set_surface_temperature_from_flux( A, Ap ); CHKERRQ(ierr);

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

static PetscScalar get_viscous_mantle_cooling_rate( Ctx *E, PetscScalar Qin )
{
    PetscErrorCode ierr;
    PetscScalar    Qout;
    PetscScalar    G0, R0, R1, R2, E0, E1, E2, Q2, fwt, phi0;
    PetscInt       ind;
    Mesh           *M = &E->mesh;
    Parameters     *P = &E->parameters;
    Solution       *S = &E->solution;

    /* enable the ability for the magma ocean to cool at a rate dictated 
       by the upper mantle cooling rate.  This helps to prevent a viscous
       lid from forming at the top */

    /* for weight of different fluxes */
    ind = 0;
    ierr = VecGetValues(S->phi,1,&ind,&phi0); CHKERRQ(ierr);
    /* TODO: SWIDTH or PHI_WIDTH most appropriate choice here? */
    fwt = tanh_weight( phi0, P->phi_critical, P->matprop_smooth_width );

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

PetscErrorCode set_core_mantle_flux( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          ix,numpts_b; // ix2 // see TODO below
    PetscScalar       fac,area1,Qin;//,temp_cmb; // see TODO below
    PetscMPIInt       rank,size;

    Mesh        const *M = &E->mesh;
    Parameters  const *P = &E->parameters;
    Solution          *S = &E->solution;

    PetscFunctionBeginUser;
    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ix  = numpts_b-1; // index of last basic node
    //ix2 = ix-1; // see TODO below

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
          // TODO (optional?): add switch
          //ierr = VecGetValues(S->temp_s,1,&ix2,&temp_cmb);
          //if(temp_cmb<1.0){
          //    fac = 1.0;
          //    Qin = isothermal_or_cooling_cmb( E, fac );
          //}
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
    Parameters const *P = &E->parameters;
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

PetscScalar tsurf_param( PetscScalar temp, const AtmosphereParameters *Ap )
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

PetscErrorCode solve_dxdts( Ctx *E )
{
    PetscErrorCode ierr;
    SNES           snes;
    Vec            x,r;
    PetscScalar    *xx;

    Atmosphere                 *A = &E->atmosphere;
    // TODO: seems like below are not used anymore
    //Parameters           const *P = &E->parameters;
    //AtmosphereParameters const *Ap = &P->atmosphere_parameters;
    //VolatileParameters   const *CO2_parameters = &Ap->CO2_parameters;
    //VolatileParameters   const *H2O_parameters = &Ap->H2O_parameters;
    Volatile                   *CO2 = &A->CO2;
    Volatile                   *H2O = &A->H2O;

    PetscFunctionBeginUser;

    ierr = SNESCreate( PETSC_COMM_WORLD, &snes );CHKERRQ(ierr);

    /* Use this to address this specific SNES (nonlinear solver) from the command
       line or options file, e.g. -atmosic_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"atmosts_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    ierr = VecSetSizes( x, PETSC_DECIDE, 2 );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    ierr = SNESSetFunction(snes,r,FormFunction2,E);CHKERRQ(ierr);

    /* initialise vector x with initial guess */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    xx[0] = 0.0;
    xx[1] = 0.0;
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* Inform the nonlinear solver to generate a finite-difference approximation
       to the Jacobian */
    ierr = PetscOptionsSetValue(NULL,"-atmosts_snes_mf",NULL);CHKERRQ(ierr);

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
    CO2->dxdt = xx[0];
    H2O->dxdt = xx[1];
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    // FIXME: from copying, but update for this function
    /* Sanity check on solution (since it's non-unique) */
    //if (CO2->x < 0.0 || H2O->x < 0.0) SETERRQ2(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
    //    "Unphysical initial volatile concentrations: CO2: %g, H2O: %g",CO2->x,H2O->x);

    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
