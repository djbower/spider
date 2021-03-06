#include "atmosphere.h"
#include "bc.h"
#include "interp.h"
#include "monitor.h"

PetscErrorCode set_surface_entropy_gradient_update( Ctx *E, Vec rhs )
{
    /* apply surface boundary condition in time-stepper.  This is the
       update for d/dt(dS/dr) at the surface */

    /* due to the non-linear and functional dependences of kappah,
       an update cannot be computed for most cases,  But formally
       this is the preferred approach, and cases below could be
       used in the future */

    PetscErrorCode ierr;
    Mesh const        *M = &E->mesh;
    Parameters const P = E->parameters; 
    Solution const    *S = &E->solution;
    AtmosphereParameters const Ap = P->atmosphere_parameters;
    PetscScalar    rhs_surf, rhs_below_surf, dSdt_s_surf, cp0, dSdxi0;
    const PetscScalar *arr_xi_b;
    const PetscInt ind0=0, ind1=1;
    PetscBool BC_SET = PETSC_FALSE;

    PetscFunctionBeginUser;

    ierr = DMDAVecGetArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    ierr = VecGetValues(S->dSdt_s,1,&ind0,&dSdt_s_surf);CHKERRQ(ierr);
    ierr = VecGetValues(S->dSdxi,1,&ind0,&dSdxi0);CHKERRQ(ierr);
    ierr = VecGetValues(S->cp,1,&ind0,&cp0);CHKERRQ(ierr);
    ierr = VecGetValues(rhs,1,&ind1,&rhs_below_surf);CHKERRQ(ierr);

    switch( Ap->SURFACE_BC ){
    case 1:
        /* SURFACE_BC = MO_ATMOSPHERE_TYPE_GREY_BODY: grey-body */
        break;
    case 2:
        /* SURFACE_BC = MO_ATMOSPHERE_TYPE_ZAHNLE: steam atmosphere */
        break;
    case 3:
        /* SURFACE_BC = MO_ATMOSPHERE_TYPE_VOLATILES: self-consistent atmosphere evolution
           using plane-parallel radiative equilibrium model of Abe and Matsui (1985) */
        break;
    case 4:
        /* MO_ATMOSPHERE_TYPE_HEAT_FLUX: heat flux (prescribed) */
        break;
    case 5:
        /* SURFACE_BC = MO_ATMOSPHERE_TYPE_ENTROPY: entropy */
        rhs_surf = -1.0 / (-0.5 * (arr_xi_b[1] - arr_xi_b[0]));
        rhs_surf *= dSdt_s_surf;
        BC_SET = PETSC_TRUE;
        break;
    default:
        SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported SURFACE_BC value %d provided",Ap->SURFACE_BC);
        break;
    }

    /* for most cases, we can only approximate the surface gradient update
       in absence of a better approach, assume the update at the surface
       is the same as the update at the node below */
    if( !BC_SET ){
        rhs_surf = rhs_below_surf;
    }

    ierr = VecSetValue(rhs,ind0,rhs_surf,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(rhs);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_cmb_entropy_gradient_update( Ctx *E, Vec rhs )
{
    /* apply core-mantle boundary condition in time-stepper.  This is the
       update for d/dt(dS/dr) at the cmb */

    PetscErrorCode    ierr;
    Mesh const        *M = &E->mesh;
    Parameters const  P = E->parameters;
    Solution const    *S = &E->solution;
    PetscInt          numpts_b, ind_cmb, ind_s_cmb;
    const PetscScalar *arr_xi_b;
    PetscScalar       Ecore, Etot_cmb, area_cmb, fac_cmb, rhs_cmb, dSdt_s_cmb, temp_cmb, cp_cmb;

    PetscFunctionBeginUser;
    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    ind_cmb = numpts_b-1; // index of last basic node (i.e., cmb)
    ind_s_cmb = ind_cmb-1; // index of last staggered node above cmb

    ierr = DMDAVecGetArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    ierr = VecGetValues(M->area_b,1,&ind_cmb,&area_cmb);CHKERRQ(ierr);
    ierr = VecGetValues(S->Etot,1,&ind_cmb,&Etot_cmb);CHKERRQ(ierr);
    ierr = VecGetValues(S->dSdt_s,1,&ind_s_cmb,&dSdt_s_cmb);CHKERRQ(ierr);
    ierr = VecGetValues(S->temp,1,&ind_cmb,&temp_cmb);CHKERRQ(ierr);
    ierr = VecGetValues(S->cp,1,&ind_cmb,&cp_cmb);CHKERRQ(ierr);

    /* isothermal */
    if( P->CORE_BC==3 ){
        Ecore = Etot_cmb;
    }
    /* prescribed flux from core */
    /* P->core_bc_value is set to zero if simple core cooling */
    else{
        Ecore = P->core_bc_value * area_cmb;
    }

    fac_cmb = cp_cmb / P->cp_core;
    fac_cmb /= temp_cmb * P->tfac_core_avg;
    /* recall factors of 4 pi are not included in SPIDER (only used for output) */
    fac_cmb /= 1.0/3.0 * PetscPowScalar(P->coresize,3.0) * PetscPowScalar(P->radius,3.0);
    fac_cmb /= P->rho_core;

    rhs_cmb = -Etot_cmb + Ecore;
    rhs_cmb *= fac_cmb;
    rhs_cmb -= dSdt_s_cmb;
    rhs_cmb *= 2.0 / (arr_xi_b[ind_cmb] - arr_xi_b[ind_cmb-1] );

    ierr = VecSetValue(rhs,ind_cmb,rhs_cmb,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(rhs);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_surface_flux_from_atmosphere( Ctx *E )
{
    PetscErrorCode ierr;
    Atmosphere *A = &E->atmosphere;
    Solution *S = &E->solution;
    PetscInt const ind0 = 0;

    PetscFunctionBeginUser;

    ierr = VecSetValue(S->Jtot,ind0,A->Fatm,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscScalar get_tsurf_using_parameterised_boundary_layer( PetscScalar temp, const AtmosphereParameters Ap )
{
    /* parameterisation for the ultra-thin thermal boundary layer at the surface
       of a magma ocean */

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

PetscScalar get_dtsurf_using_parameterised_boundary_layer( PetscScalar temp, const AtmosphereParameters Ap )
{
    /* derivative with respect to temperature of the surface boundary layer
       parameterisation */

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

PetscErrorCode set_cmb_entropy_constant( Ctx *E ) 
{
    /* set entropy at core-mantle boundary */
    /* an isothermal bc is naturally accommodated by the cmb bcs, so this
       is only really necessary to set the initial cmb entropy and gradient */

    PetscErrorCode   ierr;
    Mesh             *M = &E->mesh;
    Solution         *S = &E->solution;
    PetscScalar      *arr_S_s, *arr_dSdxi_b, *arr_xi_b, *arr_S_b;
    PetscInt         ihi_b, ilo_b, w_b; 
    Parameters const P = E->parameters;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b; 

    ierr = DMDAVecGetArray(E->da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(E->da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    arr_S_b[ihi_b-1] = P->ic_core_entropy;
    arr_dSdxi_b[ihi_b-1] = arr_S_b[ihi_b-1] - arr_S_s[ihi_b-2];
    arr_dSdxi_b[ihi_b-1] /= 0.5 * (arr_xi_b[ihi_b-1] - arr_xi_b[ihi_b-2]);

    ierr = DMDAVecRestoreArray(E->da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(E->da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(E->da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_surface_entropy_constant( Ctx *E ) 
{
    /* set entropy at surface */

    PetscErrorCode   ierr;
    Mesh             *M = &E->mesh;
    Solution         *S = &E->solution;
    PetscScalar      *arr_S_s, *arr_dSdxi_b, *arr_xi_b, *arr_S_b;
    Parameters const P = E->parameters;

    PetscFunctionBeginUser;

    ierr = DMDAVecGetArray(E->da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(E->da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    arr_S_b[0] = P->ic_surface_entropy;
    arr_dSdxi_b[0] = arr_S_b[0] - arr_S_s[0];
    arr_dSdxi_b[0] /= -0.5 * (arr_xi_b[1] - arr_xi_b[0]);

    ierr = DMDAVecRestoreArray(E->da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(E->da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(E->da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_cmb_entropy_from_cmb_gradient( Ctx *E )
{
    PetscErrorCode   ierr;
    Mesh             *M = &E->mesh;
    Solution         *S = &E->solution;
    PetscScalar      *arr_S_s, *arr_dSdxi_b, *arr_xi_b, *arr_S_b;
    PetscInt         ihi_b, ilo_b, w_b;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    ierr = DMDAVecGetArray(E->da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(E->da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(E->da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    /* use dS/dr at the cmb which is constrained by the core boundary condition */
    arr_S_b[ihi_b-1] = arr_dSdxi_b[ihi_b-1] * 0.5 * (arr_xi_b[ihi_b-1]-arr_xi_b[ihi_b-2]);
    arr_S_b[ihi_b-1] += arr_S_s[ihi_b-2];

    ierr = DMDAVecRestoreArray(E->da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(E->da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(E->da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_surface_entropy_from_surface_gradient( Ctx *E )
{
    PetscErrorCode  ierr;
    Mesh             *M = &E->mesh;
    Solution         *S = &E->solution;
    PetscScalar      *arr_S_s, *arr_dSdxi_b, *arr_xi_b, *arr_S_b;

    PetscFunctionBeginUser;

    ierr = DMDAVecGetArray(E->da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(E->da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(E->da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    arr_S_b[0] = -arr_dSdxi_b[0] * 0.5 * (arr_xi_b[1] - arr_xi_b[0]);
    arr_S_b[0] += arr_S_s[0];

    ierr = DMDAVecRestoreArray(E->da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(E->da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(E->da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(E->da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode get_dpdts_from_lookup( Ctx *E )
{
    PetscErrorCode             ierr;
    PetscInt                   i;
    Atmosphere                 *A  = &E->atmosphere;
    Parameters const           P  = E->parameters;
    AtmosphereParameters const Ap = P->atmosphere_parameters;

    PetscFunctionBeginUser;

    for (i=0; i<Ap->n_volatiles; ++i) {
        /* first, get dP/dT from lookup */
        ierr = SetInterp1dValue( Ap->volatile_parameters[i]->TP_interp, A->tsurf, NULL, &A->volatiles[i].dpdt );CHKERRQ(ierr);
        /* second, product rule dP/dt = dP/dT * dTsurf/dt */
        A->volatiles[i].dpdt *= A->dtsurfdt;
        /* note that the non-dimensionalisation does not affect the derivative */
        /* dlog10(P(Pa))/dT = dlog10(Phat)/dT + dlog10(P0)/dT */
        /* because the second term on the RHS is zero */
    }

    /* reactions are not compatible with pseudo-volatiles */
    for (i=0; i<Ap->n_reactions; ++i) {
        A->reactions[i].dmrdt = 0;
    }

    PetscFunctionReturn(0);
}

PetscErrorCode solve_dpdts( Ctx *E )
{
    PetscErrorCode             ierr;
    SNES                       snes;
    Vec                        x,r;
    PetscScalar                *xx,atol,rtol;
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
        /* initial guess of dp/dt, which we solve for */
        xx[i] = 0.0;
    }
    for (i=0; i<Ap->n_reactions; ++i) {
        xx[Ap->n_volatiles + i] = 0.0;
    }
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* Inform the nonlinear solver to generate a finite-difference approximation
       to the Jacobian */
    ierr = PetscOptionsSetValue(NULL,"-atmosts_snes_mf",NULL);CHKERRQ(ierr);
    /* get the accuracy from the time stepper and use for atmosts */
    ierr = PetscOptionsGetScalar(NULL,NULL,"-ts_sundials_atol",&atol,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(NULL,NULL,"-ts_sundials_rtol",&rtol,NULL);CHKERRQ(ierr);
    ierr = SNESSetTolerances(snes, atol, rtol, 0.0, PETSC_DEFAULT, PETSC_DEFAULT );CHKERRQ(ierr);
    /* need to adjust KSP tolerances to match above? Or allow to auto-determine? */

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
