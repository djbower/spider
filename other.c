/* extra prototype code previously written but not current used */

static PetscErrorCode SetInitialCMBdSdxiFromFlux(Ctx *E)
{
    PetscErrorCode ierr;
    SNES snes;
    Vec x, r;
    PetscScalar dSdxi;
    PetscInt numpts_b, ind0, ind_cmb;
    Solution *S = &E->solution;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_b, NULL, &numpts_b, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    CHKERRQ(ierr);
    ind_cmb = numpts_b - 1; // index of last basic node (i.e., cmb)
    ind0 = 0;

    ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
    CHKERRQ(ierr);
    ierr = SNESSetOptionsPrefix(snes, "cmbic_");
    CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, 1);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);
    CHKERRQ(ierr);
    ierr = VecDuplicate(x, &r);
    CHKERRQ(ierr);
    ierr = SNESSetFunction(snes, r, ObjectiveFunctionCMBdSdxiFromFlux, E);
    CHKERRQ(ierr);

    /* initial guess for cmb dS/dxi */
    ierr = VecSetValue(x, ind0, -1.0E-1, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(x);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    CHKERRQ(ierr);

    /* Inform the nonlinear solver to generate a finite-difference approximation
       to the Jacobian */
    ierr = PetscOptionsSetValue(NULL, "-cmbic_snes_mf", NULL);
    CHKERRQ(ierr);
    /* Turn off convergence based on step size and trust region tolerance */
    ierr = PetscOptionsSetValue(NULL, "-cmbic_snes_stol", "0");
    CHKERRQ(ierr);
    /* atol to solve within 1 W/m^2 */
    ierr = PetscOptionsSetValue(NULL, "-cmbic_snes_rtol", "0");
    CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL, "-cmbic_snes_atol", "1.0E0");
    CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL, "-cmbic_ksp_rtol", "1.0e-6");
    CHKERRQ(ierr);
    ierr = PetscOptionsSetValue(NULL, "-cmbic_ksp_atol", "1.0e-6");
    CHKERRQ(ierr);

    /* For solver analysis/debugging/tuning, activate a custom monitor with a flag */
    {
        PetscBool flg = PETSC_FALSE;

        ierr = PetscOptionsGetBool(NULL, NULL, "-cmbic_snes_verbose_monitor", &flg, NULL);
        CHKERRQ(ierr);
        if (flg)
        {
            ierr = SNESMonitorSet(snes, SNESMonitorVerbose, NULL, NULL);
            CHKERRQ(ierr);
        }
    }

    /* Solve */
    ierr = SNESSetFromOptions(snes);
    CHKERRQ(ierr); /* Picks up any additional options (note prefix) */
    ierr = SNESSolve(snes, NULL, x);
    CHKERRQ(ierr);
    {
        SNESConvergedReason reason;
        ierr = SNESGetConvergedReason(snes, &reason);
        CHKERRQ(ierr);
        if (reason < 0)
            SETERRQ1(PetscObjectComm((PetscObject)snes), PETSC_ERR_CONV_FAILED,
                     "Nonlinear solver didn't converge: %s\n", SNESConvergedReasons[reason]);
    }

    /* store solution */
    ierr = VecGetValues(x, 1, &ind0, &dSdxi);
    ierr = VecSetValue(S->dSdxi, ind_cmb, dSdxi, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->dSdxi);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->dSdxi);
    CHKERRQ(ierr);

    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&r);
    CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode ObjectiveFunctionCMBdSdxiFromFlux(SNES snes, Vec x, Vec f, void *ptr)
{
    PetscErrorCode ierr;
    PetscScalar res, dSdxi, S_cmb, S_abv, xi_cmb, xi_abv, Jcond, Jconv, Jmix, Jgrav, target;
    PetscInt numpts_b, ind0, ind_cmb, ind_abv;
    Ctx *E = (Ctx *)ptr;
    Parameters const P = E->parameters;
    ScalingConstants const SC = P->scaling_constants;
    Mesh const *M = &E->mesh;
    Solution *S = &E->solution;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_b, NULL, &numpts_b, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    CHKERRQ(ierr);
    ind_cmb = numpts_b - 1; // index of last basic node (i.e., cmb)
    ind_abv = ind_cmb - 1;
    ind0 = 0;

    /* gradient we are solving for */
    ierr = VecGetValues(x, 1, &ind0, &dSdxi);
    CHKERRQ(ierr);
    ierr = VecGetValues(S->S_s, 1, &ind_abv, &S_abv);
    CHKERRQ(ierr);
    ierr = VecGetValues(M->xi_b, 1, &ind_cmb, &xi_cmb);
    CHKERRQ(ierr);
    ierr = VecGetValues(M->xi_b, 1, &ind_abv, &xi_abv);
    CHKERRQ(ierr);

    /* CMB entropy using standard reconstruction */
    S_cmb = S_abv + dSdxi * 0.5 * (xi_cmb - xi_abv);

    /* write dSdxi and S at CMB to Vecs in S */
    ierr = VecSetValue(S->S, ind_cmb, S_cmb, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = VecSetValue(S->dSdxi, ind_cmb, dSdxi, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->S);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->S);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->dSdxi);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->dSdxi);
    CHKERRQ(ierr);

    /* update material properties at CMB */
    ierr = set_matprop_basic(E);
    CHKERRQ(ierr);

    /* compute energy flux */
    Jcond = GetConductiveHeatFlux(E, &ind_cmb);
    Jconv = GetConvectiveHeatFlux(E, &ind_cmb);
    Jmix = GetMixingHeatFlux(E, &ind_cmb);
    Jgrav = GetGravitationalHeatFlux(E, &ind_cmb);

    /* value we want to recover */
    target = P->core_bc_value;

    /* residual is difference of CMB flux to target in W/m^2 */
    /* scaled to W/m^2 to guide solver tolerance selection */
    res = Jcond * SC->FLUX + Jconv * SC->FLUX + Jmix * SC->FLUX + Jgrav * SC->FLUX - target * SC->FLUX;

    ierr = VecSetValue(f, ind0, res, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(f);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(f);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode SetCoreMantleFluxBC(Ctx *E)
{
    PetscErrorCode ierr;
    PetscInt ind_cmb, numpts_b;
    PetscMPIInt rank, size;
    Parameters const P = E->parameters;
    Solution *S = &E->solution;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_b, NULL, &numpts_b, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    CHKERRQ(ierr);
    ind_cmb = numpts_b - 1; // index of last basic node (i.e., cmb)

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
    CHKERRQ(ierr);

    /* conform heat flux at the core mantle boundary to the
       imposed boundary condition */

    /* assume that the last rank contains the last basic node */
    if (rank == size - 1)
    {

        /* conform basal mantle flux to core mantle boundary condition */
        switch (P->CORE_BC)
        {
        case 1:
            // core cooling
            /* do nothing since core cools by mantle heat flux
               as determined above */
            break;
        case 2:
            // constant heat flux
            ierr = VecSetValue(S->Jtot, ind_cmb, P->core_bc_value, INSERT_VALUES);
            CHKERRQ(ierr);
            break;
        case 3:
            // isothermal core, i.e. CMB entropy/temperature fixed
            ierr = VecSetValue(S->Jtot, ind_cmb, 0.0, INSERT_VALUES);
            CHKERRQ(ierr);
            break;
        default:
            SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unsupported CORE_BC value %d provided", P->CORE_BC);
        }

        ierr = VecAssemblyBegin(S->Jtot);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(S->Jtot);
        CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
