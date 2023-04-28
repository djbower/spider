#include "bc.h"
#include "monitor.h"
#include "util.h"

/* Helper routine to prepend the root directory to a relative path */
/* https://gcc.gnu.org/onlinedocs/gcc-4.9.0/cpp/Stringification.html */
#define STRINGIFY(x) STRINGIFY2(x)
#define STRINGIFY2(x) #x
#define SPIDER_ROOT_DIR_STR STRINGIFY(SPIDER_ROOT_DIR)
PetscErrorCode MakeRelativeToSourcePathAbsolute(char *path)
{
    PetscErrorCode ierr;
    char tmp[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    ierr = PetscStrcpy(tmp, path);
    CHKERRQ(ierr);
    ierr = PetscStrcpy(path, SPIDER_ROOT_DIR_STR);
    CHKERRQ(ierr);
    ierr = PetscStrcat(path, "/");
    CHKERRQ(ierr); /* not portable */
    ierr = PetscStrcat(path, tmp);
    CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
#undef SPIDER_ROOT_DIR_STR

PetscErrorCode set_solution_from_temperature(Ctx *E, Vec sol)
{
    /* Clone data in S->dTdxi and top staggered node value (from S->T_s)
       to Vec sol */

    PetscErrorCode ierr;
    Solution *S = &E->solution;
    PetscScalar T0;
    PetscInt ind0 = 0;
    Vec *subVecs;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(E->numFields, &subVecs);
    CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol, sol, E->numFields, NULL, subVecs);
    CHKERRQ(ierr);
    /* set dT/dxi at basic nodes to solution Vec */
    ierr = VecCopy(S->dTdxi, subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DTDXI_B]]);
    CHKERRQ(ierr);
    /* set T0 to solution Vec.  T0 is the top staggered node temperature. */
    ierr = VecGetValues(S->T_s, 1, &ind0, &T0);
    CHKERRQ(ierr);
    ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_T0]], 0, T0, INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_T0]]);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_T0]]);
    CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccessArray(E->dm_sol, sol, E->numFields, NULL, subVecs);
    CHKERRQ(ierr);
    ierr = PetscFree(subVecs);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_temperature_from_solution(Ctx *E, Vec sol)
{
    /* Set temperature-related Vecs in E to be consistent with the Vec sol, and additionally
        compute values for the outer boundaries of the basic mesh */

    PetscErrorCode ierr;
    Solution *S = &E->solution;
    PetscScalar T0;
    PetscMPIInt size;
    Vec *subVecs;

    const PetscInt ind0 = 0;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
    CHKERRQ(ierr);
    if (size > 1)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "This code has only been correctly implemented for serial runs");

    ierr = PetscMalloc1(E->numFields, &subVecs);
    CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol, sol, E->numFields, NULL, subVecs);
    CHKERRQ(ierr);
    /* set S->dT/dxi at basic nodes */
    ierr = VecCopy(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DTDXI_B]], S->dTdxi);
    CHKERRQ(ierr);
    /* get first staggered node value (store as T0) */
    ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_T0]], 1, &ind0, &T0);
    CHKERRQ(ierr);

    ierr = set_temperature_reconstruction_from_ctx(E, T0);
    CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccessArray(E->dm_sol, sol, E->numFields, NULL, subVecs);
    CHKERRQ(ierr);
    PetscFree(subVecs);

    PetscFunctionReturn(0);
}

PetscErrorCode set_temperature_reconstruction_from_ctx(Ctx *E, PetscScalar T0)
{
    PetscErrorCode ierr;
    Mesh *M = &E->mesh;
    Solution *S = &E->solution;
    PetscScalar *arr_T_b, *arr_T_s, *arr_dTdxi_b, *arr_xi_s, *arr_xi_b;
    PetscInt i, ihi_b, ilo_b, w_b;
    DM da_s = E->da_s, da_b = E->da_b;

    PetscFunctionBeginUser;

    /* for looping over basic nodes */
    ierr = DMDAGetCorners(da_b, &ilo_b, 0, 0, &w_b, 0, 0);
    CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    ierr = DMDAVecGetArray(da_b, S->T, &arr_T_b);
    CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s, S->T_s, &arr_T_s);
    CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b, S->dTdxi, &arr_dTdxi_b);
    CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b, M->xi_b, &arr_xi_b);
    CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s, M->xi_s, &arr_xi_s);
    CHKERRQ(ierr);

    /* T (absolute) at all staggered and basic internal nodes */
    arr_T_s[0] = 0.0;
    /* note plus one for start of loop */
    for (i = ilo_b + 1; i < ihi_b - 1; ++i)
    {
        /* T (absolute) at staggered nodes */
        arr_T_s[i] = arr_dTdxi_b[i] * (arr_xi_s[i] - arr_xi_s[i - 1]);
        arr_T_s[i] += arr_T_s[i - 1]; // dT relative to first staggered value
        arr_T_b[i] = arr_dTdxi_b[i] * 0.5 * (arr_xi_b[i] - arr_xi_b[i - 1]);
        arr_T_b[i] += arr_T_s[i - 1];
        arr_T_s[i - 1] += T0; // add at end to try and retain precision
        arr_T_b[i - 1] += T0; // add at end to try and retain precision
    }
    /* loop above terminates before we have added the constant offset
       to the last staggered node value.  So do that here */
    arr_T_s[ihi_b - 2] += T0; // add at end to try and retain precision
    arr_T_b[ihi_b - 2] += T0; // add at end to try and retain precision

    ierr = DMDAVecRestoreArray(da_b, S->T, &arr_T_b);
    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s, S->T_s, &arr_T_s);
    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, S->dTdxi, &arr_dTdxi_b);
    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b, M->xi_b, &arr_xi_b);
    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s, M->xi_s, &arr_xi_s);
    CHKERRQ(ierr);

    /* core-mantle boundary */
    ierr = set_cmb_temperature_from_cmb_gradient(E);
    CHKERRQ(ierr);

    /* surface boundary */
    ierr = set_surface_temperature_from_surface_gradient(E);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_partial_pressures_from_solution(Ctx *E, Vec sol)
{
    /* clone the current atmosphere solution (volatile partial
       pressures and reaction masses) into the atmos struct */

    PetscErrorCode ierr;
    Atmosphere *A = &E->atmosphere;
    Parameters const P = E->parameters;
    AtmosphereParameters const Ap = P->atmosphere_parameters;
    PetscInt i;
    PetscMPIInt size;
    Vec *subVecs;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);
    CHKERRQ(ierr);
    if (size > 1)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, "This code has only been correctly implemented for serial runs");

    ierr = PetscMalloc1(E->numFields, &subVecs);
    CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol, sol, E->numFields, NULL, subVecs);
    CHKERRQ(ierr);

    /* partial pressures */
    for (i = 0; i < Ap->n_volatiles; ++i)
    {
        ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]], 1, &i, &A->volatiles[i].p);
        CHKERRQ(ierr);
        /* pseudo-volatiles track log10(P(Pa)), so convert to actual non-dimensional pressure */
        if (Ap->PSEUDO_VOLATILES)
        {
            A->volatiles[i].p = PetscPowScalar(10.0, A->volatiles[i].p);
        }
    }

    /* mass reaction terms */
    for (i = 0; i < Ap->n_reactions; ++i)
    {
        ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]], 1, &i, &A->mass_reaction[i]);
        CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(E->dm_sol, sol, E->numFields, NULL, subVecs);
    CHKERRQ(ierr);
    PetscFree(subVecs);

    PetscFunctionReturn(0);
}

PetscErrorCode set_solution_from_partial_pressures(Ctx *E, Vec sol)
{
    /* clone the values in the atmos struct to the solution Vec */

    PetscErrorCode ierr;
    Parameters P = E->parameters;
    Atmosphere *A = &E->atmosphere;
    AtmosphereParameters const Ap = P->atmosphere_parameters;
    PetscInt i;
    Vec *subVecs;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(E->numFields, &subVecs);
    CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol, sol, E->numFields, NULL, subVecs);
    CHKERRQ(ierr);

    /* partial pressures */
    for (i = 0; i < Ap->n_volatiles; ++i)
    {
        if (Ap->PSEUDO_VOLATILES)
        {
            A->volatiles[i].p = PetscLog10Real(A->volatiles[i].p);
        }
        ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]], i, A->volatiles[i].p, INSERT_VALUES);
        CHKERRQ(ierr);
    }

    /* mass reaction terms */
    for (i = 0; i < Ap->n_reactions; ++i)
    {
        ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]], i, A->mass_reaction[i], INSERT_VALUES);
        CHKERRQ(ierr);
    }

    if (Ap->n_volatiles > 0)
    {
        ierr = VecAssemblyBegin(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]]);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]]);
        CHKERRQ(ierr);
    }
    if (Ap->n_reactions > 0)
    {
        ierr = VecAssemblyBegin(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]]);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]]);
        CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(E->dm_sol, sol, E->numFields, NULL, subVecs);
    CHKERRQ(ierr);
    PetscFree(subVecs);

    PetscFunctionReturn(0);
}

PetscScalar combine_matprop(PetscScalar weight, PetscScalar mat1, PetscScalar mat2)
{
    /* linear weighting of two quantities */

    PetscScalar out;

    out = weight * mat1 + (1.0 - weight) * mat2;

    return out;
}

PetscScalar tanh_weight(PetscScalar qty, PetscScalar threshold, PetscScalar width)
{
    /* tanh weight for viscosity profile and smoothing */

    PetscScalar fwt, z;

    z = (qty - threshold) / width;
    fwt = 0.5 * (1.0 + PetscTanhScalar(z));
    return fwt;
}

PetscScalar get_smoothing(PetscScalar smooth_width, PetscScalar gphi)
{
    /* get smoothing across phase boundaries for a two phase composite */

    PetscScalar smth;

    /* no smoothing */
    if (smooth_width == 0.0)
    {
        smth = 1.0; // mixed phase only
        if ((gphi < 0.0) || (gphi > 1.0))
        {
            smth = 0.0; // single phase only
        }
    }

    /* tanh smoothing */
    else
    {
        if (gphi > 0.5)
        {
            smth = 1.0 - tanh_weight(gphi, 1.0, smooth_width);
        }
        else
        {
            smth = tanh_weight(gphi, 0.0, smooth_width);
        }
    }

    return smth;
}

PetscErrorCode Make2DPetscScalarArray(PetscInt arraySizeX, PetscInt arraySizeY, PetscScalar ***theArray_ptr)
{
    PetscErrorCode ierr;
    PetscInt i;
    PetscScalar **theArray;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(arraySizeX, theArray_ptr);
    CHKERRQ(ierr);
    theArray = *theArray_ptr;
    for (i = 0; i < arraySizeX; i++)
    {
        ierr = PetscMalloc1(arraySizeY, &theArray[i]);
        CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode make_vec_mask(DM dm, PetscInt index, Vec *mask_ptr)
{
    PetscErrorCode ierr;
    PetscInt numpts, i, ilo, ihi, w;
    const PetscScalar one = 1.0, zero = 0.0;

    ierr = DMDAGetInfo(dm, NULL, &numpts, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    CHKERRQ(ierr);

    /* create mask vector */
    ierr = VecCreate(PETSC_COMM_WORLD, mask_ptr);
    CHKERRQ(ierr);
    ierr = VecSetSizes(*mask_ptr, PETSC_DECIDE, numpts);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(*mask_ptr);
    CHKERRQ(ierr);
    ierr = VecSetUp(*mask_ptr);
    CHKERRQ(ierr);

    ierr = VecSet(*mask_ptr, 0.0);
    CHKERRQ(ierr);

    ierr = DMDAGetCorners(dm, &ilo, 0, 0, &w, 0, 0);
    CHKERRQ(ierr);
    ihi = ilo + w;

    for (i = ilo; i < ihi; ++i)
    {
        if (i < index)
        {
            ierr = VecSetValues(*mask_ptr, 1, &i, &one, INSERT_VALUES);
            CHKERRQ(ierr);
        }
        else
        {
            ierr = VecSetValues(*mask_ptr, 1, &i, &zero, INSERT_VALUES);
            CHKERRQ(ierr);
            break;
        }
    }

    ierr = VecAssemblyBegin(*mask_ptr);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*mask_ptr);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode average_by_mass_staggered(Ctx *E, Vec in_vec, Vec *in_mask_ptr, PetscScalar *out)
{
    PetscErrorCode ierr;
    Mesh *M = &E->mesh;
    Vec qty_s, mass_s;
    PetscScalar mass;
    PetscInt numpts_s;
    Vec in_mask = *in_mask_ptr;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_s, NULL, &numpts_s, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    CHKERRQ(ierr);

    /* create work vectors */
    ierr = VecCreate(PETSC_COMM_WORLD, &qty_s);
    CHKERRQ(ierr);
    ierr = VecSetSizes(qty_s, PETSC_DECIDE, numpts_s);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(qty_s);
    CHKERRQ(ierr);
    ierr = VecSetUp(qty_s);
    CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD, &mass_s);
    CHKERRQ(ierr);
    ierr = VecSetSizes(mass_s, PETSC_DECIDE, numpts_s);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(mass_s);
    CHKERRQ(ierr);
    ierr = VecSetUp(mass_s);
    CHKERRQ(ierr);

    /* compute mass-weighted average */
    ierr = VecPointwiseMult(qty_s, in_vec, in_mask);
    CHKERRQ(ierr);
    ierr = VecPointwiseMult(mass_s, M->mass_s, in_mask);
    CHKERRQ(ierr);
    ierr = VecPointwiseMult(qty_s, qty_s, mass_s);
    CHKERRQ(ierr);
    ierr = VecSum(qty_s, out);
    CHKERRQ(ierr);
    ierr = VecSum(mass_s, &mass);
    CHKERRQ(ierr);

    *out = *out / mass;

    /* this block would prevent 'nans' from being returned, but 'nans'
       are helpful because it indicates the rheological front is not
       yet moving */
    // if (mass > 0.0){
    //     *out = *out / mass;
    // }
    // else{
    //     *out = 0.0;
    // }

    ierr = VecDestroy(&qty_s);
    CHKERRQ(ierr);
    ierr = VecDestroy(&mass_s);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode invert_vec_mask(Vec *in_vec_ptr)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecShift(*in_vec_ptr, -1);
    CHKERRQ(ierr);
    ierr = VecScale(*in_vec_ptr, -1);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* helper functions for parsing parameters */

PetscErrorCode PetscScalarCheckPositive(PetscScalar value, const char *value_string)
{
    PetscFunctionBeginUser;
    if (value < 0.0)
    {
        SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "%s must be positive (currently %f)", value_string, value);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode PetscIntCheckPositive(PetscInt value, const char *value_string)
{
    PetscFunctionBeginUser;
    if (value < 0)
    {
        SETERRQ2(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "%s must be positive (currently %d)", value_string, value);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode PetscOptionsGetPositiveScalar(const char *value_string, PetscScalar *value_ptr, PetscScalar value_default, PetscBool *set)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    *value_ptr = value_default;
    ierr = PetscOptionsGetScalar(NULL, NULL, value_string, value_ptr, set);
    CHKERRQ(ierr);
    if (set != NULL && *set)
    {
        ierr = PetscScalarCheckPositive(*value_ptr, value_string);
        CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode PetscOptionsGetPositiveInt(const char *value_string, PetscInt *value_ptr, PetscInt value_default, PetscBool *set)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    *value_ptr = value_default;
    ierr = PetscOptionsGetInt(NULL, NULL, value_string, value_ptr, set);
    CHKERRQ(ierr);
    if (set != NULL && *set)
    {
        ierr = PetscIntCheckPositive(*value_ptr, value_string);
        CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
