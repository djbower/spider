#include "eos_output.h"

PetscErrorCode JSON_add_phase_boundary(const Ctx *E, const EOS Ep, const char *name, cJSON *json)
{
    /* add a phase boundary evaluated in temperature space at the basic nodes points */

    PetscErrorCode ierr;
    ScalingConstants SC = E->parameters->scaling_constants;
    cJSON *data;
    PetscInt i, ilo_b, ihi_b, w_b;
    DM da_b = E->da_b;
    Vec pres_b, phase_b;
    DimensionalisableField DF_phase_b;
    PetscScalar Tbound, scaling, *arr_phase_b;
    const PetscScalar *arr_pres_b;
    char namestr1[80];

    PetscFunctionBeginUser;

    pres_b = E->mesh.pressure_b;
    ierr = DMDAGetCorners(da_b, &ilo_b, 0, 0, &w_b, 0, 0);
    CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    /* not sure on portability of strcat?  Petsc equivalent? */

    /* phase boundary defined by temperature */
    strcat(strcpy(namestr1, name), "_b");
    scaling = SC->TEMP;
    ierr = DimensionalisableFieldCreate(&DF_phase_b, da_b, &scaling, PETSC_FALSE);
    CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(DF_phase_b, &phase_b);
    CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetName(DF_phase_b, namestr1);
    CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(DF_phase_b, "K");
    CHKERRQ(ierr);

    /* basic nodes */
    ierr = DMDAVecGetArrayRead(da_b, pres_b, &arr_pres_b);
    CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b, phase_b, &arr_phase_b);
    CHKERRQ(ierr);

    for (i = ilo_b; i < ihi_b; ++i)
    {
        ierr = EOSGetPhaseBoundary(Ep, arr_pres_b[i], &Tbound, NULL);
        CHKERRQ(ierr);
        arr_phase_b[i] = Tbound;
    }

    ierr = DMDAVecRestoreArrayRead(da_b, pres_b, &arr_pres_b);
    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b, phase_b, &arr_phase_b);
    CHKERRQ(ierr);

    ierr = DimensionalisableFieldToJSON(DF_phase_b, &data);
    CHKERRQ(ierr);
    cJSON_AddItemToObject(json, namestr1, data);

    ierr = DimensionalisableFieldDestroy(&DF_phase_b);
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
