#include "eos_output.h"
#include "eos_composite.h" // TODO need this?

PetscErrorCode JSON_add_phase_boundary( const Ctx *E, const EosParameters Ep, const char *name, cJSON *json )
{
    /* add a phase boundary evaluated in entropy and temperature space
       at the basic nodes points */

    PetscErrorCode         ierr;
    ScalingConstants       SC = E->parameters->scaling_constants;
    cJSON                  *data;
    PetscInt               i,ilo_b,ihi_b,w_b;
    DM                     da_b=E->da_b;
    Vec                    pres_b, phase_b, phase_temp_b;
    DimensionalisableField DF_phase_b, DF_phase_temp_b;
    PetscScalar            Sbound, scaling, *arr_phase_b, *arr_phase_temp_b;
    const PetscScalar      *arr_pres_b;
    EosEval                eos_eval;
    char                   namestr1[80], namestr2[80];

    PetscFunctionBeginUser;

    pres_b = E->mesh.pressure_b;
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    /* TODO: not sure on portability of strcat?  Petsc equivalent? */

    /* phase boundary defined by entropy */
    strcat(strcpy(namestr1,name),"_b"); // phase boundary defined by entropy
    scaling = SC->ENTROPY;
    ierr = DimensionalisableFieldCreate(&DF_phase_b,da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(DF_phase_b,&phase_b);CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetName(DF_phase_b,namestr1);CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(DF_phase_b,"J kg$^{-1}$ K$^{-1}$");CHKERRQ(ierr);

    /* phase boundary defined by temperature */
    strcat(strcpy(namestr2,name),"_temp_b"); // phase boundary defined by temperature
    scaling = SC->TEMP;
    ierr = DimensionalisableFieldCreate(&DF_phase_temp_b,da_b,&scaling,PETSC_FALSE);CHKERRQ(ierr);
    ierr = DimensionalisableFieldGetGlobalVec(DF_phase_temp_b,&phase_temp_b);CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetName(DF_phase_temp_b,namestr2);CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetUnits(DF_phase_temp_b,"K");CHKERRQ(ierr);

    /* basic nodes */
    ierr = DMDAVecGetArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,phase_b,&arr_phase_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,phase_temp_b,&arr_phase_temp_b);CHKERRQ(ierr);

    for(i=ilo_b;i<ihi_b;++i){
        ierr = SetPhaseBoundary( Ep, arr_pres_b[i], &Sbound, NULL );CHKERRQ(ierr);
        arr_phase_b[i] = Sbound;
        ierr = SetEosEval( Ep, arr_pres_b[i], Sbound, &eos_eval );CHKERRQ(ierr);
        arr_phase_temp_b[i] = eos_eval.T;
    }

    ierr = DMDAVecRestoreArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,phase_b,&arr_phase_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,phase_temp_b,&arr_phase_temp_b);CHKERRQ(ierr);

    ierr = DimensionalisableFieldToJSON(DF_phase_b,&data);CHKERRQ(ierr);
    cJSON_AddItemToObject(json,namestr1,data);

    ierr = DimensionalisableFieldToJSON(DF_phase_temp_b,&data);CHKERRQ(ierr);
    cJSON_AddItemToObject(json,namestr2,data);

    ierr = DimensionalisableFieldDestroy(&DF_phase_b);CHKERRQ(ierr);
    ierr = DimensionalisableFieldDestroy(&DF_phase_temp_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
