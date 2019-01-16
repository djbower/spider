#if !defined(SCALABLEFIELD_H_)
#define SCALABLEFIELD_H_

#include <petscdm.h>
#include "cJSON.h"

struct _p_DimensionalisableField {
  DM           dm;
  Vec          vecGlobal;
  PetscInt     numDomains;
  PetscScalar  *scaling;
  PetscBool    scaled;
  const char   *name;
  const char   *units;
  const char   **slotNames;
  const char   **slotUnits;
};
typedef struct _p_DimensionalisableField *DimensionalisableField;

PetscErrorCode DimensionalisableFieldCreate(DimensionalisableField*,DM,PetscScalar*,PetscBool);
PetscErrorCode DimensionalisableFieldCreateLocalVec(DimensionalisableField,Vec*);
PetscErrorCode DimensionalisableFieldDestroy(DimensionalisableField*);
PetscErrorCode DimensionalisableFieldDuplicate(DimensionalisableField,DimensionalisableField*);
PetscErrorCode DimensionalisableFieldGetGlobalVec(DimensionalisableField,Vec*);
PetscErrorCode DimensionalisableFieldGetScaling(DimensionalisableField,PetscInt*,PetscScalar*);
PetscErrorCode DimensionalisableFieldScale(DimensionalisableField);
PetscErrorCode DimensionalisableFieldUnscale(DimensionalisableField);
PetscErrorCode DimensionalisableFieldToJSON(DimensionalisableField const,cJSON**);
PetscErrorCode DimensionalisableFieldSetName(DimensionalisableField,const char*);
PetscErrorCode DimensionalisableFieldSetUnits(DimensionalisableField,const char*);
PetscErrorCode DimensionalisableFieldSetSlotName(DimensionalisableField,PetscInt,const char*);
PetscErrorCode DimensionalisableFieldSetSlotUnits(DimensionalisableField,PetscInt,const char*);
// REMOVE ?PetscErrorCode AddSingleValueToJSONArray(DM,PetscScalar,const char *, const char *, PetscScalar const, cJSON *);
PetscErrorCode JSON_add_single_value_to_object(DM,PetscScalar,const char *, const char *, PetscScalar const, cJSON *);

#endif
