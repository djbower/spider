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
  const char   **subdomainNames;
  const char   **subdomainUnits;
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
PetscErrorCode DimensionalisableFieldSetSubdomainName(DimensionalisableField,PetscInt,const char*);
PetscErrorCode DimensionalisableFieldSetSubdomainUnits(DimensionalisableField,PetscInt,const char*);

PetscErrorCode JSON_add_single_value_to_object(DM,PetscScalar,const char *, const char *, PetscScalar const, cJSON *);

#endif
