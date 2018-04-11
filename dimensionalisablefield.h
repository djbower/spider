#if !defined(SCALABLEFIELD_H_)
#define SCALABLEFIELD_H_

#include <petscdm.h>

struct _p_DimensionalisableField {
  DM          dm;
  Vec         vecGlobal;
  PetscInt    numDomains;
  PetscScalar *scaling;
  PetscBool   scaled;
};
typedef struct _p_DimensionalisableField *DimensionalisableField;

PetscErrorCode DimensionalisableFieldCreate(DimensionalisableField*,DM,PetscScalar*,PetscBool);
PetscErrorCode DimensionalisableFieldCreateLocalVec(DimensionalisableField,Vec*);
PetscErrorCode DimensionalisableFieldDestroy(DimensionalisableField*);
PetscErrorCode DimensionalisableFieldGetGlobalVec(DimensionalisableField,Vec*);
PetscErrorCode DimensionalisableFieldGetScaling(DimensionalisableField,PetscInt*,PetscScalar*);
PetscErrorCode DimensionalisableFieldScale(DimensionalisableField);
PetscErrorCode DimensionalisableFieldUnscale(DimensionalisableField);

#endif
