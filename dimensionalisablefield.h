#if !defined(SCALABLEFIELD_H_)
#define SCALABLEFIELD_H_

#include <petscdm.h>

struct _p_DimensionalisableField{
  DM          dm;
  Vec         vecGlobal;
  PetscInt    numDomains;
  PetscScalar *scaling;
  PetscBool   isDimensionalised;
};
typedef struct _p_DimensionalisableField *DimensionalisableField;

PetscErrorCode DimensionalisableFieldCreate(DimensionalisableField*,DM,PetscScalar*,PetscBool);
PetscErrorCode DimensionalisableFieldDestroy(DimensionalisableField*);
PetscErrorCode DimensionalisableFieldSetScaled(DimensionalisableField);
PetscErrorCode DimensionalisableFieldSetUnscaled(DimensionalisableField);
PetscErrorCode DimensionalisableFieldCreateLocalVec(DimensionalisableField,Vec*);
PetscErrorCode DimensionalisableFieldGetGlobalVec(DimensionalisableField,Vec*);
PetscErrorCode DimensionalisableFieldWriteToFile(DimensionalisableField,PetscBool);
PetscErrorCode DimensionalisableFieldReadFromFile(DimensionalisableField,PetscBool);

#endif
