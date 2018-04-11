#include <petscdmcomposite.h>
#include "dimensionalisablefield.h"

PetscErrorCode DimensionalisableFieldCreate(DimensionalisableField *pf, DM dm, PetscScalar *scalings, PetscBool isScaled)
{
  PetscErrorCode         ierr;
  PetscInt               i;
  DMType                 dmtype;
  PetscBool              isComposite;
  DimensionalisableField f;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,pf);CHKERRQ(ierr);
  f = *pf;

  f->dm = dm;

  ierr = DMGetType(f->dm,&dmtype);CHKERRQ(ierr);
  ierr = PetscStrcmp(dmtype,DMCOMPOSITE,&isComposite);CHKERRQ(ierr);
  if (isComposite) {
    ierr = DMCompositeGetNumberDM(f->dm,&(f->numDomains));CHKERRQ(ierr);
  } else {
    f->numDomains = 1;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Created DimensionalisableField over DM with %d subdomains\n",f->numDomains);CHKERRQ(ierr);

  ierr = PetscMalloc1(f->numDomains,&f->scaling);CHKERRQ(ierr);
  for (i=0; i<f->numDomains; ++i) {
    f->scaling[i] = scalings[i];
  }

  f->scaled = isScaled; // Not great design - requires the user to remember (or check) this when setting

  ierr = DMCreateGlobalVector(f->dm,&f->vecGlobal);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DimensionalisableFieldDestroy(DimensionalisableField *pf)
{
  PetscErrorCode         ierr;
  DimensionalisableField f;

  PetscFunctionBeginUser;
  f = *pf;
  if (f->vecGlobal) {ierr = VecDestroy(&f->vecGlobal);CHKERRQ(ierr);}
  ierr = PetscFree(f->scaling);CHKERRQ(ierr);
  ierr = PetscFree(*pf);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Get the global vector. Note that this is created and destroyed along with the
  DimensionalisableField object, hence this just returns the Vec object (which
  is a pointer). Contrast to DimensionalisableFieldCreateLocalVector()

  Thus, you do NOT need to destroy this vector.
 */
PetscErrorCode DimensionalisableFieldGetGlobalVec(DimensionalisableField f,Vec *v)
{
  PetscFunctionBeginUser;
  *v = f->vecGlobal;
  PetscFunctionReturn(0);
}

PetscErrorCode DimensionalisableFieldGetScaling(DimensionalisableField f,PetscInt *numDomains, PetscScalar *scaling)
{
  PetscInt i;

  PetscFunctionBeginUser;
  *numDomains = f->numDomains;
  for (i=0; i<f->numDomains; ++i) {
    scaling[i] = f->scaling[i];
  }
  PetscFunctionReturn(0);
}

/*
   Create a local vector. You should destroy this yourself.
 */
PetscErrorCode DimensionalisableFieldCreateLocalVec(DimensionalisableField f,Vec *pv)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMCreateLocalVector(f->dm,pv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DimensionalisableFieldScaleOrUnscale(DimensionalisableField f,PetscBool forward)
{
  PetscErrorCode ierr;
  PetscInt       i;
  Vec            *subVecs;

  PetscFunctionBeginUser;
  if (f->scaled == forward) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: attempting to scale/unscale an already-scaled/unscaled field. Doing nothing.\n");CHKERRQ(ierr);
  } else {
    f->scaled = forward;
    if (f->numDomains == 1) {
      if (forward) {
        ierr = VecScale(f->vecGlobal,*f->scaling);CHKERRQ(ierr);
      } else {
        ierr = VecScale(f->vecGlobal,1.0/(*f->scaling));CHKERRQ(ierr);
      }
    } else {
      ierr = PetscMalloc1(f->numDomains,&subVecs);CHKERRQ(ierr);
      ierr = DMCompositeGetAccessArray(f->dm,f->vecGlobal,f->numDomains,NULL,subVecs);CHKERRQ(ierr);
      for (i=0; i<f->numDomains; ++i) {
        if (forward) {
          ierr = VecScale(subVecs[i],f->scaling[i]);CHKERRQ(ierr);
        } else {

          ierr = VecScale(subVecs[i],1.0/(f->scaling[i]));CHKERRQ(ierr);
        }
      }
      ierr = DMCompositeRestoreAccessArray(f->dm,f->vecGlobal,f->numDomains,NULL,subVecs);CHKERRQ(ierr);
      ierr = PetscFree(subVecs);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DimensionalisableFieldScale(DimensionalisableField f)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = DimensionalisableFieldScaleOrUnscale(f,PETSC_TRUE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DimensionalisableFieldUnscale(DimensionalisableField f)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  ierr = DimensionalisableFieldScaleOrUnscale(f,PETSC_FALSE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
