#include <petscdmcomposite.h>
#include "dimensionalisablefield.h"
#include "ctx.h"

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

  ierr = PetscMalloc1(f->numDomains,&f->scaling);CHKERRQ(ierr);
  for (i=0; i<f->numDomains; ++i) {
    f->scaling[i] = scalings[i];
  }

  f->scaled = isScaled;

  ierr = DMCreateGlobalVector(f->dm,&f->vecGlobal);CHKERRQ(ierr);

  f->name = "Unnamed DimensionalisableField";

  ierr = PetscMalloc1(f->numDomains,&f->subdomainNames);CHKERRQ(ierr);
  ierr = PetscMalloc1(f->numDomains,&f->subdomainUnits);CHKERRQ(ierr);
  for (i=0; i<f->numDomains; ++i) {
    f->subdomainNames[i] = "Unnamed Subdomain";
    f->subdomainUnits[i] = "Unknown Units";
  }

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
  ierr = PetscFree(f->subdomainNames);CHKERRQ(ierr);
  ierr = PetscFree(f->subdomainUnits);CHKERRQ(ierr);
  ierr = PetscFree(*pf);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DimensionalisableFieldDuplicate(DimensionalisableField f,DimensionalisableField *pfNew)
{
  PetscErrorCode         ierr;
  DimensionalisableField fNew;
  PetscInt               i;

  PetscFunctionBeginUser;
  ierr = DimensionalisableFieldCreate(pfNew,f->dm,f->scaling,f->scaled);CHKERRQ(ierr);
  fNew = *pfNew;
  ierr = DimensionalisableFieldSetName(fNew,f->name);CHKERRQ(ierr);
  for (i=0; i<fNew->numDomains; ++i) {
    ierr = DimensionalisableFieldSetSubdomainName(fNew,i,f->subdomainNames[i]);CHKERRQ(ierr);
    ierr = DimensionalisableFieldSetSubdomainUnits(fNew,i,f->subdomainUnits[i]);CHKERRQ(ierr);
  }
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

/* Helper function to convert a scalar value to an appropriate string  */
PetscErrorCode valueToString(PetscScalar val,char *str,int maxLength)
{
  PetscErrorCode ierr;
#if PETSC_USE_REAL___FLOAT128
  ierr = quadmath_snprintf(str,maxLength,"%32.32Qg",val);
  if (ierr >= maxLength || ierr < 0) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_LIB,"quadmath_snprintf() failed");
#else
  ierr = PetscSNPrintf(str,maxLength,"%16.16g",val);CHKERRQ(ierr);
#endif
  return 0;
}

#define FORMAT_STRING_SIZE 64
PetscErrorCode DimensionalisableFieldToJSON(DimensionalisableField const f,cJSON **pjson)
{
  PetscErrorCode    ierr;
  Vec               vec;
  PetscInt          vecSize,i,d;
  cJSON             *json,*values,*subdomainArray;
  const PetscScalar *arr;
  Vec               *subVecs;

  PetscFunctionBeginUser;

  ierr = DimensionalisableFieldGetGlobalVec(f,&vec);CHKERRQ(ierr);

  *pjson = cJSON_CreateObject();
  json = *pjson;
  cJSON_AddItemToObject(json,"name",cJSON_CreateString(f->name));
  cJSON_AddItemToObject(json,"subdomains",cJSON_CreateNumber(f->numDomains));
  if (f->numDomains > 1) {
    subdomainArray = cJSON_CreateArray();
    ierr = PetscMalloc1(f->numDomains,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(f->dm,vec,f->numDomains,NULL,subVecs);CHKERRQ(ierr);
  } else {
    if (f->subdomainUnits[0]) {
      cJSON_AddItemToObject(json,"units",cJSON_CreateString(f->subdomainUnits[0]));
    }
    subdomainArray = NULL;
  }
  for (d=0; d<f->numDomains; ++d){
    cJSON *curr;
    Vec   vecCurr;
    if (f->numDomains > 1) {
      curr = cJSON_CreateObject();
      vecCurr = subVecs[d];
      cJSON_AddItemToObject(curr,"subdomain",cJSON_CreateNumber(d));
      if (f->subdomainNames[d]) {
        cJSON_AddItemToObject(curr,"description",cJSON_CreateString(f->subdomainNames[d]));
      }
      if (f->subdomainUnits[d]) {
        cJSON_AddItemToObject(curr,"units",cJSON_CreateString(f->subdomainUnits[d]));
      }
    } else {
      curr = json;
      vecCurr = vec;
    }
    ierr = VecGetSize(vecCurr,&vecSize);CHKERRQ(ierr);
    cJSON_AddItemToObject(curr,"size",cJSON_CreateNumber(vecSize));
    {
      char str[FORMAT_STRING_SIZE];
      ierr = valueToString(f->scaling[d],str,FORMAT_STRING_SIZE);CHKERRQ(ierr);
      cJSON_AddItemToObject(curr,"scaling",cJSON_CreateString(str));
    }
    cJSON_AddItemToObject(curr,"scaled",f->scaled? cJSON_CreateString("true") : cJSON_CreateString("false"));

    /* Store vector values as STRINGS to ensure precision we want */
    /* Current implementation assumes single rank */
    {
      PetscMPIInt size;
      ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
      if (size != 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Not implemented in parallel");
    }
    ierr = VecGetArrayRead(vecCurr,&arr);CHKERRQ(ierr);
    values = cJSON_CreateArray();
    {
      for (i=0; i<vecSize; ++i) {
        char str[FORMAT_STRING_SIZE];
        ierr = valueToString(arr[i],str,FORMAT_STRING_SIZE);CHKERRQ(ierr);
        cJSON_AddItemToArray(values,cJSON_CreateString(str));
      }
    }
    ierr = VecRestoreArrayRead(vecCurr,&arr);CHKERRQ(ierr);

    cJSON_AddItemToObject(curr,"values",values);
    if (f->numDomains > 1) {
      cJSON_AddItemToArray(subdomainArray,curr);
    }
  }
  if (f->numDomains > 1 ) {
    cJSON_AddItemToObject(json,"subdomain data",subdomainArray);
    ierr = DMCompositeRestoreAccessArray(f->dm,vec,f->numDomains,NULL,subVecs);CHKERRQ(ierr);
    ierr = PetscFree(subVecs);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
#undef FORMAT_STRING_SIZE

PetscErrorCode DimensionalisableFieldSetName(DimensionalisableField f,const char *name)
{
  PetscFunctionBeginUser;
  f->name = name;
  PetscFunctionReturn(0);
}

PetscErrorCode DimensionalisableFieldSetUnits(DimensionalisableField f,const char *units)
{
  PetscFunctionBeginUser;
  if (f->numDomains != 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Not supported for multiple-subdomain DimensionalisableField objects. Use DimensionalisableFieldSetSubdomainUnits");
  f->subdomainUnits[0] = units;
  PetscFunctionReturn(0);
}

PetscErrorCode DimensionalisableFieldSetSubdomainName(DimensionalisableField f,PetscInt subdomain,const char *name)
{
  PetscFunctionBeginUser;
  f->subdomainNames[subdomain] = name;
  PetscFunctionReturn(0);
}

PetscErrorCode DimensionalisableFieldSetSubdomainUnits(DimensionalisableField f,PetscInt subdomain,const char *units)
{
  PetscFunctionBeginUser;
  f->subdomainUnits[subdomain] = units;
  PetscFunctionReturn(0);
}

PetscErrorCode JSON_add_single_value_to_object( DM dm, PetscScalar scaling, const char *name, const char *units, const PetscScalar value, cJSON *data )
{
  PetscErrorCode         ierr;
  cJSON                  *item;
  DimensionalisableField dfield;

  PetscFunctionBeginUser;

  ierr = DimensionalisableFieldCreate(&dfield,dm,&scaling,PETSC_FALSE);CHKERRQ(ierr);
  ierr = DimensionalisableFieldSetName(dfield,name);CHKERRQ(ierr);
  ierr = DimensionalisableFieldSetUnits(dfield,units);CHKERRQ(ierr);
  ierr = VecSetValue(dfield->vecGlobal,0,value,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(dfield->vecGlobal);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(dfield->vecGlobal);CHKERRQ(ierr);
  ierr = DimensionalisableFieldToJSON(dfield,&item);CHKERRQ(ierr);
  cJSON_AddItemToObject(data,name,item);
  ierr = DimensionalisableFieldDestroy(&dfield);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
