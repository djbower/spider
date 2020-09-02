#include "eos_lookup.h"
#include "util.h"

PetscErrorCode EOSLookupPureCreate(EOSLookupPure* p) {
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, p);CHKERRQ(ierr);
  // TODO other stuff to init? (to zero?)
  PetscFunctionReturn(0);
}

// TODO refactor this stuff moved here for now, into the class impl

/* First material model option is via lookup tables */
/* lookup material properties (default) */
/* TODO: these next set of functions can be redone as function pointers in the EosParameters struct
   when a lookup is used */
static PetscErrorCode SetLookupRho( const Lookup, PetscScalar, PetscScalar, PetscScalar *);
static PetscErrorCode SetLookupAlpha( const Lookup, PetscScalar, PetscScalar, PetscScalar *);
static PetscErrorCode SetLookupCp( const Lookup, PetscScalar, PetscScalar, PetscScalar *);
static PetscErrorCode SetLookupdTdPs( const Lookup, PetscScalar, PetscScalar, PetscScalar *);
static PetscErrorCode SetLookupTemperature( const Lookup, PetscScalar, PetscScalar, PetscScalar *);


PetscErrorCode LookupCreate( Lookup *lookup_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscMalloc1(1,lookup_ptr);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode LookupDestroy( Lookup *lookup_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscFree(*lookup_ptr);CHKERRQ(ierr);
    *lookup_ptr = NULL;
    PetscFunctionReturn(0);
}


PetscErrorCode LookupFilenameSet( const char* property, const char* prefix, char* lookup_filename, PetscBool *IS_SET )
{
    PetscErrorCode ierr;
    char           buf1[1024]; /* max size */
    char           buf2[1024]; /* max size */
    PetscBool      set_rel_to_src,set;

    PetscFunctionBeginUser;

    /* Based on input options, determine which files to load.  Options ending
       with _rel_to_src indicate a path relative to the source code. In this 
       case we prepend a string, SPIDER_ROOT_DIR_STR, and /. The corresponding
       option without this overrides. */     

    /* TODO: add default file location */

    /* check for relative path name */
    ierr = PetscSNPrintf(buf1,sizeof(buf1),"%s%s%s%s","-",prefix,property,"_filename_rel_to_src");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,buf1,lookup_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);CHKERRQ(ierr);
    ierr = MakeRelativeToSourcePathAbsolute(lookup_filename);CHKERRQ(ierr);
    /* check for absolute path name */
    ierr = PetscSNPrintf(buf2,sizeof(buf2),"%s%s%s%s","-",prefix,property,"_filename");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,buf2,lookup_filename,PETSC_MAX_PATH_LEN,&set);CHKERRQ(ierr);

    /* if IS_SET is NULL, then we require a valid lookup_filename to be returned */
    if( IS_SET==NULL ){
        /* must return a valid lookup_filename */
        if ( !set && !set_rel_to_src ){
            SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"Missing argument %s or %s",buf1,buf2);
        }
    }

    /* if IS_SET is not NULL, then a valid lookup_filename is optional */
    if( IS_SET!=NULL ){
        if( set || set_rel_to_src ){
            *IS_SET = PETSC_TRUE;
        }
        else{
            *IS_SET = PETSC_FALSE;
        }
    }

    /* absolute path name always overrides relative path name */
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"%s%s%s%s%s","Warning: ",buf1," ignored because ",buf2," provided\n");CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

static PetscErrorCode SetLookupTemperature( const Lookup lookup, PetscScalar P, PetscScalar S, PetscScalar *T)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = SetInterp2dValue( lookup->temp, P, S, T );CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode SetLookupRho( const Lookup lookup, PetscScalar P, PetscScalar S, PetscScalar *rho)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = SetInterp2dValue( lookup->rho, P, S, rho );CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode SetLookupAlpha( const Lookup lookup, PetscScalar P, PetscScalar S, PetscScalar *alpha)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = SetInterp2dValue( lookup->alpha, P, S, alpha );CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode SetLookupCp( const Lookup lookup, PetscScalar P, PetscScalar S, PetscScalar *Cp)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = SetInterp2dValue( lookup->cp, P, S, Cp );CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode SetLookupdTdPs( const Lookup lookup, PetscScalar P, PetscScalar S, PetscScalar *dTdPs)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = SetInterp2dValue( lookup->dTdPs, P, S, dTdPs );CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode SetEosEvalFromLookup( const Lookup lookup, PetscScalar P, PetscScalar S, EosEval *eos_eval )
{
    PetscFunctionBeginUser;

    eos_eval->P = P;
    eos_eval->S = S;
    SetLookupTemperature( lookup, P, S, &eos_eval->T );
    SetLookupCp( lookup, P, S, &eos_eval->Cp );
    SetLookupRho( lookup, P, S, &eos_eval->rho );
    SetLookupdTdPs( lookup, P, S, &eos_eval->dTdPs );
    SetLookupAlpha( lookup, P, S, &eos_eval->alpha );
    /* lookup does not know about these quantities, since they are not used by
       SPIDER */
    eos_eval->Cv = 0.0;
    eos_eval->V = 0.0;

    PetscFunctionReturn(0);
}
