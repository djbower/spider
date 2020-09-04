#include "eos_lookup.h"
#include "util.h"

/* Prototypes for helpers used in EOS interface functions */
static PetscErrorCode EOSLookup_FilenameSet(const char*,const char*,char*,PetscBool*);

/* EOS interface functions */
static PetscErrorCode EOSEval_Lookup(EOS eos, PetscScalar P, PetscScalar S, EosEval *eval)
{
  PetscErrorCode  ierr;
  data_EOSLookup *lookup = (data_EOSLookup*) eos->impl_data;

  PetscFunctionBegin;
  eval->P = P;
  eval->S = S;
  ierr = SetInterp2dValue( lookup->temp, P, S, &eval->T );CHKERRQ(ierr);
  ierr = SetInterp2dValue( lookup->rho, P, S, &eval->rho );CHKERRQ(ierr);
  ierr = SetInterp2dValue( lookup->cp, P, S, &eval->Cp );CHKERRQ(ierr);
  ierr = SetInterp2dValue( lookup->dTdPs, P, S, &eval->dTdPs );CHKERRQ(ierr);
  ierr = SetInterp2dValue( lookup->alpha, P, S, &eval->alpha );CHKERRQ(ierr);
  /* lookup does not know about these quantities, since they are not used by
     SPIDER */
  eval->Cv = 0.0;
  eval->V = 0.0;
  PetscFunctionReturn(0);
}

static PetscErrorCode EOSDestroy_Lookup(EOS eos)
{
  PetscErrorCode ierr;
  data_EOSLookup *lookup = (data_EOSLookup*) eos->impl_data;

  PetscFunctionBegin;
  ierr = Interp2dDestroy(&lookup->alpha); CHKERRQ(ierr);
  ierr = Interp2dDestroy(&lookup->cp); CHKERRQ(ierr);
  ierr = Interp2dDestroy(&lookup->dTdPs); CHKERRQ(ierr);
  ierr = Interp2dDestroy(&lookup->rho); CHKERRQ(ierr);
  ierr = Interp2dDestroy(&lookup->temp); CHKERRQ(ierr);
  ierr = PetscFree(eos->impl_data);CHKERRQ(ierr);
  eos->impl_data = NULL;
  PetscFunctionReturn(0);
}


PetscErrorCode EOSSetUpFromOptions_Lookup(EOS eos, const char *prefix, const FundamentalConstants FC, const ScalingConstants SC)
{
  PetscErrorCode  ierr;
  data_EOSLookup *data = (data_EOSLookup*) eos->impl_data;

  PetscFunctionBegin;
  (void) FC; // unused
  /* lookup, set filenames (does not allocate memory for Interp structs) */
  /* leading underscore is clunky, but to enable the same function to
     process a variety of input strings */
  ierr = EOSLookup_FilenameSet( "_alpha", prefix, data->alpha_filename, NULL );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( data->alpha_filename, &data->alpha, SC->PRESSURE, SC->ENTROPY, 1.0/SC->TEMP );CHKERRQ(ierr);
  ierr = EOSLookup_FilenameSet( "_cp", prefix, data->cp_filename, NULL ); CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( data->cp_filename, &data->cp, SC->PRESSURE, SC->ENTROPY, SC->ENTROPY );CHKERRQ(ierr);
  ierr = EOSLookup_FilenameSet( "_dTdPs", prefix, data->dTdPs_filename, NULL  );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( data->dTdPs_filename, &data->dTdPs, SC->PRESSURE, SC->ENTROPY, SC->DTDP );CHKERRQ(ierr);
  ierr = EOSLookup_FilenameSet( "_rho", prefix, data->rho_filename, NULL );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( data->rho_filename, &data->rho, SC->PRESSURE, SC->ENTROPY, SC->DENSITY );CHKERRQ(ierr);
  ierr = EOSLookup_FilenameSet( "_temp", prefix, data->temp_filename, NULL );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( data->temp_filename, &data->temp, SC->PRESSURE, SC->ENTROPY, SC->TEMP );CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Creation Function */
PetscErrorCode EOSCreate_Lookup(EOS eos) {
  PetscFunctionBeginUser;
  eos->eval = EOSEval_Lookup;
  eos->destroy = EOSDestroy_Lookup;
  eos->setupfromoptions = EOSSetUpFromOptions_Lookup;
  PetscFunctionReturn(0);
}

/* Helper functions */
static PetscErrorCode EOSLookup_FilenameSet( const char* property, const char* prefix, char* lookup_filename, PetscBool *IS_SET )
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

// TODO ---- below here, remove once new EOS class is finished ---------

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
