#include "eos_lookup.h"
#include "util.h"

/* Prototypes for helpers used in EOS interface functions */
static PetscErrorCode EOSLookup_FilenameSet(const char*,const char*,char*,PetscBool*);

/* EOS interface functions */
static PetscErrorCode EOSEval_Lookup(EOS eos, PetscScalar P, PetscScalar S, EOSEvalData *eval)
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
  ierr = EOSEvalSetViscosity(eos, eval);CHKERRQ(ierr);
  eval->cond = eos->cond; // conductivity constant
  eval->phase_fraction = 1.0; // by definition, since only one phase
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
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, (data_EOSLookup**) (&eos->impl_data));CHKERRQ(ierr);
  eos->eval = EOSEval_Lookup;
  eos->destroy = EOSDestroy_Lookup;
  eos->setupfromoptions = EOSSetUpFromOptions_Lookup;
  PetscFunctionReturn(0);
}

/* Helper functions */
static PetscErrorCode EOSLookup_FilenameSet( const char* property, const char* prefix, char* lookup_filename, PetscBool *IS_SET )
{
    PetscErrorCode ierr;
    char           buf1[PETSC_MAX_PATH_LEN]; /* max size */
    char           buf2[PETSC_MAX_PATH_LEN]; /* max size */
    PetscBool      set_rel_to_src,set;

    PetscFunctionBeginUser;

    /* Based on input options, determine which files to load.  Options ending
       with _rel_to_src indicate a path relative to the source code. In this 
       case we prepend a string, SPIDER_ROOT_DIR_STR, and /. The corresponding
       option without this overrides. */     

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
