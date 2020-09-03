#include "eos.h"
#include "eos_composite.h" // TODO ultimately don't want this here
#include "eos_lookup.h" // TODO ultimately don't want this here
#include "eos_rtpress.h" // TODO ultimately don't want this here
#include "monitor.h"
#include "parameters.h"
#include "util.h"
#include "interp.h" // TODO ultimately don't want this here

/* Protypes for helpers for EOS interface functions */
static PetscErrorCode EOSEval_SetViscosity(EOS,EosEval*);

/* EOS interface functions (public API) */
PetscErrorCode EOSCreate(EOS* p_eos, EOSType type)
{
  PetscErrorCode ierr;
  PetscBool      flg;
  EOS            eos;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, p_eos);CHKERRQ(ierr);
  eos = *p_eos;
  eos->is_setup = PETSC_FALSE;
  ierr = PetscStrcmp(type, SPIDER_EOS_LOOKUP, &flg);CHKERRQ(ierr);
  if (flg) {
    ierr = EOSCreate_Lookup(eos);CHKERRQ(ierr);
  }
  ierr = PetscStrcmp(type, SPIDER_EOS_RTPRESS, &flg);CHKERRQ(ierr);
  if (flg) {
    ierr = EOSCreate_RTpress(eos);CHKERRQ(ierr);
  }
  ierr = PetscStrcmp(type, SPIDER_EOS_COMPOSITE, &flg);CHKERRQ(ierr);
  if (flg) {
    ierr = EOSCreate_Composite(eos);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EOSEval(const EOS eos, PetscScalar P , PetscScalar S, EosEval* eval)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  if (!eos->is_setup) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONGSTATE,"You must cannot evaluate the EOS before setting it up");

  /* Implementation-specific logic */
  ierr = (*(eos->eval))(eos,P,S,eval);CHKERRQ(ierr);

  /* Common Logic */
  eval->cond = eos->cond; // conductivity constant
  SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Not Implemented!");
  // TODO duplicate this function
  //ierr = SetEosEvalViscosity( Ep, eval );CHKERRQ(ierr);
  eval->phase_fraction = 1.0; // by definition, since only one phase
  eval->fusion = 0.0; // meaningless for a single phase
  PetscFunctionReturn(0);
}

PetscErrorCode EOSDestroy(EOS *p_eos)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  {
    EOS eos = *p_eos;

    ierr = (*eos->destroy)(eos);CHKERRQ(ierr);
    if (eos->PHASE_BOUNDARY) {
      ierr = Interp1dDestroy(&eos->phase_boundary);CHKERRQ(ierr);
    }
  }
  ierr = PetscFree(p_eos);CHKERRQ(ierr);
  *p_eos = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode EOSSetUpFromOptions(EOS eos, const char *prefix, const FundamentalConstants FC, const ScalingConstants SC)
{
  PetscErrorCode ierr;
  char           buf[1024]; /* max size */

  PetscFunctionBegin;
  if (eos->setupfromoptions) {
    ierr = (*eos->setupfromoptions)(eos, prefix, FC, SC);CHKERRQ(ierr);
  }

  /* conductivity (w/m/K) */
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_cond");CHKERRQ(ierr);
  eos->cond = 4.0; 
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&eos->cond,NULL);CHKERRQ(ierr);
  eos->cond /= SC->COND;

  /* viscosity-related, may eventually move into their own struct */
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_log10visc");CHKERRQ(ierr);
  eos->log10visc = 21.0; // FIXME: default is for solid only
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&eos->log10visc,NULL);CHKERRQ(ierr);
  eos->log10visc -= SC->LOG10VISC;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_activation_energy");CHKERRQ(ierr);
  eos->activation_energy = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&eos->activation_energy,NULL);CHKERRQ(ierr);
  /* scalings for Ea and Va include the gas constant that appears in the
     denominator of the Arrhenius viscosity law, so we then do not have
     to pass FC->GAS to the viscosity functions */
  eos->activation_energy /= SC->ENERGY * FC->GAS;

  /* activation volume (m^3/mol) */
  /* The numerical value in units of m^3/mol is the same as that in units of J/mol/Pa */
  /* You can convince yourself of this by using the scalings for ENERGY and PRESSURE to
     see that this is true */
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_activation_volume");CHKERRQ(ierr);
  eos->activation_volume = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&eos->activation_volume,NULL);CHKERRQ(ierr);
  /* as with activation energy, include gas constant in denominator */
  eos->activation_volume *= SC->PRESSURE / (SC->ENERGY * FC->GAS);

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_activation_volume_pressure_scale");CHKERRQ(ierr);
  eos->activation_volume_pressure_scale = -1.0; /* negative is not set */
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&eos->activation_volume_pressure_scale,NULL);CHKERRQ(ierr);
  eos->activation_volume_pressure_scale /= SC->PRESSURE;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_visc_comp");CHKERRQ(ierr);
  eos->visc_comp = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&eos->visc_comp,NULL);CHKERRQ(ierr);
  /* no scaling necessary since this is a ratio */

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_visc_ref_temp");CHKERRQ(ierr);
  eos->visc_ref_temp = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&eos->visc_ref_temp,NULL);CHKERRQ(ierr);
  eos->visc_ref_temp /= SC->TEMP;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_visc_ref_pressure");CHKERRQ(ierr);
  eos->visc_ref_pressure = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&eos->visc_ref_pressure,NULL);CHKERRQ(ierr);
  eos->visc_ref_pressure /= SC->PRESSURE;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",prefix,"_visc_ref_comp");CHKERRQ(ierr);
  eos->visc_ref_comp = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&eos->visc_ref_comp,NULL);CHKERRQ(ierr);
  /* no scaling necessary since this is a ratio */

  /* phase boundary */
  ierr = LookupFilenameSet( "_phase_boundary", prefix, eos->phase_boundary_filename, &eos->PHASE_BOUNDARY );CHKERRQ(ierr);
  if( eos->PHASE_BOUNDARY ){
    ierr = Interp1dCreateAndSet( eos->phase_boundary_filename, &eos->phase_boundary, SC->PRESSURE, SC->ENTROPY );CHKERRQ(ierr);
  }

  /* Mark as set up */
  eos->is_setup = PETSC_TRUE;

  PetscFunctionReturn(0);
}

PetscErrorCode EOSGetPhaseBoundary(EOS eos, PetscScalar P, PetscScalar *boundary, PetscScalar *dboundary)
{
  PetscFunctionBeginUser;
  SetInterp1dValue(eos->phase_boundary, P, boundary, dboundary); /* entropy S and derivative dS/dP */
  PetscFunctionReturn(0);
}

/* Helper Functions */
static PetscScalar GetCompositionalViscosityPrefactor( PetscScalar Mg_Si ){

    /* These expressions were worked out by Rob Spaargaren as part
       of his MSc thesis (2018) are are explained in Spaargaren et al. (2020) */

    /* Mg_Si is molar mantle Mg/Si */

    PetscScalar fac;

    if (Mg_Si <= 1.0)
        fac = 0.5185 * (1 - Mg_Si)/0.3; //
    else if (Mg_Si <= 1.25)
        /* Earth has Mg/Si = 1.08 */
        fac = -1.4815 * (Mg_Si - 1)/0.25; // -1.4815 = log10(0.033)
    else if (Mg_Si <= 1.5)
        fac = -2 + (0.5185) * (1.5 - Mg_Si)/0.25; // 0.5185 = log10(0.033) - -2
    else
        /* Fp-rich composition (Ballmer et al. 2017) */
        fac = -2;

/* this is the original formulation used for the draft version of Rob's paper
   during the review stage, the formulation was changed to that above */
#if 0
    if(Mg_Si <= 0.5)
        /* St-rich composition (Xu et al., 2017) */
        fac = 2;
    else if (Mg_Si <= 0.7)
        fac = 2 - 1.4815 * (Mg_Si - 0.5)/0.2; // 1.4815 = 2 - log10(3.3)
    else if (Mg_Si <= 1.0)
        /* fac is zero for Mg_Si = 1.0 */
        fac = 0.5185 * (1 - Mg_Si)/0.3; // 0.5185 = log10(3.3)
    else if (Mg_Si <= 1.25)
        /* Earth has Mg/Si = 1.08 */
        fac = -1.4815 * (Mg_Si - 1)/0.25; // -1.4815 = log10(0.033)
    else if (Mg_Si <= 1.5)
        fac = -2 + (0.5185) * (1.5 - Mg_Si)/0.25; // 0.5185 = log10(0.033) - -2
    else
        /* Fp-rich composition (Ballmer et al. 2017) */
        fac = -2;
#endif

    return fac;
}

static PetscErrorCode EOSEval_SetViscosity(EOS eos, EosEval *eval)
{
    PetscScalar A, log10C, dP, dT;
    PetscScalar fac1 = 1.0, fac2 = 1.0;

    PetscFunctionBeginUser;

    /* reference viscosity */
    eval->log10visc = eos->log10visc; // i.e., log10(eta_0)

    /* temperature and pressure contribution
       A(T,P) = (E_a + V_a P) / RT - (E_a + V_a Pref)/ RTref
       eta = eta_0 * exp(A)
       log10(eta) = log10(eta0) + log10(exp(A))
       log10(eta) = P->eos2_parameters.log10visc + A/ln(10) */

    /* with Ps = activation_volume_pressure_scale:
           V_a(P) = V_a exp (-P/Ps) */

    A = 0.0;

    /* pin viscosity profile to reference values */
    if( (eos->visc_ref_pressure >= 0.0) || (eos->visc_ref_temp >= 0.0) ){
        dT = ( eos->visc_ref_temp - eval->T ) / eos->visc_ref_temp;
        if( eos->activation_volume_pressure_scale > 0.0 ){
            fac1 = PetscExpReal( -eval->P / eos->activation_volume_pressure_scale );
            fac2 = PetscExpReal( -eos->visc_ref_pressure / eos->activation_volume_pressure_scale );
        }
        /* else fac1 and fac2 retain unity scalings according to initialisation above */
        dP = fac1 * eval->P - fac2 * eos->visc_ref_pressure * (eval->T / eos->visc_ref_temp );
    }
    /* do not pin viscosity profile */
    else{
        dT = 1.0;
        dP = eval->P;
        if( eos->activation_volume_pressure_scale > 0.0 ){
            dP *= PetscExpReal( -eval->P / eos->activation_volume_pressure_scale );
        }
    }

    if( eos->activation_energy > 0.0){
        A += eos->activation_energy * dT;
    }
    if( eos->activation_volume > 0.0){
        A += eos->activation_volume * dP;
    }

    /* division by R (gas constant) was already done during the scaling of parameters */
    A *= 1.0 / eval->T;
    eval->log10visc += A / PetscLogReal(10.0);

    /* compositional (Mg/Si) contribution */
    /* always pinned to some reference given by eos->visc_ref_comp */
    if( eos->visc_comp > 0.0 ){
        log10C = GetCompositionalViscosityPrefactor( eos->visc_comp );
        log10C -= GetCompositionalViscosityPrefactor( eos->visc_ref_comp );
        eval->log10visc += log10C;
    }

    /* TODO: add viscous lid */

    /* TODO: add viscosity cutoff */

    PetscFunctionReturn(0);

}


// TODO --- below here, eventually remove, once the new EOS class is complete ---------------

/* TODO: many of these Getters should actually be Setters, since they
   set the final argument rather than returning a value */



/* evaluate viscosity */
static PetscScalar GetCompositionalViscosityPrefactor( PetscScalar );

/* two phase composite eos (for mixed phase region) */
/* TODO: you'll see that the convention is for the melt phase to be listed first */
static PetscErrorCode SetEosCompositeEvalFromTwoPhase( const EosComposite, PetscScalar, PetscScalar, EosEval *);

PetscErrorCode EosParametersCreate( EosParameters* eos_parameters_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscMalloc1(1,eos_parameters_ptr);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode EosParametersDestroy( EosParameters* eos_parameters_ptr )
{
    PetscErrorCode ierr;
    EosParameters eos_parameters;
    Lookup * lookup_ptr;
    Lookup lookup;

    PetscFunctionBeginUser;

    eos_parameters = *eos_parameters_ptr;

    switch( eos_parameters->TYPE ){
        case 1:
            lookup_ptr = &eos_parameters->lookup;
            lookup = *lookup_ptr;
            ierr = Interp2dDestroy( &lookup->alpha ); CHKERRQ(ierr);
            ierr = Interp2dDestroy( &lookup->cp ); CHKERRQ(ierr);
            ierr = Interp2dDestroy( &lookup->dTdPs ); CHKERRQ(ierr);
            ierr = Interp2dDestroy( &lookup->rho ); CHKERRQ(ierr);
            ierr = Interp2dDestroy( &lookup->temp ); CHKERRQ(ierr);
            ierr = LookupDestroy( lookup_ptr );CHKERRQ(ierr);
            break;

        case 2:
            ierr = RTpressParametersDestroy( &eos_parameters->rtpress_parameters );CHKERRQ(ierr);
            break;
    }

    if( eos_parameters->PHASE_BOUNDARY ){
        ierr = Interp1dDestroy( &eos_parameters->phase_boundary ); CHKERRQ(ierr);
    }

    ierr = PetscFree(*eos_parameters_ptr);CHKERRQ(ierr);
    *eos_parameters_ptr = NULL;

    PetscFunctionReturn(0);
}


PetscErrorCode EosParametersSetFromOptions( EosParameters Ep, const FundamentalConstants FC, const ScalingConstants SC )
{
  /* creates and sets the structs that are nested within EosParameters */

  PetscErrorCode ierr;
  char           buf[1024]; /* max size */
  PetscBool      set; 

  PetscFunctionBeginUser;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",Ep->prefix,"_TYPE");CHKERRQ(ierr);
  Ep->TYPE = 1; /* default is lookup */
  ierr = PetscOptionsGetInt(NULL,NULL,buf, &Ep->TYPE,&set);CHKERRQ(ierr);
  //if (!set) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"Missing argument %s",bur);

  /* probably a good idea to initialise pointers to NULL */
  Ep->rtpress_parameters = NULL;
  Ep->lookup = NULL;

  switch( Ep->TYPE ){
      case 1:
          /* lookup, set filenames (does not allocate memory for Interp structs) */
          /* leading underscore is clunky, but to enable the same function to
             process a variety of input strings */
          ierr = LookupCreate( &Ep->lookup );CHKERRQ(ierr);
          ierr = LookupFilenameSet( "_alpha", Ep->prefix, Ep->lookup->alpha_filename, NULL );CHKERRQ(ierr);
          ierr = Interp2dCreateAndSet( Ep->lookup->alpha_filename, &Ep->lookup->alpha, SC->PRESSURE, SC->ENTROPY, 1.0/SC->TEMP );CHKERRQ(ierr);
          ierr = LookupFilenameSet( "_cp", Ep->prefix, Ep->lookup->cp_filename, NULL ); CHKERRQ(ierr);
          ierr = Interp2dCreateAndSet( Ep->lookup->cp_filename, &Ep->lookup->cp, SC->PRESSURE, SC->ENTROPY, SC->ENTROPY );CHKERRQ(ierr);
          ierr = LookupFilenameSet( "_dTdPs", Ep->prefix, Ep->lookup->dTdPs_filename, NULL  );CHKERRQ(ierr);
          ierr = Interp2dCreateAndSet( Ep->lookup->dTdPs_filename, &Ep->lookup->dTdPs, SC->PRESSURE, SC->ENTROPY, SC->DTDP );CHKERRQ(ierr);
          ierr = LookupFilenameSet( "_rho", Ep->prefix, Ep->lookup->rho_filename, NULL );CHKERRQ(ierr);
          ierr = Interp2dCreateAndSet( Ep->lookup->rho_filename, &Ep->lookup->rho, SC->PRESSURE, SC->ENTROPY, SC->DENSITY );CHKERRQ(ierr);
          ierr = LookupFilenameSet( "_temp", Ep->prefix, Ep->lookup->temp_filename, NULL );CHKERRQ(ierr);
          ierr = Interp2dCreateAndSet( Ep->lookup->temp_filename, &Ep->lookup->temp, SC->PRESSURE, SC->ENTROPY, SC->TEMP );CHKERRQ(ierr);
          break;

      case 2:
          /* analytical RTpress */
          ierr = RTpressParametersCreateAndSet( &Ep->rtpress_parameters, FC );CHKERRQ(ierr);
          break;
  }

  /* conductivity (w/m/K) */
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",Ep->prefix,"_cond");CHKERRQ(ierr);
  Ep->cond = 4.0; 
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&Ep->cond,NULL);CHKERRQ(ierr);
  Ep->cond /= SC->COND;

  /* viscosity-related, may eventually move into their own struct */
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",Ep->prefix,"_log10visc");CHKERRQ(ierr);
  Ep->log10visc = 21.0; // FIXME: default is for solid only
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&Ep->log10visc,NULL);CHKERRQ(ierr);
  Ep->log10visc -= SC->LOG10VISC;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",Ep->prefix,"_activation_energy");CHKERRQ(ierr);
  Ep->activation_energy = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&Ep->activation_energy,NULL);CHKERRQ(ierr);
  /* scalings for Ea and Va include the gas constant that appears in the
     denominator of the Arrhenius viscosity law, so we then do not have
     to pass FC->GAS to the viscosity functions */
  Ep->activation_energy /= SC->ENERGY * FC->GAS;

  /* activation volume (m^3/mol) */
  /* The numerical value in units of m^3/mol is the same as that in units of J/mol/Pa */
  /* You can convince yourself of this by using the scalings for ENERGY and PRESSURE to
     see that this is true */
  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",Ep->prefix,"_activation_volume");CHKERRQ(ierr);
  Ep->activation_volume = 0.0;
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&Ep->activation_volume,NULL);CHKERRQ(ierr);
  /* as with activation energy, include gas constant in denominator */
  Ep->activation_volume *= SC->PRESSURE / (SC->ENERGY * FC->GAS);

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",Ep->prefix,"_activation_volume_pressure_scale");CHKERRQ(ierr);
  Ep->activation_volume_pressure_scale = -1.0; /* negative is not set */
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&Ep->activation_volume_pressure_scale,NULL);CHKERRQ(ierr);
  Ep->activation_volume_pressure_scale /= SC->PRESSURE;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",Ep->prefix,"_visc_comp");CHKERRQ(ierr);
  Ep->visc_comp = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&Ep->visc_comp,NULL);CHKERRQ(ierr);
  /* no scaling necessary since this is a ratio */

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",Ep->prefix,"_visc_ref_temp");CHKERRQ(ierr);
  Ep->visc_ref_temp = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&Ep->visc_ref_temp,NULL);CHKERRQ(ierr);
  Ep->visc_ref_temp /= SC->TEMP;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",Ep->prefix,"_visc_ref_pressure");CHKERRQ(ierr);
  Ep->visc_ref_pressure = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&Ep->visc_ref_pressure,NULL);CHKERRQ(ierr);
  Ep->visc_ref_pressure /= SC->PRESSURE;

  ierr = PetscSNPrintf(buf,sizeof(buf),"%s%s%s","-",Ep->prefix,"_visc_ref_comp");CHKERRQ(ierr);
  Ep->visc_ref_comp = -1.0; // negative is not set
  ierr = PetscOptionsGetScalar(NULL,NULL,buf,&Ep->visc_ref_comp,NULL);CHKERRQ(ierr);
  /* no scaling necessary since this is a ratio */

  /* phase boundary */
  ierr = LookupFilenameSet( "_phase_boundary", Ep->prefix, Ep->phase_boundary_filename, &Ep->PHASE_BOUNDARY );CHKERRQ(ierr);
  if( Ep->PHASE_BOUNDARY ){
      ierr = Interp1dCreateAndSet( Ep->phase_boundary_filename, &Ep->phase_boundary, SC->PRESSURE, SC->ENTROPY );CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode SetEosEval( const EosParameters Ep, PetscScalar P, PetscScalar S, EosEval *eos_eval )
{
  /* TODO: could instead use function pointers, rather than pushing every evaluation of an EOS
     through this switch case */

  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  switch( Ep->TYPE ){
      case 1:
          /* lookup */
           ierr = SetEosEvalFromLookup( Ep->lookup, P, S, eos_eval );CHKERRQ(ierr);
           break;
      case 2:
          /* analytical RTpress */
          ierr = SetEosEvalFromRTpress( Ep->rtpress_parameters, P, S, eos_eval );CHKERRQ(ierr);
          break;
  }

  eos_eval->cond = Ep->cond; // conductivity constant
  ierr = SetEosEvalViscosity( Ep, eos_eval );CHKERRQ(ierr);
  eos_eval->phase_fraction = 1.0; // by definition, since only one phase
  eos_eval->fusion = 0.0; // meaningless for a single phase

  PetscFunctionReturn(0);
}


PetscErrorCode SetEosEvalViscosity( const EosParameters Ep, EosEval *eos_eval )
{
    PetscScalar A, log10C, dP, dT;
    PetscScalar fac1 = 1.0, fac2 = 1.0;

    PetscFunctionBeginUser;

    /* reference viscosity */
    eos_eval->log10visc = Ep->log10visc; // i.e., log10(eta_0)

    /* temperature and pressure contribution
       A(T,P) = (E_a + V_a P) / RT - (E_a + V_a Pref)/ RTref
       eta = eta_0 * exp(A)
       log10(eta) = log10(eta0) + log10(exp(A))
       log10(eta) = P->eos2_parameters.log10visc + A/ln(10) */

    /* with Ps = activation_volume_pressure_scale:
           V_a(P) = V_a exp (-P/Ps) */

    A = 0.0;

    /* pin viscosity profile to reference values */
    if( (Ep->visc_ref_pressure >= 0.0) || (Ep->visc_ref_temp >= 0.0) ){
        dT = ( Ep->visc_ref_temp - eos_eval->T ) / Ep->visc_ref_temp;
        if( Ep->activation_volume_pressure_scale > 0.0 ){
            fac1 = PetscExpReal( -eos_eval->P / Ep->activation_volume_pressure_scale );
            fac2 = PetscExpReal( -Ep->visc_ref_pressure / Ep->activation_volume_pressure_scale );
        }
        /* else fac1 and fac2 retain unity scalings according to initialisation above */
        dP = fac1 * eos_eval->P - fac2 * Ep->visc_ref_pressure * (eos_eval->T / Ep->visc_ref_temp );
    }
    /* do not pin viscosity profile */
    else{
        dT = 1.0;
        dP = eos_eval->P;
        if( Ep->activation_volume_pressure_scale > 0.0 ){
            dP *= PetscExpReal( -eos_eval->P / Ep->activation_volume_pressure_scale );
        }
    }

    if( Ep->activation_energy > 0.0){
        A += Ep->activation_energy * dT;
    }
    if( Ep->activation_volume > 0.0){
        A += Ep->activation_volume * dP;
    }

    /* division by R (gas constant) was already done during the scaling of parameters */
    A *= 1.0 / eos_eval->T;
    eos_eval->log10visc += A / PetscLogReal(10.0);

    /* compositional (Mg/Si) contribution */
    /* always pinned to some reference given by Ep->visc_ref_comp */
    if( Ep->visc_comp > 0.0 ){
        log10C = GetCompositionalViscosityPrefactor( Ep->visc_comp );
        log10C -= GetCompositionalViscosityPrefactor( Ep->visc_ref_comp );
        eos_eval->log10visc += log10C;
    }

    /* TODO: add viscous lid */

    /* TODO: add viscosity cutoff */

    PetscFunctionReturn(0);

}




