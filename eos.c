#include "eos.h"
#include "eos_composite.h"
#include "eos_lookup.h"
#include "eos_rtpress.h"
#include "monitor.h"
#include "parameters.h"
#include "util.h"
#include "interp.h"

/* TODO: many of these Getters should actually be Setters, since they
   set the final argument rather than returning a value */

/* First material model option is via lookup tables */
/* lookup material properties (default) */
static PetscErrorCode LookupCreate( Lookup * );
static PetscErrorCode LookupDestroy( Lookup * );
static PetscErrorCode LookupFilenameSet( const char *, const char *, char *, PetscBool * );
/* TODO: these next set of functions can be redone as function pointers in the EosParameters struct
   when a lookup is used */
static PetscErrorCode SetLookupRho( const Lookup, PetscScalar, PetscScalar, PetscScalar *);
static PetscErrorCode SetLookupAlpha( const Lookup, PetscScalar, PetscScalar, PetscScalar *);
static PetscErrorCode SetLookupCp( const Lookup, PetscScalar, PetscScalar, PetscScalar *);
static PetscErrorCode SetLookupdTdPs( const Lookup, PetscScalar, PetscScalar, PetscScalar *);
static PetscErrorCode SetLookupTemperature( const Lookup, PetscScalar, PetscScalar, PetscScalar *);
/* we use all the above to set all applicable values in an EosEval struct */
static PetscErrorCode SetEosEvalFromLookup( const Lookup, PetscScalar, PetscScalar, EosEval * );

/* Second material model option is an analytical model, which for here is rtpress (but
   in principle, others could be added) */
/* rtpress material properties (Wolf and Bower, 2018) */
static PetscErrorCode RTpressParametersCreate( RTpressParameters * );
static PetscErrorCode RTpressParametersCreateAndSet( RTpressParameters *, const FundamentalConstants );
static PetscErrorCode RTpressParametersDestroy( RTpressParameters * );
/* compared to the lookup approach above, rtpress is derived in terms of temperature and volume, which
   means we must take the current pressure and entropy and jointly solve for T and V.  Then we can use this T
   and V to solve for the other material properties.  This partly motivates the idea of updating material properties
   once in an EosEval struct, since otherwise we have to keep resolving for V and T from P and S, which will add
   a lot of additional solves when evaluating EOS */
static PetscErrorCode SetRTpressPressure( const RTpressParameters, PetscScalar, PetscScalar, PetscScalar * );
static PetscErrorCode SetRTpressEntropy( const RTpressParameters, PetscScalar, PetscScalar, PetscScalar * );
static PetscErrorCode SetRTpressVolumeTemperature( const RTpressParameters, PetscScalar, PetscScalar, PetscScalar *, PetscScalar * );
/* arguments to these functions are V and T and not the primary solution variables of P and S */
static PetscErrorCode SetRTpressRho( const RTpressParameters, PetscScalar, PetscScalar, PetscScalar * );
static PetscErrorCode SetRTpressAlpha( const RTpressParameters, PetscScalar, PetscScalar, PetscScalar * );
static PetscErrorCode SetRTpressCv( const RTpressParameters, PetscScalar, PetscScalar, PetscScalar * );
static PetscErrorCode SetRTpressCp( const RTpressParameters, PetscScalar, PetscScalar, PetscScalar * );
static PetscErrorCode SetRTpressdTdPs( const RTpressParameters, PetscScalar, PetscScalar, PetscScalar * );
/* solve for volume and temperature from pressure and entropy */
static PetscErrorCode RTpressObjectiveFunctionVolumeTemperature( SNES, Vec, Vec, void * );
static PetscErrorCode SetEosEvalFromRTpress( const RTpressParameters, PetscScalar, PetscScalar, EosEval * );

/* evaluate viscosity */
static PetscScalar GetCompositionalViscosityPrefactor( PetscScalar );

/* two phase composite eos (for mixed phase region) */
/* TODO: you'll see that the convention is for the melt phase to be listed first */
static PetscErrorCode SetEosCompositeEvalFromTwoPhase( const EosComposite, PetscScalar, PetscScalar, EosEval *);

#if 0
/* TODO: update these to general framework for evaluating Eos */
static PetscErrorCode set_rtpress_struct_SI( PetscScalar, PetscScalar, Ctx * );
static PetscErrorCode set_rtpress_struct_non_dimensional( Ctx * );
#endif

/* helper functions */
static PetscScalar per_atom_to_specific( PetscScalar, PetscScalar, PetscScalar );
static PetscScalar specific_to_per_atom( PetscScalar, PetscScalar, PetscScalar );
static PetscScalar joule_to_eV( PetscScalar );
static PetscScalar eV_to_joule( PetscScalar );

PetscErrorCode SetPhaseBoundary( const EosParameters Ep, PetscScalar P, PetscScalar *boundary, PetscScalar *dboundary )
{

    /* TODO: this function could contain a switch, to determine phase boundary by
       a means other than lookup (currently not required) */

    PetscFunctionBeginUser;
    SetInterp1dValue( Ep->phase_boundary, P, boundary, dboundary ); /* entropy S and derivative dS/dP */
    PetscFunctionReturn(0);
}

/*
 ******************************************************************************
 * Lookup EOS from data files
 ******************************************************************************
*/

static PetscErrorCode LookupCreate( Lookup *lookup_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscMalloc1(1,lookup_ptr);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode LookupDestroy( Lookup *lookup_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscFree(*lookup_ptr);CHKERRQ(ierr);
    *lookup_ptr = NULL;
    PetscFunctionReturn(0);
}




static PetscErrorCode LookupFilenameSet( const char* property, const char* prefix, char* lookup_filename, PetscBool *IS_SET )
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

static PetscErrorCode SetEosEvalFromLookup( const Lookup lookup, PetscScalar P, PetscScalar S, EosEval *eos_eval )
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


/*
 ******************************************************************************
 * RTpress model (analytical) from Wolf and Bower (2018)
 ******************************************************************************
*/

/* TODO: this analytical model should be plugged into the general EosParameters
   workflow, although the RTpress model is currently not returning the correct
   values (need to debug units) */

static PetscErrorCode RTpressParametersCreate( RTpressParameters* rtpress_parameters_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscMalloc1(1,rtpress_parameters_ptr);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode RTpressParametersDestroy( RTpressParameters* rtpress_parameters_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    ierr = PetscFree(*rtpress_parameters_ptr);CHKERRQ(ierr);
    *rtpress_parameters_ptr = NULL;
    PetscFunctionReturn(0);
}

static PetscErrorCode RTpressParametersCreateAndSet( RTpressParameters* rtpress_parameters_ptr, const FundamentalConstants FC )
{
    /* EOS parameters for rtpress, taken from jupyter notebook and Wolf and Bower (2018).
       TODO: Simplest to keep these in dimensional form and scale the returned values as a
       last step? Ideally, would scale everything here instead, as with parameters.c */

    PetscErrorCode    ierr;
    RTpressParameters rtp;

    PetscFunctionBeginUser;

    ierr = RTpressParametersCreate( rtpress_parameters_ptr );CHKERRQ(ierr);
    rtp = *rtpress_parameters_ptr;

    /* for unit conversion */
    rtp->PV_UNIT = 160.21766208; /* GPa*Ang^3/eV */
    /* first value is KBOLTZ in units of eV/K */
    /* TODO: use BOLTZMANN from FundamentalConstants */
    rtp->KBOLTZ = 8.617333262145e-5 * rtp->PV_UNIT; /* GPa*Ang^3/K, this is energy / temperature */
    /* TODO: check with ASW: is 3000=T0 and 0.6=m? */
    /* 0.6 could instead be a reference compression? (V/V0) */
    /* still seems likely that T0 is a reference temperature scale */
    rtp->bscale = 0.6/3000 * rtp->PV_UNIT / rtp->KBOLTZ; /* K/eV, FIXME: if leading constant are unitless */

    rtp->V0 = 14.74; /* Ang^3/atom */
    rtp->T0 = 3000.0; /* K */
    rtp->S0 = 0.0; /* GPa*Ang^3/K, this is energy / temperature  */
    rtp->K0 = 9.77; /* GPa */
    rtp->KP0 = 7.42; /* non-dimensional */
    rtp->E0 = -6.85 * rtp->PV_UNIT; /* GPa*Ang^3/atom */
    rtp->gamma0 = 0.282; /* non-dimensional */
    rtp->gammaP0 = -1.35; /* non-dimensional, since gammaP0 = V0 (dgamma/dV)_0 */
    rtp->m = 0.6; /* non-dimensional */
    rtp->b0 = 1.118 * rtp->bscale; /* FIXME: units? */
    rtp->b1 = -0.05 * rtp->bscale; /* FIXME: units? */
    rtp->b2 = 2.1 * rtp->bscale; /* FIXME: units? */
    rtp->b3 = 12.9 * rtp->bscale; /* FIXME: units? */
    rtp->b4 = 15.5 * rtp->bscale; /* FIXME: units? */

    /* MgSiO3 has a molar mass of 100.39 g/mol and 5 atoms */
    /* So (average) molar mass of an atom is 20 g/mol/atom */
    rtp->mavg = 20.0E-3; /* kg/mol/atom */

    rtp->AVOGADRO_ptr = &FC->AVOGADRO;

    PetscFunctionReturn(0);
}

static PetscScalar per_atom_to_specific( PetscScalar per_atom_value, PetscScalar mavg, PetscScalar Avogadro )
{
    /* convert from per atom to specific */

    PetscScalar specific_value;
    specific_value = per_atom_value * Avogadro / mavg;
    return specific_value;
}

static PetscScalar specific_to_per_atom( PetscScalar specific_value, PetscScalar mavg, PetscScalar Avogadro )
{
    /* convert from specific to per atom */

    PetscScalar per_atom_value;
    per_atom_value = specific_value * mavg / Avogadro;
    return per_atom_value;
}

static PetscScalar joule_to_eV( PetscScalar joule )
{
    /* convert from joules to eV */
    return joule * 6.242E18;
}

static PetscScalar eV_to_joule( PetscScalar eV )
{
    /* convert from eV to joules */
    return eV * 1.60218E-19;
}

static PetscErrorCode SetRTpressPressure( const RTpressParameters rtp, PetscScalar V, PetscScalar T, PetscScalar *P )
{
    /* pressure = function( volume, temperature ) */

    PetscScalar const K0 = rtp->K0;
    PetscScalar const V0 = rtp->V0;
    PetscScalar const KP0 = rtp->KP0;
    PetscScalar const gamma0 = rtp->gamma0;
    PetscScalar const gammaP0 = rtp->gammaP0;
    PetscScalar const T0 = rtp->T0;
    PetscScalar const m = rtp->m;
    PetscScalar const b0 = rtp->b0;
    PetscScalar const b1 = rtp->b1;
    PetscScalar const b2 = rtp->b2;
    PetscScalar const b3 = rtp->b3;
    PetscScalar const b4 = rtp->b4;

    PetscFunctionBeginUser;

    /* FIXME: if I try using PetscCbrtReal then I get the following warning at compilation:
           warning: implicit declaration of function 'cbrt' is invalid in C99 [-Wimplicit-function-declaration] */
    /* TODO: should I use Real or PetscScalar functions?  i.e. PetscPowReal or PetscPowScalar? */

    *P =  -9*K0*V0*(-1.0/3.0*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpReal((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/V - 1.0/3.0*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpReal((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/V)/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)));

    PetscFunctionReturn(0);
}

static PetscErrorCode SetRTpressEntropy( const RTpressParameters rtp, PetscScalar V, PetscScalar T, PetscScalar *S )
{
    /* entropy = function( volume, temperature ) */

    PetscScalar const gamma0 = rtp->gamma0;
    PetscScalar const gammaP0 = rtp->gammaP0;
    PetscScalar const V0 = rtp->V0;
    PetscScalar const S0 = rtp->S0;
    PetscScalar const T0 = rtp->T0;
    PetscScalar const m = rtp->m;
    PetscScalar const b0 = rtp->b0;
    PetscScalar const b1 = rtp->b1;
    PetscScalar const b2 = rtp->b2;
    PetscScalar const b3 = rtp->b3;
    PetscScalar const b4 = rtp->b4;

    *S = S0 + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T)) + 0.041419469658752747*m/(T*(2*m - 2)) - 0.041419469658752747/(T*(2*m - 2))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*m*PetscLogScalar(T)/(2*m - 2) - 0.041419469658752747*m*PetscLogScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))/(2*m - 2) - 0.020709734829376374 - 0.041419469658752747*PetscLogScalar(T)/(2*m - 2) + 0.041419469658752747*PetscLogScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))/(2*m - 2) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T;

    PetscFunctionReturn(0);
}

#if 0
PetscScalar GetRTpressEntropy_test( Ctx *E )
{
    /* test function (lldb, then step through to evaluate S for
       the V, T conditions defined below.  Compare with entropy
       plot in Jupyter notebook*/

    /* seems to confirm S is correct */

    RTpressParameters rtp = E->parameters->eos_parameters[0]->rtpress_parameters;
    PetscScalar   S, T, V, volfrac;

    volfrac = 0.6;
    T = 2500.0;
    V = volfrac * rtp->V0;
    S = GetRTpressEntropy( rtp, V, T ); /* -0.025934761412030771 */
    S /= rtp->KBOLTZ; /* -1.8784471379548613 */

    volfrac = 0.8;
    T = 2500.0;
    V = volfrac * rtp->V0;
    S = GetRTpressEntropy( rtp, V, T ); /* -0.016066633649871438 */
    S /= rtp->KBOLTZ; /* -1.1637015477678558 */

    volfrac = 1.0;
    T = 2500.0;
    V = volfrac * rtp->V0;
    S = GetRTpressEntropy( rtp, V, T ); /* -0.010551508142269357 */
    S /= rtp->KBOLTZ; /* -0.76424263003857285 */

    volfrac = 0.6;
    T = 5000.0;
    V = volfrac * rtp->V0;
    S = GetRTpressEntropy( rtp, V, T ); /* 0.010228437833469547 */
    S /= rtp->KBOLTZ; /* 0.74084274263333616 */

    volfrac = 0.8;
    T = 5000.0;
    V = volfrac * rtp->V0;
    S = GetRTpressEntropy( rtp, V, T ); /* 0.021941049097879431 */
    S /= rtp->KBOLTZ; /* 1.5891837301622371 */

    volfrac = 1.0;
    T = 5000.0;
    V = volfrac * rtp->V0;
    S = GetRTpressEntropy( rtp, V, T ); /* 0.027130677516848993 */
    S /= rtp->KBOLTZ; /* 1.9650669895370627 */

    return S;
}
#endif

static PetscErrorCode SetRTpressRho( const RTpressParameters rtp, PetscScalar V, PetscScalar T, PetscScalar *rho_ptr )
{

    /* returns density in SI units, kg/m^3 */

    PetscScalar const mavg = rtp->mavg;
    PetscScalar rho = *rho_ptr;

    PetscFunctionBeginUser;

    rho = mavg / V; /* kg/Ang^3/mol */

    /* convert to SI */
    rho *= PetscPowScalar( 10.0, 30.0 ); /* Ang^3/m^3 */
    rho /= *rtp->AVOGADRO_ptr;

    PetscFunctionReturn(0);
}

static PetscErrorCode SetRTpressAlpha( const RTpressParameters rtp, PetscScalar V, PetscScalar T, PetscScalar *alpha_ptr )
{
    /* thermal expansion = function( volume, temperature ) */
    /* returns thermal expansion coefficient in SI units, 1/K */

    PetscScalar alpha = *alpha_ptr;
    PetscScalar const gamma0 = rtp->gamma0;
    PetscScalar const gammaP0 = rtp->gammaP0;
    PetscScalar const V0 = rtp->V0;
    PetscScalar const K0 = rtp->K0;
    PetscScalar const KP0 = rtp->KP0;
    PetscScalar const T0 = rtp->T0;
    PetscScalar const m = rtp->m;
    PetscScalar const b0 = rtp->b0;
    PetscScalar const b1 = rtp->b1;
    PetscScalar const b2 = rtp->b2;
    PetscScalar const b3 = rtp->b3;
    PetscScalar const b4 = rtp->b4;

    PetscFunctionBeginUser;

    alpha =  -(T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/T)/(V*(-9*K0*V0*((1.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2))/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)))));

    PetscFunctionReturn(0);
}

static PetscErrorCode SetRTpressCv( const RTpressParameters rtp, PetscScalar V, PetscScalar T, PetscScalar *Cv_ptr )
{
    /* thermal heat capacity at constant volume = function( volume, temperature ) */

    PetscScalar Cv = *Cv_ptr;
    PetscScalar const V0 = rtp->V0;
    PetscScalar const T0 = rtp->T0;
    PetscScalar const m = rtp->m;
    PetscScalar const b0 = rtp->b0;
    PetscScalar const b1 = rtp->b1;
    PetscScalar const b2 = rtp->b2;
    PetscScalar const b3 = rtp->b3;
    PetscScalar const b4 = rtp->b4;

    PetscFunctionBeginUser;

    Cv = 0.020709734829376374 + 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T;

    /* FIXME: possibly wrong */
    /* convert to SI */
    Cv = per_atom_to_specific( Cv, rtp->mavg, *rtp->AVOGADRO_ptr );
    Cv = eV_to_joule( Cv );

    PetscFunctionReturn(0);
}

static PetscErrorCode SetRTpressCp( const RTpressParameters rtp, PetscScalar V, PetscScalar T, PetscScalar *Cp_ptr )
{
    /* thermal heat capacity at constant pressure = function( volume, temperature ) */

    PetscScalar Cp = *Cp_ptr;
    PetscScalar const K0 = rtp->K0;
    PetscScalar const V0 = rtp->V0;
    PetscScalar const KP0 = rtp->KP0;
    PetscScalar const gamma0 = rtp->gamma0;
    PetscScalar const gammaP0 = rtp->gammaP0;
    PetscScalar const T0 = rtp->T0;
    PetscScalar const m = rtp->m;
    PetscScalar const b0 = rtp->b0;
    PetscScalar const b1 = rtp->b1;
    PetscScalar const b2 = rtp->b2;
    PetscScalar const b3 = rtp->b3;
    PetscScalar const b4 = rtp->b4;

    PetscFunctionBeginUser;

    Cp = (0.020709734829376374 + 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T)*(-T*PetscPowScalar(T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/T, 2)/((0.020709734829376374 + 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T)*(-9*K0*V0*((1.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2))/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))))) + 1);

    /* FIXME: possibly wrong */
    /* convert to SI */
    Cp = per_atom_to_specific( Cp, rtp->mavg, *rtp->AVOGADRO_ptr );
    Cp = eV_to_joule( Cp );

    PetscFunctionReturn(0);
}

static PetscErrorCode SetRTpressdTdPs( const RTpressParameters rtp, PetscScalar V, PetscScalar T, PetscScalar *dTdPs_ptr )
{

    PetscScalar dTdPs = *dTdPs_ptr;
    PetscScalar const K0 = rtp->K0;
    PetscScalar const V0 = rtp->V0;
    PetscScalar const KP0 = rtp->KP0;
    PetscScalar const gamma0 = rtp->gamma0;
    PetscScalar const gammaP0 = rtp->gammaP0;
    PetscScalar const T0 = rtp->T0;
    PetscScalar const m = rtp->m;
    PetscScalar const b0 = rtp->b0;
    PetscScalar const b1 = rtp->b1;
    PetscScalar const b2 = rtp->b2;
    PetscScalar const b3 = rtp->b3;
    PetscScalar const b4 = rtp->b4;

    PetscFunctionBeginUser;

    dTdPs = -T*(T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/T)/((0.020709734829376374 + 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T)*(-T*PetscPowScalar(T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/T, 2)/((0.020709734829376374 + 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T)*(-9*K0*V0*((1.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2))/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))))) + 1)*(-9*K0*V0*((1.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2))/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)))));

    PetscFunctionReturn(0);
}

#if 0
PetscScalar GetRTpressPressure_test( Ctx *E )
{
    /* test function (lldb, then step through to evaluate P for
       the V, T conditions defined below.  Compare with pressure
       plot in Jupyter notebook*/

    /* seems to confirm P is correct
       FIXME: need to truncate negative values?  Perhaps not for
       smoothness of solver during inversion? */

    RTpressParameters rtp = E->parameters->eos_parameters[0]->rtpress_parameters;
    PetscScalar   P, T, V, volfrac;

    volfrac = 0.6;
    T = 2500.0;
    V = volfrac * rtp->V0;
    P = GetRTpressPressure( V, T, rtp ); /* 26.964412351908095 */

    volfrac = 0.8;
    T = 2500.0;
    V = volfrac * rtp->V0;
    P = GetRTpressPressure( V, T, rtp ); /* 3.5394262799831875 */

    volfrac = 1.0;
    T = 2500.0;
    V = volfrac * rtp->V0;
    /* FIXME: note negative value */
    P = GetRTpressPressure( V, T, rtp ); /* -0.54571824296751803 */

    volfrac = 0.6;
    T = 5000.0;
    V = volfrac * rtp->V0;
    P = GetRTpressPressure( V, T, rtp ); /* 38.568776519196767 */

    volfrac = 0.8;
    T = 5000.0;
    V = volfrac * rtp->V0;
    P = GetRTpressPressure( V, T, rtp ); /* 10.39877886001743 */

    volfrac = 1.0;
    T = 5000.0;
    V = volfrac * rtp->V0;
    P = GetRTpressPressure( V, T, rtp ); /* 2.107013748263614 */

    return P;
}
#endif

static PetscErrorCode RTpressObjectiveFunctionVolumeTemperature( SNES snes, Vec x, Vec f, void *ptr)
{
    PetscErrorCode     ierr;
    const PetscScalar  *xx;
    PetscScalar        *ff;
    PetscScalar        V, T, P, S;
    RTpressEval        *rtpress_eval = (RTpressEval*) ptr;

    PetscFunctionBeginUser;

    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

    V = xx[0];
    T = xx[1];

    SetRTpressPressure( rtpress_eval->rtp, V, T, &P );
    SetRTpressEntropy( rtpress_eval->rtp, V, T, &S );

    /* compute residual */
    ff[0] = P - rtpress_eval->P;
    ff[1] = S - rtpress_eval->S;

    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode SetRTpressVolumeTemperature( const RTpressParameters rtp, PetscScalar P, PetscScalar S, PetscScalar *V, PetscScalar *T )
{
    PetscErrorCode ierr;
    SNES           snes;
    Vec            x,r;
    PetscScalar    *xx;
    PetscInt       i;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"GetRTpressVolumeTemperature()\n");CHKERRQ(ierr);

    ierr = SNESCreate( PETSC_COMM_WORLD, &snes );CHKERRQ(ierr);

    /* Use this to address this specific SNES (nonlinear solver) from the command
       line or options file, e.g. -atmosic_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"rtpress_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    ierr = VecSetSizes( x, PETSC_DECIDE, 2 );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    /* temporary struct to pass arguments to solver */
    RTpressEval rtpress_eval = {.P=P, .S=S, .rtp=rtp };

    ierr = SNESSetFunction(snes,r,RTpressObjectiveFunctionVolumeTemperature,&rtpress_eval);CHKERRQ(ierr);

    /* initialise vector x with initial guess */
    /* TODO: could get initial guess from reference pressure profile ? */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);

    /* initial guesses are reasonable mid-values of the model */
    xx[0] = 0.6 * rtp->V0; /* compression ratio of 0.6 */
    xx[1] = 3500.0; /* slightly higher temperature than reference profile */

    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* Inform the nonlinear solver to generate a finite-difference approximation
       to the Jacobian */
    /* TODO: using sympy we can presumably compute analytically the entries for the
       Jacobian to improve convergence and speed */
     ierr = PetscOptionsSetValue(NULL,"-rtpress_snes_mf",NULL);CHKERRQ(ierr);

    /* For solver analysis/debugging/tuning, activate a custom monitor with a flag */
    {   
      PetscBool flg = PETSC_FALSE;

      ierr = PetscOptionsGetBool(NULL,NULL,"-rtpress_snes_verbose_monitor",&flg,NULL);CHKERRQ(ierr);
      if (flg) {
        ierr = SNESMonitorSet(snes,SNESMonitorVerbose,NULL,NULL);CHKERRQ(ierr);
      }
    }

    /* Solve */
    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr); /* Picks up any additional options (note prefix) */
    ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);
    {   
      SNESConvergedReason reason;
      ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
      if (reason < 0) SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
          "Nonlinear solver didn't converge: %s\n",SNESConvergedReasons[reason]);
    }   

    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    for (i=0; i<2; ++i) {
        if( xx[i] < 0.0 ){
            /* Sanity check on solution (since it's non-unique) */
            SETERRQ2(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
                "Unphysical rtpress property: slot %d, x: %g",i,xx[i]);
        }
    }

    /* store solution in struct */
    *V = xx[0];
    *T = xx[1];

    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode SetEosEvalFromRTpress( const RTpressParameters rtp, PetscScalar P, PetscScalar S, EosEval *eos_eval )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    eos_eval->P = P;
    eos_eval->S = S;
    ierr = SetRTpressVolumeTemperature( rtp, P, S, &eos_eval->V, &eos_eval->T );CHKERRQ(ierr);
    ierr = SetRTpressdTdPs( rtp, eos_eval->V, eos_eval->T, &eos_eval->dTdPs );
    ierr = SetRTpressCp( rtp, eos_eval->V, eos_eval->T, &eos_eval->Cp );
    ierr = SetRTpressCv( rtp, eos_eval->V, eos_eval->T, &eos_eval->Cv );
    ierr = SetRTpressRho( rtp, eos_eval->V, eos_eval->T, &eos_eval->rho );
    ierr = SetRTpressAlpha( rtp, eos_eval->V, eos_eval->T, &eos_eval->alpha );

    PetscFunctionReturn(0);
}

#if 0
static PetscErrorCode set_rtpress_struct_SI( PetscScalar P, PetscScalar S, Ctx *E )
{
    /* solve for volume and temperature from pressure and entropy
       once, to avoid unnecessary computations.  This updates all the 
       eos_eval struct with the material properties */

    ScalingConstants const     SC = E->parameters->scaling_constants;
    RTpressParameters          rtp = E->parameters->eos_parameters[0]->rtpress_parameters;
    EosEval                    *eos_eval = &E->eos1_eval;

    PetscFunctionBeginUser;

    /* store values to lookup in struct, with correct units for rtpress
       analytical expressions */
    /* rtpress wants GPa */
    eos_eval->P = P;
    eos_eval->P *= SC->PRESSURE * 1.0E-9;

    /* rtpress wants eV/atom/K FIXME: K correct? */
    eos_eval->S = S;
    eos_eval->S *= SC->ENTROPY;
    eos_eval->S = specific_to_per_atom( eos_eval->S, rtp->mavg, *rtp->AVOGADRO_ptr ); 
    eos_eval->S = joule_to_eV( eos_eval->S );

    /* ensure solver is working */
    /* this test passes */
    eos_eval->P = 26.964412351908095;
    eos_eval->S = -0.025934761412030771;
    /* eos_eval->V/V0 = 0.6; */
    /* eos_eval->T = 2500; */

    /* below updates V and T in eos_eval */
    /* FIXME: cutoffs required for solver? */
    GetRTpressVolumeTemperature( E );
    /* eos_eval->V/V0 = 0.6; */
    /* eos_eval->T = 2500; */

    /* now update the struct with other material properties */
    SetRTpressRho( rtp, eos_eval ); // THIS WORKS (SI UNITS)

    SetRTpressAlpha( rtp, eos_eval ); // THIS WORKS (SI UNITS)
    /* eos_eval->alpha = 0.000028799762298150523 */ /* 1/K */

    /* FIXME: heat capacities are not scaled correctly - why? */
    SetRTpressCv( rtp, eos_eval );
    /* eos_eval->Cv = 0.056735422979028394 */ /* FIXME: units? */

    SetRTpressCp( rtp, eos_eval );
    /* eos_eval->Cp = 0.059432797887255438 */ /* FIXME: units? */

    /* FIXME: units */
    SetRTpressdTdPs( rtp, eos_eval );

    PetscFunctionReturn(0);
}

static PetscErrorCode set_rtpress_struct_non_dimensional( Ctx *E )
{
    ScalingConstants const     SC = E->parameters->scaling_constants;
    EosEval                    *eos_eval = &E->eos1_eval;

    PetscFunctionBeginUser;

    /* FIXME: for consistency should probably non-dimensionalise all
       entries in this struct, and not just a selection */

    eos_eval->rho /= SC->DENSITY;
    eos_eval->Cp /= SC->ENTROPY;
    eos_eval->Cv /= SC->ENTROPY;
    eos_eval->T /= SC->TEMP;
    eos_eval->alpha *= SC->TEMP;
    eos_eval->dTdPs /= SC->DTDP;

    PetscFunctionReturn(0);
}

PetscErrorCode set_rtpress_struct( PetscScalar P, PetscScalar S, Ctx *E )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = set_rtpress_struct_SI( P, S, E ); CHKERRQ(ierr);
    ierr = set_rtpress_struct_non_dimensional( E ); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
#endif

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

////////////// WIP

/* Base Class */

PetscErrorCode EOSCreate(EOS* p_eos, EOSType type)
{
  PetscErrorCode ierr;
  PetscBool      flg;
  EOS            eos;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1, p_eos);CHKERRQ(ierr);
  eos = *p_eos;

  ierr = PetscStrcmp(type, SPIDER_EOS_LOOKUP_PURE, &flg);CHKERRQ(ierr);
  if (flg) {
    ierr = EOSLookupPureCreate((EOSLookupPure*) eos->impl_data);CHKERRQ(ierr);
  }
  ierr = PetscStrcmp(type, SPIDER_EOS_RTPRESS_PURE, &flg);CHKERRQ(ierr);
  if (flg) {
    ierr = EOSRTpressPureCreate((EOSRTpressPure*) eos->impl_data);CHKERRQ(ierr);
  }
  ierr = PetscStrcmp(type, SPIDER_EOS_COMPOSITE, &flg);CHKERRQ(ierr);
  if (flg) {
    ierr = EOSCompositeCreate((EOSComposite*) eos->impl_data);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode EOSEval(const EOS eos, PetscScalar P , PetscScalar S, EosEval* eval)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  // TODO boilerplate eval from EosSetEval
  // TODO  call impl eval function (from EosSetEval*)
  PetscFunctionReturn(0);
}

PetscErrorCode EOSDestroy(EOS *p_eos)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscFree(p_eos);CHKERRQ(ierr);
  // TODO call impl destroy function
  *p_eos = NULL;
  PetscFunctionReturn(0);
}
