#include "eos.h"
#include "monitor.h"
#include "parameters.h"
#include "util.h"

/* lookup material properties (default) */
static PetscErrorCode Interp1dCreateAndSet( const char *, Interp1d *, PetscScalar, PetscScalar );
static PetscErrorCode Interp1dDestroy( Interp1d * );
static PetscErrorCode Interp2dCreateAndSet( const char *, Interp2d *, PetscScalar, PetscScalar, PetscScalar );
static PetscErrorCode Interp2dDestroy( Interp2d * );
static PetscErrorCode LookupSolidCreate( Lookup, const ScalingConstants );
static PetscErrorCode LookupMeltCreate( Lookup, const ScalingConstants );
static PetscErrorCode LookupLiquidusCreate( Lookup, const ScalingConstants );
static PetscErrorCode LookupSolidusCreate( Lookup, const ScalingConstants );

/* rtpress material properties (Wolf and Bower, 2018) */
static PetscErrorCode RTpressParametersCreate( RTpressParameters * );
static PetscErrorCode RTpressParametersCreateAndSet( RTpressParameters *, const FundamentalConstants );
static PetscErrorCode RTpressParametersDestroy( RTpressParameters * );
static PetscScalar get_rtpress_pressure( PetscScalar, PetscScalar, const RTpressParameters );
static PetscScalar get_rtpress_entropy( PetscScalar, PetscScalar,  const RTpressParameters );
static PetscErrorCode set_rtpress_struct_SI( PetscScalar, PetscScalar, Ctx * );
static PetscErrorCode set_rtpress_struct_non_dimensional( Ctx * );
static PetscErrorCode set_rtpress_density( RTpressParameters const, EosEval * );
static PetscErrorCode set_rtpress_thermal_expansion( RTpressParameters const, EosEval * );
static PetscErrorCode set_rtpress_heat_capacity_constant_volume( RTpressParameters const, EosEval * );
static PetscErrorCode set_rtpress_heat_capacity_constant_pressure( RTpressParameters const, EosEval * );
static PetscErrorCode set_rtpress_isentropic_gradient( RTpressParameters const, EosEval * );
/* solve for volume and temperature from pressure and entropy */
static PetscErrorCode objective_function_rtpress_volume_temperature( SNES, Vec, Vec, void * );
static PetscErrorCode solve_for_rtpress_volume_temperature( Ctx * );

/* helper functions */
static PetscScalar per_atom_to_specific( PetscScalar, PetscScalar, PetscScalar );
static PetscScalar specific_to_per_atom( PetscScalar, PetscScalar, PetscScalar );
static PetscScalar joule_to_eV( PetscScalar );
static PetscScalar eV_to_joule( PetscScalar );

/* set equation of state (EOS).  Currently for a melt and a solid phase, but
   can be extended for more phases in the future */

PetscErrorCode set_eos( Parameters P )
{
    PetscErrorCode ierr;
    ScalingConstants const SC = P->scaling_constants;
    EosParameters eosp1 = P->eos1_parameters; /* melt */
    EosParameters eosp2 = P->eos2_parameters; /* solid */

    PetscFunctionBeginUser;

    /* set solid eos */
    switch( P->SOLID_EOS ){
        case 1:
            /* lookup */
            eosp2->LOOKUP_FLAG = PETSC_TRUE;
            eosp2->RTPRESS_FLAG = PETSC_FALSE;
            ierr = LookupSolidCreate( eosp2->lookup, SC );CHKERRQ(ierr);
            break;
        default:
            SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported SOLID_EOS value %d provided",P->SOLID_EOS);
            break;
    }

    /* set melt eos */
    switch( P->MELT_EOS ){
        case 1:
            /* lookup */
            eosp1->LOOKUP_FLAG = PETSC_TRUE;
            eosp1->RTPRESS_FLAG = PETSC_FALSE;
            ierr = LookupMeltCreate( eosp1->lookup, SC );CHKERRQ(ierr);
            break;
        case 2:
            /* analytical RTpress */
            eosp1->LOOKUP_FLAG = PETSC_FALSE;
            eosp1->RTPRESS_FLAG = PETSC_TRUE;
            ierr = RTpressParametersCreateAndSet( &eosp1->rtpress_parameters, P->fundamental_constants );CHKERRQ(ierr);
            break;
        default:
            SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unsupported MELT_EOS value %d provided",P->MELT_EOS);
            break;
    }

    /* conductivity (W/m/K) */
    eosp2->cond = 4.0;
    ierr = PetscOptionsGetScalar(NULL,NULL,"-cond_sol",&eosp2->cond,NULL);CHKERRQ(ierr);
    eosp2->cond /= SC->COND;
    eosp1->cond = 4.0;
    ierr = PetscOptionsGetScalar(NULL,NULL,"-cond_mel",&eosp1->cond,NULL);CHKERRQ(ierr);
    eosp1->cond /= SC->COND;

    /* viscosity (Pa.s) */
    /* for solid, this is a prefactor if activation_energy_sol or activation_volume_sol are non-zero */
    eosp2->log10visc = 21.0;
    ierr = PetscOptionsGetScalar(NULL,NULL,"-log10visc_sol",&eosp2->log10visc,NULL);CHKERRQ(ierr);
    eosp2->log10visc -= SC->LOG10VISC;
    eosp1->log10visc = 2.0;
    ierr = PetscOptionsGetScalar(NULL,NULL,"-log10visc_mel",&eosp1->log10visc,NULL);CHKERRQ(ierr);
    eosp1->log10visc -= SC->LOG10VISC;

    /* set liquidus and solidus from lookup */
    /* TODO: could also have analytical functions for these as well by using
       a case structure as above */
    /* TODO: currently, melting curves are stored in both the solid and the
       melt Eos struct */
    ierr = LookupLiquidusCreate( eosp1->lookup, SC );CHKERRQ(ierr);
    ierr = LookupLiquidusCreate( eosp2->lookup, SC );CHKERRQ(ierr);
    ierr = LookupSolidusCreate( eosp1->lookup, SC );CHKERRQ(ierr);
    ierr = LookupSolidusCreate( eosp2->lookup, SC );CHKERRQ(ierr);

    PetscFunctionReturn(0);

}


/*
 ******************************************************************************
 * Lookup EOS from data files
 ******************************************************************************
*/

static PetscErrorCode Interp1dCreateAndSet( const char * filename, Interp1d *interp_ptr, PetscScalar xconst, PetscScalar yconst )
{
    PetscErrorCode ierr;
    Interp1d interp;
    FILE *fp;
    PetscInt i=0;
    char string[PETSC_MAX_PATH_LEN];
#if (defined PETSC_USE_REAL___FLOAT128)
    char xtemp[30], ytemp[30];
#endif
    PetscScalar x, y, xscale=0.0, yscale=0.0;
    PetscInt HEAD, NX;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(1,interp_ptr);CHKERRQ(ierr);
    interp = *interp_ptr;

    fp = fopen( filename, "r" );

    if (!fp) {
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_FILE_OPEN,"Could not open file %s",filename);
    }
    // fgets reads in string, sscanf processes it
    while(fgets(string, sizeof(string), fp) != NULL) {

        /* first header contains size of arrays to create */
        /* FIXME: if the header information is incorrect it can lead
           to various memory errors.  Add some safety checks? */
        if( i==0 ){
            /* remove # at start of line */
            memmove( string, string+1, strlen(string) );
            /* TODO: will break for 64 bit integers (but presumably
               elsewhere in the code will also break) */
            sscanf( string, "%d %d", &HEAD, &NX );
            /* make arrays based on the sizes in the header */
            ierr = PetscMalloc1( NX, &interp->xa ); CHKERRQ(ierr);
            ierr = PetscMalloc1( NX, &interp->ya ); CHKERRQ(ierr);
        }

        /* get column scalings from last line of header */
        else if( i==HEAD-1 ){
           /* remove # at start of line */
            memmove( string, string+1, strlen(string) );
#if (defined PETSC_USE_REAL___FLOAT128)
            sscanf( string, "%s %s", xtemp, ytemp );
            xscale = strtoflt128(xtemp, NULL);
            yscale = strtoflt128(ytemp, NULL);
#else
            sscanf(string, "%lf %lf", &xscale, &yscale );
#endif
        }
        else if( i>=HEAD ){
#if (defined PETSC_USE_REAL___FLOAT128)
            sscanf( string, "%s %s", xtemp, ytemp );
            x = strtoflt128(xtemp, NULL);
            y = strtoflt128(ytemp, NULL);
#else
            sscanf(string, "%lf %lf", &x, &y );
#endif
            interp->xa[i-HEAD] = x;
            interp->xa[i-HEAD] *= xscale;
            interp->xa[i-HEAD] /= xconst;
            interp->ya[i-HEAD] = y;
            interp->ya[i-HEAD] *= yscale;
            interp->ya[i-HEAD] /= yconst;
            }
        ++i;
    }

    fclose( fp );

    /* for debugging */
#if (defined DEBUGOUTPUT)
    PetscErrorCode ierr;
    PetscMPIInt rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    for (i=0; i<NX; i++ ){
        ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] %d %f\n", rank, i, (double) interp->xa[i]);CHKERRQ(ierr);
    }
    for (i=0; i<NX; i++ ){
        ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] %d %f\n", rank, i, (double) interp->ya[i]);CHKERRQ(ierr);
    }
#endif

    interp->xmin= interp->xa[0];
    interp->xmax= interp->xa[NX-1];
    /* ymin and ymax are not used at present, and it might be
       dangerous to do so since for the middle-out solidus and
       liquidus data the maximum entropy is not at the end of the
       array since the maximum occurs around mid-mantle pressures */
    interp->ymin= interp->ya[0];
    interp->ymax= interp->ya[NX-1];
    interp->NX = NX;

    PetscFunctionReturn(0);
}

static PetscErrorCode Interp2dCreateAndSet( const char * filename, Interp2d *interp_ptr, PetscScalar xconst, PetscScalar yconst, PetscScalar zconst )
{
    PetscErrorCode ierr;
    Interp2d interp;
    FILE *fp;
    PetscInt i=0, j=0, k=0;
    char string[PETSC_MAX_PATH_LEN];
#if (defined PETSC_USE_REAL___FLOAT128)
    char xtemp[30], ytemp[30], ztemp[30];
#endif
   /* without assigning values below, the compiler warns about possible
      uninitialised values for just the quadruple precision case */
    PetscScalar xscale=0.0, yscale=0.0, zscale=0.0;
    PetscScalar x, y, z;
    PetscInt HEAD, NX, NY; // replaced by entries in first header line
    PetscInt xind, yind;

    PetscFunctionBeginUser;

    // NEW - allocate memory for struct
    ierr = PetscMalloc1(1, interp_ptr);CHKERRQ(ierr);
    interp = *interp_ptr;

    fp = fopen( filename, "r" );

    if(fp==NULL) {
      SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_FILE_OPEN,"Could not open file %s",filename);
    }   

    // fgets reads in string, sscanf processes it
    while(fgets(string, sizeof(string), fp) != NULL) {

        /* first header contains size of arrays to create */
        /* FIXME: if the header information is incorrect it can lead
           to various memory errors.  Add some safety checks? */
        if( i==0 ){
            /* remove # at start of line */
            memmove( string, string+1, strlen(string) );
            /* TODO: will break for 64 bit integers (but presumably
               elsewhere in the code will also break) */
            sscanf( string, "%d %d %d", &HEAD, &NX, &NY );
            /* make arrays based on the sizes in the header */
            ierr = PetscMalloc1( NX, &interp->xa ); CHKERRQ(ierr);
            ierr = PetscMalloc1( NY, &interp->ya ); CHKERRQ(ierr);
            interp->za = Make2DPetscScalarArray( NX, NY );
        }   

        /* get column scalings from last line of header */
        else if( i==HEAD-1 ){
            /* remove # at start of line */
            memmove( string, string+1, strlen(string) );
#if (defined PETSC_USE_REAL___FLOAT128)
            sscanf( string, "%s %s %s", xtemp, ytemp, ztemp );
            xscale = strtoflt128(xtemp, NULL);
            yscale = strtoflt128(ytemp, NULL);
            zscale = strtoflt128(ztemp, NULL);
#else
            sscanf( string, "%lf %lf %lf", &xscale, &yscale, &zscale );
#endif
        }   
        else if( i>=HEAD ){
#if (defined PETSC_USE_REAL___FLOAT128)
            sscanf( string, "%s %s %s", xtemp, ytemp, ztemp );
            x = strtoflt128(xtemp, NULL);
            y = strtoflt128(ytemp, NULL);
            z = strtoflt128(ztemp, NULL);
#else
            sscanf(string, "%lf %lf %lf", &x, &y, &z );
#endif

            /* need to determine x and y indices to insert data in
               correct position in the za (2-D) array */
            xind = (i-HEAD) % NX; // e.g., repeats 0->2019, 0->2019, 0->2019, for each set
            yind = (i-HEAD) / NX; // integer division gives correct yind (I think)

            /* lookup value */
            interp->za[xind][yind] = z;
            interp->za[xind][yind] *= zscale;
            interp->za[xind][yind] /= zconst;

            /* x coordinate */
            if( i<HEAD+NX ){
                interp->xa[j] = x;
                interp->xa[j] *= xscale;
                interp->xa[j] /= xconst;
                ++j;
            }   
            /* y coordinate */
            if( (i-HEAD) % NX ==0 ){
                interp->ya[k] = y;
                interp->ya[k] *= yscale;
                interp->ya[k] /= yconst;
                ++k;
            }

        }
        ++i;
    }

    fclose( fp );

    /* for debugging */
#if (defined DEBUGOUTPUT)
    PetscErrorCode ierr;
    PetscMPIInt rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    for (i=0; i<NX; i++ ){
        ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] %d %f\n", rank, i, (double) interp->xa[i]);CHKERRQ(ierr);
    }
    for (j=0; j<NY; j++ ){
        ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] %d %f\n", rank, j, (double) interp->ya[j]);CHKERRQ(ierr);
    }
    for(i=0; i<NX; i++ ){
        for(j=0; j<NY; j++ ){
            ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] %d %d %f\n", rank, i, j, (double) interp->za[i][j]);CHKERRQ(ierr);
        }
    }
#endif

    interp->NX = NX;
    interp->NY = NY;
    interp->xmin= interp->xa[0];
    interp->xmax= interp->xa[NX-1];
    interp->ymin= interp->ya[0];
    interp->ymax= interp->ya[NY-1];
    /* if we store the x and y step, we can more quickly locate the
       relevant indices in the arrays by direct calculation, if the
       data has constant spacing */
    interp->dx = interp->xa[1]-interp->xa[0];
    interp->dy = interp->ya[1]-interp->ya[0];

    PetscFunctionReturn(0);
}

PetscScalar get_val1d( Interp1d const interp, PetscScalar x )
{   
    /* wrapper for evaluating a 1-D lookup
       linear interpolation with truncation for values
       that fall outside the data lookup range */
    
    //PetscErrorCode    ierr;
    PetscScalar       w1, result;
    PetscScalar const *xa, *ya;
    PetscScalar       xmin, xmax;
    PetscInt          ind, NX;
    
    NX = interp->NX;
    xa = interp->xa;
    xmin = interp->xmin;
    xmax = interp->xmax;
    ya = interp->ya;
    
    /* to reproduce the behaviour of scipy.interpolate.interp1d the
       code should produce a ValueError if interpolation is
       attempted on a value outside of the range of x (where
       extrapolation is necessary). Here we truncate instead. */
    
    if( x<xmin ){
      //ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val1d: x<xmin, %f<%f.  Truncating\n",(double)x,(double)xmin);CHKERRQ(ierr);
      ind = 0; // minimum index, max index is always +1
      x = xmin;
    }
    else if( x>xmax ){
      //ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val1d: x>xmax, %f>%f.  Truncating\n",(double)x,(double)xmax);CHKERRQ(ierr);
      ind = NX-2; // minimum index, max index is always +1
      x = xmax;
    }
    else{
      // loop to find minimum index
      /* trivial algorithm to find minimum index when x data
         is not evenly spaced */
      ind = 0;
      while( (x-xa[ind])>=0) {
        ind += 1;
      }
      /* loop exits when sign changes, meaning that previous index
         is the minimum index */
      ind -= 1;
    }
    
    // w1 is 0 at leftmost (minimum) x, 1 at rightmost (maximum) x
    w1 = (x-xa[ind]) / (xa[ind+1]-xa[ind]); // weighting
    
    result = ya[ind] * (1.0-w1) + ya[ind+1] * w1;
    
    return result;
}

PetscScalar get_val2d( Interp2d const interp, PetscScalar x, PetscScalar y )
{
    /* wrapper for evaluating a 2-D lookup using bilinear
       interpolation.

       Note that this assumes that x data (pressure) is evenly
       spaced so we use a faster lookup approach by computing
       indices directly rather than looping through data */

    //PetscErrorCode ierr;
    PetscScalar z1, z2, z3, z4;
    PetscScalar w1, w2, w3, w4; // weights
    PetscScalar result;
    PetscScalar const *xa, *ya;
    PetscInt NX, NY;
    PetscScalar dx;
    PetscScalar xmin, xmax, ymin, ymax;
    // below only if y data is evenly spaced
    //PetscScalar dy;
    PetscInt indx, indy;

    NX = interp->NX;
    xa = interp->xa;
    xmin = interp->xmin;
    xmax = interp->xmax;
    dx = interp->dx;
    NY = interp->NY;
    ya = interp->ya;
    ymin = interp->ymin;
    ymax = interp->ymax;
    // below only if y data is evenly spaced
    //dy = interp->dy;

    /* to reproduce the behaviour of scipy.interpolate.RectBivariateSpline
       the code should truncate if interpolation is attempted on a
       value outside of the lookup data range */

    /* for pressure (x), constant spacing assumed */
    if( x<xmin ){
      //ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val2d: x<xmin, %f<%f.  Truncating\n",(double)x,(double)xmin);CHKERRQ(ierr);
      indx = 0; // minimum index, max index is always +1
      x = xmin;
    }
    else if( x>xmax ){
      //ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val2d: x>xmax, %f>%f.  Truncating\n",(double)x,(double)xmax);CHKERRQ(ierr);
      indx = NX-2; // minimum index, max index is always +1
      x = xmax;
    }
    else{
      indx = PetscFloorReal( (x-xmin)/dx );  // minimum index
    }

    // x weights
    w1 = x-xa[indx]; // x-x1
    w2 = xa[indx+1]-x; // x2-x


    /* for entropy (y), irregular spacing assumed */
    if( y<ymin ){
      //ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val2d: y<ymin, %f<%f.  Truncating\n",(double)y,(double)ymin);CHKERRQ(ierr);
      indy = 0; // minimum index, max index is always +1
      y = ymin;
    }
    else if( y>ymax ){
      //ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val2d: y>ymax, %f>%f.  Truncating\n",(double)y,(double)ymax);CHKERRQ(ierr);
      indy = NY-2; // minimum index, max index is always +1
      y = ymax;
    }
    else{
      // loop to find minimum index
      /* trivial algorithm to find minimum index when y data
         is not evenly spaced */
      indy = 0;
      while( (y-ya[indy])>=0) {
        indy += 1;
      }
      /* loop exits when sign changes, meaning that previous index
         is the minimum index */
      indy -= 1;
    }
    // y weights
    w3 = y-ya[indy]; // y-y1
    w4 = ya[indy+1]-y; // y2-y

    // min P (x), min S (y)
    z1 = interp->za[indx][indy];
    // max P (x), min S (y)
    z2 = interp->za[indx+1][indy];
    // min P (x), max S (y)
    z3 = interp->za[indx][indy+1];
    // max P (x), max S (y)
    z4 = interp->za[indx+1][indy+1];

    // bilinear interpolation
    result = z1 * w2 * w4;
    result += z2 * w1 * w4;
    result += z3 * w2 * w3;
    result += z4 * w1 * w3;
    result /= dx; // dx
    result /= ya[indy+1]-ya[indy]; // dy

    return result;
}

static PetscErrorCode Interp1dDestroy( Interp1d *interp_ptr )
{
    PetscErrorCode ierr;
    Interp1d interp = *interp_ptr;

    PetscFunctionBeginUser;

    ierr = PetscFree( interp->xa ); CHKERRQ(ierr);
    ierr = PetscFree( interp->ya ); CHKERRQ(ierr);

    ierr = PetscFree(*interp_ptr);CHKERRQ(ierr);
    *interp_ptr = NULL;

    PetscFunctionReturn(0);

}
                                                             
static PetscErrorCode Interp2dDestroy( Interp2d *interp_ptr )
{

    PetscErrorCode ierr;
    Interp2d interp = *interp_ptr;

    PetscFunctionBeginUser;

    PetscInt i;
    PetscInt NX = interp->NX;

    for(i=0; i<NX; i++) {
        free( interp->za[i] );
    }
    free( interp->za );
    PetscFree( interp->xa );
    PetscFree( interp->ya );

    ierr = PetscFree(*interp_ptr);CHKERRQ(ierr);
    *interp_ptr = NULL;

    PetscFunctionReturn(0);

}

/* Helper routine to prepend the root directory to a relative path */
/* https://gcc.gnu.org/onlinedocs/gcc-4.9.0/cpp/Stringification.html */
#define STRINGIFY(x) STRINGIFY2(x)
#define STRINGIFY2(x) #x
#define SPIDER_ROOT_DIR_STR STRINGIFY(SPIDER_ROOT_DIR)
static PetscErrorCode MakeRelativeToSourcePathAbsolute(char* path) {
  PetscErrorCode ierr;
  char tmp[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  ierr = PetscStrcpy(tmp,path);CHKERRQ(ierr);
  ierr = PetscStrcpy(path,SPIDER_ROOT_DIR_STR);CHKERRQ(ierr);
  ierr = PetscStrcat(path,"/");CHKERRQ(ierr); /* not portable */
  ierr = PetscStrcat(path,tmp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#undef SPIDER_ROOT_DIR_STR


static PetscErrorCode LookupMeltCreate( Lookup lookup, ScalingConstants const SC )
{
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;

  /* Based on input options, determine which files to load.
     Options ending with _rel_to_src indicate a path relative
     to the source code. In this case we prepend a string, SPIDER_ROOT_DIR_STR,
     and /. The corresponding option without this overrides. */

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->alpha_filename,"lookup_data/RTmelt/thermal_exp_melt.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-alphaMel_filename_rel_to_src",lookup->alpha_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->alpha_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-alphaMel_filename",lookup->alpha_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -alphaMel_filename_rel_to_src ignored because -alphaMel_filename provided\n");CHKERRQ(ierr);
    }
  }

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->cp_filename,"lookup_data/RTmelt/heat_capacity_melt.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-cpMel_filename_rel_to_src",lookup->cp_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->cp_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-cpMel_filename",lookup->cp_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -cpMel_filename_rel_to_src ignored because -cpMel_filename provided\n");CHKERRQ(ierr);
    }
  }

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->dTdPs_filename,"lookup_data/RTmelt/adiabat_temp_grad_melt.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-dtdpsMel_filename_rel_to_src",lookup->dTdPs_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->dTdPs_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-dtdpsMel_filename",lookup->dTdPs_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -dtdpsMel_filename_rel_to_src ignored because -dtdpsMel_filename provided\n");CHKERRQ(ierr);
    }
  }

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->rho_filename,"lookup_data/RTmelt/density_melt.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-rhoMel_filename_rel_to_src",lookup->rho_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->rho_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-rhoMel_filename",lookup->rho_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -rhoMel_filename_rel_to_src ignored because -rhoMel_filename provided\n");CHKERRQ(ierr);
    }
  }

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->temp_filename,"lookup_data/RTmelt/temperature_melt.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-tempMel_filename_rel_to_src",lookup->temp_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->temp_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-tempMel_filename",lookup->temp_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -tempMel_filename_rel_to_src ignored because -tempMel_filename provided\n");CHKERRQ(ierr);
    }
  }

  /* melt lookups */
  /* 2d */
  ierr = Interp2dCreateAndSet( lookup->alpha_filename, &lookup->alpha, SC->PRESSURE, SC->ENTROPY, 1.0/SC->TEMP );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( lookup->cp_filename, &lookup->cp, SC->PRESSURE, SC->ENTROPY, SC->ENTROPY );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( lookup->dTdPs_filename, &lookup->dTdPs, SC->PRESSURE, SC->ENTROPY, SC->DTDP );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( lookup->rho_filename, &lookup->rho, SC->PRESSURE, SC->ENTROPY, SC->DENSITY );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( lookup->temp_filename, &lookup->temp, SC->PRESSURE, SC->ENTROPY, SC->TEMP );CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

static PetscErrorCode LookupSolidCreate( Lookup lookup, const ScalingConstants SC )
{
  PetscErrorCode  ierr;

  PetscFunctionBeginUser;

  /* Based on input options, determine which files to load.
     Options ending with _rel_to_src indicate a path relative
     to the source code. In this case we prepend a string, SPIDER_ROOT_DIR_STR,
     and /. The corresponding option without this overrides. */

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->alpha_filename,"lookup_data/RTmelt/thermal_exp_solid.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-alphaSol_filename_rel_to_src",lookup->alpha_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->alpha_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-alphaSol_filename",lookup->alpha_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -alphaSol_filename_rel_to_src ignored because -alphaSol_filename provided\n");CHKERRQ(ierr);
    }
  }

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->cp_filename,"lookup_data/RTmelt/heat_capacity_solid.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-cpSol_filename_rel_to_src",lookup->cp_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->cp_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-cpSol_filename",lookup->cp_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -cpSol_filename_rel_to_src ignored because -cpSol_filename provided\n");CHKERRQ(ierr);
    }
  }

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->dTdPs_filename,"lookup_data/RTmelt/adiabat_temp_grad_solid.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-dtdpsSol_filename_rel_to_src",lookup->dTdPs_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->dTdPs_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-dtdpsSol_filename",lookup->dTdPs_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -dtdpsSol_filename_rel_to_src ignored because -dtdpsSol_filename provided\n");CHKERRQ(ierr);
    }
  }

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->rho_filename,"lookup_data/RTmelt/density_solid.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-rhoSol_filename_rel_to_src",lookup->rho_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->rho_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-rhoSol_filename",lookup->rho_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -rhoSol_filename_rel_to_src ignored because -rhoSol_filename provided\n");CHKERRQ(ierr);
    }
  }

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->temp_filename,"lookup_data/RTmelt/temperature_solid.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-tempSol_filename_rel_to_src",lookup->temp_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->temp_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-tempSol_filename",lookup->temp_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -tempSol_filename_rel_to_src ignored because -tempSol_filename provided\n");CHKERRQ(ierr);
    }
  }

  ierr = Interp2dCreateAndSet( lookup->alpha_filename, &lookup->alpha, SC->PRESSURE, SC->ENTROPY, 1.0/SC->TEMP );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( lookup->cp_filename, &lookup->cp, SC->PRESSURE, SC->ENTROPY, SC->ENTROPY );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( lookup->dTdPs_filename, &lookup->dTdPs, SC->PRESSURE, SC->ENTROPY, SC->DTDP );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( lookup->rho_filename, &lookup->rho, SC->PRESSURE, SC->ENTROPY, SC->DENSITY );CHKERRQ(ierr);
  ierr = Interp2dCreateAndSet( lookup->temp_filename, &lookup->temp, SC->PRESSURE, SC->ENTROPY, SC->TEMP );CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

static PetscErrorCode LookupLiquidusCreate( Lookup lookup, ScalingConstants const SC )
{
  /* TODO: this is not ideal, since the liquidus lookup is stored to both eos structs */

  PetscErrorCode  ierr;

  PetscFunctionBeginUser;

  /* Based on input options, determine which files to load.
     Options ending with _rel_to_src indicate a path relative
     to the source code. In this case we prepend a string, SPIDER_ROOT_DIR_STR,
     and /. The corresponding option without this overrides. */

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->liquidus_filename,"lookup_data/RTmelt/liquidus_andrault2011.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-liquidus_filename_rel_to_src",lookup->liquidus_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->liquidus_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-liquidus_filename",lookup->liquidus_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -liquidus_filename_rel_to_src ignored because -liquidus_filename provided\n");CHKERRQ(ierr);
    }      
  }

  ierr = Interp1dCreateAndSet( lookup->liquidus_filename, &lookup->liquidus, SC->PRESSURE, SC->ENTROPY );CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

static PetscErrorCode LookupSolidusCreate( Lookup lookup, ScalingConstants const SC )
{
  /* TODO: this is not ideal, since the liquidus lookup is stored to both eos structs */

  PetscErrorCode  ierr;

  PetscFunctionBeginUser;

  /* Based on input options, determine which files to load.
     Options ending with _rel_to_src indicate a path relative
     to the source code. In this case we prepend a string, SPIDER_ROOT_DIR_STR,
     and /. The corresponding option without this overrides. */

  {
    PetscBool set_rel_to_src,set;
    ierr = PetscStrcpy(lookup->solidus_filename,"lookup_data/RTmelt/solidus_andrault2011.dat");CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-solidus_filename_rel_to_src",lookup->solidus_filename,PETSC_MAX_PATH_LEN,&set_rel_to_src);
    ierr = MakeRelativeToSourcePathAbsolute(lookup->solidus_filename);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-solidus_filename",lookup->solidus_filename,PETSC_MAX_PATH_LEN,&set);
    if (set && set_rel_to_src) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -solidus_filename_rel_to_src ignored because -solidus_filename provided\n");CHKERRQ(ierr);
    }
  }

  ierr = Interp1dCreateAndSet( lookup->solidus_filename, &lookup->solidus, SC->PRESSURE, SC->ENTROPY );CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

/*
 ******************************************************************************
 * RTpress model (analytical) from Wolf and Bower (2018)
 ******************************************************************************
*/

PetscErrorCode RTpressParametersCreateAndSet( RTpressParameters* rtpress_parameters_ptr, const FundamentalConstants FC )
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

PetscScalar get_rtpress_pressure_test( Ctx *E )
{
    /* test function (lldb, then step through to evaluate P for
       the V, T conditions defined below.  Compare with pressure
       plot in Jupyter notebook*/

    /* seems to confirm P is correct
       FIXME: need to truncate negative values?  Perhaps not for
       smoothness of solver during inversion? */

    RTpressParameters rtp = E->parameters->eos1_parameters->rtpress_parameters;
    PetscScalar   P, T, V, volfrac;

    volfrac = 0.6;
    T = 2500.0;
    V = volfrac * rtp->V0;
    P = get_rtpress_pressure( V, T, rtp ); /* 26.964412351908095 */

    volfrac = 0.8;
    T = 2500.0;
    V = volfrac * rtp->V0;
    P = get_rtpress_pressure( V, T, rtp ); /* 3.5394262799831875 */

    volfrac = 1.0;
    T = 2500.0;
    V = volfrac * rtp->V0;
    /* FIXME: note negative value */
    P = get_rtpress_pressure( V, T, rtp ); /* -0.54571824296751803 */

    volfrac = 0.6;
    T = 5000.0;
    V = volfrac * rtp->V0;
    P = get_rtpress_pressure( V, T, rtp ); /* 38.568776519196767 */

    volfrac = 0.8;
    T = 5000.0;
    V = volfrac * rtp->V0;
    P = get_rtpress_pressure( V, T, rtp ); /* 10.39877886001743 */

    volfrac = 1.0;
    T = 5000.0;
    V = volfrac * rtp->V0;
    P = get_rtpress_pressure( V, T, rtp ); /* 2.107013748263614 */

    return P;

}


static PetscScalar get_rtpress_pressure( PetscScalar V, PetscScalar T, RTpressParameters const rtp )
{
    /* pressure = function( volume, temperature ) */

    PetscScalar P;
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

    /* FIXME: if I try using PetscCbrtReal then I get the following warning at compilation:
           warning: implicit declaration of function 'cbrt' is invalid in C99 [-Wimplicit-function-declaration] */
    /* TODO: should I use Real or PetscScalar functions?  i.e. PetscPowReal or PetscPowScalar? */

    P =  -9*K0*V0*(-1.0/3.0*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpReal((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/V - 1.0/3.0*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpReal((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/V)/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)));

    return P;

}

PetscScalar get_rtpress_entropy_test( Ctx *E )
{
    /* test function (lldb, then step through to evaluate S for
       the V, T conditions defined below.  Compare with entropy
       plot in Jupyter notebook*/

    /* seems to confirm S is correct */

    RTpressParameters rtp = E->parameters->eos1_parameters->rtpress_parameters;
    PetscScalar   S, T, V, volfrac;

    volfrac = 0.6;
    T = 2500.0;
    V = volfrac * rtp->V0;
    S = get_rtpress_entropy( V, T, rtp ); /* -0.025934761412030771 */
    S /= rtp->KBOLTZ; /* -1.8784471379548613 */

    volfrac = 0.8;
    T = 2500.0;
    V = volfrac * rtp->V0;
    S = get_rtpress_entropy( V, T, rtp ); /* -0.016066633649871438 */
    S /= rtp->KBOLTZ; /* -1.1637015477678558 */

    volfrac = 1.0;
    T = 2500.0;
    V = volfrac * rtp->V0;
    S = get_rtpress_entropy( V, T, rtp ); /* -0.010551508142269357 */
    S /= rtp->KBOLTZ; /* -0.76424263003857285 */

    volfrac = 0.6;
    T = 5000.0;
    V = volfrac * rtp->V0;
    S = get_rtpress_entropy( V, T, rtp ); /* 0.010228437833469547 */
    S /= rtp->KBOLTZ; /* 0.74084274263333616 */

    volfrac = 0.8;
    T = 5000.0;
    V = volfrac * rtp->V0;
    S = get_rtpress_entropy( V, T, rtp ); /* 0.021941049097879431 */
    S /= rtp->KBOLTZ; /* 1.5891837301622371 */

    volfrac = 1.0;
    T = 5000.0;
    V = volfrac * rtp->V0;
    S = get_rtpress_entropy( V, T, rtp ); /* 0.027130677516848993 */
    S /= rtp->KBOLTZ; /* 1.9650669895370627 */

    return S;

}

static PetscScalar get_rtpress_entropy( PetscScalar V, PetscScalar T, RTpressParameters const rtp )
{
    /* entropy = function( volume, temperature ) */

    PetscScalar S;
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

    S = S0 + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T)) + 0.041419469658752747*m/(T*(2*m - 2)) - 0.041419469658752747/(T*(2*m - 2))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*m*PetscLogScalar(T)/(2*m - 2) - 0.041419469658752747*m*PetscLogScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))/(2*m - 2) - 0.020709734829376374 - 0.041419469658752747*PetscLogScalar(T)/(2*m - 2) + 0.041419469658752747*PetscLogScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))/(2*m - 2) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T;

    return S;

}

static PetscErrorCode set_rtpress_density( RTpressParameters const rtp, EosEval *eos_eval )
{

    /* returns density in SI units, kg/m^3 */

    PetscScalar const mavg = rtp->mavg;
    PetscScalar const V = eos_eval->V;

    PetscFunctionBeginUser;

    eos_eval->rho = mavg / V; /* kg/Ang^3/mol */

    /* convert to SI */
    eos_eval->rho *= PetscPowScalar( 10.0, 30.0 ); /* Ang^3/m^3 */
    eos_eval->rho /= *rtp->AVOGADRO_ptr;

    PetscFunctionReturn(0);

}


static PetscErrorCode set_rtpress_thermal_expansion( RTpressParameters const rtp, EosEval *eos_eval )
{
    /* thermal expansion = function( volume, temperature ) */
    /* returns thermal expansion coefficient in SI units, 1/K */

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

    PetscScalar const T = eos_eval->T;
    PetscScalar const V = eos_eval->V;

    PetscFunctionBeginUser;

    eos_eval->alpha =  -(T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/T)/(V*(-9*K0*V0*((1.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2))/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)))));

    PetscFunctionReturn(0);

}

static PetscErrorCode set_rtpress_heat_capacity_constant_volume( RTpressParameters const rtp, EosEval *eos_eval )
{
    /* thermal heat capacity at constant volume = function( volume, temperature ) */

    PetscScalar const V0 = rtp->V0;
    PetscScalar const T0 = rtp->T0;
    PetscScalar const m = rtp->m;
    PetscScalar const b0 = rtp->b0;
    PetscScalar const b1 = rtp->b1;
    PetscScalar const b2 = rtp->b2;
    PetscScalar const b3 = rtp->b3;
    PetscScalar const b4 = rtp->b4;

    PetscScalar const T = eos_eval->T;
    PetscScalar const V = eos_eval->V;

    PetscFunctionBeginUser;

    eos_eval->Cv = 0.020709734829376374 + 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T;

    /* FIXME: possibly wrong */
    /* convert to SI */
    eos_eval->Cv = per_atom_to_specific( eos_eval->Cv, rtp->mavg, *rtp->AVOGADRO_ptr );
    eos_eval->Cv = eV_to_joule( eos_eval->Cv );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_rtpress_heat_capacity_constant_pressure( RTpressParameters const rtp, EosEval *eos_eval )
{
    /* thermal heat capacity at constant pressure = function( volume, temperature ) */

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

    PetscScalar const T = eos_eval->T;
    PetscScalar const V = eos_eval->V;

    PetscFunctionBeginUser;

    eos_eval->Cp = (0.020709734829376374 + 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T)*(-T*PetscPowScalar(T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/T, 2)/((0.020709734829376374 + 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T)*(-9*K0*V0*((1.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2))/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))))) + 1);

    /* FIXME: possibly wrong */
    /* convert to SI */
    eos_eval->Cp = per_atom_to_specific( eos_eval->Cp, rtp->mavg, *rtp->AVOGADRO_ptr );
    eos_eval->Cp = eV_to_joule( eos_eval->Cp );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_rtpress_isentropic_gradient( RTpressParameters const rtp, EosEval *eos_eval )
{

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

    PetscScalar const T = eos_eval->T;
    PetscScalar const V = eos_eval->V;

    PetscFunctionBeginUser;

    eos_eval->dTdPs =  -T*(T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/T)/((0.020709734829376374 + 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T)*(-T*PetscPowScalar(T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/T, 2)/((0.020709734829376374 + 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T)*(-9*K0*V0*((1.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2))/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))))) + 1)*(-9*K0*V0*((1.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*PetscPowScalar(V/V0, 2.0/3.0)*PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2) + (2.0/9.0)*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpScalar((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/PetscPowScalar(V, 2))/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*PetscPowScalar(m, 2)*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*PetscPowScalar(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V, 2)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*m*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(2*b2/PetscPowScalar(V0, 2) + 6*b3*(V/V0 - 1)/PetscPowScalar(V0, 2) + 12*b4*PetscPowScalar(V/V0 - 1, 2)/PetscPowScalar(V0, 2))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.055225959545003665*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) - 2*T0*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/6.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 3.0/2.0) + 2*T0*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-4*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 4*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 3) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) - 0.041419469658752747*m*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(2*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V + (1.0/3.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*PetscPowScalar(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1, 2)) + 0.041419469658752747*((5.0/3.0)*gamma0*PetscPowScalar(V0/V, 2.0/3.0)/PetscPowScalar(V, 2) + (1.0/18.0)*PetscPowScalar(V0/V, 4.0/3.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2) + (5.0/18.0)*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/PetscPowScalar(V, 2))/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)))));

    PetscFunctionReturn(0);

}

static PetscErrorCode solve_for_rtpress_volume_temperature( Ctx *E )
{
    PetscErrorCode ierr;
    SNES           snes;
    Vec            x,r;
    PetscScalar    *xx;
    PetscInt       i;
    RTpressParameters  rtp = E->parameters->eos1_parameters->rtpress_parameters;
    EosEval        *eos_eval = &E->eos1_eval;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"solve_for_rtpress_volume_temperature()\n");CHKERRQ(ierr);

    ierr = SNESCreate( PETSC_COMM_WORLD, &snes );CHKERRQ(ierr);

    /* Use this to address this specific SNES (nonlinear solver) from the command
       line or options file, e.g. -atmosic_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"rtpress_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    ierr = VecSetSizes( x, PETSC_DECIDE, 2 );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    ierr = SNESSetFunction(snes,r,objective_function_rtpress_volume_temperature,E);CHKERRQ(ierr);

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
    eos_eval->V = xx[0];
    eos_eval->T = xx[1];

    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode objective_function_rtpress_volume_temperature( SNES snes, Vec x, Vec f, void *ptr)
{
    PetscErrorCode             ierr;
    const PetscScalar          *xx;
    PetscScalar                *ff;
    PetscScalar                V, T, P, S;
    Ctx                        *E = (Ctx*) ptr;
    RTpressParameters          rtp = E->parameters->eos1_parameters->rtpress_parameters;
    EosEval                    *eos_eval = &E->eos1_eval;
    PetscScalar Ptarget = eos_eval->P;
    PetscScalar Starget = eos_eval->S;

    PetscFunctionBeginUser;

    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

    V = xx[0];
    T = xx[1];

    P = get_rtpress_pressure( V, T, rtp );
    S = get_rtpress_entropy( V, T, rtp );

    /* compute residual */
    ff[0] = P - Ptarget;
    ff[1] = S - Starget;

    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_rtpress_struct_SI( PetscScalar P, PetscScalar S, Ctx *E )
{
    /* solve for volume and temperature from pressure and entropy
       once, to avoid unnecessary computations.  This updates all the 
       eos_eval struct with the material properties */

    ScalingConstants const     SC = E->parameters->scaling_constants;
    RTpressParameters          rtp = E->parameters->eos1_parameters->rtpress_parameters;
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
    solve_for_rtpress_volume_temperature( E );
    /* eos_eval->V/V0 = 0.6; */
    /* eos_eval->T = 2500; */

    /* now update the struct with other material properties */
    set_rtpress_density( rtp, eos_eval ); // THIS WORKS (SI UNITS)

    set_rtpress_thermal_expansion( rtp, eos_eval ); // THIS WORKS (SI UNITS)
    /* eos_eval->alpha = 0.000028799762298150523 */ /* 1/K */

    /* FIXME: heat capacities are not scaled correctly - why? */
    set_rtpress_heat_capacity_constant_volume( rtp, eos_eval );
    /* eos_eval->Cv = 0.056735422979028394 */ /* FIXME: units? */

    set_rtpress_heat_capacity_constant_pressure( rtp, eos_eval );
    /* eos_eval->Cp = 0.059432797887255438 */ /* FIXME: units? */

    /* FIXME: units */
    set_rtpress_isentropic_gradient( rtp, eos_eval );

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

/*
 ******************************************************************************
 * Create and Destroy functions for parameter structs
 ******************************************************************************
 */

static PetscErrorCode LookupCreate( Lookup* lookup_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(1,lookup_ptr);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode LookupDestroy( Lookup* lookup_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = PetscFree(*lookup_ptr);CHKERRQ(ierr);
    *lookup_ptr = NULL;

    PetscFunctionReturn(0);

}

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

PetscErrorCode EosParametersCreate( EosParameters* eos_parameters_ptr )
{
    PetscErrorCode ierr;
    EosParameters eos_parameters;
    Lookup * lookup_ptr;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(1,eos_parameters_ptr);CHKERRQ(ierr);
    eos_parameters = *eos_parameters_ptr;
    lookup_ptr = &eos_parameters->lookup;
    ierr = LookupCreate( lookup_ptr );

    /* the Interp1d and Interp2d structs within Lookup are created and set
       later by (for example) Interp1dCreateAndSet() */

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
    lookup_ptr = &eos_parameters->lookup;
    lookup = *lookup_ptr;

    /* melting curves are (currently) always read in and stored as a lookup */
    ierr = Interp1dDestroy( &lookup->liquidus ); CHKERRQ(ierr);
    ierr = Interp1dDestroy( &lookup->solidus ); CHKERRQ(ierr);

    if( eos_parameters->RTPRESS_FLAG ){
        ierr = RTpressParametersDestroy( &eos_parameters->rtpress_parameters );CHKERRQ(ierr);
    }

    if( eos_parameters->LOOKUP_FLAG ){
        ierr = Interp2dDestroy( &lookup->alpha ); CHKERRQ(ierr);
        ierr = Interp2dDestroy( &lookup->cp ); CHKERRQ(ierr);
        ierr = Interp2dDestroy( &lookup->dTdPs ); CHKERRQ(ierr);
        ierr = Interp2dDestroy( &lookup->rho ); CHKERRQ(ierr);
        ierr = Interp2dDestroy( &lookup->temp ); CHKERRQ(ierr);
    }

    ierr = LookupDestroy( lookup_ptr );CHKERRQ(ierr);
    ierr = PetscFree(*eos_parameters_ptr);CHKERRQ(ierr);
    *eos_parameters_ptr = NULL;

    PetscFunctionReturn(0);

}
