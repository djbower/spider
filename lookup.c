#include "lookup.h"

static PetscErrorCode set_interp2d( const char *, Interp2d * );
static PetscErrorCode set_interp1d( const char *, Interp1d *, PetscInt );

PetscErrorCode set_lookups( Ctx *E ) 
{
    /* set all 1-D and 2-D lookups */

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    {
      PetscErrorCode ierr;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"set_lookup:\n" );CHKERRQ(ierr);
    }
#endif

    /* solid lookups */
    /* 2d */
    set_interp2d( ALPHA_SOL, &E->solid_prop.alpha );
    set_interp2d( CP_SOL, &E->solid_prop.cp );
    set_interp2d( DTDPS_SOL, &E->solid_prop.dTdPs );
    set_interp2d( RHO_SOL, &E->solid_prop.rho );
    set_interp2d( TEMP_SOL, &E->solid_prop.temp );
    /* const */
    E->solid_prop.cond = COND_SOL;
    E->solid_prop.log10visc = LOG10VISC_SOL;

    /* melt lookups */
    /* 2d */
    set_interp2d( ALPHA_MEL, &E->melt_prop.alpha );
    set_interp2d( CP_MEL, &E->melt_prop.cp );
    set_interp2d( DTDPS_MEL, &E->melt_prop.dTdPs );
    set_interp2d( RHO_MEL, &E->melt_prop.rho );
    set_interp2d( TEMP_MEL, &E->melt_prop.temp );
    /* const */
    E->melt_prop.cond = COND_MEL;
    E->melt_prop.log10visc = LOG10VISC_MEL;

    /* liquidus and solidus */
    /* 1d */
    set_interp1d( LIQUIDUS, &E->solid_prop.liquidus, NLS );
    set_interp1d( LIQUIDUS, &E->melt_prop.liquidus, NLS );
    /* duplication here, but want to remain flexible for future
       approaches for a multicomponent system */
    set_interp1d( SOLIDUS, &E->solid_prop.solidus, NLS );
    set_interp1d( SOLIDUS, &E->melt_prop.solidus, NLS );

    PetscFunctionReturn(0);
}

static PetscErrorCode set_interp2d( const char * filename, Interp2d *interp )
{
    FILE *fp;
    PetscInt i=0, j=0, k=0;
    char string[100];
#if (defined PETSC_USE_REAL___FLOAT128)
    char xtemp[30], ytemp[30], ztemp[30];
#endif
    PetscScalar xscale, yscale, zscale;
    PetscScalar x, y, z;
    PetscScalar xa[NX], ya[NY], za[NX*NY];

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    {
      PetscErrorCode ierr;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"set_interp2d:\n");CHKERRQ(ierr);
    }
#endif

    fp = fopen( filename, "r" );

    if(fp==NULL) {
        perror("Error opening file.\n");
        exit(-1);
    }

    // fgets reads in string, sscanf processes it
    while(fgets(string, sizeof(string), fp) != NULL) {
        /* get column scalings from last line of header */
        if( i==HEAD-1 ){
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
        if( i>=HEAD ){
#if (defined PETSC_USE_REAL___FLOAT128)
            sscanf( string, "%s %s %s", xtemp, ytemp, ztemp );
            x = strtoflt128(xtemp, NULL);
            y = strtoflt128(ytemp, NULL);
            z = strtoflt128(ztemp, NULL);
#else
            sscanf(string, "%lf %lf %lf", &x, &y, &z );
#endif
            /* lookup value */
            za[i-HEAD] = z * zscale;
            /* x coordinate */
            if( i<HEAD+NX ){
                xa[j] = x * xscale;
                ++j;
            }
            /* y coordinate */
            if( (i-HEAD) % NX ==0 ){
                ya[k] = y * yscale;
                ++k;
            }

        }
        ++i;
    }

    fclose( fp );

    /* for debugging */
#if (defined DEBUGOUTPUT)
    PetscMPIInt rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    for (i=0; i<NX; i++ ){
        ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] %d %f\n", rank, i, (double) xa[i]);CHKERRQ(ierr);
    }
    for (j=0; j<NY; j++ ){
        ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] %d %f\n", rank, j, (double) ya[j]);
    }
#endif

    interp->xmin= xa[0];
    interp->xmax= xa[NX-1];
    interp->ymin= ya[0];
    interp->ymax= ya[NY-1];

    memmove( interp->xa, xa, sizeof interp->xa );
    memmove( interp->ya, ya, sizeof interp->ya );
    memmove( interp->za, za, sizeof interp->za );

    /* if we store the x and y step, we can more quickly locate the
       relevant indices in the arrays by direct calculation, if the
       data has constant spacing */
    interp->dx = xa[1]-xa[0];
    interp->dy = ya[1]-ya[0];

    PetscFunctionReturn(0);
}

static PetscErrorCode set_interp1d( const char * filename, Interp1d *interp, PetscInt n )
{
    FILE *fp;
    PetscInt i=0;
    char string[100];
#if (defined PETSC_USE_REAL___FLOAT128)
    char xtemp[30], ytemp[30];
#endif
    PetscScalar x, y, xscale, yscale, xa[n], ya[n];

    PetscFunctionBeginUser;

#if (defined VERBOSE)
    {
      PetscErrorCode ierr;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"set_interp1d:\n");CHKERRQ(ierr);
    }
#endif

    fp = fopen( filename, "r" );

    if (!fp) {
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_OPEN,"Could not open file");
    }
    // fgets reads in string, sscanf processes it
    while(fgets(string, sizeof(string), fp) != NULL) {
        /* get column scalings from last line of header */
        if( i==HEAD-1 ){
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
        if( i>=HEAD ){
#if (defined PETSC_USE_REAL___FLOAT128)
            sscanf( string, "%s %s", xtemp, ytemp );
            x = strtoflt128(xtemp, NULL);
            y = strtoflt128(ytemp, NULL);
#else
            sscanf(string, "%lf %lf", &x, &y );
#endif
            xa[i-HEAD] = x * xscale;
            ya[i-HEAD] = y * yscale;
            }
        ++i;
    }

    fclose( fp );

    memmove( interp->xa, xa, sizeof interp->xa );
    interp->xmin= xa[0];
    interp->xmax= xa[n-1];
    memmove( interp->ya, ya, sizeof interp->ya );
    /* ymin and ymax are not used at present, and it might be
       dangerous to do so since for the middle-out solidus and
       liquidus data the maximum entropy is not at the end of the
       array since the maximum occurs around mid-mantle pressures */
    interp->ymin= ya[0];
    interp->ymax= ya[n-1];
    interp->n = n;

    PetscFunctionReturn(0);
}

PetscScalar get_val1d( Interp1d *I, PetscScalar x )
{
    /* wrapper for evaluating a 1-D lookup
       linear interpolation with truncation for values
       that fall outside the data lookup range */

    PetscScalar w1, result;
    PetscScalar *xa, *ya, xmax, xmin;
    PetscInt ind, n;

    xa = I->xa;
    xmax = I->xmax;
    xmin = I->xmin;
    ya = I->ya;
    n = I->n;

    /* to reproduce the behaviour of scipy.interpolate.interp1d the
       code should produce a ValueError if interpolation is
       attempted on a value outside of the range of x (where
       extrapolation is necessary). Here we truncate instead. */

    if( x<=xmin ){
#if (defined VERBOSE)
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val1d: x<xmin.  Truncating\n");CHKERRQ(ierr);
#endif
      ind = 0; // minimum index, max index is always +1
      x = xmin;
    }
    else if( x>=xmax ){
#if (defined VERBOSE)
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val1d: x>xmax.  Truncating\n");CHKERRQ(ierr);
#endif
      ind = n-2; // minimum index, max index is always +1
      x = xmax;
    }
    else{
      // loop to find minimum index
      /* trivial algorithm to find minimum index when x data
         is not evenly spaced */
      ind = 0;
      while( (x-xa[ind])>0) {
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

PetscScalar get_val2d( Interp2d *I, PetscScalar x, PetscScalar y )
{
    /* wrapper for evaluating a 2-D lookup using bilinear 
       interpolation.

       Note that this assumes that x data (pressure) is evenly
       spaced so we use a faster lookup approach by computing
       indices directly rather than looping through data */

#if (defined VERBOSE)
    PetscErrorCode ierr;
#endif
    PetscScalar z1, z2, z3, z4;
    PetscScalar w1, w2, w3, w4; // weights
    PetscScalar result;
    PetscScalar dx, *xa, *ya, *za;
    PetscScalar xmin, xmax, ymin, ymax;
    // below only if y data is evenly spaced
    //PetscScalar dy;
    PetscInt indx, indy, indz1, indz2, indz3, indz4;

    dx = I->dx;
    xa = I->xa;
    ya = I->ya;
    za = I->za;
    xmin = I->xmin;
    xmax = I->xmax;
    ymin = I->ymin;
    ymax = I->ymax;
    // below only if y data is evenly spaced
    //dy = I->dy;

    /* to reproduce the behaviour of scipy.interpolate.RectBivariateSpline
       the code should truncate if interpolation is attempted on a
       value outside of the lookup data range */

    /* for pressure (x), constant spacing assumed */
    if( x<=xmin ){
#if (defined VERBOSE)
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val2d: x<xmin.  Truncating\n");CHKERRQ(ierr);
#endif
      indx = 0; // minimum index, max index is always +1
      x = xmin;
    }
    else if( x>=xmax ){
#if (defined VERBOSE)
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val2d: x>xmax.  Truncating\n");CHKERRQ(ierr);
#endif
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
    if( y<=ymin ){
#if (defined VERBOSE)
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val2d: y<ymin.  Truncating\n");CHKERRQ(ierr);
#endif
      indy = 0; // minimum index, max index is always +1
      y = ymin;
    }
    else if( y>=ymax ){
#if (defined VERBOSE)
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: get_val2d: y>ymax.  Truncating\n");CHKERRQ(ierr);
#endif
      indy = NY-2; // minimum index, max index is always +1
      y = ymax;
    }
    else{
      // loop to find minimum index
      /* trivial algorithm to find minimum index when x data
         is not evenly spaced */
      indy = 0;
      while( (y-ya[indy])>0) {
        indy += 1;
      }
      /* loop exits when sign changes, meaning that previous index
         is the minimum index */
      indy -= 1;
    }

    // y weights
    w3 = y-ya[indy]; // y-y1
    w4 = ya[indy+1]-y; // y2-y

    indz1 = indy*NX+indx; // min S (y), min P (x)
    z1 = za[indz1];
    indz2 = indz1+1; // min S (y), max P (x)
    z2 = za[indz2];
    indz3 = indz1+NX; // max S (y), min P (x)
    z3 = za[indz3];
    indz4 = indz3+1; // max S (y), max P (x)
    z4 = za[indz4];

    // bilinear interpolation
    result = z1 * w2 * w4;
    result += z2 * w1 * w4;
    result += z3 * w2 * w3;
    result += z4 * w1 * w3;
    result /= dx; // dx
    result /= ya[indy+1]-ya[indy]; // dy

    return result;
}
