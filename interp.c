#include "interp.h"
#include "util.h"

PetscErrorCode Interp1dCreateAndSet( const char *filename, Interp1d *interp_ptr, PetscScalar xconst, PetscScalar yconst )
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
            /* will break for 64 bit integers (but presumably
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

PetscErrorCode Interp2dCreateAndSet( const char * filename, Interp2d *interp_ptr, PetscScalar xconst, PetscScalar yconst, PetscScalar zconst )
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
            /* will break for 64 bit integers (but presumably
               elsewhere in the code will also break) */
            sscanf( string, "%d %d %d", &HEAD, &NX, &NY );
            /* make arrays based on the sizes in the header */
            ierr = PetscMalloc1( NX, &interp->xa ); CHKERRQ(ierr);
            ierr = PetscMalloc1( NY, &interp->ya ); CHKERRQ(ierr);
            ierr = Make2DPetscScalarArray( NX, NY, &interp->za ); CHKERRQ(ierr);
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

PetscErrorCode SetInterp1dValue( const Interp1d interp, PetscScalar x, PetscScalar *val, PetscScalar *dval )
{
    /* wrapper for evaluating a 1-D lookup
       linear interpolation with truncation for values
       that fall outside the data lookup range */

    PetscScalar const *xa, *ya;
    PetscScalar       w1, xmin, xmax;
    PetscInt          ind, NX;

    PetscFunctionBeginUser;

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

    /* interpolated value */
    if( val != NULL ){
        // w1 is 0 at leftmost (minimum) x, 1 at rightmost (maximum) x
        w1 = (x-xa[ind]) / (xa[ind+1]-xa[ind]); // weighting
        *val = ya[ind] * (1.0-w1) + ya[ind+1] * w1;
    }

    /* first-order derivative estimate */
    if( dval != NULL ){
        *dval = ( ya[ind+1] - ya[ind] ) / ( xa[ind+1] - xa[ind] );
    }

    PetscFunctionReturn(0);

}

PetscErrorCode SetInterp2dValue( const Interp2d interp, PetscScalar x, PetscScalar y, PetscScalar *val )
{
    /* wrapper for evaluating a 2-D lookup using bilinear
       interpolation.

       Note that this assumes that x data (pressure) is evenly
       spaced so we use a faster lookup approach by computing
       indices directly rather than looping through data */

    PetscScalar z1, z2, z3, z4;
    PetscScalar w1, w2, w3, w4; // weights
    PetscScalar const *xa, *ya;
    PetscInt NX, NY;
    PetscScalar dx, xmin, xmax, ymin, ymax;
    // below only if y data is evenly spaced
    //PetscScalar dy;
    PetscInt indx, indy;

    PetscFunctionBeginUser;

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
    *val = z1 * w2 * w4;
    *val += z2 * w1 * w4;
    *val += z3 * w2 * w3;
    *val += z4 * w1 * w3;
    *val /= dx; // dx
    *val /= ya[indy+1]-ya[indy]; // dy

    PetscFunctionReturn(0);
}

PetscErrorCode Interp1dDestroy( Interp1d *interp_ptr )
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

PetscErrorCode Interp2dDestroy( Interp2d *interp_ptr )
{

    PetscErrorCode ierr;
    Interp2d interp = *interp_ptr;

    PetscFunctionBeginUser;

    PetscInt i;
    PetscInt NX = interp->NX;

    for(i=0; i<NX; i++) {
        ierr = PetscFree(interp->za[i]);CHKERRQ(ierr);
    }
    ierr = PetscFree( interp->za );CHKERRQ(ierr);

    PetscFree( interp->xa );
    PetscFree( interp->ya );

    ierr = PetscFree(*interp_ptr);CHKERRQ(ierr);
    *interp_ptr = NULL;

    PetscFunctionReturn(0);

}
