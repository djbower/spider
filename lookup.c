#include "parameters.h"
#include "lookup.h"
#include "util.h"

PetscErrorCode set_interp2d( const char * filename, Interp2d *interp, PetscScalar xconst, PetscScalar yconst, PetscScalar zconst )
{
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
            interp->xa = Make1DPetscScalarArray( NX );
            interp->ya = Make1DPetscScalarArray( NY );
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

void free_interp2d( Interp2d *interp ){

    PetscInt i;
    PetscInt NX = interp->NX;

    for(i=0; i<NX; i++) {
        free( interp->za[i] );
    }
    free( interp->za );
    free( interp->xa );
    free( interp->ya );
}


PetscErrorCode set_interp1d( const char * filename, Interp1d *interp, PetscScalar xconst, PetscScalar yconst )
{
    FILE *fp;
    PetscInt i=0;
    char string[PETSC_MAX_PATH_LEN];
#if (defined PETSC_USE_REAL___FLOAT128)
    char xtemp[30], ytemp[30];
#endif
    PetscScalar x, y, xscale=0.0, yscale=0.0;
    PetscInt HEAD, NX;

    PetscFunctionBeginUser;
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
            /* FIXME: works for both double and quad? */
            sscanf( string, "%d %d", &HEAD, &NX );
            /* make arrays based on the sizes in the header */
            interp->xa = Make1DPetscScalarArray( NX );
            interp->ya = Make1DPetscScalarArray( NX );
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

PetscScalar get_val1d( Interp1d const *interp, PetscScalar x )
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

PetscScalar get_val2d( Interp2d const *interp, PetscScalar x, PetscScalar y )
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

void free_interp1d( Interp1d *interp ){

    free( interp->xa );
    free( interp->ya );

}
