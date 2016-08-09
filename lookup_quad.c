#include "lookup.h"

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
    size_t i=0, j=0, k=0;
    char string[100], xtemp[30], ytemp[30], ztemp[30];
    PetscScalar x, y, z;
    PetscScalar xa[NX], ya[NY], za[NX*NY];
    PetscScalar xscale, yscale, zscale;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    {
      PetscErrorCode ierr;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"set_interp2d:\n");CHKERRQ(ierr);
    }
#endif

    // DJB disabled for time being since my own functions support quad
    /*if (sizeof(PetscScalar) != sizeof(double)){
      perror("PetscScalar must be double to use the dataio functions here");
      exit(-1);
    }*/

    /* bilinear interpolation */
    // DJB double-precision
    //const gsl_interp2d_type *T = gsl_interp2d_bilinear;
    //gsl_spline2d *spline = gsl_spline2d_alloc( T, NX, NY );
    //gsl_interp_accel *xacc = gsl_interp_accel_alloc();
    //gsl_interp_accel *yacc = gsl_interp_accel_alloc();

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
            /* PDS: this general form of reading a string and then
               converting to a float128 seems to be the way of doing it.
               But does Petsc have its own function for doing this? */
            // DJB for quad precision only
            scanf( string, "%s %s %s", &xtemp, &ytemp, &ztemp );
            xscale = strtoflt128(xtemp, NULL);
            yscale = strtoflt128(ytemp, NULL);
            zscale = strtoflt128(ztemp, NULL);
            // DJB for double precision only
            //sscanf( string, "%lf %lf %lf", &xscale, &yscale, &zscale );
        }
        if( i>=HEAD ){
            // DJB for quad precision only
            scanf( string, "%s %s %s", &xtemp, &ytemp, &ztemp );
            x = strtoflt128(xtemp, NULL);;
            y = strtoflt128(ytemp, NULL);
            z = strtoflt128(ztemp, NULL);
            // DJB for double precision only
            //sscanf(string, "%lf %lf %lf", &x, &y, &z );
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
    // DJB double precision
    //gsl_spline2d_init( spline, xa, ya, za, NX, NY );

    /* for debugging */
#if 0
    {
      PetscMPIInt rank;

      ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
      for (i=0; i<NX; i++ ){
          ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] %d %f\n", rank, i, xa[i]);CHKERRQ(ierr);
      }

      for (j=0; j<NY; j++ ){
          ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] %d %f\n", rank, j, ya[j]);
      }
    }
#endif

    interp->xmin= xa[0];
    interp->xmax= xa[NX-1];
    interp->ymin= ya[0];
    interp->ymax= ya[NY-1];

    // DJB for quad-precision
    // PDS: sizeof Petsc/MPI compatible?
    memmove( interp->xa, xa, sizeof interp->xa );
    memmove( interp->ya, ya, sizeof interp->ya );
    memmove( interp->za, za, sizeof interp->za );

    /* if we store the x and y step, we can more quickly locate the
       relevant indices in the arrays */

    /* DJB pressure is always constantly spaced (I think) */
    interp->dx = xa[1]-xa[0]; // TODO: confirm sign

    /* DJB TODO: this does not work at present because the entropy
       coordinate in the data tables is not constant */
    //interp->dy = ya[1]-ya[0]; // TODO: confirm sign

    // DJB double-precision
    //interp->interp = spline;
    //interp->xacc = xacc;
    //interp->yacc = yacc;

    PetscFunctionReturn(0);
}

static PetscErrorCode set_interp1d( const char * filename, Interp1d *interp, PetscInt n )
{

    FILE *fp;
    size_t i=0;
    char string[100], xtemp[30], ytemp[30];
    PetscScalar x, y, xscale, yscale;
    PetscScalar xa[n], ya[n];

    PetscFunctionBeginUser;

#if (defined VERBOSE)
    {
      PetscErrorCode ierr;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"set_interp1d:\n");CHKERRQ(ierr);
    }
#endif

    // DJB disabled for time being since my own function supports quad
    /*if (sizeof(PetscScalar) != sizeof(double)){
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"PetscScalar must be double to use the dataio functions here");
    } */

    // DJB double-precision
    /* linear interpolation */
    //const gsl_interp_type *T = gsl_interp_linear;
    //gsl_interp *interpolation = gsl_interp_alloc( T, n );
    //gsl_interp_accel *acc = gsl_interp_accel_alloc();

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
            // DJB quad precision only
            scanf( string, "%s %s", &xtemp, &ytemp );
            xscale = strtoflt128(xtemp, NULL);
            yscale = strtoflt128(ytemp, NULL);
            // DJB double precision only
            //sscanf( string, "%lf %lf", &xscale, &yscale );
            }
        if( i>=HEAD ){
            // DJB quad precision only
            scanf( string, "%s %s", &xtemp, &ytemp );
            x = strtoflt128(xtemp, NULL);
            y = strtoflt128(ytemp, NULL);
            // DJB double precision only
            //sscanf(string, "%lf %lf", &x, &y );
            xa[i-HEAD] = x * xscale;
            ya[i-HEAD] = y * yscale;
            }
        ++i;
    }

    fclose( fp );

    // DJB double-precision
    //gsl_interp_init( interpolation, xa, ya, n );

    // TODO: is this the correct way of copying an array?
    memmove( interp->xa, xa, sizeof interp->xa );
    interp->xmin= xa[0];
    interp->xmax= xa[n-1];
    memmove( interp->ya, ya, sizeof interp->ya );
    interp->ymin= ya[0];
    interp->ymax= ya[n-1];

    // DJB quad-precision
    // TODO: check that spacing is constant
    interp->dx = xa[1]-xa[0];

    // DJB double-precision
    //interp->interp = interpolation;
    //interp->acc = acc;

    PetscFunctionReturn(0);
}

PetscScalar get_val1d( Interp1d *I, PetscScalar x )
{
    /* wrapper for evaluating a 1-D lookup */

    PetscScalar weight, x1, y1, y2, result;
    PetscScalar dx, *xa, *ya;
    PetscInt ind;

    dx = I->dx;
    xa = I->xa;
    ya = I->ya;

    /* simple linear interpolation */
    // PDS: PetscFloorReal correct?
    /* DJB this assumes data is evenly spaced in input data file to
       enable us to locate the index in the array directly */
    /* TODO: this only works because zero pressure is always the first
       output.  Potential for a bug here */
    ind = PetscFloorReal( x/dx );
    x1 = xa[ind];  // x (pressure) to left
    y1 = ya[ind]; // y (quantity) to left
    y2 = ya[ind+1]; // y (quantity) to right
    weight = (x-x1) / dx; // 1.0 must equal y2, 0.0 must equal y1

    result = y2 * weight + y1 * (1.0-weight);

    // DJB double-precision
    //result = gsl_interp_eval( I->interp, I->xa, I->ya, x, I->acc );

    return result;
}

PetscScalar get_val2d( Interp2d *I, PetscScalar x, PetscScalar y )
{
    /* wrapper for evaluating a 2-D lookup */
    /* bilinear interpolation that supports quad precision */

    /* DJB TODO: does not work at present because entropy coordinate
       in the data tables is not constantly spaced.  But it looks
       like it can be, in which case this function will work well */

    PetscScalar x1, x2, y1, y2, z1, z2, z3, z4;
    PetscScalar result;
    PetscScalar dx, *xa, *ya, *za, xmin; //dy,ymin
    PetscInt indx, indy, indz1, indz2, indz3, indz4;
    PetscInt indt;

    dx = I->dx;
    xa = I->xa;
    // DJB only if entropy is evenly spaced in data files
    //dy = I->dy;
    ya = I->ya;
    za = I->za;
    xmin = I->xmin;
    // DJB only if entropy is evenly spaced in data files
    //ymin = I->ymin;

    indx = PetscFloorReal( (x-xmin)/dx );
    x1 = xa[indx]; // local min x
    x2 = xa[indx+1]; // local max x

    /* DJB this requires entropy to be evenly spaced, which it is not */
    //indy = PetscFloorReal( (y-ymin)/dy );

    /* DJB more generic approach to finding minimum index when entropy
       is not evenly spaced */
    indt = 0;
    while( (ya[indt]-y)<0) {
        indt += 1;
    }
    indy = indt-1;

    y1 = ya[indy]; // local min y
    y2 = ya[indy+1]; // local max y*/

    indz1 = indy*NX+indx; // min S, min P
    z1 = za[indz1];
    indz2 = indz1+1; // min S, max P
    z2 = za[indz2];
    indz3 = indz1+NX; // max S, min P
    z3 = za[indz3];
    indz4 = indz3+1; // max S, max P
    z4 = za[indz4];

    result = z1*(x2-x)*(y2-y);
    result += z2*(x-x1)*(y2-y);
    result += z3*(x2-x)*(y-y1);
    result += z4*(x-x1)*(y-y1);
    result /= dx*(y2-y1);

    // DJB double-precision
    //result = gsl_spline2d_eval( I->interp, x, y, I->xacc, I->yacc );

    return result;
}