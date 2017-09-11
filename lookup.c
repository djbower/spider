#include "lookup.h"

static PetscErrorCode set_interp2d( const char *, Interp2d * );
static PetscErrorCode set_interp1d( const char *, Interp1d *, PetscInt );

/* Helper routine to prepend the root directory to a relative path */
static PetscErrorCode MakeRelativePathAbsolute(char* path) {
  PetscErrorCode ierr;
  char tmp[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  ierr = PetscStrcpy(tmp,path);CHKERRQ(ierr);
  ierr = PetscStrcpy(path,MAGMA_ROOT_DIR_STR);CHKERRQ(ierr); /* global_defs.h */
  ierr = PetscStrcat(path,"/");CHKERRQ(ierr); /* not portable */
  ierr = PetscStrcat(path,tmp);CHKERRQ(ierr); 
  PetscFunctionReturn(0);
}

PetscErrorCode set_lookups( Ctx *E )
{
    PetscErrorCode ierr;
    char liquidusFilename[PETSC_MAX_PATH_LEN];
    char solidusFilename[PETSC_MAX_PATH_LEN];
    char alphaSolFilename[PETSC_MAX_PATH_LEN];
    char alphaMelFilename[PETSC_MAX_PATH_LEN];
    char cpSolFilename[PETSC_MAX_PATH_LEN];
    char cpMelFilename[PETSC_MAX_PATH_LEN];
    char dtdpsSolFilename[PETSC_MAX_PATH_LEN];
    char dtdpsMelFilename[PETSC_MAX_PATH_LEN];
    char rhoSolFilename[PETSC_MAX_PATH_LEN];
    char rhoMelFilename[PETSC_MAX_PATH_LEN];
    char tempSolFilename[PETSC_MAX_PATH_LEN];
    char tempMelFilename[PETSC_MAX_PATH_LEN];

    /* set all 1-D and 2-D lookups */

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"set_lookup:\n" );CHKERRQ(ierr);
    }
#endif

    /* Based on a user-supplied flag, determine which files to load.
       Note that we prepend a string, MAGMA_ROOT_DIR_STR, and /, to make the 
       paths absolute, so the code can be run (tested) from anywhere.
       See global_defs.h for the paths */
    {
      char lscTypeString[PETSC_MAX_PATH_LEN] = "default";
      PetscBool isDefault=PETSC_FALSE,isStixrude2009=PETSC_FALSE,isAndrault2011=PETSC_FALSE;
      ierr = PetscOptionsGetString(NULL,NULL,"-curves",lscTypeString,PETSC_MAX_PATH_LEN,NULL);CHKERRQ(ierr);
      ierr = PetscStrcmp(lscTypeString,"default",&isDefault);CHKERRQ(ierr);
      ierr = PetscStrcmp(lscTypeString,"stixrude2009",&isStixrude2009);CHKERRQ(ierr);
      ierr = PetscStrcmp(lscTypeString,"andrault2011",&isAndrault2011);CHKERRQ(ierr);
      if (!(isDefault || isStixrude2009 || isAndrault2011)){
        ierr = PetscPrintf(PETSC_COMM_WORLD,"*************** WARNING ***************\n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Unrecognized -curves choice %s provided. Using defaults.\nCurrent options include:\n -curves stixrude2009\n -curves andrault2011\n",lscTypeString);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"***************************************\n");CHKERRQ(ierr);
      }
      if (isStixrude2009 || isDefault) { /* Default */
        ierr = PetscStrcpy(liquidusFilename,LIQUIDUS_STIXRUDE2009);CHKERRQ(ierr);
        ierr = MakeRelativePathAbsolute(liquidusFilename);CHKERRQ(ierr);
        ierr = PetscStrcpy(solidusFilename,SOLIDUS_STIXRUDE2009);CHKERRQ(ierr);
        ierr = MakeRelativePathAbsolute(solidusFilename);CHKERRQ(ierr);
      } else if (isAndrault2011) {
        ierr = PetscStrcpy(liquidusFilename,LIQUIDUS_ANDRAULT2011);CHKERRQ(ierr);
        ierr = MakeRelativePathAbsolute(liquidusFilename);CHKERRQ(ierr);
        ierr = PetscStrcpy(solidusFilename,SOLIDUS_ANDRAULT2011);CHKERRQ(ierr);
        ierr = MakeRelativePathAbsolute(solidusFilename);CHKERRQ(ierr);
      } else {
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"-curves logic processing failed");
      }

      ierr = PetscStrcpy(alphaSolFilename,ALPHA_SOL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(alphaSolFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(cpSolFilename,CP_SOL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(cpSolFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(dtdpsSolFilename,DTDPS_SOL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(dtdpsSolFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(rhoSolFilename,RHO_SOL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(rhoSolFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(tempSolFilename,TEMP_SOL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(tempSolFilename);CHKERRQ(ierr);

      ierr = PetscStrcpy(alphaMelFilename,ALPHA_MEL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(alphaMelFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(cpMelFilename,CP_MEL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(cpMelFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(dtdpsMelFilename,DTDPS_MEL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(dtdpsMelFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(rhoMelFilename,RHO_MEL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(rhoMelFilename);CHKERRQ(ierr);
      ierr = PetscStrcpy(tempMelFilename,TEMP_MEL_DEFAULT);CHKERRQ(ierr);
      ierr = MakeRelativePathAbsolute(tempMelFilename);CHKERRQ(ierr);
    }

    /* solid lookups */
    /* 2d */
    set_interp2d( alphaSolFilename, &E->solid_prop.alpha );
    set_interp2d( cpSolFilename, &E->solid_prop.cp );
    set_interp2d( dtdpsSolFilename, &E->solid_prop.dTdPs );
    set_interp2d( rhoSolFilename, &E->solid_prop.rho );
    set_interp2d( tempSolFilename, &E->solid_prop.temp );
    /* const */
    E->solid_prop.cond = COND_SOL;
    E->solid_prop.log10visc = LOG10VISC_SOL;

    /* melt lookups */
    /* 2d */
    set_interp2d( alphaMelFilename, &E->melt_prop.alpha );
    set_interp2d( cpMelFilename, &E->melt_prop.cp );
    set_interp2d( dtdpsMelFilename, &E->melt_prop.dTdPs );
    set_interp2d( rhoMelFilename, &E->melt_prop.rho );
    set_interp2d( tempMelFilename, &E->melt_prop.temp );
    /* const */
    E->melt_prop.cond = COND_MEL;
    E->melt_prop.log10visc = LOG10VISC_MEL;

    /* liquidus and solidus */
    /* 1d */

    set_interp1d( liquidusFilename, &E->solid_prop.liquidus, NLS );
    set_interp1d( liquidusFilename, &E->melt_prop.liquidus, NLS );
    /* duplication here, but want to remain flexible for future
       approaches for a multicomponent system */
    set_interp1d( solidusFilename, &E->solid_prop.solidus, NLS );
    set_interp1d( solidusFilename, &E->melt_prop.solidus, NLS );

    PetscFunctionReturn(0);
}

static PetscErrorCode set_interp2d( const char * filename, Interp2d *interp )
{
    FILE *fp;
    PetscInt i=0, j=0, k=0;
    char string[PETSC_MAX_PATH_LEN];
#if (defined PETSC_USE_REAL___FLOAT128)
    char xtemp[30], ytemp[30], ztemp[30];
#endif
   /* TODO: without assigning values below, the compiler warns about possible
      uninitialised values for just the quadruple precision case */
    PetscScalar xscale=0.0, yscale=0.0, zscale=0.0;
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
            za[i-HEAD] = z;
#if (defined LOOKUPDIM)
            za[i-HEAD] *= zscale;
#endif

            /* x coordinate */
            if( i<HEAD+NX ){
                xa[j] = x;
#if (defined LOOKUPDIM)
                xa[j] *= xscale;
#endif
                ++j;
            }
            /* y coordinate */
            if( (i-HEAD) % NX ==0 ){
                ya[k] = y;
#if (defined LOOKUPDIM)
                ya[k] *= yscale;
#endif
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
    char string[PETSC_MAX_PATH_LEN];
#if (defined PETSC_USE_REAL___FLOAT128)
    char xtemp[30], ytemp[30];
#endif
    PetscScalar x, y, xscale=0.0, yscale=0.0, xa[n], ya[n];

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
            xa[i-HEAD] = x;
            ya[i-HEAD] = y;
#if (defined LOOKUPDIM)
            xa[i-HEAD] *= xscale;
            ya[i-HEAD] *= yscale;
#endif
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
