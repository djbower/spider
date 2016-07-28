#include "lookup.h"

PetscErrorCode set_lookups( Ctx *E ) 
{
    PetscErrorCode ierr;
    /* set all 1-D and 2-D lookups */

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_lookup:\n" );CHKERRQ(ierr);
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

    PetscErrorCode ierr;
    FILE *fp;
    size_t i=0, j=0, k=0;
    char string[100];
    PetscScalar x, y, z;
    PetscScalar xa[NX], ya[NY], za[NX*NY];
    PetscScalar xscale, yscale, zscale;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_interp2d:\n");CHKERRQ(ierr);
#endif
    if (sizeof(PetscScalar) != sizeof(double)){
      perror("PetscScalar must be double to use the dataio functions here");
      exit(-1);
    }

    /* bilinear interpolation */
    const gsl_interp2d_type *T = gsl_interp2d_bilinear;
    gsl_spline2d *spline = gsl_spline2d_alloc( T, NX, NY );
    gsl_interp_accel *xacc = gsl_interp_accel_alloc();
    gsl_interp_accel *yacc = gsl_interp_accel_alloc();

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
            sscanf( string, "%lf %lf %lf", &xscale, &yscale, &zscale );
        }
        if( i>=HEAD ){
            sscanf(string, "%lf %lf %lf", &x, &y, &z );
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
    gsl_spline2d_init( spline, xa, ya, za, NX, NY );

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

    interp->interp = spline;
    interp->xacc = xacc;
    interp->yacc = yacc;

    PetscFunctionReturn(0);
}

static PetscErrorCode set_interp1d( const char * filename, Interp1d *interp, PetscInt n )
{

    PetscErrorCode ierr;
    FILE *fp;
    size_t i=0;
    char string[100];
    PetscScalar x, y, xscale, yscale;
    PetscScalar xa[n], ya[n];

    PetscFunctionBeginUser;

#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_interp1d:\n");CHKERRQ(ierr);
#endif
    if (sizeof(PetscScalar) != sizeof(double)){
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"PetscScalar must be double to use the dataio functions here");
    }

    /* linear interpolation */
    const gsl_interp_type *T = gsl_interp_linear;
    gsl_interp *interpolation = gsl_interp_alloc( T, n );
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

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
            sscanf( string, "%lf %lf", &xscale, &yscale );
            }
        if( i>=HEAD ){
            sscanf(string, "%lf %lf", &x, &y );
            xa[i-HEAD] = x * xscale;
            ya[i-HEAD] = y * yscale;
            }
        ++i;
    }

    fclose( fp );

    gsl_interp_init( interpolation, xa, ya, n );

    // TODO: is this the correct way of copying an array?
    memmove( interp->xa, xa, sizeof interp->xa );
    interp->xmin= xa[0];
    interp->xmax= xa[n-1];
    memmove( interp->ya, ya, sizeof interp->ya );
    interp->ymin= ya[0];
    interp->ymax= ya[n-1];

    interp->interp = interpolation;
    interp->acc = acc;

    PetscFunctionReturn(0);
}

PetscScalar get_val1d( Interp1d *I, PetscScalar x )
{
    /* wrapper for evaluating a 1-D lookup */

    PetscScalar result;

    result = gsl_interp_eval( I->interp, I->xa, I->ya, x, I->acc );

    return result;
}

PetscScalar get_val2d( Interp2d *I, PetscScalar x, PetscScalar y )
{
    /* wrapper for evaluating a 2-D lookup */

    PetscScalar result;

    result = gsl_spline2d_eval( I->interp, x, y, I->xacc, I->yacc );

    return result;
}


PetscErrorCode free_interp1d( Interp1d *I )
{
    /* free memory of Interp1d */

    PetscFunctionBeginUser;

    gsl_interp_free( I->interp );
    gsl_interp_accel_free( I->acc );

    PetscFunctionReturn(0);
}

PetscErrorCode free_interp2d( Interp2d *I )
{
    /* free memory of Interp2d */

    PetscFunctionBeginUser;

    gsl_spline2d_free( I->interp );
    gsl_interp_accel_free( I->xacc );
    gsl_interp_accel_free( I->yacc );

    PetscFunctionReturn(0);
}

PetscErrorCode free_memory_interp( Ctx *E )
{
    /* free memory allocated by interpolation functions */

    PetscFunctionBeginUser;

    /* liquidus and solidus lookups */
    free_interp1d( &E->solid_prop.liquidus );
    free_interp1d( &E->solid_prop.solidus );
    free_interp1d( &E->melt_prop.liquidus );
    free_interp1d( &E->melt_prop.solidus );

    /* solid properties lookup */
    free_interp2d( &E->solid_prop.alpha );
    free_interp2d( &E->solid_prop.cp );
    free_interp2d( &E->solid_prop.dTdPs );
    free_interp2d( &E->solid_prop.rho );
    free_interp2d( &E->solid_prop.temp );

    /* melt properties lookup */
    free_interp2d( &E->melt_prop.alpha );
    free_interp2d( &E->melt_prop.cp );
    free_interp2d( &E->melt_prop.dTdPs );
    free_interp2d( &E->melt_prop.rho );
    free_interp2d( &E->melt_prop.temp );

    PetscFunctionReturn(0);
}
