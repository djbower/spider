#include "ic.h"
static PetscErrorCode set_ic_from_perturbation( Ctx *, Vec );
//static PetscErrorCode set_ic_from_melt_fraction( Ctx *, Vec );
//static PetscErrorCode set_ic_from_file( Ctx *, Vec );

PetscErrorCode set_initial_condition(Ctx *E, Vec dSdr_in) 
{

    PetscFunctionBeginUser;

    /* ic from ic_dSdr */
    set_ic_from_perturbation( E, dSdr_in );

    /* ic from melt fraction */
    /*set_ic_from_melt_fraction( E, dSdr_in ); */

    /* read ic from file */
    //set_ic_from_file( E, dSdr_in );

    PetscFunctionReturn(0);
}

/* TODO: refresh for new approach with dSdr_b */
/* possibly broken at present */
//static PetscErrorCode set_ic_from_file( Ctx *E, Vec dSdr_in )
//{
//    PetscErrorCode    ierr;
//    FILE              *fp;
//    PetscInt          i=0;
    /* initial condition file must be called this */
//    char              filename[7] = "X_s.0";
//    char              string[100];
//#if (defined PETSC_USE_REAL___FLOAT128)
//    char              xtemp[30], ytemp[30];
//#endif
//    PetscScalar       *arr;
//    PetscScalar       x, y;
    //PetscScalar       xa[NUMPTS_S_DEFAULT];//, ya[NUMPTS_S_DEFAULT];

//    PetscFunctionBeginUser;

//    ierr = VecGetArray( dSdr_in, &arr ); CHKERRQ(ierr);

//    fp = fopen( filename, "r" );

//    if(fp==NULL) {
//        perror("Error opening file.\n");
//        exit(-1);
//    }

    // fgets reads in string, sscanf processes it
//    while(fgets(string, sizeof(string), fp) != NULL) {
//#if (defined PETSC_USE_REAL___FLOAT128)
//        sscanf( string, "%s %s", xtemp, ytemp );
//        x = strtoflt128(xtemp, NULL);
//        y = strtoflt128(ytemp, NULL);
//#else
//        sscanf(string, "%lf %lf", &x, &y );
//#endif
//        arr[i] = y;
        //ya[i] = y;
//        ++i;
//    }

//    fclose( fp );

 //   ierr = VecRestoreArray( dSdr_in, &arr ); CHKERRQ(ierr);

//    PetscFunctionReturn(0);

//}

static PetscErrorCode set_ic_from_perturbation( Ctx *E, Vec dSdr_in )
{

    PetscErrorCode     ierr;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_ic_from_perturbation:\n");CHKERRQ(ierr);
#endif

    ierr = VecSet( dSdr_in, IC_DSDR ); CHKERRQ( ierr );

    /* these next two lines are simply convenient reminders that the
       first and last values are meaningless because the fluxes here
       are controlled by boundary conditions.  But for debugging and
       clarity it is convenient to explicitly set these values to
       zero */
    ierr = VecSetValue( dSdr_in, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
    // TODO: is next line compatible with command line input? */
    ierr = VecSetValue( dSdr_in, NUMPTS_B_DEFAULT-1, 0.0, INSERT_VALUES); CHKERRQ(ierr);

    VecAssemblyBegin( dSdr_in );
    VecAssemblyEnd( dSdr_in );

    PetscFunctionReturn(0);
}
