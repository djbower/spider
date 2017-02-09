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

    ierr = VecSet( dSdr_in, -IC_DSDR ); CHKERRQ( ierr );

    /* these next two lines are simply convenient reminders that the
       first and last values are meaningless because the fluxes here
       are controlled by boundary conditions.  But for debugging and
       clarity it is convenient to explicitly set these values to
       zero */
    ierr = VecSetValue( dSdr_in, 0, 0.0, INSERT_VALUES); CHKERRQ(ierr);
    /* TODO: is next line compatible with command line input? */
    ierr = VecSetValue( dSdr_in, NUMPTS_B_DEFAULT-1, 0.0, INSERT_VALUES); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/* TODO: refresh for new approach with dSdr_b */
/* possibly broken at present */
//static PetscErrorCode set_ic_from_melt_fraction( Ctx *E, Vec S_in )
//{

//    PetscErrorCode    ierr;
//    PetscScalar       meltf, val, maxval, step;
//    PetscInt          i, ilo;
//    Vec               liq_s, dfus_s;
//    Solution          *S;
//    PetscScalar       *arr;

//    PetscFunctionBeginUser;

  //  S = &E->solution;
  //  dfus_s = S->fusion_s;
  //  liq_s = S->liquidus_s;

    /* melt fraction contour to follow */
 //   meltf = 0.99;
 //   val = 1.0 - meltf;

  //  ierr = VecCopy( liq_s, S_in ); CHKERRQ(ierr);
 //   ierr = VecAXPY( S_in, -val, dfus_s ); CHKERRQ(ierr);
  //  ierr = VecShift( S_in, -1.0 ); CHKERRQ(ierr); // (1.0 is reference)

    /* find overturn point and set everything to right to overturn value */
 //   ierr = VecMax( S_in, &ilo, &maxval );
 //   ierr = VecGetArray( S_in, &arr );
  //  for(i=ilo; i<NUMPTS_S_DEFAULT; ++i){
 //       arr[i] = maxval;
 //   }

    /* entropy drop.  Everything to left cannot drop by more than this
       value relative to the overturn value */
 //   step = 0.0075;
 //   for(i=0; i<ilo; ++i){
 //       if(arr[i] < maxval-step){
 //           arr[i] = maxval-step;
 //       }
 //   }

 //   ierr = VecRestoreArray( S_in, &arr );

    /* now add small gradient to initial perturbed value */
//    make_super_adiabatic( E, S_in );

//    PetscFunctionReturn(0);
//}
