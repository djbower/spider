#include "ctx.h"

static PetscErrorCode make_super_adiabatic( Ctx *E )
{

    PetscErrorCode ierr;
#if (defined VERBOSE)
    printf("make_super_adiabatic:\n");
#endif

    PetscInt i,ilo_s,ihi_s;
    PetscScalar fac=1.0,*arr_S_s,val,pres_b_last; // percent
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Vec pres_b,pres_s;
    Mesh *M;
    Solution *S;
    PetscMPIInt rank,size;
    const PetscInt ix = NUMPTS-1; // Dangerous if PetscInt is not int!

    PetscFunctionBeginUser;
    M = &E->mesh;
    S = &E->solution;
    pres_b = M->pressure_b;
    pres_s = M->pressure_s;
    ierr = VecGetOwnershipRange(S->S_s,&ilo_s,&ihi_s);CHKERRQ(ierr);

    /* Scatter the last value to all procs */
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (rank == size-1) { // kind of a kludge
      ierr = VecGetValues(pres_b,1,&ix,&val);
    }
    pres_b_last = 0.0;
    MPI_Scatter(&val,1,MPIU_SCALAR,&pres_b_last,1,MPIU_SCALAR,size-1,PETSC_COMM_WORLD);

    
#if (defined DEBUGOUTPUT)
  ierr = PetscPrintf(PETSC_COMM_SELF,"rank %d, value of last point is %f\n",rank,pres_b_last);CHKERRQ(ierr);
#endif

    ierr = VecGetArray(S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = VecGetArrayRead(pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = VecGetArrayRead(pres_s,&arr_pres_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        arr_S_s[i] *= 1.0 + 0.01 * fac * arr_pres_s[i]/pres_b_last;
    }
    ierr = VecRestoreArray(S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(pres_s,&arr_pres_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_initial_condition( Ctx *E )
{
    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    printf("set_initial_condition:\n");
#endif
    S = &E->solution;

    ierr = VecSet(S->S_s,SINIT);CHKERRQ(ierr);

    make_super_adiabatic( E );

    PetscFunctionReturn(0);

}
static PetscErrorCode spherical_area( Vec radius, Vec area )
{
    /* non-dimensional area.  Not sure this really matters, would it
       would seem that keeping these scalings close to 1.0 is
       preferred for numerical accuracy? */

    PetscErrorCode    ierr;
    PetscScalar       *arr_a;
    const PetscScalar *arr_r;
    PetscInt          i,ilo,ihi;

    PetscFunctionBeginUser;
    ierr = VecGetOwnershipRange(radius,&ilo,&ihi);CHKERRQ(ierr);
    ierr = VecGetArrayRead(radius,&arr_r);CHKERRQ(ierr);
    ierr = VecGetArray(area,&arr_a);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        arr_a[i] = PetscPowScalar( arr_r[i]/RADOUT, 2.0 );
    }
    ierr = VecRestoreArray(area,&arr_a);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(radius,&arr_r);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode spherical_volume( Vec radius, Vec volume )
{
    /* non-dimensional volume.  Not sure this really matters, but it
       would seem that keeping these scalings close to 1.0 is 
       preferred for numerical accuracy? */

    PetscErrorCode    ierr;
    PetscScalar       *arr_v;
    const PetscScalar *arr_r;
    PetscInt          i,ilo,ihi;

    PetscFunctionBeginUser;
    ierr = VecGetOwnershipRange(radius,&ilo,&ihi);CHKERRQ(ierr);
    ierr = VecGetArrayRead(radius,&arr_r);CHKERRQ(ierr);
    ierr = VecGetArray(volume,&arr_v);CHKERRQ(ierr);
    for(i=ilo; i<ihi-1; ++i){ // Note upper loop bound is adjusted
        arr_v[i] = PetscPowScalar(arr_r[i]/RADOUT,3.0) - PetscPowScalar(arr_r[i+1]/RADOUT,3.0);
        arr_v[i] *= 1.0 / 3.0;
    }
    ierr = VecRestoreArray(volume,&arr_v);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(radius,&arr_r);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode mixing_length( Vec radius, Vec mix )
{
    PetscErrorCode    ierr;
    PetscScalar       *arr_m;
    const PetscScalar *arr_r;
    PetscInt          i,ilo,ihi;
    PetscScalar       rad1, rad2;

    PetscFunctionBeginUser;
    ierr = VecGetOwnershipRange(radius,&ilo,&ihi);CHKERRQ(ierr);
    ierr = VecGetArrayRead(radius,&arr_r);CHKERRQ(ierr);
    ierr = VecGetArray(mix,&arr_m);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        rad1 = arr_r[i] - RADIN;
        rad2 = RADOUT - arr_r[i];
        arr_m[i] = PetscMin( rad1, rad2 );
    }
    ierr = VecRestoreArray(mix,&arr_m);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(radius,&arr_r);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
static PetscErrorCode aw_pressure( Vec radius, Vec pressure )
{
    PetscErrorCode    ierr;
    PetscScalar       dep,*arr_p;
    const PetscScalar *arr_r;
    PetscInt          i,ilo,ihi;

    PetscFunctionBeginUser;
    ierr = VecGetOwnershipRange(radius,&ilo,&ihi);CHKERRQ(ierr);
    ierr = VecGetArrayRead(radius,&arr_r);CHKERRQ(ierr);
    ierr = VecGetArray(pressure,&arr_p);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        dep = RADOUT - arr_r[i];
        arr_p[i] = -RHOS * GRAVITY / BETA;
        arr_p[i] *= PetscExpScalar( BETA * dep ) - 1.0;
    }
    ierr = VecRestoreArray(pressure,&arr_p);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(radius,&arr_r);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode aw_pressure_gradient( Vec radius, Vec grad )
{
    PetscErrorCode    ierr;
    PetscScalar       dep,*arr_g;
    const PetscScalar *arr_r;
    PetscInt          i,ilo,ihi;

    PetscFunctionBeginUser;
    ierr = VecGetOwnershipRange(radius,&ilo,&ihi);CHKERRQ(ierr);
    ierr = VecGetArray(grad,&arr_g);CHKERRQ(ierr);
    ierr = VecGetArrayRead(radius,&arr_r);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        dep = RADOUT - arr_r[i];
        arr_g[i] = RHOS * GRAVITY * PetscExpScalar( BETA * dep );
    }
    ierr = VecRestoreArray(grad,&arr_g);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(radius,&arr_r);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode set_mesh( Ctx *E)
{
#if (defined VERBOSE)
    printf( "set_mesh:\n" );
#endif

    PetscErrorCode ierr;
    PetscScalar    *arr;
    PetscInt       i,ilo_b,ihi_b,ilo_s,ihi_s;
    Mesh           *M;

    PetscFunctionBeginUser;

    M = &E->mesh;

    /* Create vectors required for the mesh */
    for (i=0;i<NUMMESHVECS;++i) {
      ierr = DMCreateGlobalVector(E->da_b,&M->meshVecs[i]);CHKERRQ(ierr);
    }
    // TODO: We should be creating ALL the vecs at the same time - do that as we roll all this into one function (create these vecs outside this function)

    M->area_b     = M->meshVecs[0]; 
    M->dPdr_b     = M->meshVecs[1];
    M->pressure_b = M->meshVecs[2];
    M->radius_b   = M->meshVecs[3];
    M->mix_b      = M->meshVecs[4];

    for (i=0;i<NUMMESHVECSS;++i) {
      ierr = DMCreateGlobalVector(E->da_s,&M->meshVecsS[i]);CHKERRQ(ierr);
    }
    M->pressure_s = M->meshVecsS[0];
    M->radius_s   = M->meshVecsS[1];
    M->volume_s   = M->meshVecsS[2];

    /* basic node spacing (negative) */
    M->dx_b = -(RADOUT-RADIN) / (NUMPTS-1);

    /* staggered node spacing (negative) */
    M->dx_s = M->dx_b;

    /* radius at basic nodes */
    ierr = VecGetOwnershipRange(M->radius_b,&ilo_b,&ihi_b);CHKERRQ(ierr);
    ierr = VecGetArray(M->radius_b,&arr);CHKERRQ(ierr);
    for (i=ilo_b; i<ihi_b; ++i){
        //M->radius_b[NUMPTS-1-i] = RADIN - i*M->dx_b;
        arr[i] = RADIN - (NUMPTS-1-i)*M->dx_b; //TODO: check
    }
    ierr = VecRestoreArray(M->radius_b,&arr);CHKERRQ(ierr);

    /* radius at staggered nodes */
    ierr = VecGetOwnershipRange(M->radius_s,&ilo_s,&ihi_s);CHKERRQ(ierr);
    ierr = VecGetArray(M->radius_s,&arr);CHKERRQ(ierr);
    for (i=ilo_s;i<ihi_s;++i){
        //M->radius_s[NUMPTS-2-i] = RADIN - 0.5*M->dx_b - i*M->dx_b;
        arr[i] = RADIN - 0.5*M->dx_b - (NUMPTSS-1-i)*M->dx_b; // TODO: check
    }
    ierr = VecRestoreArray(M->radius_s,&arr);CHKERRQ(ierr);

    /* Adams-Williamson EOS */

    /* pressure at basic nodes */
    aw_pressure( M->radius_b, M->pressure_b);

    /* dP/dr at basic nodes */
    aw_pressure_gradient( M->radius_b, M->dPdr_b);

    /* pressure at staggered nodes */
    aw_pressure( M->radius_s, M->pressure_s);

    /* surface area at basic nodes, without 4*pi term */
    /* and now non-dimensional */
    spherical_area( M->radius_b, M->area_b);

    /* volume of spherical cells, without 4*pi term */
    /* and now non-dimenisonal */
    spherical_volume( M->radius_b, M->volume_s);

    /* mixing length is minimum distance from boundary */
    mixing_length( M->radius_b, M->mix_b);

    PetscFunctionReturn(0);

}

PetscErrorCode set_lookups( Ctx *E )
{
  PetscFunctionBeginUser;
#if (defined VERBOSE)
    printf( "set_lookup:\n" );
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


void free_interp1d( Interp1d *I )
{

#if (defined VERBOSE)
    printf( "free_interp1d:\n" );
#endif

    gsl_interp_free( I->interp );
    gsl_interp_accel_free( I->acc );

}

void free_interp2d( Interp2d *I )
{

#if (defined VERBOSE)
    printf( "free_interp2d:\n" );
#endif

    gsl_spline2d_free( I->interp );
    gsl_interp_accel_free( I->xacc );
    gsl_interp_accel_free( I->yacc );
}

PetscErrorCode free_memory_interp( Ctx *E )
{
    PetscFunctionBeginUser;
    /* free memory allocated by interpolation functions */

#if (defined VERBOSE)
    printf( "free_memory_interp:\n" );
#endif

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

/* time-independent quantities */

static PetscErrorCode set_liquidus( Ctx *E )
{
#if (defined VERBOSE)
    printf("set_liquidus:\n");
#endif

    PetscErrorCode ierr;
    PetscInt       i,ilo_b,ihi_b,ilo_s,ihi_s;
    Vec            pres_b,pres_s;
    PetscScalar    z,*arr_liquidus,*arr_liquidus_rho,*arr_liquidus_temp,*arr_liquidus_s,*arr_liquidus_rho_s,*arr_liquidus_temp_s;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Interp1d       *I;
    Solution       *S;
    Interp2d       *IR, *IT;


    PetscFunctionBeginUser;
    S = &E->solution;
    I = &E->melt_prop.liquidus;
    IR = &E->melt_prop.rho;
    IT = &E->melt_prop.temp;

    pres_b = E->mesh.pressure_b;
    pres_s = E->mesh.pressure_s;

    ierr = VecGetOwnershipRange(pres_b,&ilo_b,&ihi_b);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(pres_s,&ilo_s,&ihi_s);CHKERRQ(ierr);

    /* basic nodes */
    ierr = VecGetArrayRead(pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = VecGetArray(S->liquidus,&arr_liquidus);CHKERRQ(ierr);
    ierr = VecGetArray(S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = VecGetArray(S->liquidus_temp,&arr_liquidus_temp);CHKERRQ(ierr);
    for(i=ilo_b;i<ihi_b;++i){
        z = get_val1d( I, arr_pres_b[i] );
        arr_liquidus[i] = z;
        arr_liquidus_rho[i] = get_val2d( IR, arr_pres_b[i], z );
        arr_liquidus_temp[i] = get_val2d( IT, arr_pres_b[i], z );
    }
    ierr = VecRestoreArrayRead(pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->liquidus,&arr_liquidus);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->liquidus_temp,&arr_liquidus_temp);CHKERRQ(ierr);


    /* staggered nodes */
    /* need in order to compute dfus/dr at basic internal nodes */ 
    ierr = VecGetArrayRead(pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = VecGetArray(S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = VecGetArray(S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = VecGetArray(S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        z = get_val1d( I, arr_pres_s[i] );
        arr_liquidus_s[i] = z;
        arr_liquidus_rho_s[i] = get_val2d( IR, arr_pres_s[i], z );
        arr_liquidus_temp_s[i] = get_val2d( IT, arr_pres_s[i], z );
    }
    ierr = VecRestoreArrayRead(pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_solidus( Ctx *E )
{
#if (defined VERBOSE)
    printf("set_solidus:\n");
#endif

    PetscErrorCode ierr;
    PetscInt       i,ilo_b,ihi_b,ilo_s,ihi_s;
    Vec            pres_b,pres_s;
    PetscScalar    z,*arr_solidus,*arr_solidus_rho,*arr_solidus_temp,*arr_solidus_s,*arr_solidus_rho_s,*arr_solidus_temp_s;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Interp1d       *I;
    Solution       *S;
    Interp2d       *IR, *IT;

    PetscFunctionBeginUser;
    S = &E->solution;
    I = &E->melt_prop.solidus;
    IR = &E->melt_prop.rho;
    IT = &E->melt_prop.temp;

    pres_b = E->mesh.pressure_b;
    pres_s = E->mesh.pressure_s;

    ierr = VecGetOwnershipRange(pres_b,&ilo_b,&ihi_b);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(pres_s,&ilo_s,&ihi_s);CHKERRQ(ierr);

    /* basic nodes */
    ierr = VecGetArrayRead(pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = VecGetArray(S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = VecGetArray(S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);
    ierr = VecGetArray(S->solidus_temp,&arr_solidus_temp);CHKERRQ(ierr);
    for(i=ilo_b;i<ihi_b;++i){
        z = get_val1d( I, arr_pres_b[i] );
        arr_solidus[i] = z;
        arr_solidus_rho[i] = get_val2d( IR, arr_pres_b[i], z );
        arr_solidus_temp[i] = get_val2d( IT, arr_pres_b[i], z );
    }
    ierr = VecRestoreArrayRead(pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->solidus_temp,&arr_solidus_temp);CHKERRQ(ierr);


    /* staggered nodes */
    /* need in order to compute dfus/dr at basic internal nodes */ 
    ierr = VecGetArrayRead(pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = VecGetArray(S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = VecGetArray(S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = VecGetArray(S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        z = get_val1d( I, arr_pres_s[i] );
        arr_solidus_s[i] = z;
        arr_solidus_rho_s[i] = get_val2d( IR, arr_pres_s[i], z );
        arr_solidus_temp_s[i] = get_val2d( IT, arr_pres_s[i], z );
    }
    ierr = VecRestoreArrayRead(pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

static PetscErrorCode set_fusion( Ctx *E )
{
#if (defined VERBOSE)
    printf("set_fusion:\n");
#endif

    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;
    S = &E->solution;

    /* basic nodes : fusion = liquidus - solidus*/
    ierr = VecWAXPY(S->fusion,-1.0,S->solidus,S->liquidus);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_rho,-1.0,S->solidus_rho,S->liquidus_rho);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_temp,-1.0,S->solidus_temp,S->liquidus_temp);CHKERRQ(ierr);

    /* staggered nodes */
    // TODO: Check - Intentionally no Rho here? 
    ierr = VecWAXPY(S->fusion_s,-1.0,S->solidus_s,S->liquidus_s);CHKERRQ(ierr);
    //ierr = VecWAXPY(S->fusion_rho_s,-1.0,S->solidus_rho_s,S->liquidus_rho_s);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_temp_s,-1.0,S->solidus_temp_s,S->liquidus_temp_s);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


static PetscErrorCode set_fusion_curve( Ctx *E )
{
    PetscErrorCode ierr;
    Solution *S;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    printf("set_fusion_curve:\n");
#endif
    S = &E->solution;

    /* basic nodes fusion_curve = solidus + 0.5*fusion*/
    ierr = VecWAXPY(S->fusion_curve,0.5,S->fusion,S->solidus);
    ierr = VecWAXPY(S->fusion_curve_temp,0.5,S->fusion_temp,S->solidus_temp);

    /* staggered nodes */
    ierr = VecWAXPY(S->fusion_curve_s,0.5,S->fusion_s,S->solidus_s);
    ierr = VecWAXPY(S->fusion_curve_temp_s,0.5,S->fusion_temp_s,S->solidus_temp_s);

    d_dr( E, S->fusion_curve_s, S->dfusdr );
    d_dr( E, S->fusion_curve_temp_s, S->dfusdr_temp );
    PetscFunctionReturn(0);
}

static PetscErrorCode set_mixed_phase( Ctx *E )
{
    PetscErrorCode ierr;
    Mesh *M;
    Solution *S;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    printf("set_mixed_phase:\n");
#endif

    M = &E->mesh;
    S = &E->solution;

    /* dTdrs_mix = dfusdr * (fusion_temp/fusion) + dfusdr_temp */
    ierr = VecPointwiseDivide(S->dTdrs_mix,S->fusion_temp,S->fusion);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->dTdrs_mix,S->dTdrs_mix,S->dfusdr);CHKERRQ(ierr);
    ierr = VecAXPY(S->dTdrs_mix,1.0,S->dfusdr_temp);CHKERRQ(ierr);

    ierr = VecPointwiseMult(S->dTdPs_mix,S->dTdrs_mix,M->dPdr_b);CHKERRQ(ierr);

    ierr = VecPointwiseDivide(S->cp_mix,S->fusion,S->fusion_temp);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->cp_mix,S->cp_mix,S->fusion_curve_temp);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
PetscErrorCode set_time_independent( Ctx *E )
{
  
    PetscFunctionBeginUser;

    /* linear interpolation using datafiles.  In the python code we
       fit a high order polynomial instead for smoothness, but
       perhaps this does not matter so much */
    set_liquidus( E );
    set_solidus( E );

    /* these terms all need the liquidus and solidus to be set */
    set_fusion( E );
    set_fusion_curve( E );
    set_mixed_phase( E );
    
    PetscFunctionReturn(0);
}


/* end of time-independent quantities */

//NOTE: This function MUST be called first when computing the RHS
PetscErrorCode set_capacitance( Ctx *E, Vec S_in )
{
    PetscErrorCode ierr;
    PetscInt i,ilo_s,ihi_s;
    Mesh *M;
    Solution *S;
    Vec pres_s;
    Interp2d *IRM, *ITM, *IRS, *ITS;
    const PetscScalar *arr_S_s,*arr_pres_s,*arr_liquidus_rho_s,*arr_solidus_rho_s,*arr_liquidus_temp_s,*arr_solidus_temp_s;
    PetscScalar *arr_phi_s,*arr_rho_s,*arr_temp_s;

    M = &E->mesh;
    S = &E->solution;
    pres_s = M->pressure_s;

    
    PetscFunctionBeginUser;

    IRM = &E->melt_prop.rho;
    ITM = &E->melt_prop.temp;
    IRS = &E->solid_prop.rho;
    ITS = &E->solid_prop.temp;

    ierr = VecCopy(S_in,S->S_s);CHKERRQ(ierr);

    // S->phi_s = (S->S_s - S->solidus_s)/S->fusion_s
    ierr = VecWAXPY(S->phi_s,-1.0,S->solidus_s,S->S_s);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(S->phi_s,S->phi_s,S->fusion_s);CHKERRQ(ierr);

    ierr = VecGetOwnershipRange(S->S_s,&ilo_s,&ihi_s);CHKERRQ(ierr);

    ierr = VecGetArrayRead(S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = VecGetArrayRead(pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    ierr = VecGetArray(S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = VecGetArray(S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = VecGetArray(S->temp_s,&arr_temp_s);CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){
        /* melt only */
        if( arr_phi_s[i] >= 1.0 ){
            arr_phi_s[i] = 1.0;
            arr_rho_s[i] = get_val2d( IRM, arr_pres_s[i], arr_S_s[i] );
            arr_temp_s[i] = get_val2d( ITM, arr_pres_s[i], arr_S_s[i] );
        }
        /* solid only */
        else if( arr_phi_s[i] <= 0.0 ){
            arr_phi_s[i] = 0.0;
            arr_rho_s[i] = get_val2d( IRS, arr_pres_s[i], arr_S_s[i] );
            arr_temp_s[i] = get_val2d( ITS, arr_pres_s[i], arr_S_s[i] );
        }
        /* mixed phase */
        else{
            /* density */
            arr_rho_s[i] = combine_matprop( arr_phi_s[i], 1.0/arr_liquidus_rho_s[i], 1.0/arr_solidus_rho_s[i] );
            arr_rho_s[i] = 1.0 / arr_rho_s[i];
            /* temperature */
            arr_temp_s[i] = combine_matprop( arr_phi_s[i], arr_liquidus_temp_s[i], arr_solidus_temp_s[i] );
        }

    }

    ierr = VecRestoreArrayRead(S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->temp_s,&arr_temp_s);CHKERRQ(ierr);

    // S->lhs_si = M->volume_si * S->rho_si * S->temp_si;
    ierr = VecPointwiseMult(S->lhs_s,M->volume_s,S->rho_s);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->lhs_s,S->lhs_s,S->rho_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscScalar viscosity_mix( PetscScalar meltf )
{
    PetscScalar arg, fwt, log10visc, visc;

    /* simple hyperbolic tan to get up and running */
    /* need to code up function with skew */
    arg = (1.0-meltf-F_THRESHOLD) / DF_TRANSITION;
    fwt = 0.5 * ( 1.0 + tanh(arg) );

    log10visc = fwt * LOG10VISC_SOL + (1.0 - fwt) * LOG10VISC_MEL;
    visc = PetscPowScalar( 10.0, log10visc );

    return visc;

}

PetscErrorCode set_matprop_and_flux( Ctx *E )
{
    PetscErrorCode ierr;
    PetscInt i,ilo_b,ihi_b,ilo,ihi;
    Vec pres;
    PetscScalar dr;
    PetscScalar *arr_S,*arr_phi,*arr_nu,*arr_gsuper,*arr_kappah,*arr_Etot,*arr_dSdr,*arr_dTdPs,*arr_dTdrs,*arr_dphidr,*arr_alpha,*arr_temp,*arr_cp,*arr_cond,*arr_visc,*arr_rho,*arr_Jcond,*arr_Jconv,*arr_Jtot,*arr_Jheat,*arr_Jmix,*arr_Jgrav,*arr_Jmass;
    const PetscScalar *arr_S_s,*arr_phi_s,*arr_solidus,*arr_fusion,*arr_pres,*arr_area_b,*arr_dPdr_b,*arr_liquidus_rho,*arr_solidus_rho,*arr_cp_mix,*arr_dTdrs_mix,*arr_liquidus_temp,*arr_solidus_temp,*arr_fusion_rho,*arr_fusion_temp,*arr_mix_b;
    Mesh *M;
    Solution *S;

    PetscFunctionBeginUser;
    M = &E->mesh;
    S = &E->solution;

    dr = M->dx_s;
    pres = M->pressure_b;

    /* loop over all basic internal nodes */
    ierr = VecGetOwnershipRange(M->pressure_b,&ilo_b,&ihi_b);CHKERRQ(ierr); //TODO: would be better to get these from the DM, but whatever for now
    ilo = ilo_b == 0      ? 1          : ilo_b;
    ihi = ihi_b == NUMPTS ? NUMPTS - 1 : ihi_b;

    ierr = VecGetArray(S->dphidr,&arr_dphidr);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = VecGetArray(S->phi,&arr_phi);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->fusion,&arr_fusion);CHKERRQ(ierr);
    ierr = VecGetArray(S->S,&arr_S);CHKERRQ(ierr);
    ierr = VecGetArrayRead(pres,&arr_pres);CHKERRQ(ierr);
    ierr = VecGetArray(S->dTdPs,&arr_dTdPs);CHKERRQ(ierr);
    ierr = VecGetArray(S->dTdrs,&arr_dTdrs);CHKERRQ(ierr);
    ierr = VecGetArray(S->cp,&arr_cp);CHKERRQ(ierr);
    ierr = VecGetArray(S->temp,&arr_temp);CHKERRQ(ierr);
    ierr = VecGetArray(S->alpha,&arr_alpha);CHKERRQ(ierr);
    ierr = VecGetArray(S->cond,&arr_cond);CHKERRQ(ierr);
    ierr = VecGetArray(S->visc,&arr_visc);CHKERRQ(ierr);
    ierr = VecGetArray(S->rho,&arr_rho);CHKERRQ(ierr);
    ierr = VecGetArrayRead(M->dPdr_b,&arr_dPdr_b);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->cp_mix,&arr_cp_mix);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->dTdrs_mix,&arr_dTdrs_mix);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->liquidus_temp,&arr_liquidus_temp);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->fusion_rho,&arr_fusion_rho);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->fusion_temp,&arr_fusion_temp);CHKERRQ(ierr);
    ierr = VecGetArray(S->nu,&arr_nu);CHKERRQ(ierr);
    ierr = VecGetArray(S->gsuper,&arr_gsuper);CHKERRQ(ierr);
    ierr = VecGetArray(S->dSdr,&arr_dSdr);CHKERRQ(ierr);
    ierr = VecGetArray(S->kappah,&arr_kappah);CHKERRQ(ierr);
    ierr = VecGetArray(S->Etot,&arr_Etot);CHKERRQ(ierr);
    ierr = VecGetArrayRead(M->area_b,&arr_area_b);CHKERRQ(ierr);
    ierr = VecGetArrayRead(M->mix_b,&arr_mix_b);CHKERRQ(ierr);
    ierr = VecGetArray(S->Jcond,&arr_Jcond);CHKERRQ(ierr);
    ierr = VecGetArray(S->Jconv,&arr_Jconv);CHKERRQ(ierr);
    ierr = VecGetArray(S->Jheat,&arr_Jheat);CHKERRQ(ierr);
    ierr = VecGetArray(S->Jtot,&arr_Jtot);CHKERRQ(ierr);
    ierr = VecGetArray(S->Jmix,&arr_Jmix);CHKERRQ(ierr);
    ierr = VecGetArray(S->Jgrav,&arr_Jgrav);CHKERRQ(ierr);
    ierr = VecGetArray(S->Jmass,&arr_Jmass);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->solidus_temp,&arr_solidus_temp);CHKERRQ(ierr);
    ierr = VecGetArrayRead(S->phi_s,&arr_phi_s);CHKERRQ(ierr);

    for(i=ilo; i<ihi; ++i){ // Note that these correspond to basic nodes being updated, but we also assume an ordering on the staggered nodes!

      /* these quantities are determined from staggered node
         quantities, and these are previously computed in
         set_capacitance */
      /* entropy: from simple average */
      arr_S[i] = average( arr_S_s[i], arr_S_s[i-1] ); 

      /* dSdr: central difference, 2nd order accurate */
      arr_dSdr[i] = 1.0/dr * (arr_S_s[i]-arr_S_s[i-1]);

      /* dphidr: central difference, 2nd order accurate */
      arr_dphidr[i] = 1.0/dr * (arr_phi_s[i]-arr_phi_s[i-1]);
      /* melt fraction */
      arr_phi[i] = (arr_S[i] - arr_solidus[i]) / arr_fusion[i];

      /* melt only */
      if( arr_phi[i] >= 1.0 ){
        Lookup *L = &E->melt_prop;
        arr_phi[i] = 1.0;

        /* density */
        arr_rho[i] = get_val2d( &L->rho, arr_pres[i], arr_S[i] );

        /* adiabatic temperature gradient */
        arr_dTdPs[i] = get_val2d( &L->dTdPs, arr_pres[i], arr_S[i] );
        arr_dTdrs[i] = arr_dTdPs[i] * arr_dPdr_b[i];

        /* heat capacity */
        arr_cp[i] = get_val2d( &L->cp, arr_pres[i], arr_S[i] );

        /* temperature */
        arr_temp[i] = get_val2d( &L->temp, arr_pres[i], arr_S[i] );

        /* thermal expansion coefficient */
        arr_alpha[i] = get_val2d( &L->alpha, arr_pres[i], arr_S[i] );

        /* thermal conductivity */
        arr_cond[i] = COND_MEL;

        /* viscosity */
        arr_visc[i] = PetscPowScalar( 10.0, LOG10VISC_MEL );
      }

      /* solid only */
      else if( arr_phi[i] <= 0.0 ){
        Lookup *L = &E->solid_prop;
        arr_phi[i] = 0.0;

        /* density */
        arr_rho[i] = get_val2d( &L->rho, arr_pres[i], arr_S[i] );

        /* adiabatic temperature gradient */
        arr_dTdPs[i] = get_val2d( &L->dTdPs, arr_pres[i], arr_S[i] );
        arr_dTdrs[i] = arr_dTdPs[i] * arr_dPdr_b[i];

        /* heat capacity */
        arr_cp[i] = get_val2d( &L->cp, arr_pres[i], arr_S[i] );

        /* temperature */
        arr_temp[i] = get_val2d( &L->temp, arr_pres[i], arr_S[i] );

        /* thermal expansion coefficient */
        arr_alpha[i] = get_val2d( &L->alpha, arr_pres[i], arr_S[i] );

        /* thermal conductivity */
        arr_cond[i] = COND_SOL;

        /* viscosity */
        arr_visc[i] = pow( 10, LOG10VISC_SOL );
      }

      /* mixed phase */
      else{
        /* density */

        arr_rho[i] = combine_matprop( arr_phi[i], 1.0/arr_liquidus_rho[i], 1.0/arr_solidus_rho[i] );
        arr_rho[i] = 1.0 / arr_rho[i];

        /* adiabatic temperature gradient */
        arr_dTdrs[i] = arr_dTdrs_mix[i];
        arr_dTdPs[i] = arr_dTdrs[i] * 1.0 / arr_dPdr_b[i];

        /* heat capacity */
        arr_cp[i] = arr_cp_mix[i];

        /* temperature */
        arr_temp[i] = combine_matprop( arr_phi[i], arr_liquidus_temp[i], arr_solidus_temp[i] );

        /* thermal expansion coefficient */
        arr_alpha[i] = -arr_fusion_rho[i] / arr_fusion_temp[i] * 1.0 / arr_rho[i];

        /* thermal conductivity */
        arr_cond[i] = combine_matprop( arr_phi[i], COND_MEL, COND_SOL );

        /* viscosity */
        arr_visc[i] = viscosity_mix( arr_phi[i] );

        /* other useful material properties */
        /* kinematic viscosity */
        arr_nu[i] = arr_visc[i] / arr_rho[i];

        /* gravity * super-adiabatic temperature gradient */
        arr_gsuper[i] = GRAVITY * arr_temp[i] / arr_cp[i] * arr_dSdr[i];

        /* eddy diffusivity */
        PetscScalar kh,crit;
        crit = 81.0 * PetscPowScalar(arr_nu[i],2);
        crit /= 4.0 * arr_alpha[i] * PetscPowScalar(arr_mix_b[i],4);

        if( arr_gsuper[i] < 0.0 ){
          /* no convection, subadiabatic */
          kh = 0.0;
        } else if( arr_gsuper[i] > crit ){
          /* inviscid scaling from Vitense (1953) */
          kh = 0.25 * PetscPowScalar(arr_mix_b[i],2) * PetscSqrtScalar(arr_alpha[i]*arr_gsuper[i]);
        } else{
          /* viscous scaling */
          kh = arr_alpha[i] * arr_gsuper[i] * PetscPowScalar(arr_mix_b[i],4) / (18.0*arr_nu[i]);
        }
        arr_kappah[i] = kh;

        /* can now compute energy fluxes */

        /* conductive heat flux */
        arr_Jcond[i] = arr_temp[i] / arr_cp[i] * arr_dSdr[i] + arr_dTdrs[i];
        arr_Jcond[i] *= -arr_cond[i];

        /* convective heat flux */
        arr_Jconv[i] = -arr_rho[i] * arr_temp[i] * arr_kappah[i] * arr_dSdr[i];

        /* also resets arrays since these change every time rhs
           is called */
        arr_Jheat[i] = arr_Jcond[i] + arr_Jconv[i];
        arr_Jtot[i] = arr_Jheat[i];

        //TODO: Need to clean up these declarations..
        PetscScalar dfus = arr_fusion[i];
        PetscScalar kappah = arr_kappah[i];
        PetscScalar rho = arr_rho[i];
        PetscScalar rhol = arr_liquidus_rho[i];
        PetscScalar rhos = arr_solidus_rho[i];
        PetscScalar temp = arr_temp[i];

        PetscScalar pref = temp * dfus;

        /* convective mixing */
        arr_Jmix[i] = -pref * kappah * rho * arr_dphidr[i];

        /* gravitational separation */
        // TODO: This all needs serious cleanup
        PetscScalar phi = arr_phi[i];

        PetscScalar cond1 = rhol / (11.993*rhos + rhol);
        PetscScalar cond2 = rhol / (0.29624*rhos + rhol);
        PetscScalar F;

        if(phi<cond1){
          F = 0.001*PetscPowScalar(rhos,2)*PetscPowScalar(phi,3)/(PetscPowScalar(rhol,2)*(1.0-phi));
        } else if(phi>cond2){
          F = 2.0/9.0 * phi * (1.0-phi);
        } else{
          F = 5.0/7.0*PetscPowScalar(rhos,4.5)*PetscPowScalar(phi,5.5)*(1.0-phi);
          F /= PetscPowScalar( rhol+(rhos-rhol)*phi, 4.5 );
        }

        arr_Jgrav[i] = (rhol-rhos) * rhol * rhos / rho;
        arr_Jgrav[i] *= pref * PetscPowScalar(GRAIN,2) * GRAVITY * F;
        arr_Jgrav[i] /= PetscPowScalar(10.0, LOG10VISC_MEL);

        arr_Jmass[i] = arr_Jmix[i] + arr_Jgrav[i];
        arr_Jtot[i] += arr_Jmass[i];

        arr_Etot[i] = arr_Jtot[i] * arr_area_b[i];

      }
    }
    ierr = VecRestoreArray(S->dphidr,&arr_dphidr);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->phi,&arr_phi);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->fusion,&arr_fusion);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->S,&arr_S);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(pres,&arr_pres);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->dTdPs,&arr_dTdPs);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->dTdrs,&arr_dTdrs);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->cp,&arr_cp);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->temp,&arr_temp);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->alpha,&arr_alpha);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->cond,&arr_cond);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->visc,&arr_visc);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->rho,&arr_rho);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(M->dPdr_b,&arr_dPdr_b);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->cp_mix,&arr_cp_mix);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->dTdrs_mix,&arr_dTdrs_mix);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->liquidus_temp,&arr_liquidus_temp);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->fusion_rho,&arr_fusion_rho);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->fusion_temp,&arr_fusion_temp);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->nu,&arr_nu);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->gsuper,&arr_gsuper);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->dSdr,&arr_dSdr);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->kappah,&arr_kappah);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->Etot,&arr_Etot);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(M->area_b,&arr_area_b);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(M->mix_b,&arr_mix_b);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->Jcond,&arr_Jcond);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->Jconv,&arr_Jconv);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->Jheat,&arr_Jheat);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->Jtot,&arr_Jtot);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->Jmix,&arr_Jmix);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->Jgrav,&arr_Jgrav);CHKERRQ(ierr);
    ierr = VecRestoreArray(S->Jmass,&arr_Jmass);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->solidus_temp,&arr_solidus_temp);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(S->phi_s,&arr_phi_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/* this is still needed for fusion curve derivative */
PetscErrorCode d_dr( Ctx *E, Vec in_s, Vec out_b )
{

    /* use staggered nodes to compute spatial derivative at basic
       (internal) nodes.  This is 2nd order accurate for a mesh with
       constant spacing */

    PetscErrorCode ierr;
    PetscInt i,ilo_s,ihi_s,ilo,ihi;
    PetscScalar dr,*arr_out_b;
    const PetscScalar *arr_in_s;

    PetscFunctionBeginUser;
    dr = E->mesh.dx_s;

    ierr = VecGetOwnershipRange(in_s,&ilo_s,&ihi_s);CHKERRQ(ierr);
    ilo = ilo_s ? 1 : ilo_s; // TODO: check 
    ihi = ihi_s; // TODO: check
    ierr = VecGetArrayRead(in_s,&arr_in_s);CHKERRQ(ierr);
    ierr = VecGetArray(out_b,&arr_out_b);CHKERRQ(ierr);
    for(i=ilo; i<ihi; i++){
        arr_out_b[i] = 1.0/dr * ( arr_in_s[i]-arr_in_s[i-1] );
    }
    ierr = VecRestoreArrayRead(in_s,&arr_in_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(out_b,&arr_out_b);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscScalar average( PetscScalar a, PetscScalar b )
{
    PetscScalar out;

    out = 0.5 * (a+b);

    return out;
}

PetscScalar combine_matprop( PetscScalar weight, PetscScalar mat1, PetscScalar mat2 )
{

    PetscScalar out;

    out = weight * mat1 + (1.0-weight) * mat2;

    return out;

}

void set_interp2d( const char * filename, Interp2d *interp )
{
#if (defined VERBOSE)
    printf("set_interp2d:\n");
#endif

    FILE *fp;
    size_t i=0, j=0, k=0;
    char string[100];
    PetscScalar x, y, z;
    PetscScalar xa[NX], ya[NY], za[NX*NY];
    PetscScalar xscale, yscale, zscale;

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
    /*for (i=0; i<NX; i++ ){
        printf("%d %f\n", i, xa[i]);
    }*/

    /*for (j=0; j<NY; j++ ){
        printf("%d %f\n", j, ya[j]);
    }*/

    interp->xmin= xa[0];
    interp->xmax= xa[NX-1];
    interp->ymin= ya[0];
    interp->ymax= ya[NY-1];

    interp->interp = spline;
    interp->xacc = xacc;
    interp->yacc = yacc;

}

void set_interp1d( const char * filename, Interp1d *interp, int n )
{
#if (defined VERBOSE)
    printf("set_interp1d:\n");
#endif

    FILE *fp;
    size_t i=0;
    char string[100];
    PetscScalar x, y, xscale, yscale;
    PetscScalar xa[n], ya[n];

    if (sizeof(PetscScalar) != sizeof(double)){
      perror("PetscScalar must be double to use the dataio functions here");
      exit(-1);
    }

    /* linear interpolation */
    const gsl_interp_type *T = gsl_interp_linear;
    gsl_interp *interpolation = gsl_interp_alloc( T, n );
    gsl_interp_accel *acc = gsl_interp_accel_alloc();

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

    /* for debugging */
    /*for (i=0; i<NLS; i++ ){
        printf("%lu %f %f\n", i, xa[i], ya[i]);
    }*/

    gsl_interp_init( interpolation, xa, ya, n );

    // I think this is the correct way of copying an array
    memmove( interp->xa, xa, sizeof interp->xa );
    interp->xmin= xa[0];
    interp->xmax= xa[n-1];
    memmove( interp->ya, ya, sizeof interp->ya );
    interp->ymin= ya[0];
    interp->ymax= ya[n-1];

    interp->interp = interpolation;
    interp->acc = acc;

}
