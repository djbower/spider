#include "ctx.h"

/* Set up the Context */
PetscErrorCode setup_ctx(Ctx* ctx)
{
  PetscErrorCode ierr;
  PetscInt       i;
  Vec            S_s;

  PetscFunctionBeginUser;

  /* Initialize context with parameters (most parameters are constants in global_defs.h, though) */
  set_lookups(ctx);

  /* Set up a parallel structured grid as DMComposite with two included DMDAs
     This is used to define vectors which hold the solution. The included 
     DMDAs are used to create vectors which hold various properties
     living on the same primal and staggered nodes.
  */
  const PetscInt stencilWidth = 1; //TODO: check that this is appropriate
  const PetscInt dof = 1;
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,NUMPTS,dof,stencilWidth,NULL,&ctx->da_b);CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,NUMPTSS,dof,stencilWidth,NULL,&ctx->da_s);CHKERRQ(ierr);

#if (defined DEBUGOUTPUT)
  {
    PetscInt ilo_b,w_b,ihi_b,ilo_s,w_s,ihi_s;
    PetscMPIInt rank;

    ierr = DMDAGetCorners(ctx->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ierr = DMDAGetCorners(ctx->da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] index ranges b:%d<=i<%d, s:%d<=i<%d\n",rank,ilo_b,ihi_b,ilo_s,ihi_s);CHKERRQ(ierr);
  }
#endif


  // TODO (maybe): PETSc-fy this more by getting rid of the NUMPTS and NUMPTSS parameters and instead letting the DMDAs themselves define this information (hence allowing more command-line control)

  /* Continue to initialize context with distributed data */

  for (i=0;i<NUMSOLUTIONVECS;++i){
    ierr = DMCreateGlobalVector(ctx->da_b,&(ctx->solution.solutionVecs[i]));CHKERRQ(ierr);
  }
  ctx->solution.alpha               = ctx->solution.solutionVecs[0];
  ctx->solution.alpha_mix           = ctx->solution.solutionVecs[1];
  ctx->solution.cond                = ctx->solution.solutionVecs[2];
  ctx->solution.cp                  = ctx->solution.solutionVecs[3];
  ctx->solution.cp_mix              = ctx->solution.solutionVecs[4];  // TI
  ctx->solution.dfusdr              = ctx->solution.solutionVecs[5];  // TI
  ctx->solution.dfusdr_temp         = ctx->solution.solutionVecs[6];  // TI
  ctx->solution.dphidr              = ctx->solution.solutionVecs[7];
  ctx->solution.dSdr                = ctx->solution.solutionVecs[8];  // TI
  ctx->solution.dTdPs               = ctx->solution.solutionVecs[9];
  ctx->solution.dTdPs_mix           = ctx->solution.solutionVecs[10]; // TI
  ctx->solution.dTdrs               = ctx->solution.solutionVecs[11];
  ctx->solution.dTdrs_mix           = ctx->solution.solutionVecs[12]; // TI
  ctx->solution.Etot                = ctx->solution.solutionVecs[13];
  ctx->solution.fusion              = ctx->solution.solutionVecs[14]; // TI
  ctx->solution.fusion_curve        = ctx->solution.solutionVecs[15]; // TI
  ctx->solution.fusion_curve_temp   = ctx->solution.solutionVecs[16]; // TI
  ctx->solution.fusion_rho          = ctx->solution.solutionVecs[17]; // TI
  ctx->solution.fusion_temp         = ctx->solution.solutionVecs[18]; // TI
  ctx->solution.gsuper              = ctx->solution.solutionVecs[19];
  ctx->solution.Jcond               = ctx->solution.solutionVecs[20];
  ctx->solution.Jconv               = ctx->solution.solutionVecs[21];
  ctx->solution.Jgrav               = ctx->solution.solutionVecs[22];
  ctx->solution.Jheat               = ctx->solution.solutionVecs[23];
  ctx->solution.Jmass               = ctx->solution.solutionVecs[24];
  ctx->solution.Jmix                = ctx->solution.solutionVecs[25];
  ctx->solution.Jtot                = ctx->solution.solutionVecs[26];
  ctx->solution.kappah              = ctx->solution.solutionVecs[27];
  ctx->solution.liquidus            = ctx->solution.solutionVecs[28]; // TI
  ctx->solution.liquidus_rho        = ctx->solution.solutionVecs[29]; // TI
  ctx->solution.liquidus_temp       = ctx->solution.solutionVecs[30]; // TI
  ctx->solution.nu                  = ctx->solution.solutionVecs[31];
  ctx->solution.phi                 = ctx->solution.solutionVecs[32];
  ctx->solution.rho                 = ctx->solution.solutionVecs[33];
  ctx->solution.S                   = ctx->solution.solutionVecs[34];
  ctx->solution.solidus             = ctx->solution.solutionVecs[35]; // TI
  ctx->solution.solidus_rho         = ctx->solution.solutionVecs[36]; // TI
  ctx->solution.solidus_temp        = ctx->solution.solutionVecs[37]; // TI
  ctx->solution.temp                = ctx->solution.solutionVecs[38];
  ctx->solution.visc                = ctx->solution.solutionVecs[39];

  ierr = DMCreateLocalVector(ctx->da_b,&ctx->work_local_b);CHKERRQ(ierr);

  for (i=0;i<NUMSOLUTIONVECSS;++i){
    ierr = DMCreateGlobalVector(ctx->da_s,&(ctx->solution.solutionVecsS[i]));CHKERRQ(ierr);
  }
  ctx->solution.fusion_s            = ctx->solution.solutionVecsS[0]; // TI
  ctx->solution.fusion_curve_s      = ctx->solution.solutionVecsS[1]; // TI
  ctx->solution.fusion_curve_temp_s = ctx->solution.solutionVecsS[2]; // TI
  ctx->solution.fusion_temp_s       = ctx->solution.solutionVecsS[3]; // TI
  ctx->solution.lhs_s               = ctx->solution.solutionVecsS[4];
  ctx->solution.liquidus_rho_s      = ctx->solution.solutionVecsS[5]; // TI
  ctx->solution.liquidus_s          = ctx->solution.solutionVecsS[6]; // TI
  ctx->solution.liquidus_temp_s     = ctx->solution.solutionVecsS[7]; // TI
  ctx->solution.phi_s               = ctx->solution.solutionVecsS[8];
  ctx->solution.rhs_s               = ctx->solution.solutionVecsS[9]; 
  ctx->solution.rho_s               = ctx->solution.solutionVecsS[10];
  ctx->solution.S_s                 = ctx->solution.solutionVecsS[11];
  ctx->solution.solidus_s           = ctx->solution.solutionVecsS[12]; // TI
  ctx->solution.solidus_rho_s       = ctx->solution.solutionVecsS[13]; // TI
  ctx->solution.solidus_temp_s      = ctx->solution.solutionVecsS[14]; // TI
  ctx->solution.temp_s              = ctx->solution.solutionVecsS[15];

  ierr = DMCreateLocalVector(ctx->da_s,&ctx->work_local_s);CHKERRQ(ierr);

  set_mesh(ctx);

  set_time_independent(ctx);

  /* Create a solution vector (S at the staggered nodes) and fill with an initial condition 
      Note that this is inside the Ctx 
  */
  ierr = DMCreateGlobalVector(ctx->da_s,&(ctx->solution.S_s));CHKERRQ(ierr); 
  S_s = ctx->solution.S_s;
  set_initial_condition(ctx);CHKERRQ(ierr); 

  PetscFunctionReturn(0);
}

PetscErrorCode destroy_ctx(Ctx* ctx)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBeginUser;

  /* Destroy data allocated in Ctx */
  free_memory_interp(ctx);
  ierr = DMDestroy(&ctx->da_s);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_b);CHKERRQ(ierr);
  for (i=0;i<NUMMESHVECS;++i){
    ierr = VecDestroy(&ctx->mesh.meshVecs[i]);CHKERRQ(ierr);
  }
  for (i=0;i<NUMMESHVECSS;++i){
    ierr = VecDestroy(&ctx->mesh.meshVecsS[i]);CHKERRQ(ierr);
  }
  for (i=0;i<NUMSOLUTIONVECS;++i){
    ierr = VecDestroy(&ctx->solution.solutionVecs[i]);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&ctx->work_local_b);CHKERRQ(ierr);
  for (i=0;i<NUMSOLUTIONVECSS;++i){
    ierr = VecDestroy(&ctx->solution.solutionVecsS[i]);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&ctx->work_local_s);CHKERRQ(ierr);

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
    gsl_interp_free( I->interp );
    gsl_interp_accel_free( I->acc );
}

void free_interp2d( Interp2d *I )
{
    gsl_spline2d_free( I->interp );
    gsl_interp_accel_free( I->xacc );
    gsl_interp_accel_free( I->yacc );
}

PetscErrorCode free_memory_interp( Ctx *E )
{
    PetscFunctionBeginUser;
    /* free memory allocated by interpolation functions */

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
    PetscErrorCode ierr;
    PetscInt       i,ilo_b,ihi_b,ilo_s,ihi_s,w_s,w_b;
    DM             da_s=E->da_s, da_b=E->da_b;
    Vec            pres_b,pres_s;
    PetscScalar    z,*arr_liquidus,*arr_liquidus_rho,*arr_liquidus_temp,*arr_liquidus_s,*arr_liquidus_rho_s,*arr_liquidus_temp_s;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Interp1d       *I;
    Solution       *S;
    Interp2d       *IR, *IT;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_liquidus:\n");CHKERRQ(ierr);
#endif

    S = &E->solution;
    I = &E->melt_prop.liquidus;
    IR = &E->melt_prop.rho;
    IT = &E->melt_prop.temp;

    pres_b = E->mesh.pressure_b;
    pres_s = E->mesh.pressure_s;


    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    /* basic nodes */
    ierr = DMDAVecGetArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->liquidus,&arr_liquidus);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->liquidus_temp,&arr_liquidus_temp);CHKERRQ(ierr);
    for(i=ilo_b;i<ihi_b;++i){
        z = get_val1d( I, arr_pres_b[i] );
        arr_liquidus[i] = z;
        arr_liquidus_rho[i] = get_val2d( IR, arr_pres_b[i], z );
        arr_liquidus_temp[i] = get_val2d( IT, arr_pres_b[i], z );
    }
    ierr = DMDAVecRestoreArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->liquidus,&arr_liquidus);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->liquidus_temp,&arr_liquidus_temp);CHKERRQ(ierr);


    /* staggered nodes */
    /* need in order to compute dfus/dr at basic internal nodes */ 
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        z = get_val1d( I, arr_pres_s[i] );
        arr_liquidus_s[i] = z;
        arr_liquidus_rho_s[i] = get_val2d( IR, arr_pres_s[i], z );
        arr_liquidus_temp_s[i] = get_val2d( IT, arr_pres_s[i], z );
    }
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_solidus( Ctx *E )
{
    PetscErrorCode ierr;
    PetscInt       i,ilo_b,ihi_b,w_b,ilo_s,ihi_s,w_s;
    DM             da_s=E->da_s,da_b=E->da_b;
    Vec            pres_b,pres_s;
    PetscScalar    z,*arr_solidus,*arr_solidus_rho,*arr_solidus_temp,*arr_solidus_s,*arr_solidus_rho_s,*arr_solidus_temp_s;
    const PetscScalar *arr_pres_b,*arr_pres_s;
    Interp1d       *I;
    Solution       *S;
    Interp2d       *IR, *IT;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    printf("set_solidus:\n");
#endif
    S = &E->solution;
    I = &E->solid_prop.solidus;
    IR = &E->solid_prop.rho;
    IT = &E->solid_prop.temp;

    pres_b = E->mesh.pressure_b;
    pres_s = E->mesh.pressure_s;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    /* basic nodes */
    ierr = DMDAVecGetArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,S->solidus_temp,&arr_solidus_temp);CHKERRQ(ierr);
    for(i=ilo_b;i<ihi_b;++i){
        z = get_val1d( I, arr_pres_b[i] );
        arr_solidus[i] = z;
        arr_solidus_rho[i] = get_val2d( IR, arr_pres_b[i], z );
        arr_solidus_temp[i] = get_val2d( IT, arr_pres_b[i], z );
    }
    ierr = DMDAVecRestoreArrayRead(da_b,pres_b,&arr_pres_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,S->solidus_temp,&arr_solidus_temp);CHKERRQ(ierr);


    /* staggered nodes */
    /* need in order to compute dfus/dr at basic internal nodes */ 
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        z = get_val1d( I, arr_pres_s[i] );
        arr_solidus_s[i] = z;
        arr_solidus_rho_s[i] = get_val2d( IR, arr_pres_s[i], z );
        arr_solidus_temp_s[i] = get_val2d( IT, arr_pres_s[i], z );
    }
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    
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
    ierr = VecWAXPY(S->fusion_s,-1.0,S->solidus_s,S->liquidus_s);CHKERRQ(ierr);
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
    ierr = VecWAXPY(S->fusion_curve,0.5,S->fusion,S->solidus);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_curve_temp,0.5,S->fusion_temp,S->solidus_temp);CHKERRQ(ierr);

    /* staggered nodes */
    ierr = VecWAXPY(S->fusion_curve_s,0.5,S->fusion_s,S->solidus_s);CHKERRQ(ierr);
    ierr = VecWAXPY(S->fusion_curve_temp_s,0.5,S->fusion_temp_s,S->solidus_temp_s);CHKERRQ(ierr);

    d_dr( E, S->fusion_curve_s,      S->dfusdr      );
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
    {
      PetscMPIInt rank;
      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] set_mixed_phase\n");CHKERRQ(ierr);
    }
#endif

    M = &E->mesh;
    S = &E->solution;

    /* dTdrs_mix = -dfusdr * (fusion_temp/fusion) + dfusdr_temp */
    ierr = VecPointwiseDivide(S->dTdrs_mix,S->fusion_temp,S->fusion);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->dTdrs_mix,S->dTdrs_mix,S->dfusdr);CHKERRQ(ierr);
    ierr = VecScale( S->dTdrs_mix, -1.0 );
    ierr = VecAXPY(S->dTdrs_mix,1.0,S->dfusdr_temp);CHKERRQ(ierr);

    ierr = VecPointwiseMult(S->dTdPs_mix,S->dTdrs_mix,M->dPdr_b);CHKERRQ(ierr);

    ierr = VecPointwiseDivide(S->cp_mix,S->fusion,S->fusion_temp);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->cp_mix,S->cp_mix,S->fusion_curve_temp);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode set_core_cooling( Ctx *E )
{
    PetscErrorCode    ierr;
    const PetscInt    ix0 = 0;        // index of first node
    const PetscInt    ix  = NUMPTS-1; // index of last basic node
    const PetscInt    ix2 = NUMPTS-2; // index of penultimate basic node
    PetscScalar       fac,radi,rado,vol,area1,area2;
    PetscMPIInt       rank, size;

    Mesh              *M;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_core_cooling:\n");CHKERRQ(ierr);
#endif
    M = &E->mesh;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

    /* Assume that rank 0 contains the first point, and grab the 
       radius at the first point to broadcast */
    if (!rank) {
        ierr = VecGetValues( M->radius_b,1,&ix0,&rado);CHKERRQ(ierr);
    } else {
        rado = 0.0; // just for warning suppression
    }
    ierr = MPI_Bcast(&rado,1,MPIU_SCALAR,0,PETSC_COMM_WORLD);CHKERRQ(ierr);

    /* Assume that the last rank contains the last two points */
    if (rank == size-1 ){
        ierr = VecGetValues( M->area_b,1,&ix,&area1);CHKERRQ(ierr);
        ierr = VecGetValues( M->area_b,1,&ix2,&area2);CHKERRQ(ierr);
        ierr = VecGetValues( M->radius_b,1,&ix,&radi);CHKERRQ(ierr);
        ierr = VecGetValues( M->volume_s,1,&ix2,&vol);CHKERRQ(ierr);
        fac = 4.0*M_PI*PetscPowScalar(rado, 3.0) * vol;
        fac *= RHO_CMB * CP_CMB;
        fac /= CP_CORE * MCORE * TFAC_CORE_AVG;
        fac += 1.0;
        fac = area2 / ( area1 * fac );
        E->BC_BOT_FAC = fac;
    }

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

    set_core_cooling( E );

    PetscFunctionReturn(0);
}


/* end of time-independent quantities */

//NOTE: This function MUST be called first when computing the RHS
PetscErrorCode set_capacitance( Ctx *E, Vec S_in )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    DM                da_s=E->da_s;
    Mesh              *M;
    Solution          *S;
    Vec               pres_s;
    Interp2d          *IRM, *ITM, *IRS, *ITS;
    const PetscScalar *arr_S_s,*arr_pres_s,*arr_liquidus_rho_s,*arr_solidus_rho_s,*arr_liquidus_temp_s,*arr_solidus_temp_s;
    PetscScalar       *arr_phi_s,*arr_rho_s,*arr_temp_s;

    PetscFunctionBeginUser;
    M = &E->mesh;
    S = &E->solution;
    pres_s = M->pressure_s;

    IRM = &E->melt_prop.rho;
    ITM = &E->melt_prop.temp;
    IRS = &E->solid_prop.rho;
    ITS = &E->solid_prop.temp;

    ierr = VecCopy(S_in,S->S_s);CHKERRQ(ierr);

    // S->phi_s = (S->S_s - S->solidus_s)/S->fusion_s
    ierr = VecWAXPY(S->phi_s,-1.0,S->solidus_s,S->S_s);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(S->phi_s,S->phi_s,S->fusion_s);CHKERRQ(ierr);


    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = DMDAVecGetArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);
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
    ierr = DMDAVecRestoreArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);

    // S->lhs_si = M->volume_si * S->rho_si * S->temp_si;
    ierr = VecPointwiseMult(S->lhs_s,M->volume_s,S->rho_s);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->lhs_s,S->lhs_s,S->temp_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscScalar zmap( PetscScalar z )
{
    /* for skewed viscosity profile */

    PetscScalar fac, fwt, zmap, shp;

    shp = SHAPE_TRANSITION;

    fac = PetscSqrtScalar(PetscSqr(shp*z)+1.0);

    zmap = shp*PetscSqr(z)+z*fac+1.0/shp \
        *PetscLogScalar(shp*z+fac);
       zmap *= 0.5;
    fwt = 0.5*(1.0+PetscTanhScalar( zmap ));

    return fwt;
}

static PetscScalar viscosity_mix_skew( PetscScalar meltf )
{
    /* skewed viscosity in mixed phase region */

    PetscScalar fsol, fwt, fwt2, fwt3, z;

    fsol = 1.0 - meltf;
    z = (fsol-F_THRESHOLD) / DF_TRANSITION;
    fwt = zmap( z );

    z = (1.0-F_THRESHOLD) / DF_TRANSITION;
    fwt2 = zmap( z );

    z = -F_THRESHOLD / DF_TRANSITION;
    fwt3 = zmap( z );

    fwt = -(fwt-fwt3)/(fwt3-fwt2);

    return fwt;
}

static PetscScalar viscosity_mix( PetscScalar meltf )
{
    PetscScalar fwt, log10visc, visc;

    /* simple hyperbolic tan to get up and running */
    /* need to code up function with skew */
    /*arg = (1.0-meltf-F_THRESHOLD) / DF_TRANSITION;
    fwt = 0.5 * ( 1.0 + tanh(arg) ); */

    fwt = viscosity_mix_skew( meltf );

    log10visc = fwt * LOG10VISC_SOL + (1.0 - fwt) * LOG10VISC_MEL;
    visc = PetscPowScalar( 10.0, log10visc );

    return visc;

}

PetscErrorCode set_matprop_and_flux( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_b,ihi_b,w_b,ilo,ihi;
    DM                da_s=E->da_s, da_b=E->da_b;
    Vec               pres;
    PetscScalar       dr;
    PetscScalar       *arr_S,*arr_phi,*arr_nu,*arr_gsuper,*arr_kappah,*arr_Etot,*arr_dSdr,*arr_dTdPs,*arr_dTdrs,*arr_dphidr,*arr_alpha,*arr_temp,*arr_cp,*arr_cond,*arr_visc,*arr_rho,*arr_Jcond,*arr_Jconv,*arr_Jtot,*arr_Jheat,*arr_Jmix,*arr_Jgrav,*arr_Jmass;
    const PetscScalar *arr_S_s,*arr_phi_s,*arr_solidus,*arr_fusion,*arr_pres,*arr_area_b,*arr_dPdr_b,*arr_liquidus_rho,*arr_solidus_rho,*arr_cp_mix,*arr_dTdrs_mix,*arr_liquidus_temp,*arr_solidus_temp,*arr_fusion_rho,*arr_fusion_temp,*arr_mix_b;
    Mesh              *M;
    Solution          *S;
    Vec               S_s_local,phi_s_local;

    PetscFunctionBeginUser;
    M = &E->mesh;
    S = &E->solution;

    dr = M->dx_s;
    pres = M->pressure_b;

    /* Create some local vectors for the staggered grid*/
    ierr = DMCreateLocalVector(da_s,&S_s_local);CHKERRQ(ierr);
    ierr = DMCreateLocalVector(da_s,&phi_s_local);CHKERRQ(ierr);

    /* loop over all basic internal nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ilo = ilo_b == 0      ? 1          : ilo_b;
    ihi = ihi_b == NUMPTS ? NUMPTS - 1 : ihi_b;

    ierr = DMDAVecGetArray(    da_b,S->dphidr,&arr_dphidr);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da_s,S->S_s,INSERT_VALUES,S_s_local);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_s,S->S_s,INSERT_VALUES,S_s_local);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S_s_local,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->phi,&arr_phi);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion,&arr_fusion);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->S,&arr_S);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,pres,&arr_pres);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->dTdPs,&arr_dTdPs);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->dTdrs,&arr_dTdrs);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cp,&arr_cp);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->temp,&arr_temp);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->alpha,&arr_alpha);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cond,&arr_cond);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->visc,&arr_visc);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->rho,&arr_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->dPdr_b,&arr_dPdr_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->cp_mix,&arr_cp_mix);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dTdrs_mix,&arr_dTdrs_mix);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus_temp,&arr_liquidus_temp);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion_rho,&arr_fusion_rho);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion_temp,&arr_fusion_temp);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->nu,&arr_nu);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->gsuper,&arr_gsuper);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->dSdr,&arr_dSdr);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->kappah,&arr_kappah);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Etot,&arr_Etot);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->area_b,&arr_area_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->mix_b,&arr_mix_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jcond,&arr_Jcond);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jconv,&arr_Jconv);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jheat,&arr_Jheat);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jtot,&arr_Jtot);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jmix,&arr_Jmix);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jgrav,&arr_Jgrav);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jmass,&arr_Jmass);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus_temp,&arr_solidus_temp);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(da_s,S->phi_s,INSERT_VALUES,phi_s_local);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_s,S->phi_s,INSERT_VALUES,phi_s_local);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,phi_s_local,&arr_phi_s);CHKERRQ(ierr);

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
        arr_visc[i] = PetscPowScalar( 10.0, LOG10VISC_SOL );
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
      }

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
    ierr = DMDAVecRestoreArray(    da_b,S->dphidr,&arr_dphidr);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->phi,&arr_phi);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion,&arr_fusion);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->S,&arr_S);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,pres,&arr_pres);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->dTdPs,&arr_dTdPs);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->dTdrs,&arr_dTdrs);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cp,&arr_cp);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->temp,&arr_temp);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->alpha,&arr_alpha);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cond,&arr_cond);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->visc,&arr_visc);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->rho,&arr_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->dPdr_b,&arr_dPdr_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus_rho,&arr_liquidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus_rho,&arr_solidus_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->cp_mix,&arr_cp_mix);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dTdrs_mix,&arr_dTdrs_mix);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus_temp,&arr_liquidus_temp);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus,&arr_solidus);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion_rho,&arr_fusion_rho);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion_temp,&arr_fusion_temp);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->nu,&arr_nu);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->gsuper,&arr_gsuper);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->dSdr,&arr_dSdr);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->kappah,&arr_kappah);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Etot,&arr_Etot);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->area_b,&arr_area_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->mix_b,&arr_mix_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jcond,&arr_Jcond);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jconv,&arr_Jconv);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jheat,&arr_Jheat);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jtot,&arr_Jtot);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jmix,&arr_Jmix);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jgrav,&arr_Jgrav);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jmass,&arr_Jmass);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus_temp,&arr_solidus_temp);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);

    ierr = VecDestroy(&S_s_local);CHKERRQ(ierr);
    ierr = VecDestroy(&phi_s_local);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/* this is still needed for fusion curve derivative */
PetscErrorCode d_dr( Ctx *E, Vec in_s, Vec out_b )
{

    /* use staggered nodes to compute spatial derivative at basic
       (internal) nodes.  This is 2nd order accurate for a mesh with
       constant spacing */

    PetscErrorCode     ierr;
    PetscInt           i,ilo_b,ihi_b,w_b,ilo,ihi;
    DM                 da_s=E->da_s,da_b=E->da_b;
    PetscScalar        dr,*arr_out_b;
    const PetscScalar *arr_in_s;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    {
      PetscMPIInt rank;
      MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
      ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] d_dr\n");CHKERRQ(ierr);
    }
#endif
    dr = E->mesh.dx_s;

    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ilo = ilo_b==0      ? 1        : ilo_b; 
    ihi = ihi_b==NUMPTS ? NUMPTS-1 : ihi_b;

    // TODO: here and elsewhere, we are a little glib about assuming things abou the way the DA's partition things. We should introduce checks for any function which involves both DAs at once, that the expected ranges apply.

    /* Scatter to local vectors, since we may require potentially off-processor ghost values */
    ierr = DMGlobalToLocalBegin(da_s,in_s,INSERT_VALUES,E->work_local_s);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_s,in_s,INSERT_VALUES,E->work_local_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,E->work_local_s,&arr_in_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_b,out_b,&arr_out_b);CHKERRQ(ierr);
    for(i=ilo; i<ihi; i++){
        arr_out_b[i] = 1.0/dr * ( arr_in_s[i]-arr_in_s[i-1] );
    }
    ierr = DMDAVecRestoreArrayRead(da_s,E->work_local_s,&arr_in_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_b,out_b,&arr_out_b);CHKERRQ(ierr);

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

PetscErrorCode set_interp2d( const char * filename, Interp2d *interp )
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

    PetscFunctionReturn(0);
}

PetscErrorCode set_interp1d( const char * filename, Interp1d *interp, int n )
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

    // I think this is the correct way of copying an array
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
