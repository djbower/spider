#include "ctx.h"

/* Set up the Context */
PetscErrorCode setup_ctx(Ctx* ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,numpts_b,numpts_s; 
  PetscBool      set;

  PetscFunctionBeginUser;

  /* Initialize context with parameters (most parameters are constants in global_defs.h, though) */
  set_lookups(ctx);

  /* Check for a command line option -n to set the number of *staggered* points */
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&numpts_s,&set);CHKERRQ(ierr);
  if (set){
    numpts_b = numpts_s + 1;
  } else {
    numpts_b = NUMPTS_B_DEFAULT;
    numpts_s = NUMPTS_S_DEFAULT;
  }

  /* Set up a parallel structured grid as DMComposite with two included DMDAs
     This is used to define vectors which hold the solution. The included 
     DMDAs are used to create vectors which hold various properties
     living on the same primal and staggered nodes.
  */
  const PetscInt stencilWidth = 1; //TODO: check that this is appropriate
  const PetscInt dof = 1;
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,numpts_b,dof,stencilWidth,NULL,&ctx->da_b);CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,numpts_s,dof,stencilWidth,NULL,&ctx->da_s);CHKERRQ(ierr);

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

  /* Continue to initialize context with distributed data */

  for (i=0;i<NUMSOLUTIONVECS_B;++i){
    ierr = DMCreateGlobalVector(ctx->da_b,&(ctx->solution.solutionVecs_b[i]));CHKERRQ(ierr);
  }
  ctx->solution.alpha               = ctx->solution.solutionVecs_b[0];
  ctx->solution.alpha_mix           = ctx->solution.solutionVecs_b[1];
  ctx->solution.cond                = ctx->solution.solutionVecs_b[2];
  ctx->solution.cp                  = ctx->solution.solutionVecs_b[3];
  ctx->solution.cp_mix              = ctx->solution.solutionVecs_b[4];
  ctx->solution.dfusdr              = ctx->solution.solutionVecs_b[5];
  ctx->solution.dfusdr_temp         = ctx->solution.solutionVecs_b[6];
  ctx->solution.dSdr                = ctx->solution.solutionVecs_b[7];
  ctx->solution.dSliqdr             = ctx->solution.solutionVecs_b[8];
  ctx->solution.dSsoldr             = ctx->solution.solutionVecs_b[9];
  ctx->solution.dTdrs               = ctx->solution.solutionVecs_b[10];
  ctx->solution.dTdrs_mix           = ctx->solution.solutionVecs_b[11];
  ctx->solution.Etot                = ctx->solution.solutionVecs_b[12];
  ctx->solution.fusion              = ctx->solution.solutionVecs_b[13];
  ctx->solution.fusion_curve        = ctx->solution.solutionVecs_b[14];
  ctx->solution.fusion_curve_temp   = ctx->solution.solutionVecs_b[15];
  ctx->solution.fusion_rho          = ctx->solution.solutionVecs_b[16];
  ctx->solution.fusion_temp         = ctx->solution.solutionVecs_b[17];
  ctx->solution.gsuper              = ctx->solution.solutionVecs_b[18];
  ctx->solution.Jcond               = ctx->solution.solutionVecs_b[19];
  ctx->solution.Jconv               = ctx->solution.solutionVecs_b[20];
  ctx->solution.Jgrav               = ctx->solution.solutionVecs_b[21];
  ctx->solution.Jmix                = ctx->solution.solutionVecs_b[22];
  ctx->solution.Jtot                = ctx->solution.solutionVecs_b[23];
  ctx->solution.kappah              = ctx->solution.solutionVecs_b[24];
  ctx->solution.liquidus            = ctx->solution.solutionVecs_b[25];
  ctx->solution.liquidus_rho        = ctx->solution.solutionVecs_b[26];
  ctx->solution.liquidus_temp       = ctx->solution.solutionVecs_b[27];
  ctx->solution.nu                  = ctx->solution.solutionVecs_b[28];
  ctx->solution.phi                 = ctx->solution.solutionVecs_b[29];
  ctx->solution.rho                 = ctx->solution.solutionVecs_b[30];
  ctx->solution.S                   = ctx->solution.solutionVecs_b[31];
  ctx->solution.solidus             = ctx->solution.solutionVecs_b[32];
  ctx->solution.solidus_rho         = ctx->solution.solutionVecs_b[33];
  ctx->solution.solidus_temp        = ctx->solution.solutionVecs_b[34];
  ctx->solution.temp                = ctx->solution.solutionVecs_b[35];
  ctx->solution.visc                = ctx->solution.solutionVecs_b[36];

  ierr = DMCreateLocalVector(ctx->da_b,&ctx->work_local_b);CHKERRQ(ierr);

  for (i=0;i<NUMSOLUTIONVECS_S;++i){
    ierr = DMCreateGlobalVector(ctx->da_s,&(ctx->solution.solutionVecs_s[i]));CHKERRQ(ierr);
  }
  ctx->solution.cp_s                = ctx->solution.solutionVecs_s[0];
  ctx->solution.cp_mix_s            = ctx->solution.solutionVecs_s[1];
  ctx->solution.fusion_s            = ctx->solution.solutionVecs_s[2];
  ctx->solution.fusion_curve_s      = ctx->solution.solutionVecs_s[3];
  ctx->solution.fusion_curve_temp_s = ctx->solution.solutionVecs_s[4];
  ctx->solution.fusion_temp_s       = ctx->solution.solutionVecs_s[5];
  ctx->solution.lhs_s               = ctx->solution.solutionVecs_s[6];
  ctx->solution.liquidus_rho_s      = ctx->solution.solutionVecs_s[7];
  ctx->solution.liquidus_s          = ctx->solution.solutionVecs_s[8];
  ctx->solution.liquidus_temp_s     = ctx->solution.solutionVecs_s[9];
  ctx->solution.phi_s               = ctx->solution.solutionVecs_s[10];
  ctx->solution.rho_s               = ctx->solution.solutionVecs_s[11];
  ctx->solution.S_s                 = ctx->solution.solutionVecs_s[12];
  ctx->solution.solidus_s           = ctx->solution.solutionVecs_s[13];
  ctx->solution.solidus_rho_s       = ctx->solution.solutionVecs_s[14];
  ctx->solution.solidus_temp_s      = ctx->solution.solutionVecs_s[15];
  ctx->solution.temp_s              = ctx->solution.solutionVecs_s[16];

  set_mesh(ctx);

  set_d_dr( ctx );

  set_twophase(ctx);

  /* Obtain a command-line argument for S_init */
  ctx->S_init = SINIT_DEFAULT;
  ierr = PetscOptionsGetScalar(NULL,NULL,"-sinit",&ctx->S_init,NULL);CHKERRQ(ierr);

  /* initial condition is set in main.c */

  PetscFunctionReturn(0);
}

PetscErrorCode destroy_ctx(Ctx* ctx)
{
  PetscErrorCode ierr;
  PetscInt       i;

  PetscFunctionBeginUser;

  /* Destroy data allocated in Ctx */
  for (i=0;i<NUMMESHVECS_B;++i){
    ierr = VecDestroy(&ctx->mesh.meshVecs_b[i]);CHKERRQ(ierr);
  }
  for (i=0;i<NUMMESHVECS_S;++i){
    ierr = VecDestroy(&ctx->mesh.meshVecs_s[i]);CHKERRQ(ierr);
  }
  for (i=0;i<NUMSOLUTIONVECS_B;++i){
    ierr = VecDestroy(&ctx->solution.solutionVecs_b[i]);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&ctx->work_local_b);CHKERRQ(ierr);
  for (i=0;i<NUMSOLUTIONVECS_S;++i){
    ierr = VecDestroy(&ctx->solution.solutionVecs_s[i]);CHKERRQ(ierr);
  }

  ierr = MatDestroy(&ctx->qty_at_b); CHKERRQ(ierr);
  ierr = MatDestroy(&ctx->ddr_at_b); CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_s);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx->da_b);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

//NOTE: This function MUST be called first when computing the RHS
PetscErrorCode set_capacitance( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    DM                da_s=E->da_s;
    Mesh              *M;
    Lookup            *L;
    Solution          *S;
    Vec               pres_s;
    const PetscScalar *arr_pres_s, *arr_liquidus_rho_s, *arr_solidus_rho_s, *arr_liquidus_temp_s, *arr_solidus_temp_s, *arr_S_s, *arr_liquidus_s, *arr_solidus_s, *arr_fusion_s, *arr_cp_mix_s;
    PetscScalar       *arr_phi_s, *arr_rho_s, *arr_temp_s, *arr_cp_s;
    PetscScalar       gphi, fwtl, fwts, fwtm;

    PetscFunctionBeginUser;

    M = &E->mesh;
    S = &E->solution;
    pres_s = M->pressure_s;

    ierr = VecWAXPY(S->phi_s,-1.0,S->solidus_s,S->S_s);CHKERRQ(ierr);
    ierr = VecPointwiseDivide(S->phi_s,S->phi_s,S->fusion_s);CHKERRQ(ierr);

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = DMDAVecGetArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fusion_s,&arr_fusion_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->cp_s,&arr_cp_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->cp_mix_s,&arr_cp_mix_s);CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){

      /* generalised melt fraction for smoothing */
      gphi = arr_phi_s[i];
      /* weights for each phase */
      fwtl = tanh_weight( gphi, 1.0, SWIDTH );
      fwts = 1.0 - tanh_weight( gphi, 0.0, SWIDTH );
      fwtm = (1.0-fwtl) * (1.0-fwts);

      /* melt fraction, also truncated here [0,1] */
      if (arr_phi_s[i] > 1.0){
        /* superliquidus */
        arr_phi_s[i] = 1.0;
      }
      if (arr_phi_s[i] < 0.0){
        /* subsolidus */
        arr_phi_s[i] = 0.0;
      }

      /////////////////
      /* mixed phase */
      /////////////////
      /* density */
      arr_rho_s[i] = combine_matprop( arr_phi_s[i], 1.0/arr_liquidus_rho_s[i], 1.0/arr_solidus_rho_s[i] );
      arr_rho_s[i] = 1.0 / arr_rho_s[i];
      arr_rho_s[i] *= fwtm;
      /* temperature */
      arr_temp_s[i] = combine_matprop( arr_phi_s[i], arr_liquidus_temp_s[i], arr_solidus_temp_s[i] );
      arr_temp_s[i] *= fwtm;
      /* heat capacity */
      arr_cp_s[i] = arr_cp_mix_s[i];
      arr_cp_s[i] *= fwtm;

      ////////////////
      /* melt phase */
      ////////////////
      /* get melt properties */
      L = &E->melt_prop;
      /* density */
      arr_rho_s[i]   += fwtl * get_val2d( &L->rho, arr_pres_s[i], arr_S_s[i] );
      /* temperature */
      arr_temp_s[i]  += fwtl * get_val2d( &L->temp, arr_pres_s[i], arr_S_s[i] );
      /* heat capacity */
      arr_cp_s[i]    += fwtl * get_val2d( &L->cp, arr_pres_s[i], arr_S_s[i] );

      /////////////////
      /* solid phase */
      /////////////////
      /* get solid properties */
      L = &E->solid_prop;
      /* density */
      arr_rho_s[i]   += fwts * get_val2d( &L->rho, arr_pres_s[i], arr_S_s[i] );
      /* temperature */
      arr_temp_s[i]  += fwts * get_val2d( &L->temp, arr_pres_s[i], arr_S_s[i] );
      /* heat capacity */
      arr_cp_s[i]    += fwts * get_val2d( &L->cp, arr_pres_s[i], arr_S_s[i] );

    }

    ierr = DMDAVecRestoreArrayRead(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,pres_s,&arr_pres_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_s,&arr_liquidus_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_s,&arr_solidus_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fusion_s,&arr_fusion_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_rho_s,&arr_liquidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_rho_s,&arr_solidus_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->liquidus_temp_s,&arr_liquidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->solidus_temp_s,&arr_solidus_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->cp_mix_s,&arr_cp_mix_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->rho_s,&arr_rho_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->temp_s,&arr_temp_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->cp_s,&arr_cp_s);CHKERRQ(ierr);

    // S->lhs_si = M->volume_si * S->rho_si * S->temp_si;
    ierr = VecPointwiseMult(S->lhs_s,S->temp_s,S->rho_s);CHKERRQ(ierr);
    ierr = VecPointwiseMult(S->lhs_s,S->lhs_s,M->volume_s);CHKERRQ(ierr);

    /* clean up */

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

    fwt = tanh_weight( zmap, 0.0, 1.0 );

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

static PetscScalar log_viscosity_mix( PetscScalar meltf )
{
    PetscScalar fwt, lvisc;

    fwt = viscosity_mix_skew( meltf );

    lvisc = fwt * LOG10VISC_SOL + (1.0 - fwt) * LOG10VISC_MEL;

    return lvisc;

}

PetscErrorCode set_matprop_and_flux( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_b,ihi_b,w_b,ilo,ihi,numpts_b;
    DM                da_b=E->da_b;
    PetscScalar       *arr_phi, *arr_nu, *arr_gsuper, *arr_kappah, *arr_Etot, *arr_dTdrs, *arr_alpha, *arr_temp, *arr_cp, *arr_cond, *arr_visc, *arr_rho, *arr_Jcond, *arr_Jconv, *arr_Jtot, *arr_Jmix, *arr_Jgrav;
    const PetscScalar *arr_dSdr, *arr_S_b, *arr_dSliqdr, *arr_dSsoldr, *arr_solidus, *arr_fusion, *arr_pres, *arr_area_b, *arr_dPdr_b, *arr_liquidus, *arr_liquidus_rho, *arr_solidus_rho, *arr_cp_mix, *arr_dTdrs_mix, *arr_liquidus_temp, *arr_solidus_temp, *arr_fusion_rho, *arr_fusion_temp, *arr_mix_b;
    Mesh              *M;
    Solution          *S;
    PetscScalar       gphi, fwtl, fwts;

    PetscFunctionBeginUser;

    M = &E->mesh;
    S = &E->solution;

    /* loop over all basic internal nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ierr = DMDAGetInfo(da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ilo = ilo_b == 0        ? 1            : ilo_b;
    ihi = ihi_b == numpts_b ? numpts_b - 1 : ihi_b;

    ierr = DMDAVecGetArrayRead(da_b,S->dSdr,&arr_dSdr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->S,&arr_S_b); CHKERRQ(ierr);
    /* mesh quantities */
    ierr = DMDAVecGetArrayRead(da_b,M->area_b,&arr_area_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->dPdr_b,&arr_dPdr_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->mix_b,&arr_mix_b); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->pressure_b,&arr_pres); CHKERRQ(ierr);
    /* material properties */
    ierr = DMDAVecGetArray(    da_b,S->alpha,&arr_alpha); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cond,&arr_cond); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->cp,&arr_cp); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->cp_mix,&arr_cp_mix); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dSliqdr,&arr_dSliqdr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dSsoldr,&arr_dSsoldr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->dTdrs,&arr_dTdrs); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dTdrs_mix,&arr_dTdrs_mix); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion,&arr_fusion); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion_rho,&arr_fusion_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->fusion_temp,&arr_fusion_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->gsuper,&arr_gsuper); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->kappah,&arr_kappah); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus,&arr_liquidus); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus_rho,&arr_liquidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->liquidus_temp,&arr_liquidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->nu,&arr_nu); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->phi,&arr_phi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->rho,&arr_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus,&arr_solidus); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus_rho,&arr_solidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->solidus_temp,&arr_solidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->temp,&arr_temp); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->visc,&arr_visc); CHKERRQ(ierr);
    /* energy fluxes */
    ierr = DMDAVecGetArray(    da_b,S->Etot,&arr_Etot); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jcond,&arr_Jcond); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jconv,&arr_Jconv); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jgrav,&arr_Jgrav); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jmix,&arr_Jmix); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(    da_b,S->Jtot,&arr_Jtot); CHKERRQ(ierr);

    for(i=ilo; i<ihi; ++i){ // Note that these correspond to basic nodes being updated, but we also assume an ordering on the staggered nodes!

      /* melt fraction */
      arr_phi[i] = (arr_S_b[i] - arr_solidus[i]) / arr_fusion[i];

      /* generalised melt fraction for smoothing across liquidus */
      gphi = arr_phi[i];
      /* fwtl -> 1.0 for gphi > 1.0 */
      /* fwtl -> 0.0 for gphi < 1.0 */
      fwtl = tanh_weight( gphi, 1.0, SWIDTH );
      /* fwts -> 1.0 for gphi > 0.0 */
      /* fwts -> 0.0 for gphi < 0.0 */
      fwts = tanh_weight( gphi, 0.0, SWIDTH );

      /* now truncate melt fraction */
      if (arr_phi[i] > 1.0){
        /* superliquidus */
        arr_phi[i] = 1.0;
      }
      if (arr_phi[i] < 0.0){
        /* subsolidus */
        arr_phi[i] = 0.0;
      }

      /////////////////
      /* mixed phase */
      /////////////////

      /* TODO: confirm that this cannot break for any value of i.
         by inspection, I do not think it can */
      /* density */
      arr_rho[i] = combine_matprop( arr_phi[i], 1.0/arr_liquidus_rho[i], 1.0/arr_solidus_rho[i] );
      arr_rho[i] = 1.0 / arr_rho[i];
      /* adiabatic temperature gradient */
      arr_dTdrs[i] = arr_dTdrs_mix[i];
      /* heat capacity */
      arr_cp[i] = arr_cp_mix[i];
      /* temperature */
      arr_temp[i] = combine_matprop( arr_phi[i], arr_liquidus_temp[i], arr_solidus_temp[i] );
      /* thermal expansion coefficient */
      arr_alpha[i] = -arr_fusion_rho[i] / arr_fusion_temp[i] / arr_rho[i];
      /* thermal conductivity */
      arr_cond[i] = combine_matprop( arr_phi[i], COND_MEL, COND_SOL );
      /* log viscosity */
      arr_visc[i] = log_viscosity_mix( arr_phi[i] );

      ////////////////
      /* melt phase */
      ////////////////
      if (gphi > 0.5){

        /* get melt properties */
        Lookup *L = &E->melt_prop;
        /* density */
        arr_rho[i] *= 1.0 - fwtl;
        arr_rho[i] += fwtl * get_val2d( &L->rho, arr_pres[i], arr_S_b[i] );
        /* adiabatic temperature gradient */
        arr_dTdrs[i] *= 1.0 - fwtl;
        arr_dTdrs[i] += fwtl * arr_dPdr_b[i] * get_val2d( &L->dTdPs, arr_pres[i], arr_S_b[i] );
        /* heat capacity */
        arr_cp[i] *= 1.0 - fwtl;
        arr_cp[i] += fwtl * get_val2d( &L->cp, arr_pres[i], arr_S_b[i] );
        /* temperature */
        arr_temp[i] *= 1.0 - fwtl;
        arr_temp[i] += fwtl * get_val2d( &L->temp, arr_pres[i], arr_S_b[i] );
        /* thermal expansion coefficient */
        arr_alpha[i] *= 1.0 - fwtl;
        arr_alpha[i] += fwtl * get_val2d( &L->alpha, arr_pres[i], arr_S_b[i] );
        /* thermal conductivity */
        arr_cond[i] *= 1.0 - fwtl;
        arr_cond[i] += fwtl * COND_MEL;
        /* viscosity */
        arr_visc[i] *= 1.0 - fwtl;
        arr_visc[i] += fwtl * LOG10VISC_MEL;
      }

      else if (gphi <= 0.5){

        /* get solid properties */
        Lookup *L = &E->solid_prop;
        /* density */
        arr_rho[i] *= fwts;
        arr_rho[i] += (1.0-fwts) * get_val2d( &L->rho, arr_pres[i], arr_S_b[i] );
        /* adiabatic temperature gradient */
        arr_dTdrs[i] *= fwts;
        arr_dTdrs[i] += (1.0-fwts) * arr_dPdr_b[i] * get_val2d( &L->dTdPs, arr_pres[i], arr_S_b[i] );
        /* heat capacity */
        arr_cp[i] *= fwts;
        arr_cp[i] += (1.0-fwts) * get_val2d( &L->cp, arr_pres[i], arr_S_b[i] );
        /* temperature */
        arr_temp[i] *= fwts;
        arr_temp[i] += (1.0-fwts) * get_val2d( &L->temp, arr_pres[i], arr_S_b[i] );
        /* thermal expansion coefficient */
        arr_alpha[i] *= fwts;
        arr_alpha[i] += (1.0-fwts) * get_val2d( &L->alpha, arr_pres[i], arr_S_b[i] );
        /* thermal conductivity */
        arr_cond[i] *= fwts;
        arr_cond[i] += (1.0-fwts) * COND_SOL;
        /* viscosity */
        arr_visc[i] *= fwts;
        arr_visc[i] += (1.0-fwts) * LOG10VISC_SOL;
      }

      /* compute viscosity */
      arr_visc[i] = PetscPowScalar( 10.0, arr_visc[i] );

      /* other useful material properties */
      /* kinematic viscosity */
      arr_nu[i] = arr_visc[i] / arr_rho[i];

      /* gravity * super-adiabatic temperature gradient */
      arr_gsuper[i] = -GRAVITY * arr_temp[i] / arr_cp[i] * arr_dSdr[i];

      /* eddy diffusivity */
      PetscScalar kh, crit;
      crit = 81.0 * PetscPowScalar(arr_nu[i],2);
      crit /= 4.0 * arr_alpha[i] * PetscPowScalar(arr_mix_b[i],4);

      if( arr_gsuper[i] <= 0.0 ){
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
      arr_Jconv[i] = -arr_dSdr[i] * arr_kappah[i] * arr_rho[i] * arr_temp[i];

      //TODO: Need to clean up these declarations..
      PetscScalar dfus = arr_fusion[i];
      PetscScalar rho = arr_rho[i];
      PetscScalar rhol = arr_liquidus_rho[i];
      PetscScalar rhos = arr_solidus_rho[i];
      PetscScalar temp = arr_temp[i];

      PetscScalar pref = temp * dfus;

      /* convective mixing */
      /* these first two lines give F_Jmix (in the python script) */
      arr_Jmix[i] = arr_dSdr[i] - arr_phi[i] * arr_dSliqdr[i];
      arr_Jmix[i] += (arr_phi[i]-1.0) * arr_dSsoldr[i];
      arr_Jmix[i] *= -arr_kappah[i] * arr_rho[i] * arr_temp[i];

      /* blend together convection and mixing */
      if(gphi > 0.5){
        arr_Jmix[i] *= 1.0 - fwtl;
      }
      else if (gphi <= 0.5){
        arr_Jmix[i] *= fwts;
      }

      arr_Jtot[i] = arr_Jconv[i];
      arr_Jtot[i] += arr_Jmix[i];

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
      arr_Jgrav[i] *= pref * PetscPowScalar(GRAIN,2) * -GRAVITY * F;
      arr_Jgrav[i] /= PetscPowScalar(10.0, LOG10VISC_MEL);

      arr_Jtot[i] += arr_Jcond[i] + arr_Jgrav[i];

      arr_Etot[i] = arr_Jtot[i] * arr_area_b[i];

    }
 
    ierr = DMDAVecRestoreArrayRead(da_b,S->dSdr,&arr_dSdr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->S,&arr_S_b); CHKERRQ(ierr);
    /* mesh quantities */
    ierr = DMDAVecRestoreArrayRead(da_b,M->area_b,&arr_area_b); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->dPdr_b,&arr_dPdr_b); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->mix_b,&arr_mix_b); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->pressure_b,&arr_pres); CHKERRQ(ierr);
    /* material properties */
    ierr = DMDAVecRestoreArray(    da_b,S->alpha,&arr_alpha); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cond,&arr_cond); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->cp,&arr_cp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->cp_mix,&arr_cp_mix); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dSliqdr,&arr_dSliqdr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dSsoldr,&arr_dSsoldr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->dTdrs,&arr_dTdrs); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dTdrs_mix,&arr_dTdrs_mix); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion,&arr_fusion); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion_rho,&arr_fusion_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->fusion_temp,&arr_fusion_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->gsuper,&arr_gsuper); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->kappah,&arr_kappah); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus,&arr_liquidus); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus_rho,&arr_liquidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->liquidus_temp,&arr_liquidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->nu,&arr_nu); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->phi,&arr_phi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->rho,&arr_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus,&arr_solidus); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus_rho,&arr_solidus_rho); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->solidus_temp,&arr_solidus_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->temp,&arr_temp); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->visc,&arr_visc); CHKERRQ(ierr);
    /* energy fluxes */
    ierr = DMDAVecRestoreArray(    da_b,S->Etot,&arr_Etot); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jcond,&arr_Jcond); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jconv,&arr_Jconv); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jgrav,&arr_Jgrav); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jmix,&arr_Jmix); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(    da_b,S->Jtot,&arr_Jtot); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
