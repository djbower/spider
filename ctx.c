#include "ctx.h"
#include "atmosphere.h"
#include "lookup.h"
#include "mesh.h"
#include "twophase.h"
#include "util.h"

/* Set up the Context */
PetscErrorCode setup_ctx(Ctx* ctx)
{
  PetscErrorCode   ierr;
  PetscInt         i;
  Parameters const *P = &ctx->parameters;

  PetscFunctionBeginUser;

  /*
   **All** parameters/settings should be inside ctx->parameters . Here, we
     1. Initialize these to default values
     2. Obtain any inputs from the user at the command line and update
     3. Print out all parameters
     After this point, no custom command line options should be processed,
     and ctx->parameters should never change (meaning that we can and should use
     pointers-to-constants to refer to ctx->parameters).
    */
  InitializeParameters(&ctx->parameters); /* Note we use ctx->parameters, not P! */
  SetParametersFromOptions(&ctx->parameters);
  PrintParameters(&ctx->parameters);

  /* Set up a parallel structured grid as DMComposite with two included DMDAs
     This is used to define vectors which hold the solution. The included
     DMDAs are used to create vectors which hold various properties
     living on the same primal and staggered nodes.
  */
  const PetscInt stencilWidth = 1;
  const PetscInt dof = 1;
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,P->numpts_b,dof,stencilWidth,NULL,&ctx->da_b);CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,P->numpts_s,dof,stencilWidth,NULL,&ctx->da_s);CHKERRQ(ierr);

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
  ctx->solution.fwtl                = ctx->solution.solutionVecs_b[18]; // weight for liquid
  ctx->solution.fwts                = ctx->solution.solutionVecs_b[19]; // weight for solid
  ctx->solution.gphi                = ctx->solution.solutionVecs_b[20];
  ctx->solution.gsuper              = ctx->solution.solutionVecs_b[21];
  ctx->solution.Jcond               = ctx->solution.solutionVecs_b[22];
  ctx->solution.Jconv               = ctx->solution.solutionVecs_b[23];
  ctx->solution.Jgrav               = ctx->solution.solutionVecs_b[24];
  ctx->solution.Jmix                = ctx->solution.solutionVecs_b[25];
  ctx->solution.Jtot                = ctx->solution.solutionVecs_b[26];
  ctx->solution.kappah              = ctx->solution.solutionVecs_b[27];
  ctx->solution.liquidus            = ctx->solution.solutionVecs_b[28];
  ctx->solution.liquidus_rho        = ctx->solution.solutionVecs_b[29];
  ctx->solution.liquidus_temp       = ctx->solution.solutionVecs_b[30];
  ctx->solution.nu                  = ctx->solution.solutionVecs_b[31];
  ctx->solution.phi                 = ctx->solution.solutionVecs_b[32];
  ctx->solution.rho                 = ctx->solution.solutionVecs_b[33];
  ctx->solution.S                   = ctx->solution.solutionVecs_b[34];
  ctx->solution.solidus             = ctx->solution.solutionVecs_b[35];
  ctx->solution.solidus_rho         = ctx->solution.solutionVecs_b[36];
  ctx->solution.solidus_temp        = ctx->solution.solutionVecs_b[37];
  ctx->solution.temp                = ctx->solution.solutionVecs_b[38];
  ctx->solution.visc                = ctx->solution.solutionVecs_b[39];

  ierr = DMCreateLocalVector(ctx->da_b,&ctx->work_local_b);CHKERRQ(ierr);

  for (i=0;i<NUMSOLUTIONVECS_S;++i){
    ierr = DMCreateGlobalVector(ctx->da_s,&(ctx->solution.solutionVecs_s[i]));CHKERRQ(ierr);
  }
  ctx->solution.cp_s                = ctx->solution.solutionVecs_s[0];
  ctx->solution.cp_mix_s            = ctx->solution.solutionVecs_s[1];
  ctx->solution.dSdt_s              = ctx->solution.solutionVecs_s[2];
  ctx->solution.fusion_s            = ctx->solution.solutionVecs_s[3];
  ctx->solution.fusion_curve_s      = ctx->solution.solutionVecs_s[4];
  ctx->solution.fusion_curve_temp_s = ctx->solution.solutionVecs_s[5];
  ctx->solution.fusion_temp_s       = ctx->solution.solutionVecs_s[6];
  ctx->solution.fwtl_s              = ctx->solution.solutionVecs_s[7]; // weight for liquid
  ctx->solution.fwts_s              = ctx->solution.solutionVecs_s[8]; // weight for solid
  ctx->solution.gphi_s              = ctx->solution.solutionVecs_s[9];
  ctx->solution.Hradio_s            = ctx->solution.solutionVecs_s[10];
  ctx->solution.Htidal_s            = ctx->solution.solutionVecs_s[11];
  ctx->solution.Htot_s              = ctx->solution.solutionVecs_s[12];
  ctx->solution.lhs_s               = ctx->solution.solutionVecs_s[13];
  ctx->solution.liquidus_rho_s      = ctx->solution.solutionVecs_s[14];
  ctx->solution.liquidus_s          = ctx->solution.solutionVecs_s[15];
  ctx->solution.liquidus_temp_s     = ctx->solution.solutionVecs_s[16];
  ctx->solution.phi_s               = ctx->solution.solutionVecs_s[17];
  ctx->solution.rho_s               = ctx->solution.solutionVecs_s[18];
  ctx->solution.S_s                 = ctx->solution.solutionVecs_s[19];
  ctx->solution.solidus_s           = ctx->solution.solutionVecs_s[20];
  ctx->solution.solidus_rho_s       = ctx->solution.solutionVecs_s[21];
  ctx->solution.solidus_temp_s      = ctx->solution.solutionVecs_s[22];
  ctx->solution.temp_s              = ctx->solution.solutionVecs_s[23];

  set_mesh(ctx);

  set_d_dr( ctx );

  set_twophase(ctx);

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

PetscErrorCode set_lookups( Parameters *P )
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

    Constants *C = &P->constants;

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
      if (isStixrude2009) {
        ierr = PetscStrcpy(liquidusFilename,LIQUIDUS_STIXRUDE2009);CHKERRQ(ierr);
        ierr = MakeRelativePathAbsolute(liquidusFilename);CHKERRQ(ierr);
        ierr = PetscStrcpy(solidusFilename,SOLIDUS_STIXRUDE2009);CHKERRQ(ierr);
        ierr = MakeRelativePathAbsolute(solidusFilename);CHKERRQ(ierr);
      } else if (isAndrault2011 || isDefault) { /* Default */
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
    set_interp2d( alphaSolFilename, &P->solid_prop.alpha, C->PRESSURE, C->ENTROPY, C->ALPHA );
    set_interp2d( cpSolFilename, &P->solid_prop.cp, C->PRESSURE, C->ENTROPY, C->ENTROPY );
    set_interp2d( dtdpsSolFilename, &P->solid_prop.dTdPs, C->PRESSURE, C->ENTROPY, C->DTDP );
    set_interp2d( rhoSolFilename, &P->solid_prop.rho, C->PRESSURE, C->ENTROPY, C->DENSITY );
    set_interp2d( tempSolFilename, &P->solid_prop.temp, C->PRESSURE, C->ENTROPY, C->TEMP );
    /* const */
    // FIXME: remove this redundancy
    P->solid_prop.cond = P->cond_sol;
    P->solid_prop.log10visc = P->log10visc_sol;

    /* melt lookups */
    /* 2d */
    set_interp2d( alphaMelFilename, &P->melt_prop.alpha, C->PRESSURE, C->ENTROPY, C->ALPHA );
    set_interp2d( cpMelFilename, &P->melt_prop.cp, C->PRESSURE, C->ENTROPY, C->ENTROPY );
    set_interp2d( dtdpsMelFilename, &P->melt_prop.dTdPs, C->PRESSURE, C->ENTROPY, C->DTDP );
    set_interp2d( rhoMelFilename, &P->melt_prop.rho, C->PRESSURE, C->ENTROPY, C->DENSITY );
    set_interp2d( tempMelFilename, &P->melt_prop.temp, C->PRESSURE, C->ENTROPY, C->TEMP );
    /* const */
    // FIXME remove this redundancy
    P->melt_prop.cond = P->cond_mel;
    P->melt_prop.log10visc = P->log10visc_mel;

    /* liquidus and solidus */
    /* 1d */

    set_interp1d( liquidusFilename, &P->solid_prop.liquidus, NLS, C->PRESSURE, C->ENTROPY );
    set_interp1d( liquidusFilename, &P->melt_prop.liquidus, NLS, C->PRESSURE, C->ENTROPY );
    /* duplication here, but want to remain flexible for future
       approaches for a multicomponent system */
    set_interp1d( solidusFilename, &P->solid_prop.solidus, NLS, C->PRESSURE, C->ENTROPY );
    set_interp1d( solidusFilename, &P->melt_prop.solidus, NLS, C->PRESSURE, C->ENTROPY );

    PetscFunctionReturn(0);
}
