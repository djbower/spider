static char help[] ="Parallel magma ocean timestepper";

#include "petsc.h"
#include "ctx.h" 

PetscErrorCode RHSFunction(TS,PetscReal,Vec,Vec,void*);

int main(int argc, char ** argv)
{
  PetscErrorCode ierr;
  TS             ts;      /* ODE solver object */
  Vec            S_s;     /* Solution Vector */
  Ctx            ctx;     /* Solver context */
  PetscInt       i;

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);

  // TODO: Once this runs, but all this setup into a single function

  /* Initialize context with parameters (most parameters are constants in global_defs.h, though) */
  set_lookups(&ctx);

  /* Set up a parallel structured grid as DMComposite with two included DMDAs
     This is used to define vectors which hold the solution. The included 
     DMDAs are used to create vectors which hold various properties
     living on the same primal and staggered nodes.
  */
  const PetscInt stencilWidth = 1; //TODO: check that this is appropriate
  const PetscInt dof = 1;
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,NUMPTS,dof,stencilWidth,NULL,&ctx.da_b);CHKERRQ(ierr);
  ierr = DMDACreate1d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,NUMPTSS,dof,stencilWidth,NULL,&ctx.da_s);CHKERRQ(ierr);

  // TODO (maybe): PETSc-fy this more by getting rid of the NUMPTS and NUMPTSS parameters and instead letting the DMDAs themselves define this information (hence allowing more command-line control)

  /* Continue to initialize context with distributed data */

  for (i=0;i<NUMSOLUTIONVECS;++i){
    ierr = DMCreateGlobalVector(ctx.da_b,&(ctx.solution.solutionVecs[i]));CHKERRQ(ierr);
  }
  ctx.solution.alpha               = ctx.solution.solutionVecs[0];
  ctx.solution.alpha_mix           = ctx.solution.solutionVecs[1];
  ctx.solution.cond                = ctx.solution.solutionVecs[2];
  ctx.solution.cp                  = ctx.solution.solutionVecs[3];
  ctx.solution.cp_mix              = ctx.solution.solutionVecs[4];  // TI
  ctx.solution.dfusdr              = ctx.solution.solutionVecs[5];  // TI
  ctx.solution.dfusdr_temp         = ctx.solution.solutionVecs[6];  // TI
  ctx.solution.dphidr              = ctx.solution.solutionVecs[7];
  ctx.solution.dSdr                = ctx.solution.solutionVecs[8];  // TI
  ctx.solution.dTdPs               = ctx.solution.solutionVecs[9];
  ctx.solution.dTdPs_mix           = ctx.solution.solutionVecs[10]; // TI
  ctx.solution.dTdrs               = ctx.solution.solutionVecs[11];
  ctx.solution.dTdrs_mix           = ctx.solution.solutionVecs[12]; // TI
  ctx.solution.Etot                = ctx.solution.solutionVecs[13];
  ctx.solution.fusion              = ctx.solution.solutionVecs[14]; // TI
  ctx.solution.fusion_curve        = ctx.solution.solutionVecs[15]; // TI
  ctx.solution.fusion_curve_temp   = ctx.solution.solutionVecs[16]; // TI
  ctx.solution.fusion_rho          = ctx.solution.solutionVecs[17]; // TI
  ctx.solution.fusion_temp         = ctx.solution.solutionVecs[18]; // TI
  ctx.solution.gsuper              = ctx.solution.solutionVecs[19];
  ctx.solution.Jcond               = ctx.solution.solutionVecs[20];
  ctx.solution.Jconv               = ctx.solution.solutionVecs[21];
  ctx.solution.Jgrav               = ctx.solution.solutionVecs[22];
  ctx.solution.Jheat               = ctx.solution.solutionVecs[23];
  ctx.solution.Jmass               = ctx.solution.solutionVecs[24];
  ctx.solution.Jmix                = ctx.solution.solutionVecs[25];
  ctx.solution.Jtot                = ctx.solution.solutionVecs[26];
  ctx.solution.kappah              = ctx.solution.solutionVecs[27];
  ctx.solution.liquidus            = ctx.solution.solutionVecs[28]; // TI
  ctx.solution.liquidus_rho        = ctx.solution.solutionVecs[29]; // TI
  ctx.solution.liquidus_temp       = ctx.solution.solutionVecs[30]; // TI
  ctx.solution.nu                  = ctx.solution.solutionVecs[31];
  ctx.solution.phi                 = ctx.solution.solutionVecs[32];
  ctx.solution.rho                 = ctx.solution.solutionVecs[33];
  ctx.solution.S                   = ctx.solution.solutionVecs[34];
  ctx.solution.solidus             = ctx.solution.solutionVecs[35]; // TI
  ctx.solution.solidus_rho         = ctx.solution.solutionVecs[36]; // TI
  ctx.solution.solidus_temp        = ctx.solution.solutionVecs[37]; // TI
  ctx.solution.temp                = ctx.solution.solutionVecs[38];
  ctx.solution.visc                = ctx.solution.solutionVecs[39];

  ierr = DMCreateLocalVector(ctx.da_b,&ctx.work_local_b);CHKERRQ(ierr);

  for (i=0;i<NUMSOLUTIONVECSS;++i){
    ierr = DMCreateGlobalVector(ctx.da_s,&(ctx.solution.solutionVecsS[i]));CHKERRQ(ierr);
  }
  ctx.solution.fusion_s            = ctx.solution.solutionVecsS[0]; // TI
  ctx.solution.fusion_curve_s      = ctx.solution.solutionVecsS[1]; // TI
  ctx.solution.fusion_curve_temp_s = ctx.solution.solutionVecsS[2]; // TI
  ctx.solution.fusion_temp_s       = ctx.solution.solutionVecsS[3]; // TI
  ctx.solution.lhs_s               = ctx.solution.solutionVecsS[4];
  ctx.solution.liquidus_rho_s      = ctx.solution.solutionVecsS[5]; // TI
  ctx.solution.liquidus_s          = ctx.solution.solutionVecsS[6]; // TI
  ctx.solution.liquidus_temp_s     = ctx.solution.solutionVecsS[7]; // TI
  ctx.solution.phi_s               = ctx.solution.solutionVecsS[8];
  ctx.solution.rhs_s               = ctx.solution.solutionVecsS[9]; 
  ctx.solution.rho_s               = ctx.solution.solutionVecsS[10];
  ctx.solution.S_s                 = ctx.solution.solutionVecsS[11];
  ctx.solution.solidus_s           = ctx.solution.solutionVecsS[12]; // TI
  ctx.solution.solidus_rho_s       = ctx.solution.solutionVecsS[13]; // TI
  ctx.solution.solidus_temp_s      = ctx.solution.solutionVecsS[14]; // TI
  ctx.solution.temp_s              = ctx.solution.solutionVecsS[15];

  ierr = DMCreateLocalVector(ctx.da_s,&ctx.work_local_s);CHKERRQ(ierr);

  set_mesh(&ctx);

  set_time_independent(&ctx);

  /* Create a solution vector (S at the staggered nodes) and fill with an initial condition */
  ierr = DMCreateGlobalVector(ctx.da_s,&S_s);CHKERRQ(ierr); 
  set_initial_condition(&ctx);CHKERRQ(ierr); // NOTE: as usual, it's redundant to store the solution in the context as well but for now we do so for consistency
  ierr = VecCopy(ctx.solution.S_s,S_s);CHKERRQ(ierr);

  {
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," *** Viewing S_s initial condition ***\n");CHKERRQ(ierr);
    ierr = VecView(ctx.solution.S_s,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  
  /* Set up the Jacobian function (omitted for now) */

  /* Set up timestepper */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,S_s);CHKERRQ(ierr);
 // ierr = TSSetType(ts,TSSUNDIALS);CHKERRQ(ierr);
  // TODO: More CVODE setup ...

  /* Set up the RHS Function */
  ierr = TSSetRHSFunction(ts,NULL,RHSFunction,&ctx);CHKERRQ(ierr);

  /* Set up the integration period */
  ierr = TSSetDuration(ts,1,1.e12);CHKERRQ(ierr); /* One time step, huge final time */
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);

  /* Accept command line options */
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  /* Set up some kind of monitor or viewer (omitted for now) */

  /* Solve */
  ierr = TSSolve(ts,S_s);CHKERRQ(ierr);

  /* For debugging, view some things stored in the context.
      Note that there are other viewer implementations, including 
      those that will write to a file which we should be able
      to open with python.
  */
  {
    PetscViewer viewer;
    ierr = PetscViewerCreate(PETSC_COMM_WORLD,&viewer);CHKERRQ(ierr);
    ierr = PetscViewerSetType(viewer,PETSCVIEWERASCII);CHKERRQ(ierr);
    ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," *** Viewing S_s ***\n");CHKERRQ(ierr);
    ierr = VecView(ctx.solution.S_s,viewer);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD," *** Viewing rhs_s ***\n");CHKERRQ(ierr);
    ierr = VecView(ctx.solution.rhs_s,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }

  /* Destroy data allocated in Ctx */
  free_memory_interp(&ctx);
  ierr = DMDestroy(&ctx.da_s);CHKERRQ(ierr);
  ierr = DMDestroy(&ctx.da_b);CHKERRQ(ierr);
  for (i=0;i<NUMMESHVECS;++i){
    ierr = VecDestroy(&ctx.mesh.meshVecs[i]);CHKERRQ(ierr);
  }
  for (i=0;i<NUMMESHVECSS;++i){
    ierr = VecDestroy(&ctx.mesh.meshVecsS[i]);CHKERRQ(ierr);
  }
  for (i=0;i<NUMSOLUTIONVECS;++i){
    ierr = VecDestroy(&ctx.solution.solutionVecs[i]);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&ctx.work_local_b);CHKERRQ(ierr);
  for (i=0;i<NUMSOLUTIONVECSS;++i){
    ierr = VecDestroy(&ctx.solution.solutionVecsS[i]);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&ctx.work_local_s);CHKERRQ(ierr);
  // TODO: move into a dedicated function in ctx.c

  /* Cleanup and finalize */
  ierr = VecDestroy(&S_s);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}
