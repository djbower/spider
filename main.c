static char help[] ="Parallel magma ocean timestepper";

#include "petsc.h"

#include "ctx.h"

// !!!!!!!!! Global fix required to use DMDAVecGetArray, not VecGetArray..

// TODO: Make sure to ensure this is valgrind clean

// !! This RHS function needs to be checked - it probably has errors
#undef __FUNCT__
#define __FUNCT__ "RHSFunction"
PetscErrorCode RHSFunction(TS ts,PetscReal t,Vec S_in,Vec rhs_s,void *ptr)
{
  PetscErrorCode ierr;
  Ctx               *E = (Ctx*) ptr;
  Solution          *S = &E->solution;
  PetscScalar       *arr_rhs_s;
  const PetscScalar *arr_Etot,*arr_lhs_s;
  PetscMPIInt       rank,size;
  PetscInt          i,ihi_s,ilo_s,ihi,ilo,w_s;
  DM                da_s = E->da_s, da_b=E->da_b;

  PetscFunctionBeginUser;
#if (defined VERBOSE)
    printf("rhs:\n"); // TODO: replace all printf with PetscPrintf
#endif

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);

    /* S_in is the solution array.  It's easiest to store this in
       the E struct for future access, and this step is done at the
       top of set_capacitance */

    /* loop over staggered nodes and populate E struct */
    set_capacitance( E, S_in );

    /* loop over basic (internal) nodes and populate E struct */
    set_matprop_and_flux( E );

    /* surface radiative boundary condition
       parameterised ultra-thin thermal boundary layer
       constants given by fit (python script) */

    // NOTE: here we assume that the first rank has the first point
    if (!rank) {
       PetscScalar temp_s_0,temp0,dT,val;
       PetscInt ind = 0;
       ierr = VecGetValues(S->temp_s,1,&ind,&temp_s_0);CHKERRQ(ierr);
       dT = CONSTBC * PetscPowScalar( temp_s_0, EXPBC ); 
       temp0 = temp_s_0 - dT;
       val = SIGMA*(PetscPowScalar(temp0,4.0)-PetscPowScalar(TEQM,4.0));

    /* Note - below is true because non-dim geom is exactly equal
       to one, so do not need to multiply by area of surface */
      ierr = VecSetValue(S->Etot,0,val,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(S->Jtot,0,val,INSERT_VALUES);CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Etot);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(S->Jtot);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(S->Jtot);CHKERRQ(ierr);

    // !! Here and elsewhere, are we just getting lucky, getting away with global vectors?

    /* loop over staggered nodes except last node */
    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ilo = ilo_s;
    ihi = ihi_s == NUMPTSS ? NUMPTSS -1 : ihi_s;
    ierr = DMDAVecGetArray(da_s,rhs_s,&arr_rhs_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->Etot,&arr_Etot);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->lhs_s,&arr_lhs_s);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        arr_rhs_s[i] = arr_Etot[i+1] - arr_Etot[i];
        arr_rhs_s[i] /= arr_lhs_s[i];
    }
    ierr = DMDAVecRestoreArray(da_s,rhs_s,&arr_rhs_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->Etot,&arr_Etot);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->lhs_s,&arr_lhs_s);CHKERRQ(ierr);


    /* for last point, just cool by a constant factor
       if CMBBC = 0.0 --> zero heat flux CMB
       if CMBBC = 1.0 --> isothermal CMB
       if 0.0 < CMBBC < 1.0 CMB is cooling proportional to
           time-dependent energy requirements at base of
           magma ocean */

    // NOTE: here, we somewhat dangerously assume that the last proc has the last point 
    if (rank == size-1) {
      PetscScalar val,val2;
      PetscInt ind;

      ind = NUMPTS-2;
      ierr = VecGetValues(S->Etot,1,&ind,&val);CHKERRQ(ierr);
      val *= CMBBC;
      ierr = VecSetValue(S->Etot,NUMPTS-1,val,INSERT_VALUES);CHKERRQ(ierr);

      ind = NUMPTSS-1;
      ierr = VecGetValues(S->lhs_s,1,&ind,&val2);CHKERRQ(ierr);
      val2 = ((1.0 - CMBBC)*val)/val2;
      ierr = VecSetValue(rhs_s,NUMPTSS-1,val2,INSERT_VALUES);
    }

    ierr = VecAssemblyBegin(rhs_s);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs_s);CHKERRQ(ierr);

    /* copy to the Ctx for monitoring (only). This can be removed leater TODO */
    ierr = VecCopy(rhs_s,E->solution.rhs_s);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

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
#if 0
  ierr = PetscPrintf(PETSC_COMM_WORLD," *** Viewing solution S_s ***\n");CHKERRQ(ierr);
  ierr = VecView(ctx.solution.S_s,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif

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
