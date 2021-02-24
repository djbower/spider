#include "monitor.h"
#include "util.h"

static PetscScalar GetModifiedMixingLength( PetscScalar, PetscScalar, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar GetConstantMixingLength( PetscScalar outer_radius, PetscScalar inner_radius );
static PetscScalar GetMixingLength( const Parameters, PetscScalar );

/* Helper routine to prepend the root directory to a relative path */
/* https://gcc.gnu.org/onlinedocs/gcc-4.9.0/cpp/Stringification.html */
#define STRINGIFY(x) STRINGIFY2(x)
#define STRINGIFY2(x) #x
#define SPIDER_ROOT_DIR_STR STRINGIFY(SPIDER_ROOT_DIR)
PetscErrorCode MakeRelativeToSourcePathAbsolute(char* path) {
  PetscErrorCode ierr;
  char tmp[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;
  ierr = PetscStrcpy(tmp,path);CHKERRQ(ierr);
  ierr = PetscStrcpy(path,SPIDER_ROOT_DIR_STR);CHKERRQ(ierr);
  ierr = PetscStrcat(path,"/");CHKERRQ(ierr); /* not portable */
  ierr = PetscStrcat(path,tmp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#undef SPIDER_ROOT_DIR_STR

PetscErrorCode set_solution_from_entropy( Ctx *E, Vec sol )
{
    /* Clone data in S->dSdxi and top staggered node value (from S->S_s)
       to Vec sol */

    PetscErrorCode ierr;
    Solution       *S = &E->solution;
    PetscScalar    S0;
    PetscInt       ind0 = 0;
    Vec            *subVecs

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    /* set dS/dxi at basic nodes to solution Vec */
    ierr = VecCopy( S->dSdxi, subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDXI_B]] );CHKERRQ(ierr);
    /* set S0 to solution Vec */
    ierr = VecGetValues( S->S_s,1,&ind0,&S0);CHKERRQ(ierr);
    ierr = VecSetValue( subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],0,S0,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]]);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]]);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    ierr = PetscFree(subVecs);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode set_entropy_from_solution( Ctx *E, Vec sol )
{
    /* Set entropy-related Vecs in E to be consistent with
       the Vec sol, and additionally compute values for
       outer boundaries of the basic mesh */

    PetscErrorCode ierr;
    Mesh           *M = &E->mesh;
    Solution       *S = &E->solution;
    PetscScalar    S0;
    PetscScalar    *arr_S_b, *arr_S_s, *arr_dSdxi_b, *arr_xi_s, *arr_xi_b;
    PetscInt       i, ihi_b, ilo_b, w_b;
    PetscMPIInt    size;
    DM             da_s = E->da_s, da_b=E->da_b;
    Vec            *subVecs;

    const PetscInt ind0 = 0;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    /* set S->dS/dxi at basic nodes */
    ierr = VecCopy( subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_DSDXI_B]], S->dSdxi );CHKERRQ(ierr);
    /* get first staggered node value (store as S0) */
    ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_S0]],1,&ind0,&S0);CHKERRQ(ierr);

    /* for looping over basic nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    ierr = DMDAVecGetArray(da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,M->xi_s,&arr_xi_s);CHKERRQ(ierr);

    /* S (absolute) at all staggered and basic internal nodes */
    arr_S_s[0] = 0.0;
    /* note plus one for start of loop */
    for(i=ilo_b+1; i<ihi_b-1; ++i){
      /* S (absolute) at staggered nodes */
      arr_S_s[i] = arr_dSdxi_b[i] * (arr_xi_s[i] - arr_xi_s[i-1] );
      arr_S_s[i] += arr_S_s[i-1]; // dS relative to first staggered value
      arr_S_b[i] = arr_dSdxi_b[i] * 0.5 * (arr_xi_b[i] - arr_xi_b[i-1] );
      arr_S_b[i] += arr_S_s[i-1];
      arr_S_s[i-1] += S0; // add at end to try and retain precision
      arr_S_b[i-1] += S0; // add at end to try and retain precision
    }
    /* loop above terminates before we have added the constant offset
       to the last staggered node value.  So do that here */
    arr_S_s[ihi_b-2] += S0; // add at end to try and retain precision
    arr_S_b[ihi_b-2] += S0; // add at end to try and retain precision

    /* surface entropy 
       extrapolate to surface using dS/dr below
       in practice, this is only used to give an estimate of the
       surface temperature for computing the atmosphere initial
       condition.  The surface entropy is then updated to adhere
       to the radiation balance at the surface */
    arr_S_b[0] = arr_S_s[0];
    arr_S_b[0] += -arr_dSdxi_b[1] * 0.5 * (arr_xi_b[1] - arr_xi_b[0]);

    /* cmb entropy */
    /* use dS/dr at the cmb which is constrained by the core boundary condition */
    arr_S_b[ihi_b-1] = arr_dSdxi_b[ihi_b-1] * 0.5 * (arr_xi_b[ihi_b-1]-arr_xi_b[ihi_b-2]);
    arr_S_b[ihi_b-1] += arr_S_s[ihi_b-2];

    ierr = DMDAVecRestoreArray(da_b,S->S,&arr_S_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->S_s,&arr_S_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,S->dSdxi,&arr_dSdxi_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_b,M->xi_b,&arr_xi_b);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,M->xi_s,&arr_xi_s);CHKERRQ(ierr);

    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    PetscFree(subVecs);

    PetscFunctionReturn(0);
}

PetscErrorCode set_partial_pressures_from_solution( Ctx *E, Vec sol )
{
    PetscErrorCode             ierr;
    Atmosphere                 *A = &E->atmosphere;
    Parameters                 P = E->parameters;
    AtmosphereParameters const Ap = P->atmosphere_parameters;
    PetscInt                    i;
    PetscMPIInt                 size;
    Vec                        *subVecs;

    PetscFunctionBeginUser;

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"This code has only been correctly implemented for serial runs");

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    /* partial pressures */
    for( i=0; i<Ap->n_volatiles; ++i) {
        ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],1,&i,&A->volatiles[i].p);CHKERRQ(ierr);
    }

    /* mass reaction terms */
    for( i=0; i<Ap->n_reactions; ++i) {
        ierr = VecGetValues(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]],1,&i,&A->mass_reaction[i]);CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    PetscFree(subVecs);

    PetscFunctionReturn(0);
}

PetscErrorCode set_solution_from_partial_pressures( Ctx *E, Vec sol )
{
    PetscErrorCode             ierr;
    Parameters                 P = E->parameters;
    Atmosphere                 *A = &E->atmosphere;
    AtmosphereParameters const Ap = P->atmosphere_parameters;
    PetscInt                   i;
    Vec                        *subVecs;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(E->numFields,&subVecs);CHKERRQ(ierr);
    ierr = DMCompositeGetAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);

    /* partial pressures */
    for( i=0; i<Ap->n_volatiles; ++i) {
        ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]],i,A->volatiles[i].p,INSERT_VALUES);CHKERRQ(ierr);
    }

    /* mass reaction terms */
    for( i=0; i<Ap->n_reactions; ++i) {
        ierr = VecSetValue(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]],i,A->mass_reaction[i],INSERT_VALUES);CHKERRQ(ierr);
    }

    if (Ap->n_volatiles > 0) {
      ierr = VecAssemblyBegin(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]]);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_VOLATILES]]);CHKERRQ(ierr);
    }
    if (Ap->n_reactions > 0) {
      ierr = VecAssemblyBegin(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]]);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(subVecs[E->solutionSlots[SPIDER_SOLUTION_FIELD_MO_REACTIONS]]);CHKERRQ(ierr);
    }

    ierr = DMCompositeRestoreAccessArray(E->dm_sol,sol,E->numFields,NULL,subVecs);CHKERRQ(ierr);
    PetscFree(subVecs);

    PetscFunctionReturn(0);
}

PetscScalar combine_matprop( PetscScalar weight, PetscScalar mat1, PetscScalar mat2 )
{
    /* linear weighting of two quantities */

    PetscScalar out;

    out = weight * mat1 + (1.0-weight) * mat2;

    return out;
}

PetscScalar tanh_weight( PetscScalar qty, PetscScalar threshold, PetscScalar width )
{
    /* tanh weight for viscosity profile and smoothing */

    PetscScalar    fwt, z;

    z = ( qty - threshold ) / width;
    fwt = 0.5 * ( 1.0 + PetscTanhScalar( z ) );
    return fwt;
}

PetscScalar get_smoothing( PetscScalar smooth_width, PetscScalar gphi )
{
    /* get smoothing across phase boundaries for a two phase composite */

    PetscScalar smth;

    /* no smoothing */
    if( smooth_width == 0.0 ){
        smth = 1.0; // mixed phase only
        if( (gphi < 0.0) || (gphi > 1.0) ){
            smth = 0.0; // single phase only
        }
    }

    /* tanh smoothing */
    else{
        if( gphi > 0.5 ){
            smth = 1.0 - tanh_weight( gphi, 1.0, smooth_width );
        }
        else{
            smth = tanh_weight( gphi, 0.0, smooth_width );
        }
    }

    return smth;
}

PetscErrorCode Make2DPetscScalarArray( PetscInt arraySizeX, PetscInt arraySizeY, PetscScalar ***theArray_ptr)
{
    PetscErrorCode ierr;
    PetscInt i;
    PetscScalar **theArray;

    PetscFunctionBeginUser;

    ierr = PetscMalloc1(arraySizeX, theArray_ptr);CHKERRQ(ierr);
    theArray = *theArray_ptr;
    for( i=0; i<arraySizeX; i++){
        ierr = PetscMalloc1(arraySizeY, &theArray[i]);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode make_vec_mask( DM dm, PetscInt index, Vec * mask_ptr )
{
    PetscErrorCode ierr;
    PetscInt       numpts,i,ilo,ihi,w;
    const PetscScalar one=1.0, zero=0.0;

    ierr = DMDAGetInfo(dm,NULL,&numpts,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* create mask vector */
    ierr = VecCreate( PETSC_COMM_WORLD, mask_ptr );CHKERRQ(ierr);
    ierr = VecSetSizes( *mask_ptr, PETSC_DECIDE, numpts );CHKERRQ(ierr);
    ierr = VecSetFromOptions( *mask_ptr );CHKERRQ(ierr);
    ierr = VecSetUp( *mask_ptr );CHKERRQ(ierr);

    ierr = VecSet( *mask_ptr, 0.0 );CHKERRQ(ierr);

    ierr = DMDAGetCorners(dm,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;

    for(i=ilo; i<ihi; ++i){
        if( i < index ){
            ierr = VecSetValues( *mask_ptr, 1, &i, &one, INSERT_VALUES );CHKERRQ(ierr);
        }
        else{
            ierr = VecSetValues( *mask_ptr, 1, &i, &zero, INSERT_VALUES );CHKERRQ(ierr);
            break;
        }
    }

    ierr = VecAssemblyBegin(*mask_ptr);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(*mask_ptr);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

PetscErrorCode average_by_mass_staggered( Ctx *E, Vec in_vec, Vec * in_mask_ptr, PetscScalar *out )
{
    PetscErrorCode    ierr;
    Mesh              *M = &E->mesh;
    Vec               qty_s, mass_s;
    PetscScalar       mass;
    PetscInt          numpts_s;
    Vec               in_mask = *in_mask_ptr;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* create work vectors */
    ierr = VecCreate( PETSC_COMM_WORLD, &qty_s ); CHKERRQ(ierr);
    ierr = VecSetSizes( qty_s, PETSC_DECIDE, numpts_s ); CHKERRQ(ierr);
    ierr = VecSetFromOptions( qty_s ); CHKERRQ(ierr);
    ierr = VecSetUp( qty_s ); CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &mass_s ); CHKERRQ(ierr);
    ierr = VecSetSizes( mass_s, PETSC_DECIDE, numpts_s ); CHKERRQ(ierr);
    ierr = VecSetFromOptions( mass_s ); CHKERRQ(ierr);
    ierr = VecSetUp( mass_s ); CHKERRQ(ierr);

    /* compute mass-weighted average */
    ierr = VecPointwiseMult(qty_s,in_vec,in_mask); CHKERRQ(ierr);
    ierr = VecPointwiseMult(mass_s,M->mass_s,in_mask); CHKERRQ(ierr);
    ierr = VecPointwiseMult(qty_s,qty_s,mass_s); CHKERRQ(ierr);
    ierr = VecSum(qty_s,out); CHKERRQ(ierr);
    ierr = VecSum(mass_s,&mass); CHKERRQ(ierr);

    *out = *out / mass;

    /* this block would prevent 'nans' from being returned, but 'nans'
       are helpful because it indicates the rheological front is not
       yet moving */
#if 0
    if (mass > 0.0){
        *out = *out / mass;
    }
    else{
        *out = 0.0;
    }
#endif

    ierr = VecDestroy( &qty_s ); CHKERRQ(ierr);
    ierr = VecDestroy( &mass_s ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

PetscErrorCode invert_vec_mask( Vec * in_vec_ptr )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    ierr = VecShift( *in_vec_ptr, -1 ); CHKERRQ(ierr);
    ierr = VecScale( *in_vec_ptr, -1 ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

/* helper functions for parsing parameters */

PetscErrorCode PetscScalarCheckPositive( PetscScalar value, const char * value_string )
{
    PetscFunctionBeginUser;
    if( value < 0.0 ){
        SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"%s must be positive (currently %f)",value_string,value);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode PetscIntCheckPositive( PetscInt value, const char * value_string )
{
    PetscFunctionBeginUser;
    if( value < 0 ){
        SETERRQ2(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"%s must be positive (currently %d)",value_string,value);
    }
    PetscFunctionReturn(0);
}

PetscErrorCode PetscOptionsGetPositiveScalar( const char *value_string, PetscScalar *value_ptr, PetscScalar value_default, PetscBool *set )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    *value_ptr = value_default;
    ierr = PetscOptionsGetScalar(NULL,NULL,value_string,value_ptr,set);CHKERRQ(ierr);
    if(set!=NULL && *set){
      ierr = PetscScalarCheckPositive(*value_ptr,value_string);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode PetscOptionsGetPositiveInt( const char *value_string, PetscInt *value_ptr, PetscInt value_default, PetscBool *set )
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    *value_ptr = value_default;
    ierr = PetscOptionsGetInt(NULL,NULL,value_string,value_ptr,set);CHKERRQ(ierr);
    if(set!=NULL && *set){
      ierr = PetscIntCheckPositive(*value_ptr,value_string);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

/* to avoid circular references with headers, the eddy diffusivity function and associated functions live here */
PetscErrorCode GetEddyDiffusivity( const EOSEvalData eos_eval, const Parameters P, PetscScalar radius, PetscScalar dSdxi, PetscScalar dxidr, PetscScalar *kappah_ptr, PetscScalar *kappac_ptr, PetscScalar *regime_ptr )
{
    PetscErrorCode ierr;
    PetscScalar    visc, kvisc, gsuper, kh, crit, mix, kappah, kappac, regime;

    PetscFunctionBeginUser;

    visc = eos_eval.log10visc;
    ierr = apply_log10visc_cutoff( P, &visc );CHKERRQ(ierr);
    visc = PetscPowScalar( 10.0, visc );
    kvisc = visc / eos_eval.rho; // kinematic viscosity
    gsuper = P->gravity * eos_eval.T / eos_eval.Cp * dSdxi * dxidr; // g * super adiabatic gradient

    crit = 81.0 * PetscPowScalar(kvisc,2);
    mix = GetMixingLength( P, radius);
    crit /= 4.0 * eos_eval.alpha * PetscPowScalar(mix,4);

    if( gsuper <= 0.0 ){
      /* no convection, subadiabatic */
      kh = 0.0;
      regime = 0.0;
    } else if( gsuper > crit ){
      /* inviscid scaling from Vitense (1953) */
      kh = 0.25 * PetscPowScalar(mix,2) * PetscSqrtScalar(eos_eval.alpha * gsuper);
      regime = 2.0;
    } else{
      /* viscous scaling */
      kh = eos_eval.alpha * gsuper * PetscPowScalar(mix,4) / ( 18.0 * kvisc );
      regime = 1.0;
    }   

    /* thermal eddy diffusivity */
    if (P->eddy_diffusivity_thermal > 0.0){
      /* scale */
      kappah = P->eddy_diffusivity_thermal * kh; 
    }   
    else{
      /* else set (and negate to account for sign flag) */
      kappah = -P->eddy_diffusivity_thermal;
    }   

    /* chemical eddy diffusivity */
    if (P->eddy_diffusivity_chemical > 0.0){
      /* scale */
      kappac = P->eddy_diffusivity_chemical * kh; 
    }   
    else{
      /* else set (and negate to account for sign flag) */
      kappac = -P->eddy_diffusivity_chemical;
    }   

    /* update pointers with data */
    if(kappah_ptr != NULL){
        *kappah_ptr = kappah;
    }   
    if(kappac_ptr != NULL){
        *kappac_ptr = kappac;
    }   
    if(regime_ptr != NULL){
        *regime_ptr = regime;
    }   

    PetscFunctionReturn(0);
}

static PetscScalar GetModifiedMixingLength( PetscScalar a, PetscScalar b, PetscScalar outer_radius, PetscScalar inner_radius, PetscScalar radius )
{
    /* See Kamata, 2018, JGR */
    /* conventional mixing length theory has a = b = 0.5 */
    /* a is location of peak in depth/radius space,
       b is amplitude of the peak */

    PetscScalar mix_length1, mix_length2, mix_length;

    mix_length1 = (radius - inner_radius) * b / (1.0 - a); 
    mix_length2 = (outer_radius - radius) * b / a;

    mix_length = PetscMin( mix_length1, mix_length2 );

    return mix_length;
}

static PetscScalar GetConstantMixingLength( PetscScalar outer_radius, PetscScalar inner_radius )
{
    PetscScalar mix_length;

    mix_length = 0.25 * (outer_radius - inner_radius );

    return mix_length;
}

static PetscScalar GetMixingLength( const Parameters P, PetscScalar radius )
{
    PetscScalar outer_radius, inner_radius;
    PetscScalar mix_length = 0.0;
    PetscScalar eps = 1.0E-10;

    /* for a single layer, P->layer_interface_radius = P->coresize (parameters.c),
       enabling this single expression to work for both a single and double
       layered mantle */
    /* due to floating points, sometimes the radius at the innermost boundary is slightly
       less than P->layer_interface_radius, so account for this with a small offset (eps) */
    /* if we try to resolve the ultra-thin boundary layer, eps might not be smaller enough? */

    if( radius >= P->radius * P->layer_interface_radius - eps ){
        outer_radius = P->radius;
        inner_radius = P->radius * P->layer_interface_radius;
    }   
    else{
        outer_radius = P->radius * P->layer_interface_radius;
        inner_radius = P->radius * P->coresize;
    }   

    if( P->mixing_length == 1){ 
        mix_length = GetModifiedMixingLength( P->mixing_length_a, P->mixing_length_b, outer_radius, inner_radius, radius );
    }   
    else if( P->mixing_length == 2){ 
        mix_length = GetConstantMixingLength( outer_radius, inner_radius );
    }   

    /* parameters.c ensures that P->mixing_length must be 1 or 2, so we cannot
       fall outside this if statement */

    return mix_length;
}

PetscErrorCode apply_log10visc_cutoff( Parameters const P, PetscScalar *viscosity )
{   
    PetscFunctionBeginUser;
    
    if(P->log10visc_min > 0.0){
        if(*viscosity < P->log10visc_min){
            *viscosity = P->log10visc_min;
        }
    }
    if(P->log10visc_max > 0.0){
        if(*viscosity > P->log10visc_max){
            *viscosity = P->log10visc_max;
        }
    }
    
    PetscFunctionReturn(0);
}
