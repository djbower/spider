#include "rheologicalfront.h"
#include "dimensionalisablefield.h"

static PetscErrorCode set_rheological_front_mask( DM, Vec, PetscScalar, PetscInt *, Vec );
//static PetscErrorCode set_rheological_front_mantle_properties( Ctx *, RheologicalFront *, PetscInt, Vec );
static PetscErrorCode JSON_add_rheological_front_mantle_properties( DM, Constants *, RheologicalFrontMantleProperties *, const char *, cJSON * );

static PetscErrorCode set_rheological_front_mask( DM dm, const Vec vec, const PetscScalar value, PetscInt *index, Vec mask_s )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    const PetscScalar *arr_vec_s, one=1.0, zero=0.0;
    PetscScalar       vec_value;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(dm,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    /* this simple algorithm counts up from the surface towards the
       CMB until the melt fraction at a staggered node is larger than
       the melt fraction value (rheological_front_phi) used to define
       the base of the magma ocean.  Once this value is reached, the 
       loop is exited and the index stored for later use.  There are 
       many instances when this algorithm could return nonsense 
       values, but as long as the magma ocean is generally 
       crystalising from the bottom-up, it should be OK */

    /* end-member cases:
           i = 0 if surface is below rheological transition
           i = ihi_s-1 if all mantle is above rheological transition */

    /* this loop should always return a meaningful value if the cooling
       sequence can be adequately modelled as bottom-up */

    ierr = DMDAVecGetArrayRead(dm,vec,&arr_vec_s); CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){
        vec_value = arr_vec_s[i];
        ierr = VecSetValues( mask_s, 1, &i, &one, INSERT_VALUES ); CHKERRQ(ierr);
        if( vec_value < value ){
            ierr = VecSetValues( mask_s, 1, &i, &zero, INSERT_VALUES ); CHKERRQ(ierr);
            break;
        }
    }

    ierr = VecAssemblyBegin(mask_s);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(mask_s);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayRead(dm,vec,&arr_vec_s);CHKERRQ(ierr);

    *index = i;

    PetscFunctionReturn(0);

}

#if 0
PetscErrorCode set_rheological_front_using_phi_critical( DM dm, RheologicalFront *Rf )
{
    PetscErrorCode    ierr;
    PetscInt          numpts_s, index;
    Vec               mask_s;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(dm,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    index = numpts_s;

    /* create mask vector */
    ierr = VecCreate( PETSC_COMM_WORLD, &mask_s ); CHKERRQ(ierr);
    ierr = VecSetSizes( mask_s, PETSC_DECIDE, numpts_s ); CHKERRQ(ierr);
    ierr = VecSetFromOptions( mask_s ); CHKERRQ(ierr);
    ierr = VecSetUp( mask_s ); CHKERRQ(ierr);
    ierr = VecSet( mask_s, 0.0 ); CHKERRQ(ierr);

    /* TODO: this computes the rheological front mask based on the
       melt fraction, but another option is to compute based on the
       dynamic criterion instead */
    ierr = set_rheological_front_mask_using_phi_critical( E, &index, mask_s ); CHKERRQ(ierr);

    ierr = set_rheological_front_mantle_properties( E, Rf, index, mask_s );

    ierr = VecDestroy( &mask_s ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
#endif

#if 0
static PetscErrorCode set_rheological_front_mantle_properties( Ctx *E, RheologicalFront *Rf, PetscInt const index, Vec mask_s )
{
    PetscErrorCode   ierr;
    const DM         da_s = E->da_s;
    const Mesh       *M = &E->mesh;
    const Parameters *P = &E->parameters;
    const Solution   *S = &E->solution;
    PetscScalar      phi, radius, pressure, temperature;
    PetscInt         numpts_s, index_above, index_below;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* rheological front coordinates */
    Rf->mesh_index = index;
    ierr = VecGetValues(M->radius_b,1,&index,&radius);CHKERRQ(ierr);
    Rf->depth = P->radius - radius;
    ierr = VecGetValues(M->pressure_b,1,&index,&pressure);CHKERRQ(ierr);
    Rf->pressure = pressure;

    /* mantle properties in magma ocean (above rheological front) */
    /* middle of layer */
    index_above = index/2; /* integer algebra */
    ierr = VecGetValues(S->phi,1,&index_above,&phi);CHKERRQ(ierr);
    Rf->above_middle.phi = phi;
    ierr = VecGetValues(M->radius_b,1,&index_above,&radius);CHKERRQ(ierr);
    Rf->above_middle.depth = P->radius - radius;
    ierr = VecGetValues(M->pressure_b,1,&index_above,&pressure);CHKERRQ(ierr);
    Rf->above_middle.pressure = pressure;
    ierr = VecGetValues(S->temp,1,&index_above,&temperature);CHKERRQ(ierr);
    Rf->above_middle.temperature = temperature;
    /* average by mass */
    ierr = average_by_mass_staggered( E, S->phi_s, mask_s, &Rf->above_mass_avg.phi); CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, M->radius_s, mask_s, &Rf->above_mass_avg.depth); CHKERRQ(ierr);
    Rf->above_mass_avg.depth = P->radius - Rf->above_mass_avg.depth;
    ierr = average_by_mass_staggered( E, M->pressure_s, mask_s, &Rf->above_mass_avg.pressure); CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, S->temp_s, mask_s, &Rf->above_mass_avg.temperature); CHKERRQ(ierr);

    /* only compute properties in the solid layer once the rheological front
       begins advancing through the mantle */
    if( index < numpts_s){
        /* mantle properties in he solid layer (below rheological front) */
        /* middle of layer */
        index_below = (numpts_s - index)/2 + index;
        ierr = VecGetValues(S->phi,1,&index_below,&phi);CHKERRQ(ierr);
        Rf->below_middle.phi = phi;
        ierr = VecGetValues(M->radius_b,1,&index_below,&radius);CHKERRQ(ierr);
        Rf->below_middle.depth = P->radius - radius;
        ierr = VecGetValues(M->pressure_b,1,&index_below,&pressure);CHKERRQ(ierr);
        Rf->below_middle.pressure = pressure;
        ierr = VecGetValues(S->temp,1,&index_below,&temperature);CHKERRQ(ierr);
        Rf->below_middle.temperature = temperature;
        /* average by mass */
        /* invert the mask so values are unity below the rheological front */
        ierr = invert_vec_mask( mask_s ); CHKERRQ(ierr);
        ierr = average_by_mass_staggered( E, S->phi_s, mask_s, &Rf->below_mass_avg.phi); CHKERRQ(ierr);
        ierr = average_by_mass_staggered( E, M->radius_s, mask_s, &Rf->below_mass_avg.depth); CHKERRQ(ierr);
        Rf->below_mass_avg.depth = P->radius - Rf->below_mass_avg.depth;
        ierr = average_by_mass_staggered( E, M->pressure_s, mask_s, &Rf->below_mass_avg.pressure); CHKERRQ(ierr);
        ierr = average_by_mass_staggered( E, S->temp_s, mask_s, &Rf->below_mass_avg.temperature); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);

}
#endif

#if 0
PetscErrorCode JSON_add_rheological_front_using_phi_critical( const DM dm, const Constants *C, const RheologicalFront *Rf, cJSON *json )
{
    PetscErrorCode  ierr;

    PetscFunctionBeginUser;

    ierr = JSON_add_rheological_front( dm, C, Rf, json ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
#endif

PetscErrorCode JSON_add_rheological_front( DM dm, Constants *C, RheologicalFront *Rf, const char *name, cJSON *json )
{
    PetscErrorCode  ierr;
    cJSON           *data;

    PetscFunctionBeginUser;

    data = cJSON_CreateObject();

    ierr = JSON_add_single_value_to_object(dm, 1.0, "mesh_index", "None", Rf->mesh_index, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, C->RADIUS, "depth", "m", Rf->depth, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, C->PRESSURE, "pressure", "Pa", Rf->pressure, data);CHKERRQ(ierr);

    JSON_add_rheological_front_mantle_properties( dm, C, &Rf->above_middle, "above_middle", data );
    JSON_add_rheological_front_mantle_properties( dm, C, &Rf->above_mass_avg, "above_mass_avg", data );
    JSON_add_rheological_front_mantle_properties( dm, C, &Rf->below_middle, "below_middle", data );
    JSON_add_rheological_front_mantle_properties( dm, C, &Rf->below_mass_avg, "below_mass_avg", data );

    cJSON_AddItemToObject(json,name,data);

    PetscFunctionReturn(0);

}

static PetscErrorCode JSON_add_rheological_front_mantle_properties( DM dm, Constants *C, RheologicalFrontMantleProperties *Rfmp, char const *name, cJSON *json )
{
    PetscErrorCode ierr;
    cJSON          *data;

    PetscFunctionBeginUser;

    data = cJSON_CreateObject();
    ierr = JSON_add_single_value_to_object(dm, 1.0, "phi", "None", Rfmp->phi, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, C->RADIUS, "depth", "m", Rfmp->depth, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, C->PRESSURE, "pressure", "Pa", Rfmp->pressure, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, C->TEMP, "temperature", "K", Rfmp->temperature, data);CHKERRQ(ierr);
    cJSON_AddItemToObject(json,name,data);

    PetscFunctionReturn(0);

}
