#include "rheologicalfront.h"
#include "dimensionalisablefield.h"

static PetscErrorCode JSON_add_rheological_front_mantle_properties( DM, Constants *, RheologicalFrontMantleProperties *, const char *, cJSON * );

PetscInt get_crossover_index( DM dm, const Vec vec, PetscScalar value, PetscInt offset )
{
    PetscErrorCode     ierr;
    PetscInt           i,ilo,ihi,w;
    const PetscScalar  *arr_vec;
    PetscScalar        vec_value;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(dm,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;

    /* for the regime field (basic vector), the regime is always = 0 at the
       top basic node, which messes up this detection algorithm.  So always
       start from the second node (even for the staggered) */
    ilo += 1;

    /* this simple algorithm counts up from the surface towards the
       CMB until the value at a node is larger than a user defined
       critical value.  The index is then returned.  There are 
       many instances when this algorithm could return nonsense 
       values, but as long as the magma ocean is generally 
       crystalising from the bottom-up, it should be OK */

    ierr = DMDAVecGetArrayRead(dm,vec,&arr_vec); CHKERRQ(ierr);

    for(i=ilo; i<ihi; ++i){
        vec_value = arr_vec[i];
        if( vec_value < value ){
            break;
        }
    }

    /* end-members:
         1. i = 0 if surface is below rheological transition
         2. i = ihi-1 if all mantle is above rheological transition

       set_rheological_front_mask sets (staggered) mask to 0 for
       i_staggered < index, and 1 for i_staggered >= 0.  Hence we
       need to adjust the index for end-member 2 to ensure that
       if the rheological front is not located, all of the mantle
       is assumed to be above the rheological front.  This is only
       necessary when a staggered mesh is being processed, which is
       identified by the 'offset' argument */

    /* end-member 1 automatically works */

    /* end-member 2 */
    /* only actually required for staggered mesh calculation */
    if( i == ihi-1){
        i += offset;
    }

    ierr = DMDAVecRestoreArrayRead(dm,vec,&arr_vec);CHKERRQ(ierr);

    return i;

}

// to remove soon
#if 0
PetscErrorCode set_rheological_front_mask( DM dm, PetscInt index, Vec mask )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo,ihi,w;
    const PetscScalar one=1.0, zero=0.0;

    PetscFunctionBeginUser;

    /* assume all mantle is solid */
    ierr = VecSet( mask, 0.0 ); CHKERRQ(ierr);

    ierr = DMDAGetCorners(dm,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;

    for(i=ilo; i<ihi; ++i){
        if( i < index ){
            ierr = VecSetValues( mask, 1, &i, &one, INSERT_VALUES ); CHKERRQ(ierr);
        }
        else{
            ierr = VecSetValues( mask, 1, &i, &zero, INSERT_VALUES ); CHKERRQ(ierr);
            break;
        }
    }

    ierr = VecAssemblyBegin(mask);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(mask);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
#endif

PetscErrorCode JSON_add_rheological_front( DM dm, Constants *C, RheologicalFront *Rf, const char *name, cJSON *json )
{
    PetscErrorCode  ierr;
    cJSON           *data;

    PetscFunctionBeginUser;

    data = cJSON_CreateObject();

    ierr = JSON_add_single_value_to_object(dm, 1, "mesh_index", "None", Rf->mesh_index, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, C->RADIUS, "depth", "m", Rf->depth, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, C->PRESSURE, "pressure", "Pa", Rf->pressure, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, C->TEMP, "temperature", "K", Rf->temperature, data);CHKERRQ(ierr);
    ierr = JSON_add_single_value_to_object(dm, 1.0, "phi_global", "None", Rf->phi_global, data);CHKERRQ(ierr);

    ierr = JSON_add_rheological_front_mantle_properties( dm, C, &Rf->above_middle, "above_middle", data );CHKERRQ(ierr);
    ierr = JSON_add_rheological_front_mantle_properties( dm, C, &Rf->above_mass_avg, "above_mass_avg", data );CHKERRQ(ierr);
    ierr = JSON_add_rheological_front_mantle_properties( dm, C, &Rf->below_middle, "below_middle", data );CHKERRQ(ierr);
    ierr = JSON_add_rheological_front_mantle_properties( dm, C, &Rf->below_mass_avg, "below_mass_avg", data );CHKERRQ(ierr);

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
