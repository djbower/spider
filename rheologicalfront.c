#include "rheologicalfront.h"
#include "dimensionalisablefield.h"

static PetscErrorCode JSON_add_rheological_front_mantle_properties( DM, Constants *, RheologicalFrontMantleProperties *, const char *, cJSON * );

PetscErrorCode set_rheological_front_mask( DM dm, const Vec vec, const PetscScalar value, PetscInt *index, Vec mask )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo,ihi,w;
    const PetscScalar *arr_vec, one=1.0, zero=0.0;
    PetscScalar       vec_value;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(dm,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;

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

    ierr = DMDAVecGetArrayRead(dm,vec,&arr_vec); CHKERRQ(ierr);

    //VecSetValues( mask, 1, 0, &one, INSERT_VALUES ); CHKERRQ(ierr);

    for(i=ilo; i<ihi; ++i){
        vec_value = arr_vec[i];
        ierr = VecSetValues( mask, 1, &i, &one, INSERT_VALUES ); CHKERRQ(ierr);
        if( vec_value < value ){
            ierr = VecSetValues( mask, 1, &i, &zero, INSERT_VALUES ); CHKERRQ(ierr);
            break;
        }
    }

    ierr = VecAssemblyBegin(mask);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(mask);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayRead(dm,vec,&arr_vec);CHKERRQ(ierr);

    *index = i;

    PetscFunctionReturn(0);

}

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
