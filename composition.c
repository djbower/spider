#include "composition.h"
#include "util.h"

#if 0
static PetscErrorCode set_magma_ocean_crystal_fraction( Ctx * );
static PetscErrorCode set_magma_ocean_bridgmanite_fraction( Ctx * );
static PetscErrorCode set_magma_ocean_mass_ratio( Ctx * );
#endif

#if 0
PetscErrorCode initialise_composition( Ctx *E )
{
    PetscErrorCode   ierr;
    Parameters *P = &E->parameters;
    Composition *Comp = &P->composition;

    PetscFunctionBeginUser;

    /* need to know the mass ratio at the liquidus to scale the
       pure-bridgmanite lookup data appropriately */

    // FIXME
    //Comp->mo_crystal_fraction = 0.0;
    ierr = set_magma_ocean_bridgmanite_fraction( E ); CHKERRQ(ierr);
    ierr = set_magma_ocean_mass_ratio( E ); CHKERRQ(ierr);

    /* this is static, because we always need this quantity */
    Comp->mass_ratio_liquidus = Comp->mo_mass_ratio;

    PetscFunctionReturn(0);

}

PetscErrorCode set_composition( Ctx *E )
{
    PetscErrorCode   ierr;

    PetscFunctionBeginUser;

    ierr = set_rheological_front( E ); CHKERRQ(ierr);
    ierr = set_magma_ocean_crystal_fraction( E ); CHKERRQ(ierr);
    ierr = set_magma_ocean_bridgmanite_fraction( E ); CHKERRQ(ierr);
    ierr = set_magma_ocean_mass_ratio( E ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
#endif

PetscErrorCode set_rheological_front( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    DM                da_s=E->da_s;
    Mesh              *M = &E->mesh;
    Solution          *S = &E->solution;
    Parameters        *P = &E->parameters;
    RheologicalFront  *RF = &P->rheological_front;
    const PetscScalar *arr_phi_s;
    PetscScalar       phi, radius, pressure, one=1.0, zero=0.0;
    PetscInt          numpts_s, i_above_avg, i_below_avg;
    Vec               mask_s;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    /* create mask vector */
    ierr = VecCreate( PETSC_COMM_WORLD, &mask_s ); CHKERRQ(ierr);
    ierr = VecSetSizes( mask_s, PETSC_DECIDE, numpts_s ); CHKERRQ(ierr);
    ierr = VecSetFromOptions( mask_s ); CHKERRQ(ierr);
    ierr = VecSetUp( mask_s ); CHKERRQ(ierr);
    ierr = VecSet( mask_s, 0.0 ); CHKERRQ(ierr);

    ierr = DMDAVecGetArrayRead(da_s,S->phi_s,&arr_phi_s); CHKERRQ(ierr);

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

    for(i=ilo_s; i<ihi_s; ++i){
        phi = arr_phi_s[i];
        ierr = VecSetValues( mask_s, 1, &i, &one, INSERT_VALUES ); CHKERRQ(ierr);
        if( phi < RF->phi_critical ){
            ierr = VecSetValues( mask_s, 1, &i, &zero, INSERT_VALUES ); CHKERRQ(ierr);
            break;
        }
    }

    ierr = VecAssemblyBegin(mask_s);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(mask_s);CHKERRQ(ierr);

    ierr = DMDAVecRestoreArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);

    // store
    RF->mesh_index = i;

    // set depth of the rheological front
    ierr = VecGetValues(M->radius_b,1,&i,&radius);CHKERRQ(ierr);
    RF->depth = P->radius - radius;

    // set pressure of the rheological front
    ierr = VecGetValues(M->pressure_b,1,&i,&pressure);CHKERRQ(ierr);
    RF->pressure = pressure;

    // set the mesh index in the middle of the magma ocean
    i_above_avg = i/2;
    ierr = VecGetValues(M->radius_b,1,&i_above_avg,&radius);CHKERRQ(ierr);
    RF->depth_above_avg = P->radius - radius;
    ierr = VecGetValues(M->pressure_b,1,&i_above_avg,&pressure);CHKERRQ(ierr);
    RF->pressure_above_avg = pressure;
    /* mask is one above the rheological front */
    ierr = average_by_mass_staggered( E, S->phi_s, mask_s, &RF->phi_above_avg ); CHKERRQ(ierr);
    ierr = average_by_mass_staggered( E, S->temp_s, mask_s, &RF->temperature_above_avg); CHKERRQ(ierr);

    /* only compute properties in the solid layer once the rheological front
       begins advancing through the mantle */
    if( i<numpts_s){
        // set the mesh index in the middle of the solid layer
        i_below_avg = (numpts_s - i)/2 + i;
        ierr = VecGetValues(M->radius_b,1,&i_below_avg,&radius);CHKERRQ(ierr);
        RF->depth_below_avg = P->radius - radius;
        ierr = VecGetValues(M->pressure_b,1,&i_below_avg,&pressure);CHKERRQ(ierr);
        RF->pressure_below_avg = pressure;
        ierr = invert_vec_mask( mask_s ); CHKERRQ(ierr);
        /* mask is now one below the rheological front */ 
        ierr = average_by_mass_staggered( E, S->phi_s, mask_s, &RF->phi_below_avg ); CHKERRQ(ierr);
        ierr = average_by_mass_staggered( E, S->temp_s, mask_s, &RF->temperature_below_avg); CHKERRQ(ierr);
    }

    ierr = VecDestroy( &mask_s ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

#if 0
static PetscErrorCode set_magma_ocean_crystal_fraction( Ctx *E )
{ 
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s;
    DM                da_s = E->da_s;
    Mesh              *M = &E->mesh;
    Parameters        *P = &E->parameters;
    CompositionalParameters *Comp = &P->compositional_parameters;
    Solution          *S = &E->solution;
    const PetscScalar *arr_phi_s, *arr_volume_s;
    PetscInt          index;
    PetscScalar       phi, phi_total, volume, volume_total;

    PetscFunctionBeginUser;

    /* determined by set_rheological_front */
    index = Comp->rheological_front_index;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,0,0,0);CHKERRQ(ierr);

    ierr = DMDAVecGetArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,M->volume_s,&arr_volume_s);CHKERRQ(ierr);

    phi_total = 0.0;
    volume_total = 0.0;
    /* FIXME: might break in parallel? */
    for(i=ilo_s; i<index; ++i){
        phi = arr_phi_s[i];
        volume = arr_volume_s[i];
        phi_total += phi * volume;
        volume_total += volume;
    }   
    phi_total /= volume_total;

    ierr = DMDAVecRestoreArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,M->volume_s,&arr_volume_s);CHKERRQ(ierr);

    /* FIXME: this will return 1.0 if index=0, but this might be
       confusing output to read? */
    Comp->mo_crystal_fraction = 1.0-phi_total;

    PetscFunctionReturn(0);

}

static PetscErrorCode set_magma_ocean_bridgmanite_fraction( Ctx *E )
{
    /* Derived as follows:
    *      phi * X_Brg^liq + (1-phi) * X_Brg^sol = X_Brg^0
    *      where phi is MELT fraction
    *  since we assume that Brg is the crystallising phase:
    *      X_Brg^sol = 1
    *  and rearranging, and noting that CRYSTAL fraction is
    *  1-MELT fraction, gives the desired expression for X_Brg^liq
    */

    Parameters               *P = &E->parameters;
    CompositionalParameters  *Comp = &P->compositional_parameters;
    PetscScalar              mo_crystal_fraction, XBrg;

    PetscFunctionBeginUser;

    mo_crystal_fraction = Comp->mo_crystal_fraction;
    XBrg = (Comp->X0Brg - mo_crystal_fraction) / (1.0 - mo_crystal_fraction);

    Comp->mo_bridgmanite_fraction = XBrg;

    PetscFunctionReturn(0);

}

static PetscErrorCode set_magma_ocean_mass_ratio( Ctx *E )
{
    /* Derived as follows:
    *      M_BSE = X_Brg M_Brg + X_res M_res
    *      and X_Brg + X_res = 1
    *  Therefore:
    *      M_BSE / M_Brg = X_Brg + (1-X_Brg) * (M_res / M_Brg)
    */

    Parameters               *P = &E->parameters;
    CompositionalParameters  *Comp = &P->compositional_parameters;
    PetscScalar              mo_mass_ratio;

    PetscFunctionBeginUser;

    mo_mass_ratio = (1.0-Comp->mo_bridgmanite_fraction) * Comp->muRes_muBrg;
    mo_mass_ratio += Comp->mo_bridgmanite_fraction;

    Comp->mo_mass_ratio = mo_mass_ratio;

    PetscFunctionReturn(0);

}
#endif
