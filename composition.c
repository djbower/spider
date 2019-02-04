#include "composition.h"
//#include "util.h"

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
#endif

PetscScalar get_BSE_Brg_mass_ratio( PetscScalar Brg_fraction, PetscScalar Res_Brg_mass_ratio )
{
    /* Derived as follows:
    *      M_BSE = X_Brg M_Brg + X_res M_res
    *      and X_Brg + X_res = 1
    *  Therefore:
    *      M_BSE / M_Brg = X_Brg + (1-X_Brg) * (M_res / M_Brg)
    */

    PetscScalar mass_ratio;

    PetscFunctionBeginUser;

    mass_ratio = (1.0-Brg_fraction) * Res_Brg_mass_ratio;
    mass_ratio += Brg_fraction;

    return mass_ratio;

}
