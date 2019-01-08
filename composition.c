#include "composition.h"

static PetscErrorCode set_rheological_front( Ctx * );
static PetscErrorCode set_magma_ocean_crystal_fraction( Ctx * );
static PetscErrorCode set_magma_ocean_bridgmanite_fraction( Ctx * );
static PetscErrorCode set_magma_ocean_mass_ratio( Ctx * );

PetscErrorCode initialise_composition( Ctx *E )
{
    PetscErrorCode   ierr;

    PetscFunctionBeginUser;

    ierr = set_rheological_front( E ); CHKERRQ(ierr);
    ierr = set_magma_ocean_crystal_fraction( E ); CHKERRQ(ierr);
    ierr = set_magma_ocean_bridgmanite_fraction( E ); CHKERRQ(ierr);
    ierr = set_magma_ocean_mass_ratio( E ); CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode set_rheological_front( Ctx *E )
{
    PetscErrorCode   ierr;
    PetscInt         i,ilo_s,ihi_s,w_s;
    DM               da_s=E->da_s;
    Mesh             *M = &E->mesh;
    Solution         *S = &E->solution;
    Parameters       *P = &E->parameters;
    CompositionalParameters *Comp = &P->compositional_parameters;
    const PetscScalar *arr_phi_s;
    PetscScalar phi, radius, pressure;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = DMDAVecGetArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);

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
        if( phi < Comp->rheological_front_phi )
            break;
    }

    ierr = DMDAVecRestoreArrayRead(da_s,S->phi_s,&arr_phi_s);CHKERRQ(ierr);

    // store in CompositionalParameters structure for later use
    Comp->rheological_front_index = i;

    // set depth of the rheological front
    VecGetValues(M->radius_b,1,&i,&radius);CHKERRQ(ierr);
    Comp->rheological_front_depth = P->radius - radius;

    // set pressure of the rheological front
    VecGetValues(M->pressure_b,1,&i,&pressure);CHKERRQ(ierr);
    Comp->rheological_front_pressure = pressure;

    PetscFunctionReturn(0);

}

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
    Parameters               *P = &E->parameters;
    CompositionalParameters  *Comp = &P->compositional_parameters;
    PetscScalar              mo_mass_ratio;

    PetscFunctionBeginUser;

    mo_mass_ratio = (1.0-Comp->mo_bridgmanite_fraction) * Comp->muRes_muBrg;
    mo_mass_ratio += Comp->mo_bridgmanite_fraction;

    Comp->mo_mass_ratio = mo_mass_ratio;

    PetscFunctionReturn(0);

}
