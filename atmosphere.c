#include "atmosphere.h"

static PetscScalar get_total_mass( Ctx * );
static PetscScalar get_liquid_mass( Ctx * );
static PetscScalar get_solid_mass( Ctx * );


PetscErrorCode atmosphere_test( Ctx *E )
{
    PetscScalar t1, t2, t3;

    PetscFunctionBeginUser;

    t1 = get_total_mass( E );
    t2 = get_liquid_mass( E );
    t3 = get_solid_mass( E );

    PetscFunctionReturn(0);

}


static PetscScalar get_total_mass( Ctx *E )
{
    PetscErrorCode ierr;
    Solution       *S = &E->solution;
    Mesh           *M = &E->mesh;
    Vec            mass_s;
    PetscScalar    mass_total;

    ierr = VecDuplicate( S->rho_s, &mass_s ); CHKERRQ(ierr);
    ierr = VecCopy( S->rho_s, mass_s ); CHKERRQ(ierr);

    /* remember that volume_s is scaled volume, i.e., without the
       4*PI prefactor */
    ierr = VecPointwiseMult( mass_s, mass_s, M->volume_s );
    ierr = VecSum( mass_s, &mass_total );

    /* now take account of 4*PI factor */
    mass_total *= 4.0 * PETSC_PI;

    /* clean up */
    VecDestroy( &mass_s );

    return mass_total;

}

static PetscScalar get_liquid_mass( Ctx *E )
{
    PetscErrorCode ierr;
    DM             da_s = E->da_s;
    Solution       *S = &E->solution;
    Mesh           *M = &E->mesh;
    Vec            mass_s;
    PetscScalar    *arr_mass_s;
    PetscScalar    mass_liquid;
    PetscInt       i, ilo_s, ihi_s, w_s;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = VecDuplicate( S->phi_s, &mass_s ); CHKERRQ(ierr);
    ierr = VecCopy( S->phi_s, mass_s ); CHKERRQ(ierr);

    /* create mask of 1 and 0's to denote liquid regions
       for a rough estimate liquid regions are defined as melt regions
       above the rheological transition */
    ierr = DMDAVecGetArray(da_s,mass_s,&arr_mass_s); CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        if (arr_mass_s[i] > PHI_CRITICAL){
            arr_mass_s[i] = 1.0;
        }
        else{
            arr_mass_s[i] = 0.0;
        }
    }
    ierr = DMDAVecRestoreArray(da_s,mass_s,&arr_mass_s); CHKERRQ(ierr);

    ierr = VecPointwiseMult( mass_s, mass_s, S->rho_s ); CHKERRQ(ierr);
    ierr = VecPointwiseMult( mass_s, mass_s, M->volume_s ); CHKERRQ(ierr);

    ierr = VecSum( mass_s, &mass_liquid );

    mass_liquid *= 4.0 * PETSC_PI;

    VecDestroy( &mass_s );

    return mass_liquid;

}

static PetscScalar get_solid_mass( Ctx *E )
{
    PetscErrorCode ierr;
    DM             da_s = E->da_s;
    Solution       *S = &E->solution;
    Mesh           *M = &E->mesh;
    Vec            mass_s;
    PetscScalar    *arr_mass_s;
    PetscScalar    mass_solid;
    PetscInt       i, ilo_s, ihi_s, w_s;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = VecDuplicate( S->phi_s, &mass_s ); CHKERRQ(ierr);
    ierr = VecCopy( S->phi_s, mass_s ); CHKERRQ(ierr);

    /* create mask of 1 and 0's to denote liquid regions
       for a rough estimate liquid regions are defined as melt regions
       above the rheological transition */
    ierr = DMDAVecGetArray(da_s,mass_s,&arr_mass_s); CHKERRQ(ierr);
    for(i=ilo_s; i<ihi_s; ++i){
        if (arr_mass_s[i] <= PHI_CRITICAL){
            arr_mass_s[i] = 1.0;
        }
        else{
            arr_mass_s[i] = 0.0;
        }
    }
    ierr = DMDAVecRestoreArray(da_s,mass_s,&arr_mass_s); CHKERRQ(ierr);

    ierr = VecPointwiseMult( mass_s, mass_s, S->rho_s ); CHKERRQ(ierr);
    ierr = VecPointwiseMult( mass_s, mass_s, M->volume_s ); CHKERRQ(ierr);

    ierr = VecSum( mass_s, &mass_solid );

    mass_solid *= 4.0 * PETSC_PI;

    VecDestroy( &mass_s );

    return mass_solid;

}
