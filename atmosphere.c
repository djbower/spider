#include "atmosphere.h"

static PetscScalar get_total_mass( Ctx * );
static PetscScalar get_liquid_mass( Ctx * );
#if 0
static PetscScalar get_solid_mass( Ctx * );
#endif
static PetscScalar get_partialP_water( PetscScalar );
static PetscScalar get_partialP_carbon( PetscScalar );
static PetscErrorCode set_partialP( Ctx * );
static PetscScalar get_atmosphere_mass( Volatile * );
static PetscScalar get_optical_depth( Volatile * );

PetscScalar get_emissivity( Ctx *E )
{
    PetscScalar    mass_liq, mass_tot, frac_liq, frac_sol;
    PetscScalar    M0, mass_atm, emiss, tau;
    //PetscScalar    mass_sol;
    Volatile       *W = &E->vol_water;
    Volatile       *C = &E->vol_carbon;

    /* this uses current concentration of volatiles to compute
       partial pressures */
    /* updates W->X and C->X */
    set_partialP( E );

    /* now solve mass balance of volatiles */
    /* need mass fraction of melt and solid */
    mass_liq = get_liquid_mass( E );
    //mass_sol = get_solid_mass( E );
    mass_tot = get_total_mass( E );

    /* total mass will change with time.  Unphysical, but
       hopefully won't change by much.  So for simplicitly, let's
       just compute fraction of liquid and total to total mass.  Then
       we can use the 'fixed' initial mass determined from the AW EOS
       at the start since this is unchanging.  This essentially enforces
       that the liquid and solid mass must total the total mass, and
       this total mass cannot change with time */
    frac_liq = mass_liq / mass_tot;
    frac_sol = 1.0 - frac_liq;
    M0 = W->M0; // same for all volatiles since initial/total mass of planet

    /* water */
    mass_atm = get_atmosphere_mass( W );
    W->X = W->X0M0 - mass_atm;
    W->X /= W->kdist*M0*frac_sol+M0*frac_liq;

    /* carbon */
    mass_atm = get_atmosphere_mass( C );
    C->X = C->X0M0 - mass_atm;
    C->X /= C->kdist*M0*frac_sol+M0*frac_liq;

    /* optical depth */
    tau = get_optical_depth( W );
    tau += get_optical_depth( C );

    /* emissivity */
    emiss = 2.0 / (tau + 2.0);

    return emiss;
}

static PetscScalar get_optical_depth( Volatile *V )
{
    PetscScalar tau, mass_atm;

    mass_atm = get_atmosphere_mass( V );
    tau = 3.0 * mass_atm / (8.0*PETSC_PI*PetscSqr(RADIUS));
    // note negative gravity!
    tau *= PetscSqrtScalar( V->kabs*-GRAVITY/(3.0*V->P0) );

    return tau;
}

static PetscScalar get_atmosphere_mass( Volatile *V )
{
    PetscScalar mass;

    mass = (4.0*PETSC_PI*PetscSqr(RADIUS))/GRAVITY;
    mass *= V->P;
    mass *= -1.0; // because GRAVITY is negative by definition

    return mass;
}


PetscErrorCode set_initial_water( Ctx *E )
{
    Mesh           *M = &E->mesh;
    Volatile       *W = &E->vol_water;

    PetscFunctionBeginUser;

    W->kdist = KDIST_WATER;
    W->X0 = X0_WATER_WT_PERCENT / 100.0;
    W->X = W->X0;
    /* setting the initial mass of the magma ocean according
       to the density profile used to compute the reference
       pressure profile */
    W->M0 = M->mass0;
    W->X0M0 = W->X0 * W->M0;
    W->kabs = KABS_WATER;
    W->P0 = P0;

    PetscFunctionReturn(0);
}

PetscErrorCode set_initial_carbon( Ctx *E )
{
    Mesh           *M = &E->mesh;
    Volatile       *C = &E->vol_carbon;

    PetscFunctionBeginUser;

    C->kdist = KDIST_CARBON;
    C->X0 = X0_CARBON_WT_PERCENT / 100.0;
    C->X = C->X0;
    /* setting the initial mass of the magma ocean according
       to the density profile used to compute the reference
       pressure profile */
    C->M0 = M->mass0;
    C->X0M0 = C->X0 * C->M0;
    C->kabs = KABS_CARBON;
    C->P0 = P0;

    PetscFunctionReturn(0);
}

/* TODO: since the masses of the liquid and solid parts of the mantle
   depend on the density structure, there is no guarantee that the
   mass of the overall mantle remains constant.  This will influence
   the volatile exchange, but probably only slightly. */

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

#if 0
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
#endif

/* partial pressures of volatiles
   x_vol is mass fraction of volatiles "in the magma" (Lebrun)
   according to my derivation, "in the magma" is for melt
   fractions higher than the rheological transition */

static PetscErrorCode set_partialP( Ctx *E )
{
    Volatile       *W = &E->vol_water;
    Volatile       *C = &E->vol_carbon;

    PetscFunctionBeginUser;

    W->P = get_partialP_water( W->X );
    C->P = get_partialP_carbon( C->X );

    PetscFunctionReturn(0);

}

static PetscScalar get_partialP_water( PetscScalar x_vol )
{
    PetscScalar p;

    /* Lebrun et al. (2013) eqn, 16 */
    p = x_vol / 6.8E-8;
    p = PetscPowScalar( p, 1.0/0.7 );

    return p;
}

static PetscScalar get_partialP_carbon( PetscScalar x_vol )
{
    PetscScalar p;

    /* Lebrun et al. (2013) eqn. 17 */
    p = x_vol / 4.4E-12;

    return p;
}
