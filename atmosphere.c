#include "atmosphere.h"

static PetscScalar get_liquid_mass( Ctx * );
#if 0
static PetscScalar get_solid_mass( Ctx * );
#endif
static PetscScalar get_atmosphere_mass( Ctx *, PetscScalar );
static PetscScalar get_optical_depth( PetscScalar, PetscScalar );
#if 0
static PetscScalar get_partialP_H2O( PetscScalar );
#endif
static PetscScalar get_pCO2( PetscScalar );
static PetscScalar get_dpCO2dx( PetscScalar );

/* general functions */

PetscScalar get_emissivity( Ctx *E, PetscScalar x0, PetscScalar x1 )
{
    PetscScalar m0, p0, tau0, tau1;
    PetscScalar tau, emissivity;

    /* CO2 */
    p0 = get_pCO2( x0 );
    m0 = get_atmosphere_mass( E, p0 );
    tau0 = get_optical_depth( m0, CO2_KABS );

    /* H2O */
    // FIXME: DJB ignore H2O for the time being
    //m1 = get_atmosphere_mass( E, H2O_INITIAL, x1, H2O_KDIST );
    //tau1 = get_optical_depth( m1, H2O_KABS );
    tau1 = 0.0;

    /* total */
    tau = tau0 + tau1;
    emissivity = 2.0 / (tau + 2.0);

    return emissivity;
}

static PetscScalar get_optical_depth( PetscScalar mass_atm, PetscScalar kabs )
{
    PetscScalar tau;

    tau = 3.0 * mass_atm / (8.0*PETSC_PI*PetscSqr(RADIUS));
    // note negative gravity!
    tau *= PetscSqrtScalar( kabs*-GRAVITY/(3.0*P0) );

    return tau;
}

static PetscScalar get_atmosphere_mass( Ctx *E, PetscScalar p )
{
    PetscScalar mass_atm;

    mass_atm = 4.0*PETSC_PI*PetscSqr(RADIUS) * p / -GRAVITY;

    return mass_atm;

}

static PetscScalar get_liquid_mass( Ctx *E )
{
    PetscErrorCode ierr;
    Solution       *S = &E->solution;
    Mesh           *M = &E->mesh;
    Vec            mass_s;
    PetscScalar    mass_liquid;

    // Mliq = sum[ phi*dm ]
    ierr = VecDuplicate( S->phi_s, &mass_s ); CHKERRQ(ierr);
    ierr = VecCopy( S->phi_s, mass_s ); CHKERRQ(ierr);

    ierr = VecPointwiseMult( mass_s, mass_s, M->mass_s ); CHKERRQ(ierr);
    ierr = VecSum( mass_s, &mass_liquid );

    VecDestroy( &mass_s );

    return mass_liquid;

}

PetscScalar get_dX0dt( Ctx *E, PetscScalar X0, Vec dSdt_s )
{

   /* update for CO2 content */

    PetscErrorCode ierr;
    Solution       *S = &E->solution;
    Mesh           *M = &E->mesh;
    PetscScalar    A, num, den, dX0dt, dphidtdm, mass_liq;
    Vec            dphidm_s;

    ierr = VecDuplicate( dSdt_s, &dphidm_s); CHKERRQ(ierr);
    ierr = VecCopy( dSdt_s, dphidm_s ); CHKERRQ(ierr);
    ierr = VecPointwiseDivide( dphidm_s, dphidm_s, S->fusion_s ); CHKERRQ(ierr);
    ierr = VecPointwiseMult( dphidm_s, dphidm_s, M->mass_s ); CHKERRQ(ierr);

    mass_liq = get_liquid_mass( E );
    A = 1.0/CO2_HENRY;

    ierr = VecSum( dphidm_s, &dphidtdm );

    num = X0 * (CO2_KDIST-1.0) * dphidtdm;
    den = CO2_KDIST*M->mass0 + (1.0-CO2_KDIST)*mass_liq;
    den += 4.0*PETSC_PI*PetscSqr(RADIUS) * A / -GRAVITY;

    dX0dt = num / den;

    return dX0dt;
}

#if 0
/* partial pressures of volatiles
   x_vol is mass fraction of volatiles "in the magma" (Lebrun)
   according to my derivation, "in the magma" is for melt
   fractions higher than the rheological transition */

static PetscScalar get_partialP_H2O( PetscScalar x_vol )
{
    PetscScalar p;

    /* what are the units?  mass fraction / weight percent? */
    /* Massol et al says ``with xH2O and XCO2 being respectively
       the mass fraction of water and CO2 dissolved in the melt
       respectively expressed in wt% and ppm */

    /* Lebrun et al. (2013) eqn. 16
       Massol et al. (2016) eqn. 17
       Salvador et al. (2017) eqn. C3 */
    p = x_vol / 6.8E-8;
    p = PetscPowScalar( p, 1.0/0.7 );

    return p;
}
#endif

////////////////////////////
/* CO2 specific functions */
////////////////////////////

PetscScalar get_initial_xCO2( Ctx *E )
{

    /* solve mass balance to get initial volatile content of the
       liquid */

    Mesh        *M = &E->mesh;
    PetscScalar result;

    result = 4.0 * PETSC_PI * PetscSqr(RADIUS);
    result /= -GRAVITY * M->mass0 * CO2_HENRY;
    result += 1.0;
    result = 1.0 / result;
    result *= CO2_INITIAL;

    return result;
}


static PetscScalar get_pCO2( PetscScalar x )
{
    /* partial pressure of CO2 */
    /* see magma_ocean_notes.tex */
    /* Lebrun et al. (2013) eqn. 17
       Massol et al. (2016) eqn. 18
       Salvador et al. (2017) eqn. C4 */

    PetscScalar p;

    /* where x is the concentration in the liquid in wt % */
    p = x / CO2_HENRY;

    return p;
}

static PetscScalar get_dpCO2dx( PetscScalar x )
{
    /* dpCO2/dx. i.e., derivative of partial pressure (Pa) with respect
       to concentration (wt %) */

    PetscScalar dpdx;

    dpdx = 1.0 / CO2_HENRY;

    return dpdx;
}
