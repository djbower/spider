#include "atmosphere.h"

static PetscScalar get_liquid_mass( Ctx * );
static PetscScalar get_solid_mass( Ctx * );
static PetscScalar get_atmosphere_mass( Ctx *, PetscScalar, PetscScalar, PetscScalar );
static PetscScalar get_optical_depth( PetscScalar, PetscScalar );
#if 0
static PetscScalar get_partialP_H2O( PetscScalar );
static PetscScalar get_partialP_CO2( PetscScalar );
#endif

PetscScalar get_emissivity( Ctx *E, PetscScalar X0, PetscScalar X1 )
{
    PetscScalar m0, tau0, m1, tau1, tau, emissivity;

    /* CO2 */
    m0 = get_atmosphere_mass( E, CO2_INITIAL, X0, CO2_KDIST );
    tau0 = get_optical_depth( m0, CO2_KABS );

    /* H2O */
    m1 = get_atmosphere_mass( E, H2O_INITIAL, X1, H2O_KDIST );
    tau1 = get_optical_depth( m1, H2O_KABS );

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

static PetscScalar get_atmosphere_mass( Ctx *E, PetscScalar Xinit, PetscScalar Xvol, PetscScalar Kdist )
{
    Mesh M = E->mesh;
    PetscScalar mass_sol, mass_liq, mass_atm;

    mass_sol = get_solid_mass( E );
    mass_liq = get_liquid_mass( E );

    mass_atm = Xinit*M.mass0; // initial total mass
    mass_atm -= Xvol*(Kdist*mass_sol+mass_liq); // minus content in magma ocean (solid and liquid)

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

static PetscScalar get_solid_mass( Ctx *E )
{
    PetscErrorCode ierr;
    Solution       *S = &E->solution;
    Mesh           *M = &E->mesh;
    Vec            mass_s;
    PetscScalar    mass_solid;

    // Msol = sum[ (1-phi)*dm ]
    ierr = VecDuplicate( S->phi_s, &mass_s ); CHKERRQ(ierr);
    ierr = VecCopy( S->phi_s, mass_s ); CHKERRQ(ierr);
    ierr = VecScale( mass_s, -1.0); CHKERRQ(ierr);
    ierr = VecShift( mass_s, 1.0 ); CHKERRQ(ierr);

    ierr = VecPointwiseMult( mass_s, mass_s, M->mass_s ); CHKERRQ(ierr);
    ierr = VecSum( mass_s, &mass_solid );

    VecDestroy( &mass_s );

    return mass_solid;

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
    A = 4.4E-12;

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

static PetscScalar get_partialP_CO2( PetscScalar x_vol )
{
    PetscScalar p;

    /* see magma_ocean_notes.tex */

    /* Lebrun et al. (2013) eqn. 17
       Massol et al. (2016) eqn. 18
       Salvador et al. (2017) eqn. C4 */

    /* where x_vol is mass fraction */
    p = x_vol / 4.4E-12;

    return p;
}
#endif
