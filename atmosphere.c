#include "atmosphere.h"

static PetscErrorCode set_Mliq( Ctx * );
static PetscErrorCode set_Msol( Ctx * );
static PetscErrorCode set_dMliqdt( Ctx * );
static PetscScalar get_atmosphere_mass( PetscScalar );
static PetscScalar get_optical_depth( PetscScalar, PetscScalar );
#if 0
static PetscScalar get_partialP_H2O( PetscScalar );
#endif
static PetscErrorCode set_pCO2( Atmosphere * );
static PetscErrorCode set_dpCO2dx( Atmosphere * );

/* general functions */

PetscErrorCode set_emissivity( Ctx *E )
{
    Atmosphere *A = &E->atmosphere;

    PetscFunctionBeginUser;

    /* CO2 */
    set_pCO2( A );
    A->m0 = get_atmosphere_mass( A->p0 );
    A->tau0 = get_optical_depth( A->m0, CO2_KABS );

    /* H2O */
    // FIXME: DJB ignore H2O for the time being
    //set_pH2O( A );
    //A->m1 = get_atmosphere_mass( A->p1 );
    //A->tau1 = get_optical_depth( A->m1, H2O_KABS );
    A->tau1 = 0.0;

    /* total */
    A->tau = A->tau0 + A->tau1;
    A->emissivity = 2.0 / (A->tau + 2.0);

    /* for the case with no coupled atmospheric growth (H2O_INITIAL
       and CO2_INITIAL both zero), we can scale by the EMISSIVITY
       to give a grey-body model with constant emissivity.  Since
       above result will give emissivity=1.0 with no atmospheric
       growth (but see below, now over-written to avoid negative
       emissivity which could arise if the user enters a negative
       initial volatile content, presumably in error!) */
    /* this is really a check, since it makes no sense to have a
       negative initial volatile content, and therefore always
       switch to a constant emissivity grey-body if the volatile
       content is either zero or negative */
    if (CO2_INITIAL <= 0.0 && H2O_INITIAL <= 0.0){
        A->emissivity = EMISSIVITY;
    }

    PetscFunctionReturn(0);

}

static PetscScalar get_optical_depth( PetscScalar mass_atm, PetscScalar kabs )
{
    PetscScalar tau;

    tau = 3.0 * mass_atm / (8.0*PETSC_PI*PetscSqr(RADIUS));
    // note negative gravity!
    tau *= PetscSqrtScalar( kabs*-GRAVITY/(3.0*P0) );

    return tau;
}

static PetscScalar get_atmosphere_mass( PetscScalar p )
{
    PetscScalar mass_atm;

    mass_atm = 4.0*PETSC_PI*PetscSqr(RADIUS) * p / -GRAVITY;

    return mass_atm;

}

static PetscErrorCode set_Mliq( Ctx *E )
{
    PetscErrorCode ierr;
    Atmosphere     *A = &E->atmosphere;
    Solution       *S = &E->solution;
    Mesh           *M = &E->mesh;
    Vec            mass_s;

    PetscFunctionBeginUser;

    // Mliq = sum[ phi*dm ]
    ierr = VecDuplicate( S->phi_s, &mass_s ); CHKERRQ(ierr);
    ierr = VecCopy( S->phi_s, mass_s ); CHKERRQ(ierr);

    ierr = VecPointwiseMult( mass_s, mass_s, M->mass_s ); CHKERRQ(ierr);
    ierr = VecSum( mass_s, &A->Mliq );

    VecDestroy( &mass_s );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_Msol( Ctx *E )
{
    PetscErrorCode ierr;
    Atmosphere     *A = &E->atmosphere;
    Solution       *S = &E->solution;
    Mesh           *M = &E->mesh;
    Vec            mass_s;

    PetscFunctionBeginUser;

    // Msol = sum[ (1-phi)*dm ]
    ierr = VecDuplicate( S->phi_s, &mass_s ); CHKERRQ(ierr);
    ierr = VecCopy( S->phi_s, mass_s ); CHKERRQ(ierr);
    ierr = VecScale( mass_s, -1.0 ); CHKERRQ(ierr);
    ierr = VecShift( mass_s, 1.0 ); CHKERRQ(ierr);

    ierr = VecPointwiseMult( mass_s, mass_s, M->mass_s ); CHKERRQ(ierr);
    ierr = VecSum( mass_s, &A->Msol );

    VecDestroy( &mass_s );

    PetscFunctionReturn(0);

}

static PetscErrorCode set_dMliqdt( Ctx *E )
{
    PetscErrorCode    ierr;
    PetscInt          i,ilo_s,ihi_s,w_s;
    DM                da_s = E->da_s;
    Atmosphere        *A = &E->atmosphere;
    Solution          *S = &E->solution;
    Mesh              *M = &E->mesh;
    Vec               result_s;
    PetscScalar       *arr_result_s;
    const PetscScalar *arr_dSdt_s, *arr_fusion_s, *arr_fwtl_s, *arr_fwts_s, *arr_phi_s, *arr_mass_s;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = VecDuplicate( S->dSdt_s, &result_s); CHKERRQ(ierr);
    ierr = VecCopy( S->dSdt_s, result_s ); CHKERRQ(ierr);

    ierr = DMDAVecGetArrayRead(da_s,M->mass_s,&arr_mass_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,S->dSdt_s,&arr_dSdt_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fusion_s,&arr_fusion_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fwtl_s,&arr_fwtl_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->fwts_s,&arr_fwts_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_s,S->phi_s,&arr_phi_s); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,result_s,&arr_result_s); CHKERRQ(ierr);

    for(i=ilo_s; i<ihi_s; ++i){

        arr_result_s[i] = arr_dSdt_s[i] * arr_mass_s[i];
        arr_result_s[i] /= arr_fusion_s[i];

        if (arr_phi_s[i] > 0.5){
            arr_result_s[i] *= ( 1.0 - arr_fwtl_s[i] );
        }

        else if (arr_phi_s[i] <=0.5){
            arr_result_s[i] *= arr_fwts_s[i];
        }
    }

    ierr = DMDAVecRestoreArrayRead(da_s,M->mass_s,&arr_mass_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,S->dSdt_s, &arr_dSdt_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fusion_s,&arr_fusion_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fwtl_s,&arr_fwtl_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->fwts_s,&arr_fwts_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da_s,S->phi_s,&arr_phi_s); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,result_s, &arr_result_s); CHKERRQ(ierr);

    ierr = VecSum( result_s, &A->dMliqdt ); CHKERRQ(ierr);

    VecDestroy( &result_s );

    PetscFunctionReturn(0);

}


PetscErrorCode set_dx0dt( Ctx *E )
{
   /* update for dissolved CO2 content in the magma ocean */
    Atmosphere     *A = &E->atmosphere;
    Mesh           *M = &E->mesh;
    PetscScalar    num, den;

    PetscFunctionBeginUser;

    set_dpCO2dx( A );
    set_Mliq( E );
    set_Msol( E );
    set_dMliqdt( E );

    num = A->x0 * (CO2_KDIST-1.0) * A->dMliqdt;
    den = CO2_KDIST * M->mass0 + (1.0-CO2_KDIST) * A->Mliq;
    den += ( 4.0*PETSC_PI*PetscSqr(RADIUS) / -GRAVITY ) * A->dp0dx;

    A->dx0dt = num / den;

    PetscFunctionReturn(0);
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
    Atmosphere  *A = &E->atmosphere;
    Mesh        *M = &E->mesh;

    PetscScalar x0;

    x0 = 4.0 * PETSC_PI * PetscSqr(RADIUS);
    x0 /= -GRAVITY * M->mass0 * CO2_HENRY;
    x0 += 1.0;
    x0 = 1.0 / x0;
    x0 *= CO2_INITIAL;

    // update struct
    A->x0 = x0;
    A->x0init = CO2_INITIAL;

    // need to return to add to augmented vector
    return x0;
}

static PetscErrorCode set_pCO2( Atmosphere *A )
{
    /* partial pressure of CO2 */
    /* see magma_ocean_notes.tex */
    /* Lebrun et al. (2013) eqn. 17
       Massol et al. (2016) eqn. 18
       Salvador et al. (2017) eqn. C4 */

    PetscScalar xmf = A->x0 / 100.0 // x as mass fraction

    PetscFunctionBeginUser;

    /* x must be mass fraction here, otherwise the units of A are
       not correct */
    A->p0 = xmf / CO2_HENRY;

    PetscFunctionReturn(0);

}

static PetscErrorCode set_dpCO2dx( Atmosphere *A )
{
    /* dpCO2/dx. i.e., derivative of partial pressure (Pa) with respect
       to concentration (wt %) */

    PetscFunctionBeginUser;

    A->dp0dx = 1.0 / CO2_HENRY;

    PetscFunctionReturn(0);

}
