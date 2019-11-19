#include "reaction.h"

static PetscScalar get_psurf_exponent( const ReactionParameters * );
static PetscScalar get_reaction_quotient( const ReactionParameters *, const Atmosphere *, PetscInt );
static PetscScalar get_reaction_quotient_time_derivative( const ReactionParameters *, const Atmosphere *, const AtmosphereParameters *, PetscInt );
static PetscScalar get_log10_equilibrium_constant( const ReactionParameters *, PetscScalar, const Constants * );
static PetscScalar get_dlog10KdT( const ReactionParameters *, PetscScalar, const Constants * );

/* Note: this could logically be included in parameters.c, but that file was getting crowded */

/* A named reaction */
PetscErrorCode ReactionParametersCreateMethane1(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "methane1";
  reaction_parameters->n_volatiles = 3;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // CO2
  reaction_parameters->stoichiometry[1] = -2.0;  // H2
  reaction_parameters->stoichiometry[2] = 1.0; // CH4
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = -16276;
  reaction_parameters->Keq_coeffs[1] = -5.4738;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 1.0;

  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"CO2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[1] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"CH4",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[2] = v;
  }
  for (i=0; i<reaction_parameters->n_volatiles; ++i) if (reaction_parameters->volatiles[i] == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Didn't find required volatiles for reaction %s",reaction_parameters->type);
  PetscFunctionReturn(0);
}

/* A named reaction */
PetscErrorCode ReactionParametersCreateAmmonia1(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "ammonia1";
  reaction_parameters->n_volatiles = 3;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // N2
  reaction_parameters->stoichiometry[1] = -3.0;  // H2
  reaction_parameters->stoichiometry[2] = 2.0; // NH3
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 5331.9;
  reaction_parameters->Keq_coeffs[1] = -11.884;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0;

  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"N2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[1] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"NH3",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[2] = v;
  }
  for (i=0; i<reaction_parameters->n_volatiles; ++i) if (reaction_parameters->volatiles[i] == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Didn't find required volatiles for reaction %s",reaction_parameters->type);
  PetscFunctionReturn(0);
}

/* A named reaction */
/* This one is used for testing, since it assumes constant fO2 and K */
/* K = pH2/pH2O = 10**2.0 = 100 */
/* simplified version used in thesis work of Taylor Kutra (2019) */
PetscErrorCode ReactionParametersCreateSimpleWater1(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "simplewater1";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // H2O
  reaction_parameters->stoichiometry[1] = 1.0;  // H2
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 0;
  reaction_parameters->Keq_coeffs[1] = 2.0;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.0;

  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2O",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[1] = v;
  }
  for (i=0; i<reaction_parameters->n_volatiles; ++i) if (reaction_parameters->volatiles[i] == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Didn't find required volatiles for reaction %s",reaction_parameters->type);
  PetscFunctionReturn(0);
}

/* A named reaction */
/* This one is used for testing, since it assumes constant fO2 and K */
/* K = pH2/pH2O = 10**1.0 = 10 */
PetscErrorCode ReactionParametersCreateSimpleWater2(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "simplewater2";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // H2O
  reaction_parameters->stoichiometry[1] = 1.0;  // H2
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 0;
  reaction_parameters->Keq_coeffs[1] = 1.0;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.0;

  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2O",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[1] = v;
  }
  for (i=0; i<reaction_parameters->n_volatiles; ++i) if (reaction_parameters->volatiles[i] == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Didn't find required volatiles for reaction %s",reaction_parameters->type);
  PetscFunctionReturn(0);
}

/* A named reaction */
/* This one is used for testing, since it assumes constant fO2 and K */
/* K = pH2/pH2O = 10**0.0 = 1 */
PetscErrorCode ReactionParametersCreateSimpleWater3(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "simplewater3";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // H2O
  reaction_parameters->stoichiometry[1] = 1.0;  // H2
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = 0;
  reaction_parameters->Keq_coeffs[1] = 0.0;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.0;

  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2O",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[1] = v;
  }
  for (i=0; i<reaction_parameters->n_volatiles; ++i) if (reaction_parameters->volatiles[i] == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Didn't find required volatiles for reaction %s",reaction_parameters->type);
  PetscFunctionReturn(0);
}

/* A named reaction */
PetscErrorCode ReactionParametersCreateWater1(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "water1";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // H2O
  reaction_parameters->stoichiometry[1] = 1.0;  // H2
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = -12794;
  reaction_parameters->Keq_coeffs[1] = 2.7768;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.5;

  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2O",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"H2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[1] = v;
  }
  for (i=0; i<reaction_parameters->n_volatiles; ++i) if (reaction_parameters->volatiles[i] == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Didn't find required volatiles for reaction %s",reaction_parameters->type);
  PetscFunctionReturn(0);
}

/* A named reaction */
PetscErrorCode ReactionParametersCreateCarbonDioxide1(ReactionParameters* reaction_parameters_ptr, const AtmosphereParameters *Ap)
{
  PetscErrorCode     ierr;
  PetscInt           i,v;
  PetscBool          flg;
  ReactionParameters reaction_parameters;

  PetscFunctionBeginUser;
  ierr = PetscMalloc1(1,reaction_parameters_ptr);CHKERRQ(ierr);
  reaction_parameters = *reaction_parameters_ptr;
  reaction_parameters->type = "carbondioxide1";
  reaction_parameters->n_volatiles = 2;
  ierr = PetscMalloc3(2,&reaction_parameters->Keq_coeffs,reaction_parameters->n_volatiles,&reaction_parameters->stoichiometry,reaction_parameters->n_volatiles,&reaction_parameters->volatiles);CHKERRQ(ierr);
  reaction_parameters->stoichiometry[0] = -1.0;  // CO2
  reaction_parameters->stoichiometry[1] = 1.0;  // CO
  /* equilibrium constant coefficients */
  reaction_parameters->Keq_coeffs[0] = -14787;
  reaction_parameters->Keq_coeffs[1] = 4.5472;
  /* fO2 stoichiometry */
  reaction_parameters->fO2_stoichiometry = 0.5;

  for (i=0; i<reaction_parameters->n_volatiles; ++i) reaction_parameters->volatiles[i] = -1; /* error value */
  for (v=0; v<Ap->n_volatiles; ++v) {
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"CO2",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[0] = v;
    ierr = PetscStrcmp(Ap->volatile_parameters[v].prefix,"CO",&flg);CHKERRQ(ierr);
    if (flg) reaction_parameters->volatiles[1] = v;
  }
  for (i=0; i<reaction_parameters->n_volatiles; ++i) if (reaction_parameters->volatiles[i] == -1) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Didn't find required volatiles for reaction %s",reaction_parameters->type);
  PetscFunctionReturn(0);
}

PetscErrorCode ReactionParametersDestroy(ReactionParameters* reaction_parameters_ptr)
{
  PetscErrorCode     ierr;

  PetscFunctionBeginUser;
  ierr = PetscFree3((*reaction_parameters_ptr)->Keq_coeffs,(*reaction_parameters_ptr)->stoichiometry,(*reaction_parameters_ptr)->volatiles);CHKERRQ(ierr);
  ierr = PetscFree(*reaction_parameters_ptr);CHKERRQ(ierr); /* must be last */
  *reaction_parameters_ptr = NULL;
  PetscFunctionReturn(0);
}

/* TODO: it will probably be helpful for debugging and output to store
   the equilibrium constant and its derivative, but currently we do not
   have an obvious place to store them since they are not parameters
   but rather time-dependent quantities */

/* Compute equilibrium constant */
static PetscScalar get_log10_equilibrium_constant( const ReactionParameters * reaction_parameters_ptr, PetscScalar temp, const Constants *C )
{
    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscScalar        log10Keq;

    temp *= C->TEMP;

    /* log10Keq = a/T + b is a standard form for computing an equilibrium constant */
    /* e.g., Schaefer and Fegley (2017) */
    log10Keq = reaction_parameters->Keq_coeffs[0] / temp + reaction_parameters->Keq_coeffs[1];

    return log10Keq;

}

/* Derivative of log10 equilibrium constant with respect to temperature */
PetscScalar get_dlog10KdT( const ReactionParameters * reaction_parameters_ptr, PetscScalar temp, const Constants *C )
{
    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscScalar        dlog10KdT;

    temp *= C->TEMP;

    dlog10KdT = -reaction_parameters->Keq_coeffs[0] / PetscPowScalar( temp, 2.0 );

    /* TODO: check, must non-dimensionalise, since we used a scaled
       (dimensional) temperature */
    dlog10KdT *= C->TEMP;

    return dlog10KdT;

}

/* Compute modified equilibrium constant */
/* This includes fO2, which helps numerically since the total quantity is better scaled */
PetscScalar get_log10_modified_equilibrium_constant( const ReactionParameters * reaction_parameters_ptr, PetscScalar temp, const Constants *C, const Atmosphere *A )
{
    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscScalar        log10G, log10K; 

    log10K = get_log10_equilibrium_constant( reaction_parameters_ptr, temp, C );

    log10G = log10K - reaction_parameters->fO2_stoichiometry * A->log10fO2;

    return log10G;

}

/* Derivative of log10 modified equilibrium constant with respect to temperature */
PetscScalar get_dlog10GdT( const  ReactionParameters * reaction_parameters_ptr, PetscScalar temp, const Constants *C, const Atmosphere *A )
{
    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscScalar        dlog10KdT, dlog10GdT;

    dlog10KdT = get_dlog10KdT( reaction_parameters_ptr, A->tsurf, C );

    dlog10GdT = dlog10KdT - reaction_parameters->fO2_stoichiometry * A->dlog10fO2dT;

    return dlog10GdT;

}

/* Exponent of extra factor of psurf to ensure reaction quotient is non-dimensional */
static PetscScalar get_psurf_exponent( const ReactionParameters * reaction_parameters_ptr )
{
    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscInt    j;
    PetscScalar expon = 0;

    for( j=0; j<reaction_parameters->n_volatiles; ++j) {
        expon += reaction_parameters->stoichiometry[j];
    }

    return expon;

}

/* Compute reaction quotient (products, numerator) */

PetscScalar get_reaction_quotient_products( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A )
{

    return get_reaction_quotient( reaction_parameters_ptr, A, 1 );

}

PetscScalar get_reaction_quotient_reactants( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A )
{

    return get_reaction_quotient( reaction_parameters_ptr, A, -1 );

}

static PetscScalar get_reaction_quotient( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A, PetscInt SIGN )
{
    /* returns numerator for SIGN=1 (products) and denominator for SIGN=-1 (reactants) */
    /* note: excludes fO2 since this is badly scaled.  fO2 is dealt with at the same time
       as the equilibrium constant */

    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscInt           j;
    PetscScalar        Q = 1;
    PetscScalar        expon;

    for (j=0; j<reaction_parameters->n_volatiles; ++j) {

        const PetscInt v = reaction_parameters->volatiles[j];

        if( SIGN * reaction_parameters->stoichiometry[j] > 0.0 ){
            /* NOTE: introduced scaling by A->psurf to improve scaling for numerical solver (FD Jacobian) */
            /* TODO: if this works, swap out A->volatiles[v0].p/A->psurf for the volume mixing ratio? */
            Q *= PetscPowScalar( (A->volatiles[v].p/A->psurf), SIGN * reaction_parameters->stoichiometry[j] );
        }

    }

    /* since we are using the volume mixing ratio (i.e. scaling by A->psurf),
       we must account for the extra factors of A->psurf to ensure that
       Q=Qp/Qr is non-dimensional */
    /* TODO: keep here for the time being, but could maybe move this elsewhere */
    expon = get_psurf_exponent( reaction_parameters_ptr );
    if( SIGN * expon > 0.0 ){
        Q *= PetscPowScalar( A->psurf, SIGN * expon );
    }

    return Q;

}

PetscScalar get_reaction_quotient_products_time_derivative( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A, const AtmosphereParameters *Ap )
{

    return get_reaction_quotient_time_derivative( reaction_parameters_ptr, A, Ap, 1 );

}

PetscScalar get_reaction_quotient_reactants_time_derivative( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A, const AtmosphereParameters *Ap )
{

    return  get_reaction_quotient_time_derivative( reaction_parameters_ptr, A, Ap, -1 );

}

/* compute reaction quotient (products, numerator) derivative with respect to time t */
static PetscScalar get_reaction_quotient_time_derivative( const ReactionParameters * reaction_parameters_ptr, const Atmosphere *A, const AtmosphereParameters *Ap, PetscInt SIGN )
{
    /* returns dQp/dt for SIGN=1 (products) and dQr/dt for SIGN=-1 (reactants) */

    ReactionParameters reaction_parameters = *reaction_parameters_ptr;
    PetscInt           j,k;
    PetscScalar        dQdt=0, dvdt;
    PetscScalar        expon, prefactor;
    PetscBool          INCLUDE_PSURF = PETSC_FALSE;

    /* this is a bit ugly, because I decide whether to include the scaling of surface pressure in the
       numerator if the stoichiometry is positive (otherwise it is included in the denominator */
    expon = get_psurf_exponent( reaction_parameters_ptr );
    if( SIGN * expon > 0.0 ){
        INCLUDE_PSURF = PETSC_TRUE;
    }

    /* dpsurf/dt previously updated in get_dxdt */

    /* what follows is the chain rule */

    /* first, account for volatiles in this reaction that are explicitly tracked
       by the code */
    for (j=0; j<reaction_parameters->n_volatiles; ++j) {

        if( SIGN * reaction_parameters->stoichiometry[j] > 0.0 ){

            /* compute derivative of this volatile */
            const PetscInt v = reaction_parameters->volatiles[j];
            prefactor = 1.0;

            /* find all other volatiles with positive stoichiometry that define the prefactor */
            for (k=0; k<reaction_parameters->n_volatiles; ++k) {
                if ( SIGN * reaction_parameters->stoichiometry[k] > 0.0 && k!=j ){
                    const PetscInt v2 = reaction_parameters->volatiles[k];
                    prefactor *= PetscPowScalar( (A->volatiles[v2].p/A->psurf), SIGN * reaction_parameters->stoichiometry[k] );
                }
            }
            /* should the prefactor include Psurf? */
            if( INCLUDE_PSURF ){
                prefactor *= PetscPowScalar( A->psurf, SIGN * expon );
            }
            /* compute the derivative of this volatile */
            /* TODO: this is also a generic formula for the time derivative of the
               volume mixing ratio and appears elsewhere in the code */
            dvdt = (1.0/A->psurf) * A->volatiles[v].dxdt * A->volatiles[v].dpdx;
            dvdt -= (A->volatiles[v].p / PetscPowScalar( A->psurf, 2.0 )) * A->dpsurfdt;

            /* so the contribution of this derivative is as follows */
            dvdt *= SIGN * reaction_parameters->stoichiometry[j] * PetscPowScalar( (A->volatiles[v].p/A->psurf), SIGN * reaction_parameters->stoichiometry[j] - 1.0 );
            dQdt += prefactor * dvdt;
        }
    }

    /* now account for dpsurf/dt */
    if( INCLUDE_PSURF ) {

        /* find all other volatiles with positive stoichiometry that define the prefactor */
        prefactor = 1.0;
        for (k=0; k<reaction_parameters->n_volatiles; ++k) {
            if ( SIGN * reaction_parameters->stoichiometry[k] > 0.0 ){
                const PetscInt v2 = reaction_parameters->volatiles[k];
                prefactor *= PetscPowScalar( (A->volatiles[v2].p/A->psurf), SIGN * reaction_parameters->stoichiometry[k] );
            }
        }

        /* now compute the time derivative of psurf */
        dvdt = A->dpsurfdt;

        /* so the contribution is */
        dvdt *= SIGN * expon * PetscPowScalar( A->psurf, SIGN * expon - 1.0 );
        dQdt += prefactor * dvdt;
    }

    return dQdt;

}
