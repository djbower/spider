#include "rtpress.h"
#include "monitor.h"

static PetscScalar get_rtpress_pressure( PetscScalar, PetscScalar, Eos * );
static PetscScalar get_rtpress_entropy( PetscScalar, PetscScalar, Eos * );
static PetscErrorCode objective_function_rtpress_volume_temperature( SNES, Vec, Vec, void * );
static PetscErrorCode solve_for_rtpress_volume_temperature( Ctx * );

PetscErrorCode set_rtpress_parameters( Eos *rtp )
{
    /* eos parameters for rtpress, taken from jupyter notebook
       TODO: check units */

    PetscFunctionBeginUser;

    rtp->V0 = 14.74;
    rtp->T0 = 3000.0;
    rtp->S0 = 0.0;
    rtp->K0 = 9.77;
    rtp->KP0 = 7.42;
    rtp->E0 = -1097.490985248;
    rtp->gamma0 = 0.282;
    rtp->gammaP0 = -1.35;
    rtp->m = 0.6;
    rtp->b0 = 2.5947702519786517;
    rtp->b1 = -0.11604518121550321;
    rtp->b2 = 4.8738976110511345;
    rtp->b3 = 29.939656753599827;
    rtp->b4 = 35.97400617680599;

    PetscFunctionReturn(0);

}

static PetscScalar get_rtpress_pressure( PetscScalar V, PetscScalar T, Eos *rtp )
{
    /* pressure = function( volume, temperature ) */

    PetscScalar P;
    PetscScalar const K0 = rtp->K0;
    PetscScalar const V0 = rtp->V0;
    PetscScalar const KP0 = rtp->KP0;
    PetscScalar const gamma0 = rtp->gamma0;
    PetscScalar const gammaP0 = rtp->gammaP0;
    PetscScalar const T0 = rtp->T0;
    PetscScalar const m = rtp->m;
    PetscScalar const b0 = rtp->b0;
    PetscScalar const b1 = rtp->b1;
    PetscScalar const b2 = rtp->b2;
    PetscScalar const b3 = rtp->b3;
    PetscScalar const b4 = rtp->b4;

    /* FIXME: if I try using PetscCbrtReal then I get the following warning at compilation:
           warning: implicit declaration of function 'cbrt' is invalid in C99 [-Wimplicit-function-declaration] */
    /* TODO: should I use Real or PetscScalar functions?  i.e. PetscPowReal or PetscPowScalar? */

    P =  -9*K0*V0*(-1.0/3.0*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0) - 1)*PetscExpReal((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/V - 1.0/3.0*cbrt(V/V0)*((3.0/2.0)*KP0 - 3.0/2.0)*PetscExpReal((1 - cbrt(V/V0))*((3.0/2.0)*KP0 - 3.0/2.0))/V)/PetscPowScalar((3.0/2.0)*KP0 - 3.0/2.0, 2) + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/m + 0.013806489886250916*T0*PetscPowScalar(T0, m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/m - T0*(0.027612979772501833*T0*PetscPowScalar(T0, m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m - 2*T0) - 0.027612979772501833*T0*m*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/((2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b1/V0 + 2*b2*(V/V0 - 1)/V0 + 3*b3*PetscPowScalar(V/V0 - 1, 2)/V0 + 4*b4*PetscPowScalar(V/V0 - 1, 3)/V0)*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(-2*T0*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) + 2*T0*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), 2) - 0.041419469658752747*m*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*(-gamma0*PetscPowScalar(V0/V, 2.0/3.0)/V - 1.0/6.0*PetscPowScalar(V0/V, 2.0/3.0)*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0)/V)/((2*m - 2)*(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)));

    return P;

}

static PetscScalar get_rtpress_entropy( PetscScalar V, PetscScalar T, Eos *rtp )
{
    /* entropy = function( volume, temperature ) */

    PetscScalar S;
    PetscScalar const gamma0 = rtp->gamma0;
    PetscScalar const gammaP0 = rtp->gammaP0;
    PetscScalar const V0 = rtp->V0;
    PetscScalar const S0 = rtp->S0;
    PetscScalar const T0 = rtp->T0;
    PetscScalar const m = rtp->m;
    PetscScalar const b0 = rtp->b0;
    PetscScalar const b1 = rtp->b1;
    PetscScalar const b2 = rtp->b2;
    PetscScalar const b3 = rtp->b3;
    PetscScalar const b4 = rtp->b4;

    S = S0 + T*(0.027612979772501833*PetscPowScalar(T, m)*T0*(2 - 2*m)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/PetscPowScalar(2*T*m - 2*T, 2) + 0.027612979772501833*PetscPowScalar(T, m)*T0*m*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(T*(2*T*m - 2*T)) + 0.041419469658752747*m/(T*(2*m - 2)) - 0.041419469658752747/(T*(2*m - 2))) + 0.027612979772501833*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(2*T*m - 2*T) - 0.027612979772501833*T0*PetscPowScalar(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1), m)*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/(2*T0*m*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1) - 2*T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1)) + 0.041419469658752747*m*log(T)/(2*m - 2) - 0.041419469658752747*m*log(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))/(2*m - 2) - 0.020709734829376374 - 0.041419469658752747*log(T)/(2*m - 2) + 0.041419469658752747*log(T0*PetscSqrtReal(6*gamma0*((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0) + (1.0/2.0)*PetscPowScalar((1.0/2.0)*PetscPowScalar(V0/V, 2.0/3.0) - 1.0/2.0, 2)*(36*PetscPowScalar(gamma0, 2) - 12*gamma0 - 18*gammaP0) + 1))/(2*m - 2) - 0.013806489886250916*PetscPowScalar(T, m)*T0*(b0 + b1*(V/V0 - 1) + b2*PetscPowScalar(V/V0 - 1, 2) + b3*PetscPowScalar(V/V0 - 1, 3) + b4*PetscPowScalar(V/V0 - 1, 4))*PetscPowScalar(1.0/T0, m)/T;

    return S;

}

static PetscErrorCode solve_for_rtpress_volume_temperature( Ctx *E )
{
    PetscErrorCode ierr;
    SNES           snes;
    Vec            x,r;
    PetscScalar    *xx;
    PetscInt       i;
    EosEval        *eos_eval = &E->eos_eval;

    PetscFunctionBeginUser;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"solve_for_rtpress_volume_temperature()\n");CHKERRQ(ierr);

    ierr = SNESCreate( PETSC_COMM_WORLD, &snes );CHKERRQ(ierr);

    /* Use this to address this specific SNES (nonlinear solver) from the command
       line or options file, e.g. -atmosic_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"rtpress_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    ierr = VecSetSizes( x, PETSC_DECIDE, 2 );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    ierr = SNESSetFunction(snes,r,objective_function_rtpress_volume_temperature,E);CHKERRQ(ierr);

    /* initialise vector x with initial guess */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    for (i=0; i<2; ++i) {
        /* FIXME: improve initial guesses */
        xx[i] = 100.0;
    }
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* Inform the nonlinear solver to generate a finite-difference approximation
       to the Jacobian */
     ierr = PetscOptionsSetValue(NULL,"-rtpress_snes_mf",NULL);CHKERRQ(ierr);

    /* For solver analysis/debugging/tuning, activate a custom monitor with a flag */
    {   
      PetscBool flg = PETSC_FALSE;

      ierr = PetscOptionsGetBool(NULL,NULL,"-rtpress_snes_verbose_monitor",&flg,NULL);CHKERRQ(ierr);
      if (flg) {
        ierr = SNESMonitorSet(snes,SNESMonitorVerbose,NULL,NULL);CHKERRQ(ierr);
      }
    }

    /* Solve */
    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr); /* Picks up any additional options (note prefix) */
    ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);
    {   
      SNESConvergedReason reason;
      ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);
      if (reason < 0) SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
          "Nonlinear solver didn't converge: %s\n",SNESConvergedReasons[reason]);
    }   

    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    for (i=0; i<2; ++i) {
        if( xx[i] < 0.0 ){
            /* Sanity check on solution (since it's non-unique) */
            SETERRQ2(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
                "Unphysical rtpress property: slot %d, x: %g",i,xx[i]);
        }
    }

    /* store solution in struct */
    eos_eval->V = xx[0];
    eos_eval->T = xx[1];

    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

static PetscErrorCode objective_function_rtpress_volume_temperature( SNES snes, Vec x, Vec f, void *ptr)
{
    PetscErrorCode             ierr;
    const PetscScalar          *xx;
    PetscScalar                *ff;
    PetscScalar                V, T, P, S;
    Ctx                        *E = (Ctx*) ptr;
    Eos                        *rtp = &E->parameters.rtpress;
    EosEval                    *eos_eval = &E->eos_eval;
    PetscScalar Ptarget = eos_eval->P;
    PetscScalar Starget = eos_eval->S;

    PetscFunctionBeginUser;

    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

    V = xx[0];
    T = xx[1];

    P = get_rtpress_pressure( V, T, rtp );
    S = get_rtpress_entropy( V, T, rtp );

    /* compute residual */
    ff[0] = P - Ptarget;
    ff[1] = S - Starget;

    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

PetscScalar get_rtpress_temperature( PetscScalar P, PetscScalar S, Ctx *E )
{
    Constants const            *C = &E->parameters.constants;
    EosEval                    *eos_eval = &E->eos_eval;

    /* store values to lookup in struct */
    eos_eval->P = P;
    eos_eval->S = S;

    /* FIXME: I assume that the rtpress expressions are dimensional */
    eos_eval->P *= C->PRESSURE;
    eos_eval->S *= C->ENTROPY;

    /* below updates V and T in eos_eval */
    solve_for_rtpress_volume_temperature( E );

    return eos_eval->T;

}
