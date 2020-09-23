#include "mesh.h"
#include "monitor.h"

static PetscErrorCode regular_mesh( Ctx * );
//static PetscErrorCode geometric_mesh( Ctx * );
static PetscErrorCode spherical_area( DM, Vec, Vec );
static PetscErrorCode spherical_volume( Ctx *, Vec, Vec );
// FIXME: TODO: REMOVE
//static PetscScalar get_layer( DM, Vec, Vec, Parameters const );
static PetscErrorCode aw_density( DM, Vec, Vec, Parameters const, PetscScalar * );
static PetscErrorCode aw_pressure( DM, Vec, Vec, Parameters const );
static PetscErrorCode aw_pressure_gradient( DM, Vec, Vec, Parameters const );
static PetscErrorCode aw_mass( Mesh * );
// below kept for testing purposes
static PetscErrorCode set_xi_from_radius( DM, Vec, Vec, Vec, Parameters const, PetscScalar );
static PetscErrorCode aw_mantle_density( const Parameters, PetscScalar * );
static PetscErrorCode aw_radius_from_xi( Ctx * );

PetscErrorCode set_mesh( Ctx *E)
{
    PetscErrorCode ierr;
    Mesh           *M = &E->mesh;
    DM             da_b=E->da_b, da_s=E->da_s;
    Parameters     P = E->parameters;

    PetscFunctionBeginUser;

    /* for regular mesh (mass coordinates) */
    ierr = regular_mesh( E );CHKERRQ(ierr);

    /* FIXME: broken for mass coordinates */
    //geometric_mesh( E );

    /* assume we know planetary size (core radius, and total radius)
       this allows us to compute the rho0 that will recover the exact
       desired bounds of the mesh (xi from 0 to P->radius, and r from
       r_cmb to P->radius).  These functions are all tied to the assumed
       relationship between rho and r, given by the Adams-Williamson
       EOS */
    if(1){

        /* Adams-Williamson EOS */

        /* determine reference density from AW EOS */
        ierr = aw_mantle_density( P, &M->mantle_density );CHKERRQ(ierr);

        /* with rho0 and the form of rho known (rho(r) from AW EOS),
           we can solve an inverse problem to determine the physical
           mesh coordinates for both the basic and staggered nodes */
        ierr = aw_radius_from_xi( E );CHKERRQ(ierr);

        /* for testing, do the inverse calculation */
        ierr = set_xi_from_radius( da_b, M->radius_b, M->xi_b, M->dxidr_b, P, M->mantle_density );CHKERRQ(ierr);
        ierr = set_xi_from_radius( da_s, M->radius_s, M->xi_s, NULL, P, M->mantle_density );CHKERRQ(ierr);

        /* with radius known, now can update other quantities such
           as pressure, using AW EOS */

        /* pressure at basic nodes */
        ierr = aw_pressure( da_b, M->radius_b, M->pressure_b, P);CHKERRQ(ierr);

        /* dP/dr at basic nodes */
        ierr = aw_pressure_gradient( da_b, M->radius_b, M->dPdr_b, P);CHKERRQ(ierr);

        /* pressure at staggered nodes */
        ierr = aw_pressure( da_s, M->radius_s, M->pressure_s, P);CHKERRQ(ierr);

        /* dP/dr at staggered nodes */
        ierr = aw_pressure_gradient( da_s, M->radius_s, M->dPdr_s, P );CHKERRQ(ierr);

    }

    /* alternatively, we could solve static structure equations to
       recover relationship between xi, radius, pressure, rho */
    else {

        /* solve static structure equations to determine planetary
           size at this time step (would need to move into time
           loop) */
        ;
    }


    /* surface area at basic nodes, without 4*pi term */
    ierr = spherical_area( da_b, M->radius_b, M->area_b);CHKERRQ(ierr);

    /* surface area at staggered nodes, without 4*pi term */
    ierr = spherical_area( da_s, M->radius_s, M->area_s );CHKERRQ(ierr);

    /* volume of spherical cells, without 4*pi term */
    ierr = spherical_volume( E, M->radius_b, M->volume_s);CHKERRQ(ierr);

    /* FIXME: this is to be REMOVED, once migrated into the new EOS object approach */
    /* layer id.  0 everywhere for single layer (as determined by
       P->mixing_length), and 0 for upper and 1 for lower layer
       when P->mixing_length==3 */
    //get_layer( da_b, M->radius_b, M->layer_b, P );

    /* density at staggered nodes */
    ierr = aw_density( da_s, M->radius_s, M->rho_s, P, &M->mantle_density );CHKERRQ(ierr);

    /* mass at staggered nodes */
    ierr = aw_mass( M );CHKERRQ(ierr);

    /* mantle mass also needed for atmosphere calculations */
    P->atmosphere_parameters->mantle_mass_ptr = &M->mantle_mass;

    PetscFunctionReturn(0);
}

static PetscErrorCode regular_mesh( Ctx *E )
{

    PetscErrorCode ierr;
    PetscScalar    *arr;
    PetscInt       i,ilo_b,ihi_b,ilo_s,ihi_s,w_b,w_s,numpts_b,numpts_s;
    Mesh           *M = &E->mesh;
    Parameters     P = E->parameters;
    DM             da_b=E->da_b, da_s=E->da_s;
    PetscScalar    dx_b;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* basic node spacing (negative) */
    /* mass coordinate enforced to go from 0 to P->radius */
    /* note that with appropriate choice of scaling RADIUS, this can be
       made to go between 0 and unity (i.e., P->radius is a scaled quantity) */
    dx_b = -P->radius / (numpts_b-1);

    /* radius at basic nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ierr = DMDAVecGetArray(da_b,M->xi_b,&arr);CHKERRQ(ierr);
    for (i=ilo_b; i<ihi_b; ++i){
        arr[i] = -(numpts_b-1-i)*dx_b;
    }
    ierr = DMDAVecRestoreArray(da_b,M->xi_b,&arr);CHKERRQ(ierr);

    /* radius at staggered nodes */
    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAVecGetArray(da_s,M->xi_s,&arr);CHKERRQ(ierr);
    for (i=ilo_s;i<ihi_s;++i){
        arr[i] = -0.5*dx_b - (numpts_s-1-i)*dx_b;
    }
    ierr = DMDAVecRestoreArray(da_s,M->xi_s,&arr);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#if 0
/* FIXME: below broken for mass coordinates */
static PetscErrorCode geometric_mesh( Ctx *E )
{

    PetscErrorCode ierr;
    PetscScalar    *arr_b, *arr_s;
    PetscInt       i,ilo_b,ihi_b,ilo_s,ihi_s,w_b,w_s,numpts_b,numpts_s;
    Mesh           *M;
    DM             da_b=E->da_b, da_s=E->da_s;
    // TODO: should not be hard-coded
    PetscScalar    dr_min = 1.5696123057604772e-10; // 1 mm
    // 45928 basic nodes
    //PetscScalar dr_max = 9.810076911002982e-06; // 62.5 m
    // 22996 basic nodes
    //PetscScalar dr_max = 1.9620153822005963e-05; // 125 m
    // 11532 basic nodes
    //PetscScalar dr_max = 3.924030764401193e-05; // 250 m
    // 5802 basic nodes
    //PetscScalar dr_max = 7.848061528802385e-05; // 500 m
    // 2939 basic nodes
    //PetscScalar    dr_max = 0.0001569612305760477; // 1 km
    // 656 basic nodes
    //PetscScalar    dr_max = 0.0007848061528802385; // 5 km
    // 372 basic nodes
    //PetscScalar    dr_max = 0.001569612305760477; // 10 km
    // 278 basic nodes
    PetscScalar    dr_max = 0.0023544184586407153; // 15 km
    PetscScalar    geom_fac = 1.2;
    PetscScalar    dr, dr_prev;

    PetscFunctionBeginUser;

    M = &E->mesh;

    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* TODO: dynamically determine number of mesh points, and
       automatically ensure that arr[0] = 1.0 and
       arr[numpts_b-1] = RADIN */

    /* FIXME: update to work in parallel */

    /* radius at basic nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ierr = DMDAVecGetArray(da_b,M->radius_b,&arr_b);CHKERRQ(ierr);
    arr_b[numpts_b-1] = RADIN;
    dr_prev = dr_min / geom_fac;
    for (i=ilo_b; i<ihi_b-1; ++i){ /* Note upper bound */
        dr = dr_prev * geom_fac;
        if( dr > dr_max){
            dr = dr_max;
        }
        dr_prev = dr;
        arr_b[numpts_b-i-2] = arr_b[numpts_b-i-1] + dr;
    }
    // TODO: below is hacky.
    arr_b[0] = 1.0;

    ierr = DMDAVecRestoreArray(da_b,M->radius_b,&arr_b);CHKERRQ(ierr);

    /* radius at staggered nodes */
    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAVecGetArray(da_s,M->radius_s,&arr_s);CHKERRQ(ierr);
    for (i=ilo_s;i<ihi_s;++i){
        arr_s[i] = 0.5 * (arr_b[i]+arr_b[i+1]);
    }
    ierr = DMDAVecRestoreArray(da_s,M->radius_s,&arr_s);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}
#endif

static PetscErrorCode spherical_area(DM da, Vec radius, Vec area )
{

    PetscErrorCode    ierr;
    PetscScalar       *arr_area;
    const PetscScalar *arr_radius;
    PetscInt          i,ilo,ihi,w;

    PetscFunctionBeginUser;
    ierr = DMDAGetCorners(da,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w; 
    ierr = DMDAVecGetArrayRead(da,radius,&arr_radius);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,area,&arr_area);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        /* excludes 4*pi prefactor */
        arr_area[i] = PetscPowScalar( arr_radius[i], 2.0 );
    }    
    ierr = DMDAVecRestoreArrayRead(da,radius,&arr_radius);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,area,&arr_area);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode spherical_volume(Ctx * E, Vec radius, Vec volume )
{

    PetscErrorCode    ierr;
    PetscScalar       *arr_volume;
    const PetscScalar *arr_radius;
    PetscInt          i,ilo_s,ihi_s,w_s,ilo,ihi;
    DM                da_b=E->da_b,da_s=E->da_s;
    Vec               radius_local=E->work_local_b;

    PetscFunctionBeginUser;
    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ilo = ilo_s;
    ihi = ihi_s;
    ierr = DMGlobalToLocalBegin(da_b,radius,INSERT_VALUES,radius_local);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da_b,radius,INSERT_VALUES,radius_local);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da_b,radius_local,&arr_radius);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da_s,volume,&arr_volume);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        arr_volume[i] = PetscPowScalar(arr_radius[i],3.0) - PetscPowScalar(arr_radius[i+1],3.0);
        /* note excludes 4*pi prefactor */
        arr_volume[i] *= 1.0 / 3.0;
    }
    // TODO: here and elsewhere, it's very dangerous to use the same indice to refer to two DAs without checking that the local ranges are valid. 
    ierr = DMDAVecRestoreArrayRead(da_b,radius_local,&arr_radius);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,volume,&arr_volume);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/* FIXME: TODO: REMOVE EVENTUALLY */
#if 0
/* TODO: need to have some notion of a layer ID, but should it be within mesh.c? */
static PetscScalar get_layer( DM da, Vec radius, Vec layer, const Parameters P )
{
    PetscErrorCode ierr;
    PetscScalar *arr_layer;
    const PetscScalar *arr_r;
    PetscInt i,ilo,ihi,w;

    PetscFunctionBeginUser;
    ierr = DMDAGetCorners(da,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;
    ierr = DMDAVecGetArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,layer,&arr_layer);CHKERRQ(ierr);

    for(i=ilo; i<ihi; ++i){
        if( P->mixing_length==3){
            if( arr_r[i] > P->mixing_length_layer_radius ){
                arr_layer[i] = 0;
            }
            else{
                arr_layer[i] = 1;
            }
        }
        else{
            arr_layer[i] = 0;
        }
    }
    ierr = DMDAVecRestoreArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,layer,&arr_layer);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
#endif

// kept for testing
# if 1
static PetscErrorCode set_xi_from_radius( DM da, Vec radius, Vec xi, Vec dxidr, const Parameters P, PetscScalar mantle_density )
{

    /* set mass coordinate from radius.  For the simplest case of a prescribed hydrostatic
       equation of state (Adams-Williamson), density is a function of r and the mapping can
       be computed directly from radius --> xi (mass coordinate) */

    PetscErrorCode ierr;
    PetscScalar    dep,*arr_xi,*arr_dxidr;
    const PetscScalar *arr_r;
    PetscInt       i,ilo,ihi,w;

    PetscFunctionBeginUser;

    ierr = DMDAGetCorners(da,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;
    ierr = DMDAVecGetArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,xi,&arr_xi);CHKERRQ(ierr);
    if( dxidr != NULL ){
        ierr = DMDAVecGetArray(da,dxidr,&arr_dxidr);CHKERRQ(ierr);
    }

    for(i=ilo; i<ihi; ++i){
        dep = P->radius - arr_r[i];
        /* below is AW density, and could instead use existing function aw_density() */
        /* evaluate integral at r */
        /* this is mass contained with shell of radius r */
        arr_xi[i] = -2/PetscPowScalar(P->beta,3) - PetscPowScalar(arr_r[i],2)/P->beta - 2*arr_r[i]/PetscPowScalar(P->beta,2);
        arr_xi[i] *= P->rhos * PetscExpScalar( P->beta * dep );
        /* minus integral at radius = 0 */
        //arr_xi[i] -= -2/PetscPowScalar(P->beta,3) * P->rhos * PetscExpScalar( P->beta * P->radius );
        /* seems better to do this instead, since the mantle does not extend to r=0 */
        /* for mass coordinates, Abe 1995 says that the choice of reference density is perfectly arbitrary, so this should be OK */
        /* minus integral at radius = core-mantle boundary */
        /* this will set the CMB radius to zero mass coordinate */
        arr_xi[i] -= (-2/PetscPowScalar(P->beta,3) - PetscPowScalar(P->radius*P->coresize,2)/P->beta - 2*P->radius*P->coresize/PetscPowScalar(P->beta,2)) * P->rhos * PetscExpScalar( P->beta * P->radius * (1.0-P->coresize) );
        /* include other prefactors, according to the formulation from radius to mass coordinate */
        arr_xi[i] *= 3 / mantle_density;
        arr_xi[i] = PetscPowScalar( arr_xi[i], 1.0/3.0 );

        /* set dxi/dr for derivative mapping */
        /* last basic node returns dxidr infty, since xi=0 */
        if( dxidr != NULL ){
            arr_dxidr[i] = P->rhos * PetscExpScalar( P->beta * dep ) / mantle_density;
            arr_dxidr[i] *= PetscPowScalar( arr_r[i] / arr_xi[i], 2.0 );
        }

    }

    ierr = DMDAVecRestoreArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,xi,&arr_xi);CHKERRQ(ierr);
    if( dxidr != NULL ){
        ierr = DMDAVecRestoreArray(da,dxidr,&arr_dxidr);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
#endif

static PetscErrorCode aw_pressure( DM da, Vec radius, Vec pressure, const Parameters P )
{
    PetscErrorCode    ierr;
    PetscScalar       dep,*arr_p;
    const PetscScalar *arr_r;
    PetscInt          i,ilo,ihi,w;

    PetscFunctionBeginUser;
    ierr = DMDAGetCorners(da,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;
    ierr = DMDAVecGetArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,pressure,&arr_p);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        dep = P->radius - arr_r[i];
        arr_p[i] = P->rhos * -P->gravity / P->beta;
        arr_p[i] *= PetscExpScalar( P->beta * dep ) - 1.0;
    }
    ierr = DMDAVecRestoreArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,pressure,&arr_p);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode aw_density( DM da, Vec radius, Vec density, const Parameters P, PetscScalar *mantle_density )
{
    PetscErrorCode    ierr;
    PetscScalar       dep, *arr_density;
    const PetscScalar *arr_r;
    PetscInt          i,ilo,ihi,w;

    PetscFunctionBeginUser;
    ierr = DMDAGetCorners(da,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;
    ierr = DMDAVecGetArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,density,&arr_density);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        dep = P->radius - arr_r[i];
        arr_density[i] = P->rhos * PetscExpScalar( P->beta * dep );
    }
    ierr = DMDAVecRestoreArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,density,&arr_density);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode aw_mantle_density( const Parameters P, PetscScalar *mantle_density )
{
    PetscScalar dep;

    /* average mantle density */
    /* integral at planetary radius P->radius */
    dep = 0.0;
    *mantle_density = -2/PetscPowScalar(P->beta,3) - PetscPowScalar(P->radius,2)/P->beta - 2*P->radius/PetscPowScalar(P->beta,2);
    *mantle_density *= P->rhos; // next is unity since dep=0: * PetscExpScalar( P->beta * dep );
    /* minus integral at r = 0 */
    //*mantle_density -= -2/PetscPowScalar(P->beta,3) * PetscExpScalar( P->beta * P->radius );
    /* seems better to do this instead, since the mantle does not extend to r=0 */
    /* for mass coordinates, Abe 1995 says that the choice of reference density is perfectly arbitrary, so this should be OK */
    /* minus integral at core-mantle boundary P->radius * P->coresize */
    dep = P->radius * (1.0-P->coresize);
    *mantle_density -= (-2/PetscPowScalar(P->beta,3) - PetscPowScalar(P->radius*P->coresize,2)/P->beta - 2*P->radius*P->coresize/PetscPowScalar(P->beta,2)) * P->rhos * PetscExpScalar( P->beta * dep );
    /* above is integrated mass from core-mantle boundary to surface radius.  Now divide by mantle volume to get density */
    //*mantle_density /= 1.0/3.0 * ( PetscPowScalar(P->radius,3.0) - PetscPowScalar(P->coresize*P->radius,3.0) );
    /* new below follows definition of mass coordinates to tie the average density to ensure that
       the outermost mass coordinate xi = P->radius (innermost mass coordinate xi = 0 at r = rcmb) */
    *mantle_density *= 3.0 / PetscPowScalar( P->radius, 3.0 );

    PetscFunctionReturn(0);
}

static PetscErrorCode aw_mass( Mesh *M )
{
    PetscErrorCode ierr;
 
    PetscFunctionBeginUser;

    ierr = VecCopy( M->rho_s, M->mass_s ); CHKERRQ(ierr);

    ierr = VecPointwiseMult( M->mass_s, M->mass_s, M->volume_s );
    /* excludes 4*pi prefactor */
    ierr = VecSum( M->mass_s, &M->mantle_mass );

    /* also determine volume, also excluding 4*pi prefactor */
    ierr = VecSum( M->volume_s, &M->mantle_volume );

    /* FIXME: M->mantle_density is computed analytically to give exact value */
    /* this here would over-ride */
    //M->mantle_density = M->mantle_mass / M->mantle_volume;

    PetscFunctionReturn(0);

}

static PetscErrorCode aw_pressure_gradient( DM da, Vec radius, Vec grad, Parameters const P )
{
    PetscErrorCode    ierr;
    PetscScalar       dep,*arr_g;
    const PetscScalar *arr_r;
    PetscInt          i,ilo,ihi,w;

    PetscFunctionBeginUser;
    ierr = DMDAGetCorners(da,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;
    ierr = DMDAVecGetArray(da,grad,&arr_g);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        dep = P->radius - arr_r[i];
        arr_g[i] = P->rhos * P->gravity * PetscExpScalar( P->beta * dep );
    }
    ierr = DMDAVecRestoreArray(da,grad,&arr_g);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode objective_function_radius( SNES snes, Vec x, Vec f, void *ptr )
{
    PetscErrorCode    ierr;
    const PetscScalar *xx, *xi_b, *xi_s;
    PetscScalar       *ff, dep, target_xi;
    Ctx               *E = (Ctx*) ptr;
    Parameters  const P = E->parameters;
    Mesh        const *M = &E->mesh;
    PetscInt          i,ilo_b,ihi_b,numpts_b,w_b,ilo_s,ihi_s,numpts_s,w_s;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAGetCorners(E->da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;

    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAGetCorners(E->da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;

    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr); /* initial guess of r */
    ierr = VecGetArrayRead(M->xi_b,&xi_b);CHKERRQ(ierr); /* target mass coordinate */
    ierr = VecGetArrayRead(M->xi_s,&xi_s);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr); /* residual function */

    /* FIXME: broken for parallel */
    for(i=0; i< numpts_b+numpts_s; ++i){
        dep = P->radius - xx[i];
        /* evaluate integral at r: mass contained with shell of radius r */
        ff[i] = -2/PetscPowScalar(P->beta,3) - PetscPowScalar(xx[i],2)/P->beta - 2*xx[i]/PetscPowScalar(P->beta,2);
        ff[i] *= P->rhos * PetscExpScalar( P->beta * dep );
        /* minus integral at radius = core-mantle boundary */
        ff[i] -= (-2/PetscPowScalar(P->beta,3) - PetscPowScalar(P->radius*P->coresize,2)/P->beta - 2*P->radius*P->coresize/PetscPowScalar(P->beta,2)) * P->rhos * PetscExpScalar( P->beta * P->radius * (1.0-P->coresize) );
        /* minus predicted mass from mass coordinate */
        /* first part relates to basic nodes */
        if (i<numpts_b){
            target_xi = xi_b[i];
        }
        /* second part relates to staggered nodes */
        else{
            target_xi = xi_s[i-numpts_b];
        }

#if 1
        // currently best here
        ff[i] -= (M->mantle_density / 3) * PetscPowScalar(target_xi,3.0);
#endif

#if 0
        // testing here
        ff[i] *= 3.0 / M->mantle_density;
        ff[i] = PetscPowScalar( ff[i], 1.0/3.0 );
        ff[i] -= target_xi;
#endif


        //ff[i] /= (M->mantle_density / 3) * PetscPowScalar(P->radius,3.0);

        /* include other prefactors, according to the formulation from radius to mass coordinate */
        //ff[i] *= 3 / M->mantle_density;

        // Dan trying to improve solver behaviour by skipping this cube root
        //ff[i] = PetscPowScalar( ff[i], 1.0/3.0 );
        //ff[i] -= PetscPowScalar(xi[i],1.0);

        /* residual is simple difference of this to desired */
        //ff[i] -= PetscPowScalar(xi[i],3.0);
        //ff[i] = PetscPowScalar(ff[i],1.0/3.0);
    }

    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(M->xi_b,&xi_b);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(M->xi_s,&xi_s);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

static PetscErrorCode aw_radius_from_xi( Ctx *E )
{
    PetscErrorCode ierr;
    SNES           snes;
    Vec            x,r;
    PetscScalar    *xx, *radius, dx;
    PetscInt       i,numpts_b,numpts_s;
    Mesh           *M = &E->mesh;
    Parameters const P = E->parameters;

    PetscFunctionBeginUser;

    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"aw_radius_from_xi()\n");CHKERRQ(ierr);

    ierr = SNESCreate( PETSC_COMM_WORLD, &snes );CHKERRQ(ierr);

    /* Use this to address this specific SNES (nonlinear solver) from the command
       line or options file, e.g. -atmosic_snes_view */
    ierr = SNESSetOptionsPrefix(snes,"xi_");CHKERRQ(ierr);

    ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);
    /* here I am combining the basic and staggered nodes, but a DMComposite construction
       is preferred */
    ierr = VecSetSizes( x, PETSC_DECIDE, numpts_b+numpts_s );CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    ierr = SNESSetFunction(snes,r,objective_function_radius,E);CHKERRQ(ierr);

    /* initialise vector x with initial guess */
    ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
    dx = P->radius * (1.0-P->coresize) / (numpts_b-1);
    /* basic nodes */
    for (i=0; i<numpts_b; ++i) {
        /* best initial guess is evenly space from surface to cmb */
        xx[i] = P->radius - i * dx;
    }
    /* staggered nodes */
    for (i=numpts_b; i<numpts_b+numpts_s; ++i) {
        /* best initial guess is evenly space from surface to cmb */
        xx[i] = P->radius - 0.5 * dx - (i-numpts_b) * dx;
    }
    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* Inform the nonlinear solver to generate a finite-difference approximation
       to the Jacobian */
    ierr = PetscOptionsSetValue(NULL,"-xi_snes_mf",NULL);CHKERRQ(ierr);

    /* Turn off convergence based on step size */
    ierr = PetscOptionsSetValue(NULL,"-xi_snes_stol","0");CHKERRQ(ierr);

    /* Turn off convergenced based on trust region tolerance */
    ierr = PetscOptionsSetValue(NULL,"-xi_snes_trtol","0");CHKERRQ(ierr);

    /* For solver analysis/debugging/tuning, activate a custom monitor with a flag */
    {
      PetscBool flg = PETSC_FALSE;

      ierr = PetscOptionsGetBool(NULL,NULL,"-xi_snes_verbose_monitor",&flg,NULL);CHKERRQ(ierr);
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

    /* extract solution for basic radius from solution vec */
    ierr = DMDAVecGetArray(E->da_b,M->radius_b,&radius);CHKERRQ(ierr);
    for (i=0; i<numpts_b; ++i) {
        if( xx[i] < 0.0 ){
            /* Sanity check on solution */
            SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
                "Unphysical radius coordinate, x: %g",xx[i]);
        }   
        else{
            radius[i] = xx[i];
        }   
    }   
    ierr = DMDAVecRestoreArray(E->da_b,M->radius_b,radius);CHKERRQ(ierr);

    /* extract solution for staggered radius from solution vec */
    ierr = DMDAVecGetArray(E->da_s,M->radius_s,&radius);CHKERRQ(ierr);
    for (i=numpts_b; i<numpts_b+numpts_s; ++i) {
        if( xx[i] < 0.0 ){
            /* Sanity check on solution */
            SETERRQ1(PetscObjectComm((PetscObject)snes),PETSC_ERR_CONV_FAILED,
                "Unphysical radius coordinate, x: %g",xx[i]);
        }   
        else{
            radius[i-numpts_b] = xx[i];
        }   
    }   
    ierr = DMDAVecRestoreArray(E->da_s,M->radius_s,radius);CHKERRQ(ierr);

    ierr = VecRestoreArray(x,&xx);CHKERRQ(ierr);

    /* now compute dxi/dr once all radius and xi are known */
    ierr = VecCopy( M->radius_b, M->dxidr_b ); CHKERRQ(ierr);
    ierr = VecScale( M->dxidr_b, -1.0 );CHKERRQ(ierr);
    ierr = VecShift( M->dxidr_b, P->radius );CHKERRQ(ierr);
    ierr = VecScale( M->dxidr_b, P->beta );CHKERRQ(ierr);
    ierr = VecExp( M->dxidr_b );CHKERRQ(ierr);
    ierr = VecScale( M->dxidr_b, P->rhos );CHKERRQ(ierr);
    ierr = VecScale( M->dxidr_b, 1.0 / M->mantle_density );CHKERRQ(ierr);
    /* multiply radius squared */
    ierr = VecPointwiseMult( M->dxidr_b, M->dxidr_b, M->radius_b );CHKERRQ(ierr);
    ierr = VecPointwiseMult( M->dxidr_b, M->dxidr_b, M->radius_b );CHKERRQ(ierr);
    /* divide xi squared */
    ierr = VecPointwiseDivide( M->dxidr_b, M->dxidr_b, M->xi_b );CHKERRQ(ierr);
    ierr = VecPointwiseDivide( M->dxidr_b, M->dxidr_b, M->xi_b );CHKERRQ(ierr);

    /* TODO: last value of M->dxidr_b is inf, and values nearer the base of the
       mantle increase substantially.  Perhaps an argument for shifting the cmb
       away from zero mass coordinate */
    
    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

