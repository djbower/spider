#include "mesh.h"

static PetscErrorCode regular_mesh( Ctx * );
//static PetscErrorCode geometric_mesh( Ctx * );
static PetscErrorCode spherical_area( DM, Vec, Vec );
static PetscErrorCode spherical_volume( Ctx *, Vec, Vec );
static PetscErrorCode mixing_length( DM, Vec, Vec, Parameters const );
static PetscScalar mixing_length_conventional( Parameters const, PetscScalar );
static PetscScalar mixing_length_average( Parameters const, PetscScalar );
static PetscScalar mixing_length_layer( Parameters const, PetscScalar );
static PetscScalar get_layer( DM, Vec, Vec, Parameters const );
static PetscErrorCode aw_density( DM, Vec, Vec, Parameters const );
static PetscErrorCode aw_pressure( DM, Vec, Vec, Parameters const );
static PetscErrorCode aw_pressure_gradient( DM, Vec, Vec, Parameters const );
static PetscErrorCode aw_mass( Mesh * );

PetscErrorCode set_mesh( Ctx *E)
{

    Mesh           *M = &E->mesh;
    DM             da_b=E->da_b, da_s=E->da_s;
    Parameters     P = E->parameters;

    PetscFunctionBeginUser;

    /* for regular mesh, although without resolving the ultra-thin
       thermal boundary layer at the base of the mantle this likely
       gives wrong results */
    regular_mesh( E );

    /* need to use geometric mesh to resolve ultra-thin thermal
       boundary layer at the base of the mantle */
    //geometric_mesh( E );

    /* Adams-Williamson EOS */

    /* pressure at basic nodes */
    aw_pressure( da_b, M->radius_b, M->pressure_b, P);

    /* dP/dr at basic nodes */
    aw_pressure_gradient( da_b, M->radius_b, M->dPdr_b, P);

    /* pressure at staggered nodes */
    aw_pressure( da_s, M->radius_s, M->pressure_s, P);

    /* dP/dr at staggered nodes */
    aw_pressure_gradient( da_s, M->radius_s, M->dPdr_s, P );

    /* surface area at basic nodes, without 4*pi term */
    spherical_area( da_b, M->radius_b, M->area_b);

    /* surface area at staggered nodes, without 4*pi term */
    spherical_area( da_s, M->radius_s, M->area_s );

    /* volume of spherical cells, without 4*pi term */
    spherical_volume( E, M->radius_b, M->volume_s);

    /* mixing length is minimum distance from boundary */
    mixing_length( da_b, M->radius_b, M->mix_b, P);

    /* layer id.  0 everywhere for single layer (as determined by
       P->mixing_length), and 0 for upper and 1 for lower layer
       when P->mixing_length==3 */
    get_layer( da_b, M->radius_b, M->layer_b, P );

    /* density at staggered nodes */
    aw_density( da_s, M->radius_s, M->rho_s, P );

    /* mass at staggered nodes */
    aw_mass( M );

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
    dx_b = -P->radius*(1.0-P->coresize) / (numpts_b-1);

    /* radius at basic nodes */
    ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ierr = DMDAVecGetArray(da_b,M->radius_b,&arr);CHKERRQ(ierr);
    for (i=ilo_b; i<ihi_b; ++i){
        arr[i] = P->radius*P->coresize - (numpts_b-1-i)*dx_b;
    }
    ierr = DMDAVecRestoreArray(da_b,M->radius_b,&arr);CHKERRQ(ierr);

    /* radius at staggered nodes */
    ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAVecGetArray(da_s,M->radius_s,&arr);CHKERRQ(ierr);
    for (i=ilo_s;i<ihi_s;++i){
        arr[i] = P->radius*P->coresize - 0.5*dx_b - (numpts_s-1-i)*dx_b;
    }
    ierr = DMDAVecRestoreArray(da_s,M->radius_s,&arr);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}

#if 0
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

static PetscErrorCode mixing_length( DM da, Vec radius, Vec mix, const Parameters P )
{
    PetscErrorCode    ierr;
    PetscScalar       *arr_m;
    const PetscScalar *arr_r;
    PetscInt          i,ilo,ihi,w;

    PetscFunctionBeginUser;
    ierr = DMDAGetCorners(da,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;
    ierr = DMDAVecGetArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,mix,&arr_m);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
      // conventional mixing length
      if(P->mixing_length == 1){
        arr_m[i] = mixing_length_conventional( P, arr_r[i] );
      }
      // avg mixing length from conventional theory, i.e. 1/4*depth */
      else if(P->mixing_length == 2){
        arr_m[i] = mixing_length_average( P, arr_r[i] );
      }
      // two layers
      else if(P->mixing_length == 3){
        arr_m[i] = mixing_length_layer( P, arr_r[i] );
      }
    }
    ierr = DMDAVecRestoreArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,mix,&arr_m);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


static PetscScalar mixing_length_conventional( const Parameters P, const PetscScalar rad_in )
{
    /* distance to nearest top or bottom boundary */

    PetscScalar rad1, rad2, out;

    rad1 = rad_in - P->radius * P->coresize;
    rad2 = P->radius - rad_in;
    out = PetscMin( rad1, rad2 );

    return out;   
}


static PetscScalar mixing_length_average( const Parameters P, const PetscScalar rad_in )
{
    /* constant, based on average of distance to nearest top or bottom
       boundary */

    PetscScalar out;

    out = 0.25 * P->radius * (1.0-P->coresize);

    return out;
}

static PetscScalar mixing_length_layer( const Parameters P, const PetscScalar rad_in )
{
    /* account for an additional interface separating two convecting
       layers */

    PetscScalar out, dist3;

    out = mixing_length_conventional( P, rad_in );
    dist3 = rad_in - P->mixing_length_layer_radius;
    dist3 = PetscAbsReal( dist3 );
    out = PetscMin( out, dist3 );

    return out;
}

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

static PetscErrorCode aw_density( DM da, Vec radius, Vec density, const Parameters P )
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

static PetscErrorCode aw_mass( Mesh *M )
{
    PetscErrorCode ierr;
 
    PetscFunctionBeginUser;

    ierr = VecCopy( M->rho_s, M->mass_s ); CHKERRQ(ierr);

    ierr = VecPointwiseMult( M->mass_s, M->mass_s, M->volume_s );
    /* excludes 4*pi prefactor */
    ierr = VecSum( M->mass_s, &M->mantle_mass );

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
