#include "mesh.h"

//static PetscErrorCode regular_mesh( Ctx * );
static PetscErrorCode geometric_mesh( Ctx * );
static PetscErrorCode spherical_area( DM, Vec, Vec );
static PetscErrorCode spherical_volume( Ctx *, Vec, Vec );
static PetscErrorCode mixing_length( DM, Vec, Vec );
static PetscErrorCode aw_pressure( DM, Vec, Vec );
static PetscErrorCode aw_pressure_gradient( DM, Vec, Vec );

PetscErrorCode set_mesh( Ctx *E)
{

    PetscErrorCode ierr;
    PetscInt       i;
    Mesh           *M;
    DM             da_b=E->da_b, da_s=E->da_s;

    PetscFunctionBeginUser;
#if (defined VERBOSE)
    ierr = PetscPrintf(PETSC_COMM_WORLD,"set_mesh:\n");CHKERRQ(ierr);
#endif

    M = &E->mesh;

    /* Create vectors required for the mesh */
    for (i=0;i<NUMMESHVECS_B;++i) {
      ierr = DMCreateGlobalVector(E->da_b,&M->meshVecs_b[i]);CHKERRQ(ierr);
    }
    // TODO: We should be creating ALL the vecs at the same time - do that as we roll all this into one function (create these vecs outside this function)

    M->area_b     = M->meshVecs_b[0];
    M->dPdr_b     = M->meshVecs_b[1];
    M->pressure_b = M->meshVecs_b[2];
    M->radius_b   = M->meshVecs_b[3];
    M->mix_b      = M->meshVecs_b[4];

    for (i=0;i<NUMMESHVECS_S;++i) {
      ierr = DMCreateGlobalVector(E->da_s,&M->meshVecs_s[i]);CHKERRQ(ierr);
    }
    M->pressure_s = M->meshVecs_s[0];
    M->radius_s   = M->meshVecs_s[1];
    M->volume_s   = M->meshVecs_s[2];
    M->dPdr_s     = M->meshVecs_s[3];
    M->area_s     = M->meshVecs_s[4];

    /* for regular mesh, although without resolving the ultra-thin
       thermal boundary layer at the base of the mantle this likely
       gives wrong results */
    //regular_mesh( E );

    /* need to use geometric mesh to resolve ultra-thin thermal
       boundary layer at the base of the mantle */
    geometric_mesh( E );

    /* Adams-Williamson EOS */

    /* pressure at basic nodes */
    aw_pressure( da_b, M->radius_b, M->pressure_b);

    /* dP/dr at basic nodes */
    aw_pressure_gradient( da_b, M->radius_b, M->dPdr_b);

    /* pressure at staggered nodes */
    aw_pressure( da_s, M->radius_s, M->pressure_s);

    /* dP/dr at staggered nodes */
    aw_pressure_gradient( da_s, M->radius_s, M->dPdr_s );

    /* surface area at basic nodes, without 4*pi term */
    /* and now non-dimensional */
    spherical_area( da_b, M->radius_b, M->area_b);

    /* surface area at staggered nodes, without 4*pi term */
    /* and now non-dimensional */
    spherical_area( da_s, M->radius_s, M->area_s );

    /* volume of spherical cells, without 4*pi term */
    /* and now non-dimensional */
    spherical_volume( E, M->radius_b, M->volume_s);

    /* mixing length is minimum distance from boundary */
    mixing_length( da_b, M->radius_b, M->mix_b);

    PetscFunctionReturn(0);
}

/*static PetscErrorCode regular_mesh( Ctx *E )
{

    PetscErrorCode ierr;
    PetscScalar    *arr;
    PetscInt       i,ilo_b,ihi_b,ilo_s,ihi_s,w_b,w_s,numpts_b,numpts_s;
    Mesh           *M;
    DM             da_b=E->da_b, da_s=E->da_s;

    PetscFunctionBeginUser;

    M = &E->mesh;

    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);*/

    /* basic node spacing (negative) */
    //M->dx_b = -(RADOUT-RADIN) / (numpts_b-1);

    /* staggered node spacing (negative) */
    //M->dx_s = M->dx_b;

    /* radius at basic nodes */
    /*ierr = DMDAGetCorners(da_b,&ilo_b,0,0,&w_b,0,0);CHKERRQ(ierr);
    ihi_b = ilo_b + w_b;
    ierr = DMDAVecGetArray(da_b,M->radius_b,&arr);CHKERRQ(ierr);
    for (i=ilo_b; i<ihi_b; ++i){
        arr[i] = RADIN - (numpts_b-1-i)*M->dx_b;
    }
    ierr = DMDAVecRestoreArray(da_b,M->radius_b,&arr);CHKERRQ(ierr);*/

    /* radius at staggered nodes */
    /*ierr = DMDAGetCorners(da_s,&ilo_s,0,0,&w_s,0,0);CHKERRQ(ierr);
    ihi_s = ilo_s + w_s;
    ierr = DMDAVecGetArray(da_s,M->radius_s,&arr);CHKERRQ(ierr);
    for (i=ilo_s;i<ihi_s;++i){
        arr[i] = RADIN - 0.5*M->dx_b - (numpts_s-1-i)*M->dx_b;
    }
    ierr = DMDAVecRestoreArray(da_s,M->radius_s,&arr);CHKERRQ(ierr);

    PetscFunctionReturn(0);

}*/

static PetscErrorCode geometric_mesh( Ctx *E )
{

    PetscErrorCode ierr;
    PetscScalar    *arr_b, *arr_s;
    PetscInt       i,ilo_b,ihi_b,ilo_s,ihi_s,w_b,w_s,numpts_b,numpts_s;
    Mesh           *M;
    DM             da_b=E->da_b, da_s=E->da_s;
    // TODO: should not be hard-coded
    PetscScalar    dr_min = 1.5696123057604772e-10; // 1 mm
    // 2939 basic nodes
    //PetscScalar    dr_max = 0.0001569612305760477; // 1 km
    // 656 basic nodes
    //PetscScalar    dr_max = 0.0007848061528802385; // 5 km
    // 372 basic nodes
    //PetscScalar    dr_max = 0.001569612305760477; // 10 km
    // 278 basic nodes
    //PetscScalar    dr_max = 0.0023544184586407153; // 15 km
    PetscScalar    geom_fac = 1.2;
    PetscScalar    dr, dr_prev;

    PetscFunctionBeginUser;

    M = &E->mesh;

    ierr = DMDAGetInfo(E->da_b,NULL,&numpts_b,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = DMDAGetInfo(E->da_s,NULL,&numpts_s,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);CHKERRQ(ierr);

    /* TODO: dynamically determine number of mesh points, and
       automatically ensure that arr[0] = RADOUT and
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
    arr_b[0] = RADOUT;

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
        arr_volume[i] *= 1.0 / 3.0;
    }
    // TODO: here and elsewhere, it's very dangerous to use the same indice to refer to two DAs without checking that the local ranges are valid. 
    ierr = DMDAVecRestoreArrayRead(da_b,radius_local,&arr_radius);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da_s,volume,&arr_volume);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode mixing_length( DM da, Vec radius, Vec mix )
{
    PetscErrorCode    ierr;
    PetscScalar       *arr_m;
    const PetscScalar *arr_r;
    PetscInt          i,ilo,ihi,w;
    PetscScalar       rad1, rad2;

    PetscFunctionBeginUser;
    ierr = DMDAGetCorners(da,&ilo,0,0,&w,0,0);CHKERRQ(ierr);
    ihi = ilo + w;
    ierr = DMDAVecGetArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da,mix,&arr_m);CHKERRQ(ierr);
    for(i=ilo; i<ihi; ++i){
        rad1 = arr_r[i] - RADIN;
        rad2 = RADOUT - arr_r[i];
        arr_m[i] = PetscMin( rad1, rad2 );
    }
    ierr = DMDAVecRestoreArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,mix,&arr_m);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode aw_pressure( DM da, Vec radius, Vec pressure )
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
        dep = RADOUT - arr_r[i];
        arr_p[i] = RHOS * GRAVITY / BETA;
        arr_p[i] *= PetscExpScalar( BETA * dep ) - 1.0;
    }
    ierr = DMDAVecRestoreArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,pressure,&arr_p);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

static PetscErrorCode aw_pressure_gradient( DM da, Vec radius, Vec grad )
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
        dep = RADOUT - arr_r[i];
        arr_g[i] = -RHOS * GRAVITY * PetscExpScalar( BETA * dep );
    }
    ierr = DMDAVecRestoreArray(da,grad,&arr_g);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da,radius,&arr_r);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

