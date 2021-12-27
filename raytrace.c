/* Ray tracing interface modified to NIP tomography. */
/*
 Copyright (C) 2004 University of Texas at Austin
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <rsf.h>
#include <stdbool.h>
#include "raytrace.h"
#include "grid2.h"
#include "atela.h"
#include "dynamic.h"
#include "interface2d.h"
#include <math.h>
/*^*/

#ifndef _raytrace_h

#define DEG2RAD SF_PI/180. /* Degrees to radians conversion */
#define ORDER 4 /* Ray tracing interpolation order */
/*^*/

typedef struct RayTrace* raytrace;
/* abstract data type */
/*^*/

#endif

struct RayTrace {
    bool sym;
    int dim, nt;
    float dt, z0;
    grid2 grd2;
};
/* concrete data type */

static void iso_rhs(void* par, float* y, float* f){}

float getvelocity(void *par, float **traj, int i)
/* Get velocity from a point in the ray trajectory
Note: The i variable is the sample index of the ray trajectory traj. This
function returns the velocity from the grid in the sample location
(z,x)=(traj[i][0],traj[i][1]).
*/
{
	raytrace rt;
	float *x;

	rt = (raytrace) par;

	x=sf_floatalloc(2);

	x[0]=traj[i][0];
	x[1]=traj[i][1];
	
	return sqrtf(1./grid2_vel(rt->grd2,x));
}

static void dyn_iso_rhs(void* par, float dvdn, float *x, float* y, float* f)
/* right-hand side for isotropic raytracing */
{    
	raytrace rt;
	float v;
	
	rt = (raytrace) par;

	v = sqrtf(1./grid2_vel(rt->grd2,x));
	
	f[0]   = v*v*y[1]; /* v^2 p */
	f[1] = (-1./v)*dvdn*y[0]; /* -1/v dv/dn q */
}

float second_derivative(void *par, float *n, float *x, float v)
/*< Second derivative >*/
{
	raytrace rt;
	float vpdx, vmdx;
	float dx=0.001;
	float *tmp;

	tmp = sf_floatalloc(2);

	rt = (raytrace) par;

	n[0]*=dx;
	n[1]*=dx;

	// ppdx
	tmp[0] = x[0]+n[1];
	tmp[1] = x[1]-n[0];
	vpdx = sqrtf(1./grid2_vel(rt->grd2,tmp));

	// pmdx
	tmp[0]=x[0]-n[1];
	tmp[1]=x[1]+n[0];
	vmdx = sqrtf(1./grid2_vel(rt->grd2,tmp));

	return (vpdx-2.*v+vmdx)/(dx*dx);
}

float calculateRNIPWithDynamicRayTracing(
					  void *par,
					  float dt,
					  float nt,
					  float **traj,
					  float v0
)
/*< Second derivative >*/
{

	float v;
	int it;
	raytrace rt;
	float *x;
	float *n;
	float *dvdn;
	float mod;
	float rnip;

    	rt = (raytrace) par;
	x = sf_floatalloc(2);
	n = sf_floatalloc(2);
	dvdn = sf_floatalloc(nt-2);

	for(it=0;it<nt-2;it++){
		x[0]=traj[it][0];
		x[1]=traj[it][1];
		n[0] = (traj[it+1][0]-traj[it][0]);
		n[1] = (traj[it+1][1]-traj[it][1]);
		mod = sqrtf(n[0]*n[0]+n[1]*n[1]);
		n[0] /= mod;
		n[1] /= mod;

		v = sqrtf(1./grid2_vel(rt->grd2,x));
		dvdn[it]=second_derivative(rt,n,x,v);
	}

	x[0]=0.; // q=0
	x[1]=1.; // p=1
	sf_dynamic_runge_init(2,nt-2,2*dt);
	rnip = sf_dynamic_runge_step(x,rt,dyn_iso_rhs,dvdn,traj,v0);
	sf_dynamic_runge_close();

	return rnip;
}

float snellslaw(float ei, float vi, float vt)
/*< Snells Law to get transmission angle >*/
{
	return asinf((vt/vi)*sinf(ei));
}

float calculateeiAngle(float **traj, int ir)
/*< Calculate incident angle >*/
{
	float *x;
	float ei;

	x = sf_floatalloc(2);
	x[0] = traj[ir][0]-traj[ir-1][0];
	x[1] = traj[ir][1]-traj[ir-1][1];

	ei = sqrtf(x[0]*x[0]+x[1]*x[1]);
	ei = acos(-x[0]/ei);

	return ei;
}

void hubralTransmissionLaw(float *rnip, float vt, float vi, float ei)
/*< Calculate transmited RNIP using Hubral's transmission law >*/
{
	float et;
	float ri;

	et = snellslaw(ei,vi,vt);

	ri = *rnip;
	ri = (1./ri);
	ri *= (cosf(ei)*cosf(ei))/(cosf(et)*cosf(et));
	ri = (vt/vi)*ri;
	*rnip = 1./ri;
}

void transmitedRNIPThroughInterface(
					void *par, /* Raytrace struct */
					void *interface,
					int *ir, /* ray sample index */
					float **traj, /* ray trajectory */
					float *rnip /* RNIP parameter */)
/*< Calculate transmited RNIP parameter through interface using Hubral laws
Note:
>*/
{
	float vi, vt=0.;
	raytrace rt;
	itf2d it2;
	int i;
	float zi;
	float ei;
	int pass=0; // Ray passed through interface?
        
	rt = (raytrace) par;
	it2 = (itf2d) interface;

	vi = getvelocity(rt,traj,*ir);

	vt = getvelocity(rt,traj,++(*ir));
	
	while(vi!=vt){
		*rnip+=2*vt*rt->dt;

		vi = vt;

		vt = getvelocity(rt,traj,++(*ir));
		zi = getZCoordinateOfInterface(it2,traj[*ir][1]);
		if(zi>traj[*ir][0] && pass!=1){
			ei = calculateeiAngle(traj,*ir-1);
			hubralTransmissionLaw(rnip,vt,vi,ei);
			pass = 1;
			vt = getvelocity(rt,traj,++(*ir));
		}
	}

}

float calculateRNIPWithHubralLaws(
				  void *par, /* Raytrace struct */
				  float** traj, /* Normal ray trajectory */
				  int nt, /* Normal ray times samples */
				  float *v, /* layers velocities */
				  float t0, /* Normal ray traveltime */
				  int itf, /* Interface index */
				  float *sz,
				  int nsz,
				  float osz,
				  float dsz,
				  int nx)
/*< Calculate RNIP parameter using Hubral's transmission and propagation laws

Note: 
>*/
{

	int i, j; // loop counter
	float rnip=0.; // RNIP parameter
	raytrace rt; // raytrace struct
	float vt, vi; // velocities
	float *x; // (z,x) position vector
	float* szz; // Z coordinates of interface being inverted
	itf2d interface;

	szz = sf_floatalloc(nsz/2);

	//interface = itf2d_init(szz,nsz/2,osz,dsz);

	rt = (raytrace) par;
	x = sf_floatalloc(2);

	x[0] = traj[0][0];
	x[1] = traj[0][1];
	vi = sqrtf(1./grid2_vel(rt->grd2,x));
	vt = vi;
	rnip+=2*vt*rt->dt;

	for(i=1;i<nt;i++){
		x[0] = traj[i][0];
		x[1] = traj[i][1];
		vt = sqrtf(1./grid2_vel(rt->grd2,x));

		/* If the ray reaches interface use transmission law */
		if(vt!=vi){
			for(j=0;j<nsz/2;j++){
				szz[j]=sz[j+((itf-1)*nsz/2)];
			}
			interface = itf2d_init(szz,nsz/2,osz,dsz);
			//itf2d_setZNodepoints(interface,szz);
			transmitedRNIPThroughInterface(rt,interface,&i,traj,&rnip);
			vi=vt;
		}

		/* Propagation law */
		rnip+=2*vt*rt->dt;
	}

	return rnip;
}

static int term(void* par, float* y)
/* grid termination */
{
    raytrace rt;
    
    rt = (raytrace) par;
	
    switch (rt->dim) {
		case 2:
			return grid2_term(rt->grd2,y);
		default:
			sf_error("%s: Cannot raytrace with dim=%d",__FILE__,rt->dim);
			return 0;
    }
}

raytrace raytrace_init(int dim            /* dimensionality (2 or 3) */, 
					   bool sym,          /* if symplectic */
					   int nt             /* number of ray tracing steps */, 
					   float dt           /* ray tracing step (in time) */,
					   int* n             /* slowness dimensions [dim] */, 
					   float* o, float* d /* slowness grid [dim] */,
					   float* slow2       /* slowness squared [n3*n2*n1] */, 
					   int order          /* interpolation order */)
/*< Initialize ray tracing object. 
 * Increasing order increases accuracy but
 decreases efficiency. Recommended values: 3 or 4.
 * slow2 can be changed or deallocated after
 raytrace_init.
 >*/
{
    raytrace rt;
    
    rt = (raytrace) sf_alloc (1,sizeof(*rt));
    
    rt->dim = dim;
    rt->sym = sym;
    rt->nt = nt;
    rt->dt = dt;
    rt->z0 = o[0];
    
    switch (dim) {
		case 2:
			rt->grd2 = grid2_init (n[0], o[0], d[0], 
								   n[1], o[1], d[1],
								   slow2, order);
			break;
		default:
			sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
    }
	
    return rt;
}

void raytrace_close (raytrace rt)
/*< Free internal storage >*/
{
    switch (rt->dim) {
		case 2:
			grid2_close (rt->grd2);
			break;
    }
    free (rt);
}

int trace_ray (raytrace rt  /* ray tracing object */, 
			   float* x     /* point location {z,y,x} [dim] */, 
			   float* p     /* ray parameter vector [dim] */, 
			   float** traj /* output ray trajectory [nt+1,dim] */)
/*< Trace a ray.
 * Values of x and p are changed inside the function.
 * The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D
 * Vector p points in the direction of the ray. 
 The length of the vector is not important.
 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*cos(a); p[2] = sin(b)*sin(a) in 3-D
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians
 * The output code for it = trace_ray(...)
 it=0 - ray traced to the end without leaving the grid
 it>0 - ray exited at the top of the grid
 it<0 - ray exited at the side or bottom of the grid
 * The total traveltime along the ray is 
 nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int i, dim, it=0, nt;
    float y[6], s2;
	
    dim = rt->dim;
    nt = rt->nt;
	
    if (!rt->sym) {
		switch (dim) {
			case 2:
				s2 = grid2_vel(rt->grd2,x);
				break;
			default:
				s2 = 0.;
				sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
		}
		
		for (i=0; i < dim; i++) {
			y[i] = x[i];
			y[i+dim] = p[i]*sqrtf(s2);
		}
		
		sf_runge_init(2*dim, nt, rt->dt);
		it = sf_ode23_step (y, rt,iso_rhs,term,traj);
		sf_runge_close();
		
		for (i=0; i < dim; i++) {
			x[i] = y[i];
			p[i] = y[i+dim];
		}
    } else {
		switch (dim) {
			case 2:
				it = atela_step (dim, nt, rt->dt, true, x, p, 
								 rt->grd2, 
								 grid2_vgrad, grid2_vel, grid2_term, traj);
				break;
			default:
				sf_error("%s: cannot handle %d dimensions",__FILE__,rt->dim);
				break;
		}
    }
    
    if (it > 0 && x[0] > rt->z0) {
		return (-it); /* exit through the side or bottom */
    } else {
		return it;
    }
}

/* 	$Id$	 */
