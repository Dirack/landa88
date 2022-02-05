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

float getVelocityForRaySampleLocation(
		  void *par, /* Raytrace struct */
		  float **traj, /* Ray trajectory */
		  int i /* Sample index in ray trajectory */)
/*< Get velocity from a point in the ray trajectory
Note: The i variable is the sample index of the ray trajectory traj. This
function returns the velocity from the grid in the sample location
(z,x)=(traj[i][0],traj[i][1]).
>*/
{
	raytrace rt;
	float x[2];

	rt = (raytrace) par;

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
/*< Calculate the second derivative in the direction normal to ray trajectory >*/
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
					  void *par, /* Raytrace struct */
					  float dt, /* Time sampling */
					  float nt, /* Number of time samples */
					  float **traj, /* Ray trajectory (z,x) */
					  float v0 /* Near surface velocity */
)
/*< Calculate RNIP with dynamic ray tracing >*/
{
	float v; // velocity
	int it; // Loop counter
	raytrace rt; // Raytrace struct
	float *x; // Sample coordinate (z,x)
	float *n; // Normal vector
	float *dvdn; // Derivative normal to ray direction
	float mod; // tmp variable
	float rnip; // RNIP parameter

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

		/* Calculate derivative for each ray sample */
		dvdn[it]=second_derivative(rt,n,x,v);

	} // Loop over ray samples

	/* Initial conditions for a point source */
	x[0]=0.; // q=0
	x[1]=1.; // p=1

	/* Fourth order Runge-Kutta dynamic ray tracing */
	sf_dynamic_runge_init(2,nt-2,2*dt);
	rnip = sf_dynamic_runge_step(x,rt,dyn_iso_rhs,dvdn,traj,v0);
	sf_dynamic_runge_close();

	return rnip;
}

float getTransmissionAngleSnellsLaw(
		float ei, /* Incident angle in radians */
		float vi, /* Velocity in incident ray layer */
		float vt /* Velocity in transmited ray layer */)
/*< Snells Law to get transmission angle >*/
{
	return asinf((vt/vi)*sinf(ei));
}

float calculateIncidentAngle(
			     void *par, /* Interface struct */
			     float **traj, /* Ray trajectory */
			     int ir /* Ray sample index */)
/*< Calculate incident angle 
Note: This function calculates a vector using the diffence between
traj[ir] and traj[ir-1]. The angle between this vector and the normal
to the interface will be the incident angle returned
>*/
{
	float x[2];
	float n[2];
	float ei;
	itf2d it2;

	it2 = (itf2d) par;

	n[0] = getZCoordinateOfInterface(it2,traj[ir][1]+0.01)
	-getZCoordinateOfInterface(it2,traj[ir][1]);
	n[1] = 0.01;

	ei = sqrtf(n[0]*n[0]+n[1]*n[1]);
	n[0]/=ei;
	n[1]/=ei;

	ei = n[0];
	n[0]=-n[1];
	n[1]=ei;

	x[0] = traj[ir][0]-traj[ir-1][0];
	x[1] = traj[ir][1]-traj[ir-1][1];

	ei = sqrtf(x[0]*x[0]+x[1]*x[1]);
	ei = acos((x[0]*n[0]+x[1]*n[1])/ei);

	return ei;
}

float calculateInterfaceCurvature(
				  void *par, /* Interface struct */
				  float x /* Coordinate (z,x) */)
/*< Calculate interface curvature using curvature formula for 2D function >*/
{

	float kf; // Interface curvature
	float f1, f2; // tmp variables
	int is; // Spline index
	itf2d it2; // Interface struct
	float coef[4]; // Cubic spline coefficients matrix
	float a, b, c; // Spline coefficients

	it2 = (itf2d) par;
	is = (x-itf2d_o(it2))/itf2d_d(it2);
	itf2d_getSplineCoefficients(it2,coef,is);

	a = coef[0];
	b = coef[1];
	c = coef[2];

	f1 = 3*a*x*x+2*b*x+c;
	f2 = fabs(6*a*x+2*b);
	f1 = 1+f1*f1;
	f1 = f1*f1*f1;
	kf = f2/(sqrtf(f1));

	return kf;
}

void getTransmitedRNIPHubralTransmissionLaw(
			   float *rnip, /* RNIP in incident ray layer */
			   float vt, /* Velocity in transmited ray layer */
			   float vi, /* Velocity in incident ray layer */
			   float ei, /* Ray incident angle*/
			   float kf /* Interface curvature */)
/*< Calculate transmited RNIP through interface using Hubral's transmission law
Note: RNIP parameter is modified inside the function to transmited RNIP value
>*/
{
	float et;
	float ri;

	et = getTransmissionAngleSnellsLaw(ei,vi,vt);

	ri = *rnip;
	ri = (1./ri);
	ri *= (cosf(ei)*cosf(ei))/(cosf(et)*cosf(et));
	ri = (vt/vi)*ri;

	ri = ri+ (1./(cosf(et)*cosf(et))*((vt/vi)*cosf(ei)-cosf(et)))*kf;

	*rnip = 1./ri;
}

void transmitedRNIPThroughInterface(
					void *par, /* Raytrace struct */
					void *interface, /* Interface struct */
					int *ir, /* ray sample index */
					float **traj, /* ray trajectory */
					float *rnip /* RNIP parameter */)
/*< Calculate transmited RNIP parameter through interface using Hubral laws
Note: Velocity model interpolation makes a transition velocity zone close to
the interface, in this zone is applied transmission law where the ray passes
through interface
>*/
{
	float vi; // Velocity - incident ray layer
	float vt=0.; // Velocity - transmited ray layer
	raytrace rt; // Raytrace struct
	itf2d it2; // Interface struct
	float zi; // Interface z coordinate
	float ei; // Incident angle (radians)
	int pass=false; // Ray passed through interface?
	float kf; // Interface curvature
        
	rt = (raytrace) par;
	it2 = (itf2d) interface;

	vi = getVelocityForRaySampleLocation(rt,traj,*ir);

	vt = getVelocityForRaySampleLocation(rt,traj,++(*ir));
	
	while(vi!=vt){
		*rnip+=2*vt*rt->dt;

		vi = vt;

		vt = getVelocityForRaySampleLocation(rt,traj,++(*ir));
		zi = getZCoordinateOfInterface(it2,traj[*ir][1]);
		if(zi>traj[*ir][0] && pass==false){
			ei = calculateIncidentAngle(it2,traj,*ir-1);
			kf = calculateInterfaceCurvature(it2,traj[*ir][1]);
			getTransmitedRNIPHubralTransmissionLaw(rnip,vt,vi,ei,kf);
			pass = true;
			vt = getVelocityForRaySampleLocation(rt,traj,++(*ir));
		}
	}

}

float calculateRNIPWithHubralLaws(
				  void *par, /* Raytrace struct */
				  float** traj, /* Normal ray trajectory */
				  int nt, /* Normal ray times samples */
				  float *v, /* layers velocities */
				  int nv, /* Number of layers */
				  float t0, /* Normal ray traveltime */
				  int itf, /* Interface index */
				  float *sz, /* Splines nodepoints */
				  int nsz, /* Number of nodepoints */
				  float osz, /* Nodepoints origin */
				  float dsz /* Nodepoints sampling */)
/*< Calculate RNIP parameter using Hubral's transmission and propagation laws >*/
{

	int i, j; // loop counter
	float rnip=0.; // RNIP parameter
	raytrace rt; // raytrace struct
	float vt, vi; // velocities
	float *x; // (z,x) position vector
	float* szz; // Z coordinates of interface being inverted
	itf2d interface; // Interface struct

	szz = sf_floatalloc(nsz/(nv-1));

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
			for(j=0;j<nsz/3;j++){
				szz[j]=sz[j+((itf-1)*nsz/3)];
			}
			interface = itf2d_init(szz,nsz/3,osz,dsz);
			transmitedRNIPThroughInterface(rt,interface,&i,traj,&rnip);
			vi=vt;
		}

		/* Propagation law */
		rnip+=2*vt*rt->dt;
	}

	return rnip;
}

void sortingXinAscendingOrder(
				float *x, /* x vector to sort */
				float *z, /* z(x) vector */
				int n /* Vectors dimension */)
/*< x vector sorting in ascending order >*/
{
	int i; // Loop counter
	float tmpx, tmpz; // Temporary variables
	int k; // Sorting key (number of changes)

	do{
		k=0;
		for(i=1;i<n;i++){
			if(x[i-1]>x[i]){
				tmpx=x[i-1]; tmpz=z[i-1];
				x[i-1]=x[i]; z[i-1]=z[i];
				x[i]=tmpx; z[i]=tmpz;
				k++;
			}
		} // Loop vector samples
	}while(k!=0);
}

void interfaceInterpolationFromNipSources(float **s, /* NIP sources */
					  int ns, /* Number of NIP sources */
					  float *sz, /* Spline nodepoints */
					  int nsz, /* Number of nodepoints */
					  float osz, /* Nodepoints origin */
					  float dsz, /* Nodepoints sampling */
					  int nsv /* Number of layers */)
/*< Use NIP sources location to draw interfaces 
Note: If the velocity model is correct the NIP sources location coincides with interfaces. So, they can be used to draw interface through cubic spline interpolation.
>*/
{
	int nxs=ns/(nsv-1); // Number of NIP sources for each interface
	int nxsz=nsz/(nsv-1); // Number of nodepoints for each interface
	int i, im; // Loop counter
	float *tsx, *tsz; // Temporary spline vector
	float *coef; // Coefficients matrix
	float xx, xs; // X coordinate
	float oxs; // Spline's origin
	int l; // Spline index

	tsx = sf_floatalloc(nxs);
	tsz = sf_floatalloc(nxs);
	coef = sf_floatalloc(4*(nxs-1));

	for(i=0;i<2;i++){
		l=0;
		for(im=0;im<nxs;im++){
			tsz[im]=s[i*nxs+im][0];
			tsx[im]=s[i*nxs+im][1];
		}
		sortingXinAscendingOrder(tsx,tsz,nxs);
		calculateSplineCoeficients(nxs,tsx,tsz,coef);
		oxs=tsx[0];
		for(im=0;im<nxsz;im++){
			xx=im*dsz+osz;
			if(xx<tsx[0] || xx>tsx[nxs-1]){
				sz[i*nxsz+im]=tsz[l];
			}else{
				while(xx>tsx[l+1]){
					l++;
					oxs=tsx[l];
				}
				xs=xx-oxs;
				sz[i*nxsz+im]=coef[l*4+0]*xs*xs*xs+coef[l*4+1]*xs*xs+coef[l*4+2]*xs+coef[l*4+3];
			}
		}
	}
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
