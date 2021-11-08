/*
	 tomography.c (c)
	 
	 Purpose: 'Mvfsacrsnh.c' library for raytracing and traveltime
	 calculation.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 19/09/2019

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <rsf.h>
#include "raytrace.h"
#include "tomography.h"

#ifndef GDB_DEBUG
	#define DSLOW 0.04
	#define DANGLE 1.0
#else
	#define DSSLOW 0.04
	#define DANGLE 1.0
#endif
/*^*/

float creTimeApproximation(float h, // Half-offset
			 float m, // CMP
			 float v0, // Near surface velocity
			 float t0, // Normal ray traveltime
			 float m0, // Central CMP
			 float RNIP, // CRE parameter
			 float BETA, // CRE parameter
			 bool cds // Use CDS condition?
			 )
/*< Calculate CRE traveltime approximation t(m,h)
Note: If cds parameter is false, it uses the CRE formula to calculate traveltime.
If cds parameter is true, it uses the non-hyperbolic CRS formula with CDS condition (RN=RNIP)
to calculate traveltime.
>*/
{ 
	float alpha; // Asymmetry parameter
	float d = m-m0; // Distance to central CMP m0
	float c1; // CRE coefficient
	float c2; // CRE coefficient
	float a1, a2, b2, b1, Fd, Fd1, Fd2; // Non-hyperbolic CRS coefficient
	float t; // traveltime t(m,h)

	if(!cds){
		c1 = (d+h)/RNIP;
		c2 = (d-h)/RNIP;
		alpha = sin(BETA)/RNIP;
		t = (t0-2*RNIP/v0)+(RNIP/v0)*sqrt(1-2*alpha*(d+h)+c1*c1)+(RNIP/v0)*sqrt(1-2*alpha*(d-h)+c2*c2);
	}else{

		a1=(2*sin(BETA))/(v0);
		a2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RNIP);
		b2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RNIP);
		b1=2*b2+a1*a1-a2;
		Fd=(t0+a1*d)*(t0+a1*d)+a2*d*d;
		Fd2=(t0+a1*(d-h))*(t0+a1*(d-h))+a2*(d-h)*(d-h);
		Fd1=(t0+a1*(d+h))*(t0+a1*(d+h))+a2*(d+h)*(d+h);
		t=sqrt((Fd+b1*h*h+sqrt(Fd2*Fd1))*0.5);
	}
	return t;
}

void rayEndpointError(float *x,float *p,float **s,float t)
/*< Output error message if ray get to the model side or bottom >*/
{
	/* TODO to correct the way you treat side rays */
	sf_warning("Ray endpoint => x=%f y=%f p[0]=%f p[1]=%f",x[1],x[0],p[0],p[1]);
	sf_warning("Ray starting point=> x=%f y=%f",s[1],s[0]);
	sf_warning("Ray traveltime => t=%f",t);
	sf_error("Bad ray angle, ray get to the model side/bottom");
}

float calculateTimeMisfit(float** s, /* NIP sources matrix (z,x) pairs */
			   float v0, /* Near surface velocity */
			   float* t0, /* Normal ray traveltime for each NIP source */
			   float* m0, /* Central CMP for each NIP source */
			   float* RNIP, /* RNIP parameter for each NIP source */
			   float* BETA, /* BETA parameter for each NIP source */
			   int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *slow, /* Slowness velociy model */
			   float *a, /* Normal ray angle for each NIP source (degrees) */
			   int ns, /* Number of NIP sources */
			   int itf, /* Interface being inverted */
			   float ***data,
			   int *data_n,
			   float *data_o,
			   float *data_d)
/*< Return L2 norm of the time misfit: The time misfit is the difference
between the traveltime calculated using raytracing and the traveltime calculated
with the CRE traveltime formula 

Values of x and p are changed inside the function.
The trajectory traj is stored as follows: {z0,y0,z1,y1,z2,y2,...} in 2-D

Note: This function traces nr reflection rays from each NIP source
(a depth point coordinate) to acquisition surface. NIP sources coordinates
are passed through s matrix. This function returns the L2 norm of the difference
between the traveltime of the reflection rays and the calculated traveltime using CRE
traveltime approximation. This difference is time misfit.

To simulate a reflection ray, this function traces a ray from the NIP source to the
source location in the acquisition surface and it stores its traveltime ts. And this
function traces a ray from the NIP source to the receiver location in the acquisition
surface and it stores the traveltime tr. The total reflection ray traveltime will be the
sum of t=ts+tr.
 >*/
{

	int is, it, i; // loop counter
	float p[2]; // slowness vector
	float t=0.; // Ray traveltime
	float normalRayAngleRad; // Normal ray angle in radians
	int nt=10000; // number of time samples in each ray
	float dt=0.001; // time sampling of rays
	raytrace rt; // raytrace struct
	float** traj; // Ray trajectory (z,x)
	float m; // CMP
	float h; // half-offset
	float tmis=0; // time misfit
	float xs=0.; // Source position
	float xr=0.; // Receiver position
	float tr=0.; // NIP to receiver ray traveltime
	float ts=0.; // NIP to source ray traveltime
	float *x; // Source position (z,x)
	float *nrnip; // Calculate normal ray rnips
	float *nbeta; // Calculate normal ray betas

	x = sf_floatalloc(2);
	nrnip = sf_floatalloc(ns);
	nbeta = sf_floatalloc(ns);

	for(is=(itf*ns);is<(itf*ns+ns);is++){

		x[0]=s[is][0];
		x[1]=s[is][1];
		
		normalRayAngleRad = a[is]*DEG2RAD;
		p[0] = -cosf(normalRayAngleRad);
		p[1] = sinf(normalRayAngleRad);

		/* initialize ray tracing object */
		rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);

		traj = sf_floatalloc2(2,nt+1);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		if(it>0){ // Ray endpoint at acquisition surface
			t = it*dt;
			ts=t;
			tr=t;
			xs=x[1];
			xr=x[1];
			i = it >= 2 ? it - 2 : it - 1;
                        /* Escape vector */
                        x[0]=traj[it][0];
                        x[1]=traj[it][1];
                        x[0]-=traj[i][0];
                        x[1]-=traj[i][1];
			/* Dot product with unit vector pointing upward */
                        t = sqrt(x[0]*x[0]+x[1]*x[1]); /* Length */
			t = acos(x[0]/t)-SF_PI;
			if(x[1]<0) t = -t;

			/* Calculate RNIP */
			nrnip[is]=sqrt((x[0]-s[is][0])*(x[0]-s[is][0])+(x[1]-s[is][1])*(x[1]-s[is][1]));	
			/* Keep BETA */
			nbeta[is]=t;
		}else if(it == 0){ // Ray endpoint inside model
			t = abs(nt)*dt;
			rayEndpointError(x,p,s,t);
		}else{ // Side or bottom ray
			rayEndpointError(x,p,s,t);
		}

		/* Raytrace close */
		raytrace_close(rt);
		free(traj);

		m = (xr+xs)/2.;
		h = (xr-xs)/2.;
		t = creTimeApproximation(h,m,v0,t0[is],m0[is],RNIP[is],BETA[is],false);
		tmis += fabs((ts+tr)-t);
		//sf_warning("rnip[%d]=%f %f",is,nrnip[is],RNIP[is]);
		//sf_warning("beta[%d]=%f %f",is,nbeta[is]*180/3.1415,BETA[is]*180/3.1415);
		//sf_warning("t0[%d]=%f %f",is,ts+tr,t0[is]);

	} /* Loop over NIP sources */

	/* L2 norm to evaluate the time misfit */
	tmis = sqrt(tmis*tmis);
	return tmis;
}

