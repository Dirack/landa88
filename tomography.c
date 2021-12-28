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
			   float ***data, /* Seismic data cube A(m,h,t) */
			   int *data_n, /* Data number of samples */
			   float *data_o, /* Data axis origin */
			   float *data_d, /* Data sampling */
			   float *sz,
			   int nsz,
			   float osz,
			   float dsz)
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
	float *x; // Source position (z,x)
	float *nrnip; // Calculate normal ray rnips
	float *nbeta; // Calculate normal ray betas
	float sumAmplitudes=0.; // Amplitudes sum
	float sumAmplitudes2=0.; // Amplitudes sum squared
	float cm0; // central CMP m0
	float ct0; // normal ray traveltime t0
	float alpha; // Asymmetry paramter
	int ih; // half-offset index
	int im; // CMP index
	int tetai; // time index
	int numSamples=0; // samples counter
	float *vv;

	vv = sf_floatalloc(3);
	vv[0]=1.5; vv[1]=1.7; vv[2]=2.0;

	x = sf_floatalloc(2);
	nrnip = sf_floatalloc(ns);
	nbeta = sf_floatalloc(ns);


	/* initialize ray tracing object */
	rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
	traj = sf_floatalloc2(2,nt+1);

	for(is=(itf*ns);is<(itf*ns+ns);is++){

		x[0]=s[is][0];
		x[1]=s[is][1];
		
		normalRayAngleRad = a[is]*DEG2RAD;
		p[0] = -cosf(normalRayAngleRad);
		p[1] = sinf(normalRayAngleRad);


		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);
		if(it>0){ // Ray endpoint at acquisition surface

			cm0 = x[1];
			ct0 = 2*it*dt;
                        x[0]=traj[it][0];
                        x[1]=traj[it][1];

                        /* Calculate RNIP */
			nrnip[is-(itf*ns)] = calculateRNIPWithHubralLaws(rt,traj,it,vv,ct0,itf,sz,nsz,osz,dsz,ns);
			//nrnip[is-(itf*ns)] = calculateRNIPWithDynamicRayTracing(rt,dt,it,traj,v0);

			//sf_warning("rnipc=%f RNIP=%f",v0*ct0,nrnip[is-(itf*ns)]);
			if(nrnip[is-(itf*ns)]<0.0) sf_warning("rnipc=%f RNIP=%f",v0*ct0,nrnip[is-(itf*ns)]);

			/* Escape vector */ 
			i = it >= 2 ? it - 2 : it - 1;
                        x[0]-=traj[i][0];
                        x[1]-=traj[i][1];

			/* Calculate BETA */
			/* Dot product with unit vector pointing upward */
                        t = sqrt(x[0]*x[0]+x[1]*x[1]); /* Length */
			t = acos(x[0]/t)-SF_PI; /* Teta */
			if(x[1]<0) t = -t;

			nbeta[is-(itf*ns)]=t;

			alpha = sinf(nbeta[is-(itf*ns)])/nrnip[is-(itf*ns)];

			/* CRE STACKING */
			sumAmplitudes = 0;

			/* CRE trajectory calculation 
			   Semblance over CRE traveltime curve */
			for(ih=0; ih < 25; ih++){

				h = ih*data_d[1]+data_o[1];

				if(alpha <= 0.001 && alpha >= -0.001){
                                        m = cm0;
                        	}else{
                                        m = cm0 + (1/(2*alpha)) * (1 - sqrt(1 + 4 * alpha * alpha * h * h));
                                }

				im = (int) (m/data_d[2]);

				tetai = (int) ((double) creTimeApproximation(h,m,v0,ct0,cm0,nrnip[is-(itf*ns)],nbeta[is-(itf*ns)],false)/data_d[0]);

				if(tetai > data_n[0] || tetai < 0){
					sumAmplitudes += 0.;
				}else{
					sumAmplitudes += data[im][ih][tetai];
				}

				sumAmplitudes2 += (sumAmplitudes*sumAmplitudes);
				numSamples++;

			} /* loop over h*/

		}else if(it == 0){ // Ray endpoint inside model
			t = abs(nt)*dt;
			rayEndpointError(x,p,s,t);
		}else{ // Side or bottom ray
			rayEndpointError(x,p,s,t);
		}

	} /* Loop over NIP sources */

	raytrace_close(rt);
	free(traj);

	/* L2 norm to evaluate the time misfit */
	// TODO: Choose the best object function criteria
	//tmis = (sumAmplitudes*sumAmplitudes)/(numSamples*sumAmplitudes2);
	tmis = sumAmplitudes;
	return tmis;
}

