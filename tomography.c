/*
	 tomography.c (c)
	 
	 Purpose: 'Mvfsacrsnh.c' library for raytracing and traveltime
	 calculation.
	 	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 19/09/2019

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <rsf.h>
#include "model2d.h"
#include "raytrace.h"
#include "tomography.h"
/*^*/

#ifndef GDB_DEBUG
	#define DSLOW 0.04
	#define DANGLE 1.0
#else
	#define DSSLOW 0.04
	#define DANGLE 1.0
#endif
/*^*/

#define OFFSET_APERTURE 50
#define CMP_APERTURE 10
#define SEMBLANCE
#define CRE_TIME_CURVE
#define DT 0.001
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
If cds parameter is true, it uses the non-hyperbolic CRS formula with CDS condition (RN=RNIP) to calculate traveltime.
>*/
{ 
	float alpha; // Asymmetry parameter
	float d = m-m0; // Distance to central CMP m0
	float c1; // CRE coefficient
	float c2; // CRE coefficient
	float a1, a2, b2, b1, Fd, Fd1, Fd2; // Non-hyperbolic CRS coefficient
	float t; // traveltime t(m,h)

	if(cds){
		a1=(2*sin(BETA))/(v0);
		a2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RNIP);
		b2=(2*cos(BETA)*cos(BETA)*t0)/(v0*RNIP);
		b1=2*b2+a1*a1-a2;
		Fd=(t0+a1*d)*(t0+a1*d)+a2*d*d;
		Fd2=(t0+a1*(d-h))*(t0+a1*(d-h))+a2*(d-h)*(d-h);
		Fd1=(t0+a1*(d+h))*(t0+a1*(d+h))+a2*(d+h)*(d+h);
		t=sqrt((Fd+b1*h*h+sqrt(Fd2*Fd1))*0.5);
	}else{
		c1 = (d+h)/RNIP;
		c2 = (d-h)/RNIP;
		alpha = sin(BETA)/RNIP;
		t = (t0-2*RNIP/v0)+(RNIP/v0)*sqrt(1-2*alpha*(d+h)+c1*c1)+(RNIP/v0)*sqrt(1-2*alpha*(d-h)+c2*c2);
	}
	return t;
}

void rayEndpointError(float *x,float *p,float **traj,float t)
/*< Output error message if ray get to the model side or bottom >*/
{
	/* TODO to correct the way you treat side rays */
	sf_warning("Ray endpoint => x=%f y=%f p[0]=%f p[1]=%f",x[1],x[0],p[0],p[1]);
	sf_warning("Ray starting point=> x=%f y=%f",traj[0][1],traj[0][0]);
	sf_warning("Ray traveltime => t=%f",t);
	sf_error("Bad ray angle, ray get to the model side/bottom");
}

void setInitialRayPointAndRayVector(float **s, /* NIP sources Matrix */
				float *x, /* Initial ray point (z,x) */
				float *p, /* Initial ray vector */
				int is, /* Source index */
				float a /* Initial ray angle (radians) */)
/*< Set the initial ray point and ray vector >*/
{
		x[0]=s[is][0];
		x[1]=s[is][1];
		p[0] = -cosf(a);
		p[1] = sinf(a);
}

void calculateEscapeVector(
			   float *x, /* Ray endpoint */
			   float **traj, /* Ray trajectory */
			   int it /* Endpoint index */)
/*< Calculate Escape vector from ray trajectory, x is changed inside the function >*/
{
	int i;
	i = it >= 2 ? it - 2 : it - 1;
	x[0]=traj[it][0];
	x[1]=traj[it][1];
	x[0]-=traj[i][0];
	x[1]-=traj[i][1];
}

float calculateBetaWithRayTrajectory(
				     float *x, /* Ray endpoint */
				     float **traj /* Ray trajectory */,
				     int it /* Endpoint index */)
/*< Calculate BETA parameter using dot product with unit vector pointing upward
Note: x is changed inside the function
>*/
{
	float xx;

	calculateEscapeVector(x,traj,it);
	xx=sqrt(x[0]*x[0]+x[1]*x[1]);
	xx=acos(-x[0]/xx);
	if(x[1]<0) xx = -xx;
	return xx;
}

int stackOverCRETimeCurve(
			   float RNIP, /* RNIP parameter */
			   float BETA, /* BETA parameter */
			   float m0, /* Central CMP */
			   float t0, /* Normal ray traveltime */
			   float v0, /* Near surface velocity */
			   float *sumAmplitudes, /* Samples sum */
			   float *sumAmplitudes2, /* Samples sum squared */
			   float ***data, /* Seismic data cube */
			   int *n, /* Data number of samples (n1,n2,n3) */
			   float *o /* Data axis origins (o1,o2,o3) */,
			   float *d /* Data samplings (d1,d2,d3) */)
/*< Calculate CRE trajectory calculation and stack over CRE traveltime curve
Note: sumAmplitudes and sumAmplitudes2 variables are changed inside function
>*/
{
	float alpha; // Asymetry parameter
	int ih, im; // Loop counter
	float h; // Half-offset
	float m; // CMP
	int tetai; // Time sample index
	int numSamples=1; // Number of time samples to stack
	float sa=0.; // Samples sum
	float sa2=0.; // samples sum squared

	alpha = sinf(BETA)/RNIP;
	//sf_warning("%d",__LINE__);

	for(ih=0; ih < OFFSET_APERTURE; ih++){

		h = ih*d[1]+o[1];

		if(alpha <= 0.001 && alpha >= -0.001){
			m = m0;
		}else{
			m = m0 + (1/(2*alpha)) * (1 - sqrt(1 + 4 * alpha * alpha * h * h));
		}

	//sf_warning("%f",m);
		im = (int) (m/d[2]);

	//sf_warning("%d",__LINE__);
		tetai = (int) round((double) creTimeApproximation(h,m,v0,t0,m0,RNIP,BETA,false)/d[0]);

	//sf_warning("im=%d ih=%d tetai=%d",im,ih,tetai);
		if(tetai > n[0] || tetai < 0 || im < 0 || im > n[2]){
			sa += 0.;
		}else{
			sa += data[im][ih][tetai];
		}

	//sf_warning("%d",__LINE__);
		sa2 += (sa*sa);
		numSamples++;

	} /* loop over half-offset */

	*sumAmplitudes = sa;
	*sumAmplitudes2 = sa2;

	return numSamples;
}

int stackOverCRETimeSurface(
			   float RNIP, /* RNIP parameter */
			   float BETA, /* BETA parameter */
			   float m0, /* Central CMP */
			   float t0, /* Normal ray traveltime */
			   float v0, /* Near surface velocity */
			   float *sumAmplitudes, /* Samples sum */
			   float *sumAmplitudes2, /* Samples sum squared */
			   float ***data, /* Seismic data cube */
			   int *n, /* Data number of samples (n1,n2,n3) */
			   float *o /* Data axis origins (o1,o2,o3) */,
			   float *d /* Data samplings (d1,d2,d3) */)
/*< Calculate CRE trajectory calculation and stack over CRE traveltime surface
Note: sumAmplitudes and sumAmplitudes2 variables are changed inside function.
This function returns the number of samples in stack
>*/
{
	int ih, im; // Loop counter
	float h; // Half-offset
	float m; // CMP
	int tetai; // Time sample index
	int numSamples=1; // Number of time samples to stack
	float sa=0.; // Samples sum
	float sa2=0.; // samples sum squared
	int im0; // Loop counter

	im0 = (int) (m0/d[2]);

	for(im=im0-CMP_APERTURE;im<(im0+CMP_APERTURE+1);im++){

		m = im*d[2]+o[2];

		for(ih=0; ih < OFFSET_APERTURE; ih++){

			h = ih*d[1]+o[1];

			tetai = (int) ((double) creTimeApproximation(h,m,v0,t0,m0,RNIP,BETA,false)/d[0]);

			if(tetai > n[0] || tetai < 0){
				sa += 0.;
			}else{
				sa += data[im][ih][tetai];
			}

			sa2 += (sa*sa);
			numSamples++;

		} /* loop over half-offset */
	} /* loop over CMP */

	*sumAmplitudes = sa;
	*sumAmplitudes2 = sa2;

	return numSamples;
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
			   float *sz, /* Interface nodepoints */
			   int nsz, /* Number of nodepoints */
			   float osz, /* Interfaces axis origin */
			   float dsz, /* Nodepoints sampling */
			   float *vv, /* Layers velocity */
			   float nv /* Number of layers */)
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

	int is, it; // loop counter
	float p[2]; // slowness vector
	float t=0.; // Ray traveltime
	float normalRayAngleRad; // Normal ray angle in radians
	int nt=10000; // number of time samples in each ray
	float dt=0.001; // time sampling of rays
	raytrace rt; // raytrace struct
	float** traj; // Ray trajectory (z,x)
	float tmis=0; // time misfit
	float *x; // Source position (z,x)
	float sumAmplitudes=0.; // Amplitudes sum
	float sumAmplitudes2=0.; // Amplitudes sum squared
	int numSamples=1; // Number of samples
	float t0p;

	x = sf_floatalloc(2);

	/* initialize ray tracing object */
	rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
	traj = sf_floatalloc2(2,nt+1);

	for(is=(itf*ns); is<(itf*ns+ns); is++){
		t0p = t0[is];

		normalRayAngleRad = a[is]*DEG2RAD;

		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(s,x,p,is,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		if(it>0){ // Ray endpoint at acquisition surface

			m0[is]= x[1];
			t0[is]= 2*it*dt;

                        /* Calculate RNIP */
			//RNIP[is] = calculateRNIPWithHubralLaws(rt,traj,it,vv,nv,t0[is],itf,sz,nsz,osz,dsz);
			RNIP[is] = calculateRNIPWithDynamicRayTracing(rt,dt,it,traj,v0);

			if(RNIP[is]<0.0){
				sf_warning("ERROR: RNIP=%f",RNIP[is]);
				sf_warning("%f",a[is]);
				sf_warning("NIP=%d teta=%f",is,a[is]);
	                        sf_warning("From: x=%f z=%f",s[is][1],s[is][0]);
                        	sf_warning("To: x=%f z=%f",traj[it][1],traj[it][0]);
			}

			/* Calculate BETA */
			BETA[is] = calculateBetaWithRayTrajectory(x,traj,it);
			a[is] = BETA[is]*180./SF_PI;

			if(a[is] < -90.){
				sf_warning("ERROR: BETA=%f",BETA[is]);
				sf_warning("%f",a[is]);
				sf_warning("NIP=%d teta=%f",is,a[is]);
	                        sf_warning("From: x=%f z=%f",s[is][1],s[is][0]);
                        	sf_warning("To: x=%f z=%f",traj[it][1],traj[it][0]);

			}

			/* STACKING */
			//sumAmplitudes = 0.;
			//sumAmplitudes2 = 0.;
			#ifdef CRE_TIME_CURVE
			numSamples = stackOverCRETimeCurve(RNIP[is],BETA[is],m0[is],t0[is],v0,&sumAmplitudes,&sumAmplitudes2,data,data_n,data_o,data_d);
			if(sumAmplitudes2<0.00001 || t0[is]>t0p+0.4 || t0[is] < t0p-0.4 || RNIP[is] < 0.){
				tmis += 0.;
			}else{
				tmis += (sumAmplitudes*sumAmplitudes)/(numSamples*sumAmplitudes2);
			}
			//tmis += sumAmplitudes;
			#else
			numSamples = stackOverCRETimeSurface(RNIP[is],BETA[is],m0[is],t0[is],v0,&sumAmplitudes,&sumAmplitudes2,data,data_n,data_o,data_d);
			tmis += (sumAmplitudes*sumAmplitudes)/(numSamples*sumAmplitudes2);
			#endif

		}else if(it == 0){ // Ray endpoint inside model
			t = abs(nt)*dt;
			rayEndpointError(x,p,traj,t);
		}else{ // Side or bottom ray
			rayEndpointError(x,p,traj,t);
		}

	} /* Loop over NIP sources */

	free(traj);

	/* L2 norm to evaluate the time misfit */
	// TODO: Choose the best object function criteria
	#ifdef SEMBLANCE
	//tmis = (sumAmplitudes*sumAmplitudes)/(numSamples*sumAmplitudes2);
	tmis /= ns;
	#else
	tmis = sumAmplitudes;
	#endif
	return tmis;
}

void interfaceSetup(
		    float **s, /* NIP sources (z,x) */
		    int ns, /* NIP sources per interface */
		    float *m0, /* m0's for each NIP */
		    float *t0, /* t0's for each NIP */
		    float *BETA, /* BETA for each NIP */
		    int itf, /* Interface index */
		    int *n, /* Model dimension */
		    float *d, /* Model sampling */
		    float *o, /* Model axis origin */
		    float *slow /* Slowness model (1/v^2) */)
/*< Interface setup using NIP sources
NOTE: This function launches normal rays from acquisition surface into model to determine NIP sources
position stored in s matrix. These are used to draw model interfaces. The normal ray traveltime is
determined by t0 parameter and starting angle is BETA rotated PI radians (180 degrees). Starting
ray position is (x=m0,z=0) at acquisition surface.
>*/
{

	float x[2];
	float p[2];
	float t;
	raytrace rt;
	float **traj;
	int nt;
	int it;
	int i;
	int is;
	float a;

	for(is=(itf*ns); is<(itf*ns+ns); is++){

		/* initialize ray tracing object */
		nt = (int) (t0[is]/(2*DT));
		rt = raytrace_init(2,true,nt,DT,n,o,d,slow,ORDER);

		/* Ray tracing */
		traj = sf_floatalloc2(2,nt+1);
		
		/* initialize position */
		x[0] = 0.; 
		x[1] = m0[is];

		/* initialize direction */
		a=(BETA[is]+SF_PI);
		p[0] = -cosf(a);
		p[1] = sinf(a);

		it = trace_ray (rt, x, p, traj);

		/* write ray end points */
		s[is][0]=traj[nt-1][0];
		s[is][1]=traj[nt-1][1];

		/* Treat ray tracing errors */
		/* it=0, Ray stoped inside model
		   it<0, Side/Bottom rays
		   it>0, Top rays */
		if(it!=0){
                        sf_warning("BAD RAY ANGLE IN NIP MODEL SETUP");
			sf_warning("NIP=%d teta=%f it=%d",is,BETA[is]-180,it);
                        sf_warning("From: x=%f z=0.",m0[is]);
                        sf_warning("To: x=%f z=%f",s[is][1],s[is][0]);
                        sf_warning("Starting angle: %f",BETA[is]);
                        sf_warning("Escape angle: %f",BETA[is]);
			sf_error("%s: %d",__FILE__,__LINE__);
		}

		/* Free raytrace storage */
		raytrace_close(rt);
		free(traj);
	}

}

float forwardModeling(
			   float** s, /* NIP sources matrix (z,x) pairs */
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
			   float *sz, /* Interface nodepoints */
			   int nsz, /* Number of nodepoints */
			   float osz, /* Interfaces axis origin */
			   float dsz, /* Nodepoints sampling */
			   float *vv, /* Layers velocity */
			   float nv, /* Number of layers */
				mod2d mod)
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

	int is, it; // loop counter
	float p[2]; // slowness vector
	float t=0.; // Ray traveltime
	float normalRayAngleRad; // Normal ray angle in radians
	int nt=10000; // number of time samples in each ray
	float dt=0.001; // time sampling of rays
	raytrace rt; // raytrace struct
	float** traj; // Ray trajectory (z,x)
	float tmis=0; // time misfit
	float *x; // Source position (z,x)
	float sumAmplitudes=0.; // Amplitudes sum
	float sumAmplitudes2=0.; // Amplitudes sum squared
	int numSamples=1; // Number of samples
	float **nip;
	float beta;
	float rnip;

	x = sf_floatalloc(2);
	nip = sf_floatalloc2(2,1);

	/* initialize ray tracing object */
	rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
	traj = sf_floatalloc2(2,nt+1);

	for(is=0; is<10; is++){

		nip[0][1]=is*0.2+1.5;
		nip[0][0]=mod2d_getZCoordinateOfInterface(mod,itf,nip[0][1]);
		normalRayAngleRad = mod2d_getNipAngleForX(mod,itf,nip[0][1]);

		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(nip,x,p,0,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);
		//sf_warning("sz=%f",sz[is]);
		//sf_warning("x=%f z=%f",x[1],x[0]);

		if(it>0){ // Ray endpoint at acquisition surface

			//m0[is]= x[1];
			//t0p = t0[is+itf*ns];
			//sf_warning("itf=%d aqui=>%f",itf,2*it*dt);

                        /* Calculate RNIP */
			//RNIP[is] = calculateRNIPWithHubralLaws(rt,traj,it,vv,nv,t0[is],itf,sz,nsz,osz,dsz);
			rnip = calculateRNIPWithDynamicRayTracing(rt,dt,it,traj,v0);

			/*if(RNIP[is]<0.0){
				sf_warning("ERROR: RNIP=%f",RNIP[is]);
				sf_warning("%f",a[is]);
				sf_warning("NIP=%d teta=%f",is,a[is]);
	                        sf_warning("From: x=%f z=%f",s[is][1],s[is][0]);
                        	sf_warning("To: x=%f z=%f",traj[it][1],traj[it][0]);
			}*/

			/* Calculate BETA */
			beta = calculateBetaWithRayTrajectory(x,traj,it);
			//a[is] = BETA[is]*180./SF_PI;

			/*if(a[is] < -90.){
				sf_warning("ERROR: BETA=%f",BETA[is]);
				sf_warning("%f",a[is]);
				sf_warning("NIP=%d teta=%f",is,a[is]);
	                        sf_warning("From: x=%f z=%f",s[is][1],s[is][0]);
                        	sf_warning("To: x=%f z=%f",traj[it][1],traj[it][0]);

			}*/

			/* STACKING */
			sumAmplitudes = 0.;
			sumAmplitudes2 = 0.;
			//sf_warning("oi");
			numSamples = stackOverCRETimeCurve(rnip,beta,x[1],2*it*dt,v0,&sumAmplitudes,&sumAmplitudes2,data,data_n,data_o,data_d);
			//sf_warning("capa");
			if(sumAmplitudes2<0.00001 || rnip < 0. ){
				tmis += 0.;
			}else{
				tmis += (sumAmplitudes*sumAmplitudes)/(numSamples*sumAmplitudes2);
			}
			
		}else if(it == 0){ // Ray endpoint inside model
			t = abs(nt)*dt;
			rayEndpointError(x,p,traj,t);
		}else{ // Side or bottom ray
			//rayEndpointError(x,p,traj,t);
		}

	} /* Loop over NIP sources */

	free(traj);
	free(x);
	free(nip);

	/* L2 norm to evaluate the time misfit */
	// TODO: Choose the best object function criteria
	#ifdef SEMBLANCE
	//tmis = (sumAmplitudes*sumAmplitudes)/(numSamples*sumAmplitudes2);
	tmis /= ns;
	#else
	tmis = sumAmplitudes;
	#endif
	return tmis;
}

