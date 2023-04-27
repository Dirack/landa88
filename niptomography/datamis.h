#include "grid2.h"
#include "atela.h"
#include "dynamic.h"
#include "interface2d.h"
#define DT 0.001
#define OFFSET_APERTURE 101

#ifndef __DATAMIS_H__
#define __DATAMIS_H__

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
		if(alpha <= 0.001 && alpha >= -0.001){
			t = (t0-2*RNIP/v0)+(2*RNIP/v0)*sqrt(1+h*h/(RNIP*RNIP));
		}else{
			t = (t0-2*RNIP/v0)+(RNIP/v0)*sqrt(1-2*alpha*(d+h)+c1*c1)+(RNIP/v0)*sqrt(1-2*alpha*(d-h)+c2*c2);
		}

	}
	return t;
}

float nmo(float t0, float v0, float d, float h, float BETA, float RNIP){
	float t;
	t = (t0*t0+(2*sin(BETA)*d/v0))*(t0*t0+(2*sin(BETA)*d/v0))+(2*t0*cos(BETA)*cos(BETA)/(RNIP*v0))*(d*d+h*h);
	return sqrt(t);
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

	for(ih=0; ih < OFFSET_APERTURE; ih++){

		h = ih*d[1]+o[1]; h/=2;

		if(alpha <= 0.001 && alpha >= -0.001){
			m = m0;
		}else{
			m = m0 + (1/(2*alpha)) * (1 - sqrt(1 + 4 * alpha * alpha * h * h));
		}

		im = (int) round((m-o[2])/d[2]);

		tetai = (int) round((double) creTimeApproximation(h,m,v0,t0,m0,RNIP,BETA,false)/d[0]);
		printf("%f\n",creTimeApproximation(h,m,v0,t0,m0,RNIP,BETA,false));
		//tetai = (int) round((double) (nmo(t0,v0,m-m0,h,BETA,RNIP)-o[0])/d[0]);
		//printf("%f\n",nmo(t0,v0,m-m0,h,BETA,RNIP));

		if(tetai > n[0] || tetai < 0){
			sa += 0.;
		}else{
			sa += data[im][ih][tetai];
		}

		sa2 += (sa*sa);
		numSamples++;

	} /* loop over half-offset */

	*sumAmplitudes = sa;
	*sumAmplitudes2 = sa2;

	return numSamples-1;
}

#endif
