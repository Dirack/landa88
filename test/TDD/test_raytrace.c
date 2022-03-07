/*
* test_raytrace.c (C)
* 
* Purpose: Unit tests of raytrace.c lib.
* 
* Site: https://dirack.github.io
* 
* Version 1.0
* 
* Programmer: Rodolfo A C Neves 26/12/2021
* 
* Email: rodolfo_profissional@hotmail.com
* 
* License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.
*/

#include "Unity/unity.h"
#include "raytrace.h"
#include <stdio.h>
#include <rsf.h>

void setUp(){};

void tearDown(){};

void test_getTransmissionAngleSnellsLaw()
/*< Function to get transmission angle using snell's law >*/
{
	float ei=SF_PI/6., vt=1.5, vi=1.7;
	TEST_ASSERT_FLOAT_WITHIN(0.001,0.456,getTransmissionAngleSnellsLaw(ei,vi,vt));
}

void test_getTransmitedRNIPHubralTransmissionLaw()
/*< Test function to get Transmited RNIP parameter through interface using Hubral's transmission law >*/
{
	float rnip=1.0;
	float vt=1.5, vi=1.7;
	float ei=SF_PI/6.;
	itf2d it2;
	int i;
	float *sz;
	float kf;

	sz = sf_floatalloc(5);

	for(i=0;i<5;i++)
		sz[i]=1.0;

	it2 = itf2d_init(sz,5,0.,0.5);

	kf = calculateInterfaceCurvature(it2,1.2);
	getTransmitedRNIPHubralTransmissionLaw(&rnip,vt,vi,ei,kf);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.218,rnip);
}

void test_calculateIncidentAngle()
/*< Test function to get incident angle for a given ray sample coordinate in ray trajectory >*/
{
	float t[4]={1.22,1.05,1.2,1.0};
	float **traj;
	itf2d it2;
	float *sz;
	int i;

	sz = sf_floatalloc(5);
	for(i=0;i<5;i++)
		sz[i]=1.0;

	it2 = itf2d_init(sz,5,0.,0.5);

	traj = sf_floatalloc2(2,2);

	traj[0][0]=t[0];
	traj[0][1]=t[1];
	traj[1][0]=t[2];
	traj[1][1]=t[3];

	TEST_ASSERT_FLOAT_WITHIN(0.01,1.19,calculateIncidentAngle(it2,traj,1));
}

void test_getVelocityForRaySampleLocation()
/*< Test of the function to get velocity from grid for a ray sample location given >*/
{
	float *slow;
	int im;
	raytrace rt;
	int n[2]={10,10};
	float o[2]={0.,-2.};
	float d[2]={0.01,0.01};
	float t[2]={0.05,-1.95};
	float **traj;

        slow =  sf_floatalloc(100);
	traj = sf_floatalloc2(2,1);
	traj[0][0] = t[0];
	traj[0][1] = t[1];

	for(im=0;im<100;im++)
		slow[im] = 1./(1.5*1.5);

	rt = raytrace_init(2,true,10,1.0,n,o,d,slow,4);

	TEST_ASSERT_FLOAT_WITHIN(0.01,1.5,getVelocityForRaySampleLocation(rt,traj,0));
}

float cubicInterface(float x)
/*< Cubic function for tests >*/
{
	return x*x*x-3*x*x+4;
}

void test_firstDerivativeFunction()
/*< Test first derivative numerical calculation >*/
{
	float fxph;
	float fxmh;
	float x[4]={0.,1.,2.,3.};
	float dydx[4]={0.,-3.,0.,9.};
	int i;

	for(i=0;i<4;i++){
		fxph = cubicInterface(x[i]+0.001);
		fxmh = cubicInterface(x[i]-0.001);

		TEST_ASSERT_FLOAT_WITHIN(0.01,dydx[i],first_deriv(0.001,fxph,fxmh));
	}
}

void test_secondDerivativeFunction()
/*< Test second derivative numerical calculation >*/
{
	float fxm7h, fxm6h, fxm5h, fxm4h, fxm3h, fxm2h, fxmh, fx;
	float x[4]={0.,1.,2.,3.};
	float dydx[4]={6.,0.,6.,12.};
	int i;

	for(i=0;i<4;i++){
		fxm7h = cubicInterface(x[i]-7*0.01);
		fxm6h = cubicInterface(x[i]-6*0.01);
		fxm5h = cubicInterface(x[i]-5*0.01);
		fxm4h = cubicInterface(x[i]-4*0.01);
		fxm3h = cubicInterface(x[i]-3*0.01);
		fxm2h = cubicInterface(x[i]-2*0.01);
		fxmh = cubicInterface(x[i]-0.01);
		fx = cubicInterface(x[i]);

		//TEST_ASSERT_FLOAT_WITHIN(0.01,dydx[i],second_deriv(0.001,fxp2h,fxph,fx));
		printf("%f %f\n",dydx[i],second_deriv(0.01,fxm7h,fxm6h,fxm5h,fxm4h,fxm3h,fxm2h,fxmh,fx));
	}
}


void test_calculateInterfaceCurvature()
/*< TODO >*/
{
	float capa[4]={6.,0.,6.,0.016};
	itf2d it2;
	int i;
	float x[5]={-1,0,1,2,3};
	float y[10]={0.,3.125,4.0,3.375,2.,0.625,0.,0.875,4.,10.125};

	it2 = itf2d_init(y,10,-1,0.5);

	//for(i=0;i<11;i++)
	//	printf("%f ",cubicInterface(i*0.5-1.));
	for(i=0;i<4;i++)
		printf("k[%d]=%f\n",i,calculateInterfaceCurvature(it2,i));
		//TEST_ASSERT_FLOAT_WITHIN(0.01,capa[i],calculateInterfaceCurvature(it2,i));
}

int main(void){

	UNITY_BEGIN();
	/*RUN_TEST(test_getTransmissionAngleSnellsLaw);
	RUN_TEST(test_getTransmitedRNIPHubralTransmissionLaw);
	RUN_TEST(test_calculateIncidentAngle);
	RUN_TEST(test_getVelocityForRaySampleLocation);*/
	//RUN_TEST(test_firstDerivativeFunction);
	//RUN_TEST(test_secondDerivativeFunction);
	RUN_TEST(test_calculateInterfaceCurvature);
	return UNITY_END();
}
