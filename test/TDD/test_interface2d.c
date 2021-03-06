/*
* test_interface2d.c (C)
* 
* Purpose: Unit tests of interface2d struct.
* 
* Site: https://dirack.github.io
* 
* Version 1.0
* 
* Programmer: Rodolfo A C Neves 25/12/2021
* 
* Email: rodolfo_profissional@hotmail.com
* 
* License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.
*/

#include "Unity/unity.h"
#include "interface2d.h"
#include <stdio.h>
#include <rsf.h>

itf2d itf;
float *sz;

void setUp(){};

void tearDown(){};

void initialize_test_interface(){
/*< Function to initialize interface struct to be used in tests >*/

	int i;

	sz = sf_floatalloc(5);
	for(i=0;i<5;i++)
		sz[i]=1.0;

	itf = itf2d_init(sz,5,-2.0,1.0);
}

void test_getZNodepoints(){
/*< Test if the nodepoints are corrected initialized >*/
	// All nodes are 1.0
	int i;
	float *s;
	s = sf_floatalloc(5);

	itf2d_getZNodepoints(itf,s);
	TEST_ASSERT_EQUAL_FLOAT_ARRAY(s,sz,5);	
}

void test_setZNodepoints(){
/*< Test nodepoints setting z(x) >*/

	int i;
	float *s1, *s2;
	itf2d it2;

	s1 = sf_floatalloc(5);
	s2 = sf_floatalloc(5);

	for(i=0;i<5;i++)
		s1[i] = 2.5;

	it2 = itf2d_init(sz,5,-2.0,1.0);
	itf2d_setZNodepoints(itf,s1);
	itf2d_getZNodepoints(itf,s2);
	TEST_ASSERT_EQUAL_FLOAT_ARRAY(s2,s1,5);	
}


void test_getNODFunctions(){
/*< Test functions to get n, o, d from interface struct >*/
	TEST_ASSERT_EQUAL(itf2d_n(itf),5);
	TEST_ASSERT_FLOAT_WITHIN(0.01,itf2d_o(itf),-2.0);
	TEST_ASSERT_FLOAT_WITHIN(0.01,itf2d_d(itf),1.0);
}

void test_interfaceInterpolation(){
/*< Test cubic spline interpolation of interfaces nodepoints with calculated values >*/

	itf2d s0;
	float z[7]={2.0,4.0,2.77,1.0,1.65,3.0,3.0};

	s0 = itf2d_init(z,7,1.0,1.0);
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.55,getZCoordinateOfInterface(s0,1.2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.99,getZCoordinateOfInterface(s0,2.9));
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.95,getZCoordinateOfInterface(s0,5.2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,3.1,getZCoordinateOfInterface(s0,6.7));
}

void test_cubicSplineCoefficientsCalculation(){
/*< Test cubic spline coefficients calculation function >*/

	float x[5]={1.0,2.0,4.0,6.0,7.0};
	float z[5]={2.0,4.0,1.0,3.0,3.0};
	float ss[4][4] = {{-0.783,0.,2.783,2.},{0.692,-2.35,0.433,4},{-0.483,1.8,-0.666,1.},{0.366,-1.1,0.733,3.}};
	float *coef;
	int n=5;
	int i, j;

	coef = sf_floatalloc(4*(n-1));
	calculateSplineCoeficients(5,x,z,coef);
	
	for(i=0;i<4;i++){
		TEST_ASSERT_FLOAT_WITHIN(0.01,coef[0+i*4],ss[i][0]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,coef[1+i*4],ss[i][1]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,coef[2+i*4],ss[i][2]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,coef[3+i*4],ss[i][3]);
	}
}

void test_cubicSplineCoefficientsCalculationPassedAsVector(){
/*< Test cubic spline coefficients calculation function passing coef matrix as vector component >*/

	float x[5]={1.0,2.0,4.0,6.0,7.0};
	float z[5]={2.0,4.0,1.0,3.0,3.0};
	float ss[4][4] = {{-0.783,0.,2.783,2.},{0.692,-2.35,0.433,4},{-0.483,1.8,-0.666,1.},{0.366,-1.1,0.733,3.}};
	float **coef;
	int n=5;
	int i, j;

	coef = sf_floatalloc2(4*(n-1),2);
	for(i=0;i<2;i++){
                calculateSplineCoeficients(n,x,z,coef[i]);
        }

	for(j=0;j<2;j++){
		for(i=0;i<4;i++){
			TEST_ASSERT_FLOAT_WITHIN(0.01,coef[j][0+i*4],ss[i][0]);
			TEST_ASSERT_FLOAT_WITHIN(0.01,coef[j][1+i*4],ss[i][1]);
			TEST_ASSERT_FLOAT_WITHIN(0.01,coef[j][2+i*4],ss[i][2]);
			TEST_ASSERT_FLOAT_WITHIN(0.01,coef[j][3+i*4],ss[i][3]);
		}
	}
}

void test_calcInterfacesZcoord(){
/*< Test interfaces Z(x) calculation >*/

	float x[5]={1.0,2.0,4.0,6.0,7.0};
	float z1[5]={1.0,1.0,1.0,1.0,1.0};
	float z2[5]={1.85,1.85,1.85,1.85,1.85};
	float **coef;
	float *zi;
	int n=5;
	int i, j;

	coef = sf_floatalloc2(4*(n-1),2);
	zi = sf_floatalloc(2);

	calculateSplineCoeficients(n,x,z1,coef[0]);
	calculateSplineCoeficients(n,x,z2,coef[1]);

	// Spline 1
	calcInterfacesZcoord(zi,2,1.5-1.0,0,coef);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.0,zi[0]);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.85,zi[1]);

	// Spline 2
	calcInterfacesZcoord(zi,2,2.5-2.0,1,coef);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.0,zi[0]);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.85,zi[1]);

	// Spline 3
	calcInterfacesZcoord(zi,2,4.5-4.0,2,coef);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.0,zi[0]);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.85,zi[1]);

	// Spline 4
	calcInterfacesZcoord(zi,2,6.5-6.0,3,coef);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.0,zi[0]);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.85,zi[1]);
}

int main(void){

	initialize_test_interface();

	UNITY_BEGIN();
	RUN_TEST(test_getZNodepoints);
	RUN_TEST(test_setZNodepoints);
	RUN_TEST(test_getNODFunctions);
	RUN_TEST(test_interfaceInterpolation);
	RUN_TEST(test_cubicSplineCoefficientsCalculation);
	RUN_TEST(test_cubicSplineCoefficientsCalculationPassedAsVector);
	RUN_TEST(test_calcInterfacesZcoord);
	return UNITY_END();
}
