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

int main(void){

	initialize_test_interface();

	UNITY_BEGIN();
	RUN_TEST(test_getZNodepoints);
	RUN_TEST(test_setZNodepoints);
	RUN_TEST(test_getNODFunctions);
	RUN_TEST(test_interfaceInterpolation);
	return UNITY_END();
}
