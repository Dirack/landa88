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

void test_snellslaw(){
	float ei=SF_PI/6., vt=1.5, vi=1.7;
	TEST_ASSERT_FLOAT_WITHIN(0.001,0.456,snellslaw(ei,vi,vt));
}

void test_hubralTransmissionLaw(){
	float rnip=1.0;
	float vt=1.5, vi=1.7;
	float ei=SF_PI/6.;

	hubralTransmissionLaw(&rnip,vt,vi,ei);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.218,rnip);
}

void test_calculateeiAngle(){
	float t[4]={1.22,1.05,1.2,1.0};
	float ei;
	float **traj;

	traj = sf_floatalloc2(2,2);

	traj[0][0]=t[0];
	traj[0][1]=t[1];
	traj[1][0]=t[2];
	traj[1][1]=t[3];

	ei = calculateeiAngle(traj,1);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.19,ei);
}

int main(void){

	UNITY_BEGIN();
	RUN_TEST(test_snellslaw);
	RUN_TEST(test_hubralTransmissionLaw);
	RUN_TEST(test_calculateeiAngle);
	return UNITY_END();
}
