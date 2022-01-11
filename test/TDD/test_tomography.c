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
#include "tomography.h"
#include <stdbool.h>
#include <stdio.h>
#include <rsf.h>

void setUp(){};

void tearDown(){};

void test_setInitialRayPointAndRayVector()
/*< Test function to set up ray initial point and initial slowness vector >*/
{
	float x[2];
	float a=SF_PI/6.;
	float p[2];
	float **s;
	float p0=-cosf(SF_PI/6), p1=sinf(SF_PI/6);

	s = sf_floatalloc2(2,1);
	s[0][0]=1.2;
	s[0][1]=1.5;

	setInitialRayPointAndRayVector(s,x,p,0,a);
	TEST_ASSERT_FLOAT_WITHIN(0.01,x[0],s[0][0]);
	TEST_ASSERT_FLOAT_WITHIN(0.01,x[1],s[0][1]);
	TEST_ASSERT_FLOAT_WITHIN(0.01,p[0],p0);
	TEST_ASSERT_FLOAT_WITHIN(0.01,p[1],p1);
}

void test_calculateEscapeVector()
/*< Test function to calculate escape vector from ray trajectory >*/
{
	float **traj;
	float x[2];

	traj = sf_floatalloc2(2,2);
	traj[0][0] = 1.0;
	traj[0][1] = 1.0;
	traj[1][0] = 2.0;
	traj[1][1] = 2.0;

	calculateEscapeVector(x,traj,1);

	TEST_ASSERT_FLOAT_WITHIN(0.01,x[0],1.0);
	TEST_ASSERT_FLOAT_WITHIN(0.01,x[1],1.0);
}

void test_calculateBetaWithRayTrajectory()
/*< Test function to calculate BETA parameter from ray trajectory >*/
{
	float **traj;
	float x[2]={0.0,0.0};

	traj = sf_floatalloc2(3,3);
	traj[0][0] = 2.0;
	traj[0][1] = 2.0;
	traj[1][0] = 1.0;
	traj[1][1] = 1.0;
	traj[2][0] = 0.0;
	traj[2][1] = 0.0;

	TEST_ASSERT_FLOAT_WITHIN(0.01,SF_PI/4.,calculateBetaWithRayTrajectory(x,traj,2));
}


int main(void){

	UNITY_BEGIN();
	RUN_TEST(test_setInitialRayPointAndRayVector);
	RUN_TEST(test_calculateEscapeVector);
	RUN_TEST(test_calculateBetaWithRayTrajectory);
	return UNITY_END();
}
