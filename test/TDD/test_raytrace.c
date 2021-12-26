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

int main(void){

	UNITY_BEGIN();
	RUN_TEST(test_snellslaw);
	return UNITY_END();
}
