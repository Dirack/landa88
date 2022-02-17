/*
* test_layer2d.c (C)
* 
* Purpose: Unit tests of layer2d struct.
* 
* Site: https://dirack.github.io
* 
* Version 1.0
* 
* Programmer: Rodolfo A C Neves 16/02/2022
* 
* Email: rodolfo_profissional@hotmail.com
* 
* License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.
*/

#include "Unity/unity.h"
#include "interface2d.h"
#include "layer2d.h"
#include <stdio.h>
#include <rsf.h>

lay2d l2;

void setUp(){};

void tearDown(){};

void initialize_layer2d_for_tests()
/*< Initialize layer2d struct for tests >*/
{
	l2 = lay2d_init(1.5,1.8,1.2);
}

void test_layer2dInitialization()
/*< Test layer2d struct initialization >*/
{
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.5,lay2d_getvel(l2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.8,lay2d_getvmax(l2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.2,lay2d_getvmin(l2));
}

void test_setLayer2dVelocities()
/*< Test layer2d velocities settings >*/
{
	lay2d_setvel(l2,2.5);
	lay2d_setvmax(l2,2.8);
	lay2d_setvmin(l2,2.2);

	TEST_ASSERT_FLOAT_WITHIN(0.01,2.5,lay2d_getvel(l2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.8,lay2d_getvmax(l2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.2,lay2d_getvmin(l2));
}

int main(void){

	initialize_layer2d_for_tests();

	UNITY_BEGIN();
	RUN_TEST(test_layer2dInitialization);
	RUN_TEST(test_setLayer2dVelocities);
	return UNITY_END();
}
