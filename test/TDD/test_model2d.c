/*
* test_model2d.c (C)
* 
* Purpose: Unit tests of model2d struct.
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
#include "model2d.h"
#include <stdio.h>
#include <rsf.h>

mod2d m2;

void setUp(){};

void tearDown(){};

void test_mod2dInitialization()
/*< Test model2d initialization >*/
{
	int nlay=4;
	float sv[4]={1.5,1.65,1.7,2.0};
	float minsv[4]={1.2,1.6,1.69,1.8};
	float maxsv[4]={1.6,1.69,2.0,2.5};
	int nsz=33;
	float osz=-2.;
	float dsz=1.;
	float sz[33]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.45,1.45,1.45,1.45,1.45,1.45,1.45,1.45,1.45,1.45,1.45,1.85,1.85,1.85,1.85,1.85,1.85,1.85,1.85,1.85,1.85,1.85,};

	m2 = mod2d_init(nlay,sv,minsv,maxsv,nsz,osz,dsz,sz);

	TEST_ASSERT_EQUAL(4,mod2d_getnumlayers(m2));
	TEST_ASSERT_EQUAL(11,mod2d_getnuminterfacesnodes(m2));

	TEST_ASSERT_FLOAT_WITHIN(0.01,1.2,mod2d_getlayervmin(m2,0));
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.6,mod2d_getlayervmin(m2,1));
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.69,mod2d_getlayervmin(m2,2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.8,mod2d_getlayervmin(m2,3));

	TEST_ASSERT_FLOAT_WITHIN(0.01,1.6,mod2d_getlayervmax(m2,0));
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.69,mod2d_getlayervmax(m2,1));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.0,mod2d_getlayervmax(m2,2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.5,mod2d_getlayervmax(m2,3));

	TEST_ASSERT_FLOAT_WITHIN(0.01,1.5,mod2d_getlayervel(m2,0));
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.65,mod2d_getlayervel(m2,1));
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.7,mod2d_getlayervel(m2,2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.0,mod2d_getlayervel(m2,3));

}

void test_setModel2dVelocities()
/*< Test model2d velocities settings >*/
{
	mod2d_setlayervel(m2,0,2.5);
	mod2d_setlayervel(m2,1,2.65);
	mod2d_setlayervel(m2,2,2.7);
	mod2d_setlayervel(m2,3,3.0);

	TEST_ASSERT_FLOAT_WITHIN(0.01,2.5,mod2d_getlayervel(m2,0));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.65,mod2d_getlayervel(m2,1));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.7,mod2d_getlayervel(m2,2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,3.0,mod2d_getlayervel(m2,3));

	mod2d_setlayervmin(m2,0,2.3);
	mod2d_setlayervmin(m2,1,2.5);
	mod2d_setlayervmin(m2,2,2.6);
	mod2d_setlayervmin(m2,3,2.9);

	TEST_ASSERT_FLOAT_WITHIN(0.01,2.3,mod2d_getlayervmin(m2,0));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.5,mod2d_getlayervmin(m2,1));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.6,mod2d_getlayervmin(m2,2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.9,mod2d_getlayervmin(m2,3));

	mod2d_setlayervmax(m2,0,2.7);
	mod2d_setlayervmax(m2,1,2.8);
	mod2d_setlayervmax(m2,2,2.9);
	mod2d_setlayervmax(m2,3,3.0);

	TEST_ASSERT_FLOAT_WITHIN(0.01,2.7,mod2d_getlayervmax(m2,0));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.8,mod2d_getlayervmax(m2,1));
	TEST_ASSERT_FLOAT_WITHIN(0.01,2.9,mod2d_getlayervmax(m2,2));
	TEST_ASSERT_FLOAT_WITHIN(0.01,3.0,mod2d_getlayervmax(m2,3));
}

void test_getInterfacesNodesFunction()
/*< Test get interfaces nodes function from model2d >*/
{
	int i;
	float z[11];

	mod2d_getinterfacesnodes(m2,0,z);

	for(i=0;i<11;i++)
		TEST_ASSERT_FLOAT_WITHIN(0.01,1.0,z[i]);

	mod2d_getinterfacesnodes(m2,1,z);

	for(i=0;i<11;i++)
		TEST_ASSERT_FLOAT_WITHIN(0.01,1.45,z[i]);

	mod2d_getinterfacesnodes(m2,2,z);

	for(i=0;i<11;i++)
		TEST_ASSERT_FLOAT_WITHIN(0.01,1.85,z[i]);
}

void test_getAllInterfacesNodesFunction()
/*< Test get all interfaces nodes function from model2d >*/
{
	float z[33];
	int i, j;
	float zz[3]={1.0,1.45,1.85};

	mod2d_getallinterfacesnodes(m2,z);

	for(i=0;i<3;i++){
		for(j=0;j<11;j++){
			TEST_ASSERT_FLOAT_WITHIN(0.01,zz[i],z[i*11+j]);
		}
	}
}

int main(void){

	UNITY_BEGIN();
	RUN_TEST(test_mod2dInitialization);
	RUN_TEST(test_setModel2dVelocities);
	RUN_TEST(test_getInterfacesNodesFunction);
	RUN_TEST(test_getAllInterfacesNodesFunction);
	return UNITY_END();
}
