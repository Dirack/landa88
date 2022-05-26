#include "Unity/unity.h"
#include "raytrace.h"
#include "setup.h"
#include <stdio.h>
#include <rsf.h>

int n[2]={301,1001};
float d[2]={0.01,0.01};
float o[2]={0.,-2.};
float v0=1.5;
float *slow;

void init(){
	int i;
	float ss=1./(v0*v0);
	slow = sf_floatalloc(n[0]*n[1]);

	for(i=0;i<n[0]*n[1];i++)
		slow[i]=ss;
}

void setUp(){};

void tearDown(){};

void test_nipModelSetupInConstantVelocityModel()
/*< Test model setup in a constant velocity model: NIP sources position and NIP angles >*/
{
	float **s;
	int ns=4;
	float m0[4]={5.,5.,5.,5.};
	float t0[4]={1.,1.,2.,2.};
	float a[4]={45.,-45.,30.,-30.};
	float dist;
	int i;
	float *teta;

	s = sf_floatalloc2(2,ns);
	teta = sf_floatalloc(ns);

	for(i=0;i<ns;i++)
		teta[i]=a[i];

	nipModelSetup(s,ns,m0,t0,a,n,d,o,slow);

	/* Test NIP angles */
	for(i=0;i<ns;i++)
		TEST_ASSERT_FLOAT_WITHIN(0.02,teta[i],a[i]);

	/* Test NIP sources position */
	for(i=0;i<ns;i++){
		dist = v0 * t0[i]/2.;
		TEST_ASSERT_FLOAT_WITHIN(0.01,dist*cosf(a[i]*SF_PI/180.),s[i][0]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,m0[i]-dist*sinf(a[i]*SF_PI/180.),s[i][1]);
	}

	free(s);
	free(teta);
}

int main(void){

	init();
        UNITY_BEGIN();
        RUN_TEST(test_nipModelSetupInConstantVelocityModel);
        return UNITY_END();
}

