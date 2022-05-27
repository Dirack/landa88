#include "Unity/unity.h"
#include "raytrace.h"
#include "forward.h"
#include <stdio.h>
#include <rsf.h>

int n[2]={301,1001};
float d[2]={0.01,0.01};
float o[2]={0.,-2.};
float v0=1.5;
float *slow;
float *slow2;

void init(){
	int i, j;
	float z;
	float ss=1./(v0*v0);
	slow = sf_floatalloc(n[0]*n[1]);

	for(i=0;i<n[0]*n[1];i++)
		slow[i]=ss;
}

void init2(){
	int i, j;
	float z;
	slow2 = sf_floatalloc(n[0]*n[1]);

	for(j=0; j<n[1]; j++){
		for(i=0; i<n[0]; i++){
			z=i*d[0];
			if(z<1.){
				slow2[i+j*n[0]] = 2.;
			}else{
				slow2[i+j*n[0]] = 2.5;
			}
		}
	}
	for(i=0;i<n[0]*n[1];i++){
		slow2[i]=1./(slow2[i]*slow2[i]);
	}
}

void setUp(){};

void tearDown(){};

void test_forwardModelingInConstantVelocityModel()
/*< Test model setup in a constant velocity model: NIP sources position and NIP angles >*/
{
	float **s;
	int ns=4;
	float m0[4]={0.,0.,0.,0.};
	float t0[4]={0.,0.,0.,0.};
	float a[4]={45.,-45,30,-30};
	float dist;
	int i;
	float *BETA;
	float *RNIP;

	s = sf_floatalloc2(2,ns);
	BETA = sf_floatalloc(ns);
	RNIP = sf_floatalloc(ns);

	for(i=0;i<ns;i++)
		BETA[i]=a[i];

	s[0][0]=1.;
	s[0][1]=3.;

	s[1][0]=2.;
	s[1][1]=3.5;

	s[2][0]=1.5;
	s[2][1]=4.5;

	s[3][0]=0.5;
	s[3][1]=2.5;

	forwardModeling(s,ns,m0,t0,a,n,d,o,slow,BETA,RNIP);

	/* Test for t0, m0 */
	for(i=0;i<ns;i++){
		dist = s[i][0]/cosf(a[i]*SF_PI/180.);
		TEST_ASSERT_FLOAT_WITHIN(0.01,2*dist/v0,t0[i]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,s[i][1]+dist*sinf(a[i]*SF_PI/180.),m0[i]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,a[i]*SF_PI/180.,BETA[i]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,2*dist,RNIP[i]);
	}

	free(BETA);
	free(RNIP);
	free(s);
}

void test_forwardModelingInTwoLayersModel()
/*< Test model setup in a constant velocity model: NIP sources position and NIP angles >*/
{
	float **s;
	int ns=4;
	float m0[4]={0.,0.,0.,0.};
	float t0[4]={0.,0.,0.,0.};
	float a[4]={0.,0.,0.,0.};
	float dist;
	int i;
	float *BETA;
	float *RNIP;

	s = sf_floatalloc2(2,ns);
	BETA = sf_floatalloc(ns);
	RNIP = sf_floatalloc(ns);

	for(i=0;i<ns;i++)
		BETA[i]=a[i];

	s[0][0]=1.5;
	s[0][1]=3.;

	s[1][0]=1.5;
	s[1][1]=3.5;

	s[2][0]=1.5;
	s[2][1]=4.5;

	s[3][0]=1.5;
	s[3][1]=2.5;

	forwardModeling(s,ns,m0,t0,a,n,d,o,slow2,BETA,RNIP);

	/* Test for t0, m0 */
	for(i=0;i<ns;i++){
		TEST_ASSERT_FLOAT_WITHIN(0.01,2*((1./2.)+(0.5/2.5)),t0[i]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,s[i][1],m0[i]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,a[i]*SF_PI/180.,BETA[i]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,2*1.5,RNIP[i]);
	}

	free(BETA);
	free(RNIP);
	free(s);
}


int main(void){

	init();
	init2();
        UNITY_BEGIN();
        RUN_TEST(test_forwardModelingInConstantVelocityModel);
        RUN_TEST(test_forwardModelingInTwoLayersModel);
        return UNITY_END();
}

