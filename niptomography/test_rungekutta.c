#include "Unity/unity.h"
#include "raytrace.h"
#include "forward.h"
#include "rungekutta.h"
#include <stdio.h>
#include <rsf.h>

int n[2]={301,1001};
float d[2]={0.01,0.01};
float o[2]={0.,-2.};
float v0=1.5;
float *slow;

void init(){
	int i, j;
	float z;
	float ss=1./(v0*v0);
	slow = sf_floatalloc(n[0]*n[1]);

	for(j=0; j<n[1]; j++){
		for(i=0; i<n[0]; i++){
			z=i*d[0];
			if(z<1.){
				slow[i+j*n[0]] = 1.508;
			}else{
				slow[i+j*n[0]] = 2.0;
			}
		}
	}
	for(i=0;i<n[0]*n[1];i++){
		slow[i]=1./(slow[i]*slow[i]);
	}

}

void setUp(){}

void tearDown(){}

void test_getRNIPUsingDynamicRayTracingInConstantVelocityModel(){
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        float normalRayAngleRad;
	float a[1]={0.,};
        int nt=10000;
        float dt=0.001;
        int it;
        float t;
        int is;
	float **s;
	float rnip=0.;
	int i;

	s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
	normalRayAngleRad = a[0]*DEG2RAD;

	for(i=0;i<5;i++){
		s[0][0] = 1.;
		s[0][1] = i*0.5+2.5;
		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		rnip = calculateRNIPWithDynamicRayTracing(rt,dt,it,traj,v0);

		TEST_ASSERT_FLOAT_WITHIN(0.1,2.0,rnip);
	}

	raytrace_close(rt);
        free(traj);

}

void test_getRNIPUsingDynamicRayTracingInTwoLayersVelocityModel(){
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        float normalRayAngleRad;
	float a[1]={0.,};
        int nt=10000;
        float dt=0.001;
        int it;
        float t;
        int is;
	float **s;
	float rnip=0.;
	int i;

	s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
	normalRayAngleRad = a[0]*DEG2RAD;

	for(i=0;i<5;i++){
		s[0][0] = 1.85;
		s[0][1] = i*0.5+2.5;
		/* Set initial ray point and ray vector */
		setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

		/* Ray tracing */
		it = trace_ray (rt, x, p, traj);

		//printf("it=%d\n",it);
		rnip = calculateRNIPWithDynamicRayTracing(rt,dt,it,traj,v0);

		TEST_ASSERT_FLOAT_WITHIN(0.1,4.25,rnip);
	}
	raytrace_close(rt);
        free(traj);

}

void test_getRNIPUsingDynamicRayTracingInTwoLayersVelocityModelForBendingRay(){
        float x[2];
        float p[2];
        raytrace rt;
        float **traj;
        float normalRayAngleRad;
	float a[1]={20,};
        int nt=10000;
        float dt=0.001;
        int it;
        float t;
        int is;
	float **s;
	float rnip=0.;
	int i;
	float tetal=SF_PI/2.-a[0]*DEG2RAD;
	float ei, et;
	float vi=2.0, vt=1.508;
	float rnipi, rnipt, rnip_hubral;

	ei = SF_PI/2. - tetal;
	et = asinf((vt/vi)*sinf(ei));

	rnipi = 2.*(0.85/sinf(tetal));

	rnipt = ((vi/vt)*((cosf(et)*cosf(et))/(cosf(ei)*cosf(ei)))*rnipi);

	rnip_hubral = fabs(2*1./cosf(et) + rnipt);
	printf("tetal=%f ei=%f et=%f ri=%f rt=%f\n",tetal*180/SF_PI,ei*180/SF_PI,et*180/SF_PI,rnipi,rnipt);
	printf("rnip_hubral=%f\n",rnip_hubral);

	s = sf_floatalloc2(2,1);

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);
	normalRayAngleRad = a[0]*DEG2RAD;

	s[0][0] = 1.85;
	s[0][1] = 5.;

	/* Set initial ray point and ray vector */
	setInitialRayPointAndRayVector(s,x,p,0,normalRayAngleRad);

	/* Ray tracing */
	it = trace_ray (rt, x, p, traj);

	printf("it=%d traj=%f\n",it,traj[it][0]);
	rnip = calculateRNIPWithDynamicRayTracing(rt,dt,it,traj,v0);
	printf("rnip=%f\n",rnip);

	TEST_ASSERT_FLOAT_WITHIN(0.3,rnip_hubral,rnip);

	raytrace_close(rt);
        free(traj);

}


int main(void){

	init();
        UNITY_BEGIN();
        //RUN_TEST(test_getRNIPUsingDynamicRayTracingInConstantVelocityModel);
        //RUN_TEST(test_getRNIPUsingDynamicRayTracingInTwoLayersVelocityModel);
        RUN_TEST(test_getRNIPUsingDynamicRayTracingInTwoLayersVelocityModelForBendingRay);
        return UNITY_END();
}

