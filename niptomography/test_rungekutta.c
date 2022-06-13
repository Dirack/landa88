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
float *slow2;

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

void init_modelx2y2()
/*< Velocity model: v(z,x)=z^2+x^2
This model is used to test partial derivative calculation in the direction normal to ray
trajectory. The velocity gradient will be grad(v)=(2z,2x) and the second derivative in
z and x direction will be equal 2.
>*/
{
	int i, j;
	float z,x;

	slow2 = sf_floatalloc(n[0]*n[1]);

	for(j=0; j<n[1]; j++){
		x=j*d[1]+o[1];
		for(i=0; i<n[0]; i++){
			z=i*d[0]+o[0];
			slow2[i+j*n[0]]=z*z+x*x;
		}
	}
	// Get slowness s=1/v^2
	for(i=0;i<n[0]*n[1];i++){
		slow2[i]=1./(slow2[i]*slow2[i]);
		//printf("s=%f\n",slow2[i]);
	}

	//sf_error("oi");
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
	float a[1]={31.,};
        int nt=10000;
        float dt=0.001;
        int it;
        float t;
        int is;
	float **s;
	float rnip=0.;
	int i;

	//TEST_IGNORE_MESSAGE("TODO");
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

	TEST_IGNORE_MESSAGE("TODO");

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

void test_getVectorFromTwoPoints()
/*< Given two points return a vector from their difference.
Example: Given (0,0) and (1,1) return (1,1), because (1,1) = (1-0,1-0)
>*/
{
	float x1[2];
	float x2[2];
	float v[2];

	// Return (1,1) for (0,0) and (1,1)
	x1[0]=0; x1[1]=0; x2[0]=1; x2[1]=1;
	getVectorFromTwoPoints(x1,x2,v);
	TEST_ASSERT_FLOAT_WITHIN(0.1,1.,v[0]);
	TEST_ASSERT_FLOAT_WITHIN(0.1,1.,v[1]);

	// Return (1,1) for (2,2) and (1,1)
	x1[0]=1; x1[1]=1; x2[0]=2; x2[1]=2;
	getVectorFromTwoPoints(x1,x2,v);
	TEST_ASSERT_FLOAT_WITHIN(0.1,1.,v[0]);
	TEST_ASSERT_FLOAT_WITHIN(0.1,1.,v[1]);
}

void test_getUnitVectorForAVector()
/*< For a vector given, calculate the unit vector in the vector direction
Example: Given (2,2) return (1/sqrt(2),1/sqrt(2)) or (0.7,0.7)
>*/
{
	float v[2];
	float ve[2];
	float no[2];

	// Return (0.7,0.7) for (2,2)
	v[0]=2; v[1]=2; ve[0]=0.7; ve[1]=0.7;
	getUnitVectorForAVector(v,no);
	TEST_ASSERT_FLOAT_WITHIN(0.1,ve[0],no[0]);
	TEST_ASSERT_FLOAT_WITHIN(0.1,ve[1],no[1]);
}

void test_rotateAVector90DegreesRight()
/*< Rotate a vector 90 degrees to right to obtain a normal vector
Example: Given (1,1) return (0,1)
>*/
{
	float v[2], ve[2], no[2];

	// Return (0,1) for (1,0)
	v[0]=1; v[1]=0; ve[0]=0; ve[1]=1;
	rotateAVector90DegreesRight(v,no);
	TEST_ASSERT_EQUAL_FLOAT_ARRAY(ve,no,2);

	// Return (-1,1) for (1,1)
	v[0]=1; v[1]=1; ve[0]=-1; ve[1]=1;
	rotateAVector90DegreesRight(v,no);
	TEST_ASSERT_EQUAL_FLOAT_ARRAY(ve,no,2);
}

void test_rotateAVector90DegreesLeft()
/*< Rotate a vector 90 degrees to right to obtain a normal vector
Example: Given (1,1) return (0,1)
>*/
{
	float v[2], ve[2], no[2];

	// Return (1,0) for (0,-1)
	v[0]=1; v[1]=0; ve[0]=0; ve[1]=-1;
	rotateAVector90DegreesLeft(v,no);
	TEST_ASSERT_EQUAL_FLOAT_ARRAY(ve,no,2);

	// Return (1,-1) for (1,1)
	v[0]=1; v[1]=1; ve[0]=1; ve[1]=-1;
	rotateAVector90DegreesLeft(v,no);
	TEST_ASSERT_EQUAL_FLOAT_ARRAY(ve,no,2);
}

void test_secondVelocityPartialDerivativeNormalToRayDirection()
/*<
Get velocity partial derivative in the direction normal to ray trajectory. This
test can be done using a known velocity field and its partial derivatives in x and y.
Example: If v(y,x)=x^2+y^2, d^2v/dx^2 and d^2v/dy^2 are equal 2. So, to vertical and
horizontal rays trajectories derivatives are known and equal 2.
>*/
{
	raytrace rt;
	float x[2];
	float no[2];
	float v;
	float der;

        rt = raytrace_init(2,true,nt,dt,n,o,d,slow2,ORDER);

	// Get d^2v/dx^2 (should be equal 2)
	x[0]=3.; x[1]=3.; no[0]=1; no[1]=0;
	v=x[0]*x[0]+x[1]*x[1];
	der=secondVelocityPartialDerivativeNormalToRayDirection(rt,no,x,v);
	TEST_ASSERT_FLOAT_WITHIN(0.1,2.,der);
	x[0]=3.; x[1]=3.; no[0]=-1; no[1]=0;
	v=x[0]*x[0]+x[1]*x[1];
	der=secondVelocityPartialDerivativeNormalToRayDirection(rt,no,x,v);
	TEST_ASSERT_FLOAT_WITHIN(0.1,2.,der);

	// Get d^2v/dy^2 (should be equal 2)
	x[0]=3.; x[1]=3.; no[0]=0; no[1]=1;
	v=x[0]*x[0]+x[1]*x[1];
	der=secondVelocityPartialDerivativeNormalToRayDirection(rt,no,x,v);
	TEST_ASSERT_FLOAT_WITHIN(0.1,2.,der);
	x[0]=3.; x[1]=3.; no[0]=0; no[1]=-1;
	v=x[0]*x[0]+x[1]*x[1];
	der=secondVelocityPartialDerivativeNormalToRayDirection(rt,no,x,v);
	TEST_ASSERT_FLOAT_WITHIN(0.1,2.,der);


	raytrace_close(rt);
}

void test_secondVelocityPartialDerivativeNormalToRayDirectionAllRayTrajectory()
/*<
*** This tests is done in all ray trajectory ***
Get velocity partial derivative in the direction normal to ray trajectory. This
test can be done using a known velocity field and its partial derivatives in x and y.
Example: If v(y,x)=x^2+y^2, d^2v/dx^2 and d^2v/dy^2 are equal 2. So, to vertical and
horizontal rays trajectories derivatives are known and equal 2.

!IMPORTANT: Velocity field bends ray, so derivative is not exact!
>*/
{
	raytrace rt;
	float x[2];
	float p[2];
	int i, it;
	float *dvdn;
	float **traj;
	float normalRayAngleRad;
	int nt=10000;
	float dt=0.001;

        rt = raytrace_init(2,true,nt,dt,n,o,d,slow2,ORDER);
        traj = sf_floatalloc2(2,nt+1);

	// Vertical ray from bottom to top
	x[0]=2.; x[1]=3.;
	normalRayAngleRad = 0.;
	p[0] = -1;//-cosf(normalRayAngleRad);
	p[1] = 0.;//sinf(normalRayAngleRad);
	
	/* Ray tracing */
	it = trace_ray (rt, x, p, traj);
	//sf_warning("%d",it);

	dvdn = sf_floatalloc(it);

	secondVelocityPartialDerivativeNormalToRayDirectionAllRayTrajectory(rt,traj,dvdn,it);
	for(i=0;i<it;i++){
		TEST_ASSERT_FLOAT_WITHIN(0.15,2.,dvdn[i]);
	}

	// Horizontal ray from left to right
	x[0]=2.; x[1]=3.;
	normalRayAngleRad = 0.;
	p[0] = 0;//-cosf(normalRayAngleRad);
	p[1] = 1;//sinf(normalRayAngleRad);
	
	/* Ray tracing */
	it = trace_ray (rt, x, p, traj);

	dvdn = sf_floatalloc(it);

	secondVelocityPartialDerivativeNormalToRayDirectionAllRayTrajectory(rt,traj,dvdn,it);
	for(i=0;i<it;i++){
		TEST_ASSERT_FLOAT_WITHIN(0.16,2.,dvdn[i]);
	}

	free(traj);
	free(dvdn);
	raytrace_close(rt);
}


int main(int argc, char* argv[]){

	init();
	init_modelx2y2();
        UNITY_BEGIN();
        RUN_TEST(test_getRNIPUsingDynamicRayTracingInConstantVelocityModel);
        RUN_TEST(test_getRNIPUsingDynamicRayTracingInTwoLayersVelocityModel);
        RUN_TEST(test_getRNIPUsingDynamicRayTracingInTwoLayersVelocityModelForBendingRay);
	RUN_TEST(test_getVectorFromTwoPoints);
	RUN_TEST(test_getUnitVectorForAVector);
	RUN_TEST(test_rotateAVector90DegreesRight);
	RUN_TEST(test_rotateAVector90DegreesLeft);
	RUN_TEST(test_secondVelocityPartialDerivativeNormalToRayDirection);
	RUN_TEST(test_secondVelocityPartialDerivativeNormalToRayDirectionAllRayTrajectory);
        return UNITY_END();
}

