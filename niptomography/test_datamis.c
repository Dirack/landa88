#include "Unity/unity.h"
#include "raytrace.h"
#include "forward.h"
#include "datamis.h"
#include <stdio.h>
#include <rsf.h>

int n[2]={301,1001};
float d[2]={0.01,0.01};
float o[2]={0.,-2.};
float v0=1.508;
float *slow;
float ***data;
int data_n[3]={1001,161,964};
float data_o[3]={0.,0.,0.};
float data_d[3]={0.004,0.025,0.00625};
sf_file in; // datacube RSF file

void init(){
	int i;
	float ss=1./(v0*v0);
	slow = sf_floatalloc(n[0]*n[1]);

	for(i=0;i<n[0]*n[1];i++)
		slow[i]=ss;
}

void setUp(){};

void tearDown(){};

void test_dataMisfitInConstantVelocityModel()
/*< Test model setup in a constant velocity model: NIP sources position and NIP angles >*/
{
	float sumAmplitudes=0., sumAmplitudes2=0.;
	int numSamples;
	float tmis;
	float **s;
	int ns=2;
	float m0[2]={0.,0.};
	float t0[2]={1.3262,1.3262};
	float a[2]={0.,0.};
	float dist;
	int i,j,k;
	float *BETA;
	float *RNIP;
	float v;
	float dv=0.1;
	float ov=1.0;
	float diff=10.;
	float otm=ov;

	s = sf_floatalloc2(2,ns);
	BETA = sf_floatalloc(ns);
	RNIP = sf_floatalloc(ns);

	for(i=0;i<ns;i++)
		BETA[i]=a[i];

	s[0][0]=1.;
	s[0][1]=5.;

	s[1][0]=1.;
	s[1][1]=3.;

	
	for(i=0;i<10;i++){
		v=i*dv+ov;
		for(j=0;j<n[0]*n[1];j++)
			slow[j]=1./(v*v);
		tmis=0.;
		sumAmplitudes=0., sumAmplitudes2=0.;
		forwardModeling(s,1,m0,t0,a,n,d,o,slow,BETA,RNIP);
		numSamples = stackOverCRETimeCurve(RNIP[0],BETA[0],m0[0],t0[0],v0,&sumAmplitudes,&sumAmplitudes2,data,data_n,data_o,data_d);
		tmis += (sumAmplitudes*sumAmplitudes)/(numSamples*sumAmplitudes2);
		if(fabs(t0[0]-1.3262)<diff){
			diff = fabs(t0[0]-1.3262);
			otm = v;
		}
	}

	TEST_ASSERT_FLOAT_WITHIN(0.01,v0,otm);

	free(BETA);
	free(RNIP);
	free(s);

}

int main(int argc, char* argv[]){

        /* Redirect the stdin to datacube file */
        freopen("data/interpolatedDataCube2.rsf","r",stdin);

        sf_init(argc,argv);
        in = sf_input("in");

        /* Read seismic data cube */
        data=sf_floatalloc3(data_n[0],data_n[1],data_n[2]);
        sf_floatread(data[0][0],data_n[0]*data_n[1]*data_n[2],in);

	init();
        UNITY_BEGIN();
        RUN_TEST(test_dataMisfitInConstantVelocityModel);
        return UNITY_END();
}

