#include "Unity/unity.h"
#include "raytrace.h"
#include "forward.h"
#include "datamis.h"
#include <time.h>
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
FILE *SEMB;
FILE *RNIP;
FILE *VEL;

void init(float v){
       int i, j;
        float z;

        for(j=0; j<n[1]; j++){
                for(i=0; i<n[0]; i++){
                        z=i*d[0];
                        if(z<1.){
                                slow[i+j*n[0]] = 1.508;
                        }else{
                                slow[i+j*n[0]] = v;
                        }
                }
        }
        for(i=0;i<n[0]*n[1];i++){
                slow[i]=1./(slow[i]*slow[i]);
        }

}

void setUp(){};

void tearDown(){};

void test_semblanceTwoLayersModel()
/*< Test model setup in a constant velocity model: NIP sources position and NIP angles >*/
{
	float **s;
        int ns=1;
        float m0[1]={0.};
        float t0[1]={2.183};
        float a[1]={0.};
	float med[9];
	float average=0.;;
        float BETA;
        float rnip;
        float v;
	int nv=1;
	float ov=1.508;
	float dv=0.01;
	int i, j;
	float lambda=1.;
	float vm=1.9;
	float semb=0.;
	float sumAmplitudes, sumAmplitudes2;
	int numSamples=1;

        s = sf_floatalloc2(1,ns);

	s[0][0]=1.85;
	s[0][1]=3.;
	m0[0]=s[0][1];
	t0[0]=1.33;
	BETA=0.;

	for(i=0;i<nv;i++){
		v=dv*i+ov;
		//init(v);
		//rnip=2*0.85*v/1.508+1.+1;
		rnip=2*v*t0[0];
		sumAmplitudes=0., sumAmplitudes2=0.;
                numSamples = stackOverCRETimeCurve(rnip,BETA,m0[0],t0[0],v0,&sumAmplitudes,&sumAmplitudes2,data,data_n,data_o,data_d);

		semb = (sumAmplitudes*sumAmplitudes)/(0.004*numSamples*sumAmplitudes2);
		fprintf(VEL,"%f\n",v);
		fprintf(SEMB,"%f\n",semb);
		fprintf(RNIP,"%f\n",rnip);
	}

}


int main(int argc, char* argv[]){

        /* Redirect the stdin to datacube file */
        //freopen("data/interpolatedDataCube2.rsf","r",stdin);
        freopen("data/teste.rsf","r",stdin);
	VEL = fopen("vel.txt","w");
	SEMB = fopen("semb.txt","w");
	RNIP = fopen("rnip.txt","w");

        sf_init(argc,argv);
        in = sf_input("in");

        /* Read seismic data cube */
        data=sf_floatalloc3(data_n[0],data_n[1],data_n[2]);
        sf_floatread(data[0][0],data_n[0]*data_n[1]*data_n[2],in);

        slow = sf_floatalloc(n[0]*n[1]);
	srand(time(NULL));
	init(1.508);
        //UNITY_BEGIN();
        //RUN_TEST(test_vfsaOptimization);
	//RUN_TEST(test_semblanceTwoLayersModel);
	test_semblanceTwoLayersModel();
	fclose(VEL);
	fclose(SEMB);
	fclose(RNIP);
        //return UNITY_END();
}
