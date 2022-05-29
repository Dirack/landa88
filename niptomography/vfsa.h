#include "forward.h"
#include "setup.h"
#include "datamis.h"
#define signal(s) ((s<0)?(-1.):(1.))
/*< Signal function >*/

void sortinginAscendingOrder(
                                float *x, /* x vector to sort */
                                int n /* Vectors dimension */)
/*< x vector sorting in ascending order >*/
{
        int i; // Loop counter
        float tmpx; // Temporary variables
        int k; // Sorting key (number of changes)

        do{
                k=0;
                for(i=1;i<n;i++){
                        if(x[i-1]>x[i]){
                                tmpx=x[i-1];
                                x[i-1]=x[i];
                                x[i]=tmpx;
                                k++;
                        }
                } // Loop vector samples
        }while(k!=0);
}

float getRandomNumberBetween0and1(){
/*< Function to get a random number between 0 and 1 >*/

        return (float)(rand()%1000)/1000;
}

float vfsaOptimization(float **s, int ns, float *m0, float *t0, float *a, int *n, float *o, float *d, float *slow,
float *BETA, float *RNIP, float v0, float ***data, int *data_n, float *data_o, float * data_d, float t0p, float lambda, float vm)
{
	float mis;
	float mis0;
	int numSamples;
	float temp;
	int q;
	float u;
	float disturbance;
	float minvel=vm-lambda, maxvel=vm+lambda;
	float vnew=1.;
	float Em0, deltaE, PM;
	float sumAmplitudes=0., sumAmplitudes2=0.;
	int i, j;
	float votm=vm+lambda;
	float v=vm+lambda;
	float semb;
	float semb0;
	float vv[2]={v,2.0};

	for(j=0;j<n[0]*n[1];j++)
		slow[j]=1./(v*v);

	nipModelSetup(s,ns,m0,t0,a,n,d,o,slow);

	forwardModeling(s,ns,m0,t0,a,n,d,o,slow,BETA,RNIP,vv,2,1);

	numSamples = stackOverCRETimeCurve(RNIP[0],BETA[0],m0[0],t0[0],v0,&sumAmplitudes,&sumAmplitudes2,data,data_n,data_o,data_d);
        mis = fabs(RNIP[0]-2.);
	mis0 = mis;
	semb0 = 0.;

        /* Very Fast Simulated Annealing (VFSA) algorithm */
        for (q=0; q<100; q++){

	       /* calculate VFSA temperature for this iteration */
                temp=1.5*expf(-0.75*pow(q,0.33));

                /* parameter disturbance */
		u=getRandomNumberBetween0and1();

		disturbance = signal(u - 0.5) * temp * (pow( (1+temp),fabs(2*u-1) )-1);

		vnew = v + (disturbance);

		if (vnew >= maxvel)
			vnew = maxvel - (maxvel-minvel) * getRandomNumberBetween0and1();

		if (vnew <= minvel)
			vnew = (maxvel-minvel) * getRandomNumberBetween0and1() + minvel;

		for(j=0;j<n[0]*n[1];j++)
                        slow[j]=1./(vnew*vnew);

		vv[0]=vnew;

                mis=0;

		nipModelSetup(s,ns,m0,t0,a,n,d,o,slow);

		forwardModeling(s,ns,m0,t0,a,n,d,o,slow,BETA,RNIP,vv,2,1);

		numSamples = stackOverCRETimeCurve(RNIP[0],BETA[0],m0[0],t0[0],v0,&sumAmplitudes,&sumAmplitudes2,data,data_n,data_o,data_d);
		semb = (sumAmplitudes*sumAmplitudes)/(numSamples*sumAmplitudes2);
		mis = fabs(RNIP[0]-2.)*fabs(RNIP[0]-2.);

                if(fabs(mis) < fabs(mis0) && semb > semb0){
                        mis0 = fabs(mis);
			votm = vnew;
			semb0 = semb;
                }

                /* VFSA parameters update condition */
                deltaE = fabs(mis) - Em0;

                /* Metrópolis criteria */
                PM = expf(-deltaE/temp);

              	if (deltaE<=0){
			v = vnew;
                        Em0 = fabs(mis);
                } else {
                        u=getRandomNumberBetween0and1();
                        if (PM > u){
				v = vnew;
                                Em0 = fabs(mis);
                        }
                }

	}

	return votm;
}

float vfsaOptimizationTwoLayersModel(float **s, int ns, float *m0, float *t0, float *a, int *n, float *o, float *d, float *slow,
float *BETA, float *RNIP, float v0, float ***data, int *data_n, float *data_o, float * data_d, float t0p, float lambda, float vm)
{
	float mis;
	float mis0;
	int numSamples;
	float temp;
	int q;
	float u;
	float disturbance;
	float minvel=vm-lambda, maxvel=vm+lambda;
	float vnew=1.;
	float Em0, deltaE, PM;
	float sumAmplitudes=0., sumAmplitudes2=0.;
	int i, j;
	float votm=vm+lambda;
	float v=vm+lambda;
	float semb;
	float semb0;
	float vv[2]={1.528,v};
	float z;

	for(j=0; j<n[1]; j++){
		for(i=0; i<n[0]; i++){
			z=i*d[0];
			if(z<1.){
				slow[i+j*n[0]] = 1.528;
			}else{
				slow[i+j*n[0]] = v;
			}
		}
	}
	for(i=0;i<n[0]*n[1];i++){
		slow[i]=1./(slow[i]*slow[i]);
	}


	nipModelSetup(s,ns,m0,t0,a,n,d,o,slow);

	forwardModeling(s,ns,m0,t0,a,n,d,o,slow,BETA,RNIP,vv,2,1);

	numSamples = stackOverCRETimeCurve(RNIP[0],BETA[0],m0[0],t0[0],v0,&sumAmplitudes,&sumAmplitudes2,data,data_n,data_o,data_d);
        mis = fabs(RNIP[0]-4.2546);
	mis0 = mis;
	semb0 = 0.;

        /* Very Fast Simulated Annealing (VFSA) algorithm */
        for (q=0; q<100; q++){

	       /* calculate VFSA temperature for this iteration */
                temp=1.5*expf(-0.75*pow(q,0.33));

                /* parameter disturbance */
		u=getRandomNumberBetween0and1();

		disturbance = signal(u - 0.5) * temp * (pow( (1+temp),fabs(2*u-1) )-1);

		vnew = v + (disturbance);

		if (vnew >= maxvel)
			vnew = maxvel - (maxvel-minvel) * getRandomNumberBetween0and1();

		if (vnew <= minvel)
			vnew = (maxvel-minvel) * getRandomNumberBetween0and1() + minvel;

		for(j=0; j<n[1]; j++){
			for(i=0; i<n[0]; i++){
				z=i*d[0];
				if(z<1.){
					slow[i+j*n[0]] = 1.528;
				}else{
					slow[i+j*n[0]] = vnew;
				}
			}
		}
		for(i=0;i<n[0]*n[1];i++){
			slow[i]=1./(slow[i]*slow[i]);
		}

		vv[1] = vnew;
                mis=0;

		nipModelSetup(s,ns,m0,t0,a,n,d,o,slow);

		forwardModeling(s,ns,m0,t0,a,n,d,o,slow,BETA,RNIP,vv,2,1);

		numSamples = stackOverCRETimeCurve(RNIP[0],BETA[0],m0[0],t0[0],v0,&sumAmplitudes,&sumAmplitudes2,data,data_n,data_o,data_d);
		semb = (sumAmplitudes*sumAmplitudes)/(numSamples*sumAmplitudes2);
		if(sumAmplitudes2<0.0001) semb=0.;
		mis = fabs(RNIP[0]-4.2546)*fabs(RNIP[0]-4.2546);

                if(fabs(mis) < fabs(mis0) && semb > semb0){
                        mis0 = fabs(mis);
			votm = vnew;
			semb0 = semb;
                }

                /* VFSA parameters update condition */
                deltaE = fabs(mis) - Em0;

                /* Metrópolis criteria */
                PM = expf(-deltaE/temp);

              	if (deltaE<=0){
			v = vnew;
                        Em0 = fabs(mis);
                } else {
                        u=getRandomNumberBetween0and1();
                        if (PM > u){
				v = vnew;
                                Em0 = fabs(mis);
                        }
                }

	}

	return votm;
}
