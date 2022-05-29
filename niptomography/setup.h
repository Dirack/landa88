
#include "raytrace.h"
#define DT 0.001

void nipModelSetup(float **s,int ns, float *m0, float *t0, float* a, int *n, float *d, float *o, float *slow){

	float x[2];
	float p[2];
        float t;
        raytrace rt;
        float **traj;
        int nt;
        int it;
        int i;
        int ir;

        for(ir=0; ir<ns; ir++){

                /* initialize ray tracing object */
                nt = (int) round(t0[ir]/(2*DT));
                rt = raytrace_init(2,true,nt,DT,n,o,d,slow,ORDER);

                /* Ray tracing */
                traj = sf_floatalloc2(2,nt+1);

                /* initialize position */
                x[0] = 0.;
                x[1] = m0[ir];

                /* initialize direction */
                a[ir]=(a[ir]+180)*DEG2RAD;
                p[0] = -cosf(a[ir]);
                p[1] = sinf(a[ir]);
              	it = trace_ray (rt, x, p, traj);

                /* write ray end points */
                s[ir][0]=traj[nt-1][0];
                s[ir][1]=traj[nt-1][1];

                /* write escape angles */
                if(it!=0){
                	sf_warning("BAD RAY ANGLE IN NIP MODEL SETUP");
			sf_warning("From: x=0. z=%f",m0[ir]);
			sf_warning("To: x=%f z=%f",s[ir][0],s[ir][1]);
			sf_warning("Starting angle: %f",a[ir]*180./SF_PI-180);
			sf_warning("Escape angle: %f",a[ir]*180./SF_PI);
                }else{
                        /* Escape vector */
                        it=nt-1;
                        i=it-2;
                        x[0]=traj[it][0];
                        x[1]=traj[it][1];
                        x[0]-=traj[i][0];
                        x[1]-=traj[i][1];

                        /* Dot product with unit vector pointing upward */
                        t = sqrt(x[0]*x[0]+x[1]*x[1]); /* Length */
                        t = acos(x[0]/t);
                        if(x[1]>0) t = -t;

                        a[ir] = t*180./SF_PI;
                }

                /* Raytrace close */
                raytrace_close(rt);
                free(traj);
        }
}
