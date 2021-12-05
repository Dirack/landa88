/* Runge-Kutta ODE solvers for dynamic ray tracing system in 2D. */

#include <rsf.h>
#include <stdlib.h>
#include <math.h>

#include "dynamic.h"

static int dim, nt;
static float dt, **k, *yk;

void sf_dynamic_runge_init(int dim1 /* dimensionality */, 
		   	   int n1   /* number of ray tracing steps */, 
		   	   float d1 /* step in time */)
/*< initialize >*/
{
    dim = dim1;
    nt = n1;
    dt = d1;
    
    yk = sf_floatalloc(dim);
    k = sf_floatalloc2(dim,4);
}

void sf_dynamic_runge_close(void)
/*< free allocated storage >*/
{
    free(yk);
    free(*k);
    free(k);
}

float sf_dynamic_runge_step (float* y    /* [dim] solution */, 
		   	   void* par   /* parameters for function evaluation */,
		   	   void (*rhs)(void*,float,float*,float*,float*) 
		   	   /* RHS function */, 
			   float *dvdn,
		   	   float** traj /* [nt+1][dim] - ray trajectory (output) */) 
/*< ODE solver for dy/dt = f where f comes from rhs(par,y,f)
  Note:
  >*/
{
    int it, i;
	float *x;
	float rnip;

	x = sf_floatalloc(2);
 
    for (it = 0; it < nt; it++) {

	x[0]=traj[it][0];
	x[1]=traj[it][1];
	rhs (par, dvdn[it],x, y, k[0]); 

	for (i=0; i < dim; i++) {
	    yk[i] = y[i] + 0.5*dt*k[0][i];
	}
      
	x[0]=traj[it+1][0];
	x[1]=traj[it+1][1];
	rhs (par, dvdn[it+1],x,yk, k[1]); 

	for (i=0; i < dim; i++) {
	    yk[i] = y[i] + 0.5*dt*k[1][i];
	}
      
	x[0]=traj[it+1][0];
	x[1]=traj[it+1][1];
	rhs (par, dvdn[it+1],x,yk, k[2]); 
	for (i=0; i < dim; i++) {
	    yk[i] = y[i] + dt*k[2][i];
	}

	x[0]=traj[it+2][0];
	x[1]=traj[it+2][1];
	rhs (par, dvdn[it+2],x,yk, k[3]); 

	for (i=0; i < dim; i++) {
	    y[i] += dt*(k[0][i]+2.*k[1][i]+2.*k[2][i]+k[3][i])/6.0;
	}
	
    }
  
	rnip = 1.5*(y[1]/y[0]);
	rnip =  1./rnip;
	return rnip;
}

