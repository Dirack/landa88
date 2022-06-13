#include "forward.h"
#include "grid2.h"
#include "atela.h"
#include "dynamic.h"
#include "interface2d.h"
#define DT 0.001

#ifndef __RUNGEKUTTA_H__
#define __RUNGEKUTTA_H__

/* Runge-Kutta */
static int dim; /* Dimension - 2D */
static int nt; /* Number of time samples in a ray */
static float dt; /* Ray time sampling */
static float **k; /* Runge-Kutta coefficient K1, k2, k3, k4 */
static float *yk; /* Runge-kutta right hand side for dy/dt=f */

void sf_dynamic_runge_init(int dim1 /* dimensionality */,
                           int n1   /* number of ray tracing steps */,
                           float d1 /* step in time */)
/*< initialize Runge-Kutta >*/
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
                           void (*rhs)(void*,float,float*,float*,float*) /* RHS function */,
                           float *dvdn,
                           /* Second derivative of velocity in normal direction */
                           float** traj, /* [nt+1][dim] - ray trajectory (output) */
                           float v0 /* Near surface velosity */)
/*< ODE solver for dy/dt = f where f comes from rhs(par,y,f)

Note: It calculates p and q functions using dynamic ray tracing system
and it returns RNIP=v0(p(u0)/q(u0)). The derivative dvdn is evaluated in
each point of a given ray in normal direction of the ray trajectory
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

        rnip = v0*(y[1]/y[0]);
        rnip =  1./rnip;
        return rnip;
}

static void dyn_iso_rhs(void* par, float dvdn, float *x, float* y, float* f)
/* right-hand side for isotropic raytracing */
{
        raytrace rt;
        float v;

        rt = (raytrace) par;

        v = sqrtf(1./grid2_vel(rt->grd2,x));

        f[0]   = v*v*y[1]; /* v^2 p */
        f[1] = (-1./v)*dvdn*y[0]; /* -1/v dv/dn q */
}

void getVectorFromTwoPoints(float *x1, float *x2, float *v)
/*< TODO >*/
{
	v[0]=x2[0]-x1[0];
	v[1]=x2[1]-x1[1];
}

void getUnitVectorForAVector(float *v, float *n)
/*< TODO >*/
{
	float mod;
	mod = sqrt(v[0]*v[0]+v[1]*v[1]);
	n[0] = v[0]/mod;
	n[1] = v[1]/mod;
}

void rotateAVector90DegreesRight(float *v, float *n)
/*< TODO >*/
{
	n[0]=-v[1];
	n[1]=v[0];
}

void rotateAVector90DegreesLeft(float *v, float *n)
/*< TODO >*/
{
	n[0]=v[1];
	n[1]=-v[0];
}

float secondVelocityPartialDerivativeNormalToRayDirection(void *par, float *n, float *x, float v)
/*< Calculate the second derivative in the direction normal to ray trajectory >*/
{
        raytrace rt;
        float vpdx, vmdx;
        float dx=0.01;
        float tmp[2];

        rt = (raytrace) par;

        n[0]*=dx;
        n[1]*=dx;

        // ppdx
	rotateAVector90DegreesLeft(n,tmp);
        tmp[0] += x[0];
        tmp[1] += x[1];
        vpdx = sqrtf(1./grid2_vel(rt->grd2,tmp));

        // pmdx
	rotateAVector90DegreesRight(n,tmp);
        tmp[0]+=x[0];
        tmp[1]+=x[1];
        vmdx = sqrtf(1./grid2_vel(rt->grd2,tmp));

        return (vpdx-2.*v+vmdx)/(dx*dx);
}

void secondVelocityPartialDerivativeNormalToRayDirectionAllRayTrajectory(void *par, float **traj, float *dvdn, int nt)
/*< TODO >*/
{
	raytrace rt;
	float x[2];
	float n[2];
	int it;
	float v;

        rt = (raytrace) par;
	
        for(it=0;it<nt;it++){
                x[0]=traj[it][0];
                x[1]=traj[it][1];
                n[0]=traj[it+1][0];
                n[1]=traj[it+1][1];
                getVectorFromTwoPoints(x,n,n);
                getUnitVectorForAVector(n,n);
                v = sqrtf(1./grid2_vel(rt->grd2,x));

                /* Calculate derivative for each ray sample */
                dvdn[it]=secondVelocityPartialDerivativeNormalToRayDirection(rt,n,x,v);

        } // Loop over ray samples

}

float calculateRNIPWithDynamicRayTracing(
                                          void *par, /* Raytrace struct */
                                          float dt, /* Time sampling */
                                          float nt, /* Number of time samples */
                                          float **traj, /* Ray trajectory (z,x) */
                                          float v0 /* Near surface velocity */
)
/*< Calculate RNIP with dynamic ray tracing >*/
{
        raytrace rt; // Raytrace struct
        float x[2]; // Sample coordinate (z,x)
        float *dvdn; // Derivative normal to ray direction
        float mod; // tmp variable
        float rnip; // RNIP parameter

        rt = (raytrace) par;
        dvdn = sf_floatalloc(nt);

	secondVelocityPartialDerivativeNormalToRayDirectionAllRayTrajectory(rt,traj,dvdn,nt);

        /* Initial conditions for a point source */
        x[0]=0.; // q=0
        x[1]=1.; // p=1

        /* Fourth order Runge-Kutta dynamic ray tracing */
        sf_dynamic_runge_init(2,nt,2*dt);
        rnip = sf_dynamic_runge_step(x,rt,dyn_iso_rhs,dvdn,traj,v0);
        sf_dynamic_runge_close();

	free(dvdn);

        return rnip;
}

#endif
