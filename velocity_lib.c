#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>
#include "velocity_lib.h"

#ifndef GDB_DEBUG
void (*updateVelocityModel)(int*, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float*, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float*, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float*, /* Velocity model disturbance */
			   float*, /* Depth coordinates of sv vector */
			   float*, /* Velocity model */
			   int, /* sv n1 dimwnsion */
			   int, /* sv n2 dimension */
			   float, /* Near surface velocity */
			   float/* Depth velocity gradient */);
#endif
/*^*/

void calculateSplineCoeficients(int n, /* Vectors (x,y) dimension */
				float* x, /* x coordinates */
				float* y, /* y coordinates */
				float** coef, /* Spline coeficients */
				int n_stripes)
/*< Function to calculate natural cubic spline coeficients

Note: It Receives n points and two vectors x and y with n dimension.
It returns a coeficients vector with 4 coeficients for each of the
n-1 natural cubic splines, coef[n-1)*4].

IMPORTANT: The number of points must be equal or major than 3 (n>3)
and x vector must be in crescent order.

>*/
{

	float s2[n]; // Second derivatives matrix
	int i, ip1, ip2, im1, m, k; // Loop counter
	float hb, ha, deltaa, deltab, t; // temporary variables
	float e[n-2]; // hi's vector
	float dp[n-2]; // main diagonal

	/* Vectors dimension must be major than 3 */
	if(n<3){
		fprintf(stderr,"Erro, n<3\n");
		exit(-1);
	}

	/* x vector must be in crescent order */
	for(i=1;i<n;i++){
		if(x[i-1]>x[i]){
			fprintf(stderr,"Erro, vetor x deve possuir ordem crescente\n");
			exit(-2);
		}
	}

	for(k=0;k<n_stripes;k++){
		
		/* Simetric tridiagonal linear system build */
		ha = x[1]-x[0]; deltaa = (y[k*n+1]-y[k*n+0])/ha; m=n-2;
		for(i=0;i<m;i++){
			ip1 = i+1; ip2 = i+2;
			hb = x[ip2]-x[ip1];
			deltab = (y[k*n+ip2]-y[k*n+ip1])/hb;
			e[i] = hb; dp[i] = 2*(ha+hb);
			s2[ip1] = 6*(deltab-deltaa);
			ha=hb; deltaa=deltab;
		}

		/* Gauss elimination */
		for(i=1;i<m;i++){
			ip1=i+1; im1=i-1;
			t = e[im1]/dp[im1];
			dp[i] = dp[i]-t*e[im1];
			s2[ip1] = s2[ip1]-t*s2[i];
		}

		/* Retroactive substitutive solution */
		s2[m]=s2[m]/dp[m-1];
		for(i=m-1;i>0;i--){
			ip1=i+1; im1=i-1;
			s2[i]=(s2[i]-e[im1]*s2[ip1])/dp[im1];
		}
		s2[0]=0; s2[n-1]=0;

		/* Calculate spline coeficients */
		for(i=0;i<n-1;i++){
			ha = x[i+1]-x[i];
			coef[k][0+i*4] = (s2[i+1]-s2[i])/(6*ha);
			coef[k][1+i*4] = s2[i]/2;
			coef[k][2+i*4] = (y[k*n+i+1]-y[k*n+i])/ha-(s2[i+1]+2*s2[i])*(ha/6);
			coef[k][3+i*4] = y[k*n+i];
		}
	}
}

void updateCubicSplineVelModel( float* slow, /* Slowness vector */
		    		int* n, /* n[0]=n1 n2=n[1] */
		    		float* o, /* o[0]=o1 o[1]=o2 */
		    		float* d, /* d[0]=d1 d[1]=d2 */
			    	int dim, /* Dimension of (z,vz) vectors */
				float* sz, /* Spline not Depth coordinates */
				float* sv, /* Spline not Velocity coordinates */
				float gzbg, /* Background velocity gradient in z */
				float v0, /* Near surface velocity */
				int n_stripes)
/*< Funcion to update spline cubic velocity model:
Note:
Make a velocity varying with depth model using the spline cubic interpolation
for a set of points (z,vz) given. TODO

 >*/
{
	int i, j=0, k;
	//int ic;
	float z=0.0;
	//float** coef;
	float v[n_stripes][n[0]];
	int app, app_len=n[1]/n_stripes;

	//coef = sf_floatalloc2(4*(dim-1),n_stripes);

	/* Calculate spline coeficients */
	//calculateSplineCoeficients(dim,sz,sv,coef,n_stripes);

	/* Calculate vel(city function */
	for(k=0;k<n_stripes;k++){

		z = o[0];
		j = 0;

		for(i=1;i<dim;i++){
			
			//ic = (i-1)*4;

			while(z<=sz[i]){
				z = (j*d[0]+o[0])-sz[i-1];
				if(j>=n[0]) break;
				//v[k][j] = coef[k][0+ic]*z*z*z+coef[k][1+ic]*z*z+coef[k][2+ic]*z+coef[k][3+ic];
				v[k][j] = v0+gzbg*z+sv[(k*dim)+i-1];
				//v[k][j] = sv[(k*dim)+i-1];
				j++;
			}
		}
	}

	/* Update slowness model */
	for(k=0;k<n_stripes;k++){

		app = (k*app_len);

		for(i=0;i<n[0];i++){

			for(j=app;j<((k+1)*app_len);j++){
				/* TODO Use 2D eno interpolation to obtain velocity model*/
				slow[j*n[0]+i]=1./(v[k][i]*v[k][i]);
				//#ifdef GDB_DEBUG
				//sf_warning("[%d][%d][%d] %f %f ",i,j,k,v[k][i],slow[j*n[0]+i]);
				//#endif
			} /* Loop over distance */
		} /* Loop over depth */
	} /* Loop over cubic spline functions */
	//#ifdef GDB_DEBUG
	//sf_error("fim");
	//#endif
}

void calcInterfacesZcoord(	float *zi,
				float *sz,
				int nsz,
				int nint){
	int i;
	int npi=nsz/nint;

	for(i=0;i<nint;i++){
		zi[i] = sz[i*npi];
	}
}

void interpolateVelModel(  int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *sv, /* Velocity model disturbance */
			   int nsv, /* sv n2 dimension */
			   float *sz, /* Depth coordinates of sv vector */
			   int nsz, /* sv n1 dimwnsion */
			   float *vel, /* Velocity model */
			   int nvel)
/*< Velocity model interpolation
Note: This function uses a sv control points grid to obtain the complete
velocity model matrix through eno 2D interpolation. The sv vector is the
velocity disturbance of a constant velocity depth gradient model, that
velocity increases linearly with depth for gzbg gradient given.
 >*/
{

	int i, j;
	int k;
	float z;
	float *zi;

	zi = sf_floatalloc(nsv);
	zi[nsv-1] = (n[0]-1)*d[0]+o[0];

	/* Calculate velocity function */
        for(j=0;j<n[1];j++){

		// Calculate interfaces z coordinates
		calcInterfacesZcoord(zi,sz,nsz,nsv-1);
		k=0;
                for(i=0;i<n[0];i++){
			z = i*d[0]+o[0];
			vel[(n[0]*j)+i] = sv[k];
			if(z>zi[k]) k++;
                } /* Loop over depth */

	} /* Loop over distance */
}

void interpolateSlowModel( int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *sv, /* Velociy disturbance */
			   int nsv,
			   float *sz, /* Depth coordinate of disturbance */
			   int nsz,
			   float *slow, /* Slowness model */
			   int nslow /* n1 dimension of sv */)
/*< Slowness model interpolation
Note: This function uses a sv control points grid to obtain the complete
slowness model matrix through eno 2D interpolation. The sv vector is the
velocity disturbance of a constant velocity depth gradient model, that
velocity increases linearly with depth for gzbg gradient given.
 >*/
{

	int i, nm; // Loop counters and indexes

	nm =n[0]*n[1];
	interpolateVelModel(n,o,d,sv,nsv,sz,nsz,slow,nm);

	/* transform velocity to slowness */
	for(i=0;i<nm;i++){
			slow[i] = 1.0/(slow[i]*slow[i]);
	}
}

void enoInterpolation2d(int *n, /* Interpolated vector dimension n1=n[0] n2=n[1] */
			float *o, /* Interpolated vector  axis origin o1=o[0] o2=o[1] */
			float *d, /* Interpolated vector sampling d1=d[0] d2=d[1] */
			float *ov, /* Original vector to interpolate */
			float *iv, /* Interpolated vector */
			int nov1, /* Original vector n1 dimension */
			int nov2 /* Orignanl vector n2 dimension */)
/*< Eno interpolation 2D function
Note: This function interpolates a vector increasing the number of
samples in the interpolated vector using eno interpolation. This vector
is a 2D matrix stored in a vector ov by columns (ov[j*n1+i]) and the new
vector iv will be the interpolated vector.
 >*/
{

	sf_eno2 map;
	float f2[2];
	int i, j, i1, i2;
	float x, y;

	map = sf_eno2_init(3,nov1,nov2);

	sf_eno2_set1(map,ov);

        for(i2=0;i2<n[1];i2++){

                for(i1=0;i1<n[0];i1++){
                        x = i1*d[0]+o[0]; i=x; x -= i;
                        y = i2*d[1]+o[1]; j=y; y -= j;
                        sf_eno2_apply(map,i,j,x,y,&iv[i2*n[0]+i1],f2,FUNC);
                }
        }
        sf_eno2_close(map);
}
