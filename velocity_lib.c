/*
	 velocity_lib.c (c)
	 
	 Purpose: Functions to update velocity model.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 19/09/2021

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>
#include "interface2d.h"
#include "velocity_lib.h"

float calculateLocationMisfit( float **s, /* NIP sources location */
			   	float *sz, /* Depth coordinates of interfaces */
			   	int nsz, /* NIP sources number for each interface */
			   	float osz, /* sz origin */
			   	float dsz, /* sz sampling */
				int nshot, /* Dimension of the sz vector */
				int itf /* Interface to invert */)
/*< Calculate misfit between NIP sources and interfaces
Note: This function calculates L2 norm of the distances between NIP
sources location in Z and interfaces. The assumption is that NIP sources
are located in the interfaces, so the best velocity model minimize the distance
between then

 >*/
{

	int i; // loop counter
	int l=0; // Splines index
	float *zi; // Temporary vector to store depth coordinates
	float *x; // X coordinates of interface being inverted
	float** coef; // Cubic splines coefficients
	float misfit = 0.; // Misfit sum
	float* szz; // Z coordinates of interface being inverted

	x = sf_floatalloc(nsz);
	szz = sf_floatalloc(nsz);

	for(i=0;i<nsz;i++){
		x[i] = i*dsz+osz;
		szz[i]=sz[i+(itf*nsz)];
	}

	/* Calculate coefficients matrix (interfaces interpolation) */
	coef = sf_floatalloc2(4*(nsz-1),1);
	calculateSplineCoeficients(nsz,x,szz,coef[0]);

	zi = sf_floatalloc(1);

	/* Calculate interfaces z coordinates and misfit */
	for(i=0;i<nshot;i++){

		l = (int) (s[i][1]-osz)/dsz;

		calcInterfacesZcoord(zi,1,s[i][1]-x[l],l,coef);

		misfit += fabs(zi[0]-s[i][0]);
	}

	return sqrt(misfit*misfit);
}

void updateVelocityModel(  int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *sv, /* Velocity model disturbance */
			   int nsv, /* Dimension of sv the vector */
			   float *sz, /* Depth coordinates of interfaces */
			   int nsz, /* Dimension sz the vector */
			   float osz, /* sz origin */
			   float dsz, /* sz sampling */
			   float *vel, /* Velocity model */
			   int nvel /* Dimension of the vel vector */)
/*< Velocity model update
Note: This function uses a sv (layers velocity) vector and sz (depth interfaces
coordinates) vector to build the depth velocity model. There is nsv constant
velocity layers in the model and nsv-1 interfaces separating them.
These interfaces are described with nsz control points (nodes) in the sz vector and
they are interpolated using natural cubic spline interpolation.
 >*/
{

	int i, j; // Loop counters
	int k; // Layers index
	int l=0; // Splines index
	float z; // Depth coordinate
	float *zi; // Temporary vector to store depth coordinates
	int nx=nsz/(nsv-1); // Number of points for each interface
	float *x; // X coordinates nodes
	float** coef; // Cubic spline coefficients
	float xx; // X coordinates in velocity model
	float *szz;

	x = sf_floatalloc(nx);
	szz = sf_floatalloc(nx);
//dsz=0.5;
	for(i=0;i<nx;i++)
		x[i] = i*dsz+osz;

	//for(i=0;i<nx;i++) sf_warning("x[%d]=%f %f %f",i,x[i],dsz,osz);
	//sf_error("deu");

	/* Calculate coefficients matrix (interfaces interpolation) */
	coef = sf_floatalloc2(4*(nx-1),nsv-1);
	for(i=0;i<nsv-1;i++){
		for(j=0;j<nx;j++) szz[j]=sz[(i*nx)+j];
		calculateSplineCoeficients(nx,x,szz,coef[i]);
	}

	zi = sf_floatalloc(nsv);
	zi[nsv-1] = (n[0]-1)*d[0]+o[0];

	/* Calculate velocity function */
        for(j=0;j<n[1];j++){

		xx = d[1]*j+o[1];
		if(xx>x[l+1]) l++;
		/* Calculate interfaces z coordinates */
		calcInterfacesZcoord(zi,nsv-1,xx-x[l],l,coef);
		k=0;
                for(i=0;i<n[0];i++){
			z = i*d[0]+o[0];
			if(z>zi[k]) k++;
			vel[n[0]*j+i] = sv[k];
			//sf_warning("v[%d]=%f xx=%f xl=%f xl1=%f l=%d",i,vel[n[0]*j+i],xx,x[l],x[l+1],l);
			//if(l>=20) sf_error("capa");
                } /* Loop over depth */
		//sf_error("capa");

	} /* Loop over distance */

	free(x);
}

void buildSlownessModelFromVelocityModel(int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			 		 float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
					 float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
					 float *sv, /* Velociy disturbance */
					 int nsv, /* Dimension of sv vector */
					 float *sz, /* Depth coordinates of interfaces */
					 int nsz, /* Dimension of sz vector */
					 float osz, /* sz origin */
					 float dsz, /* sz sampling */
					 float *vel, /* Velocity model */
					 int nslow /* Dimension of vel vector */)
/*< Slowness model build from velocity model
Note: This function is a function wrapper to updateVelocityModel function.
It calls that function to update the velocity model and build the slowness
model matrix using the slowness definition slow=(1.0/(v*v)). 
 >*/
{

	int i, nm; // Loop counters and indexes

	nm =n[0]*n[1];
	updateVelocityModel(n,o,d,sv,nsv,sz,nsz,osz,dsz,vel,nm);

	/* transform velocity to slowness */
	for(i=0;i<nm;i++){
			vel[i] = 1.0/(vel[i]*vel[i]);
	}
}

