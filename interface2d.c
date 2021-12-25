/*
	 interface2d.c (c)
	 
	 Purpose: Model interface 2D interpolated with cubic splines.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 25/12/2021

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/
#include <rsf.h>
#include "interface2d.h"
#include "cubicsplineint.h"

#ifndef _interface2d_h_

typedef struct Interface2d *itf2d;
/* abstract data type */
/*^*/

#endif

struct Interface2d{
	float *coef; // Cubic spline coefficients Matrix
	float *z; // Nodepoints z(x)
	int n; // Number of nodes
	float o; // Axis origin
	float d; // Nodes sampling
};
/* Concrete data type */

itf2d itf2d_init(float *sz, /* Interface nodepoints z(x) */
		 int n1, /* Number of nodes */
		 float o1, /* Axis origin */
		 float d1 /* nodes sampling */)
/*< Initialize interface struct >*/
{
	itf2d itf;
	int i;
	float *x;

	itf = (itf2d) sf_alloc(1,sizeof(struct Interface2d));

	itf->n = n1;
	itf->o = o1;
	itf->d = d1;

	x = sf_floatalloc(n1);
	itf->z = sf_floatalloc(n1);

	for(i=0;i<n1;i++){
		x[i] = i*d1+o1;
		itf->z[i] = sz[i];
	}

	/* Calculate coefficients matrix (interfaces interpolation) */
	itf->coef = sf_floatalloc(4*(n1-1));
	calculateSplineCoeficients(n1,x,sz,&itf->coef,1);

	return itf;
}

float getZCoordinateOfInterface(
				itf2d itf, /* Interface */
				float x /* x coordinate */)
/*< Get z coordinate of the interface given x coordinate >*/
{
	float z;
	int is;
	float xs;

	is = (x-itf->o)/itf->d;

	xs = x - ((is*itf->d)+itf->o);

	z = itf->coef[is*4+0]*xs*xs*xs+itf->coef[is*4+1]*xs*xs+itf->coef[is*4+2]*xs+itf->coef[is*4+3];

	return z;
}

int itf2d_n(itf2d itf)
/*< Get number of interface nodepoints >*/
{return itf->n;}

float itf2d_o(itf2d itf)
/*< Get interface axis origin >*/
{return itf->o;}

float itf2d_d(itf2d itf)
/*< Get interface nodepoints sampling >*/
{return itf->d;}

void itf2d_setz(
		itf2d itf, /* Interface */
		float *z /* Nodepoints z(x) */)
/*< Update interface nodepoints z(x) >*/
{
	int i;
	for(i=0;i<itf->n;i++)
		itf->z[i]=z[i];
}
