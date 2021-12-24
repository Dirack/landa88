#include <rsf.h>
#include "interface2d.h"
#include "cubicsplineint.h"

#ifndef _interface2d_h_

typedef struct Interface2d *itf2d;
/* abstract data type */
/*^*/

#endif

struct Interface2d{
	float *coef;
	int n;
	float o;
	float d;
};
/* Concrete data type */

itf2d itf2d_init(float *sz,
		 int n1,
		 float o1,
		 float d1)
/*< TODO >*/
{
	itf2d itf;
	int i;
	float *x;

	itf = (itf2d) sf_alloc(1,sizeof(struct Interface2d));

	itf->n = n1;
	itf->o = o1;
	itf->d = d1;

	x = sf_floatalloc(n1);

	for(i=0;i<n1;i++){
		x[i] = i*d1+o1;
	}

	/* Calculate coefficients matrix (interfaces interpolation) */
	itf->coef = sf_floatalloc(4*(n1-1));
	calculateSplineCoeficients(n1,x,sz,&itf->coef,1);

	return itf;
}

float getZCoordinateOfInterface(
				itf2d itf,
				float x)
/*< TODO >*/
{
	float z;
	int is;
	float xs;

	is = (x-itf->o)/itf->d;

	xs = x - ((is*itf->d)+itf->o);

	z = itf->coef[is*4+0]*xs*xs*xs+itf->coef[is*4+1]*xs*xs+itf->coef[is*4+2]*xs+itf->coef[is*4+3];

	return z;
}
