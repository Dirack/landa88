/*
	 layer2d.c (c)
	 
	 Purpose: Model layer 2D.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 10/02/2022

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#include <rsf.h>
#include "layer2d.h"

#ifndef _layer2d_h_

typedef struct Layer2d *lay2d;
/* abstract data type */
/*^*/

#endif

struct Layer2d{
	float v;
	float vmax;
	float vmin;
};
/* Concrete data type */

lay2d lay2d_init(float v, /* Layer velocity */
		 float vmax, /* Max layer velocity */
		 float vmin /* Min layer velocity */)
/*< Initialize layer struct >*/
{
	lay2d layer;

	layer = (lay2d) sf_alloc(1,sizeof(struct Layer2d));

	layer->v = v;
	layer->vmax = vmax;
	layer->vmin = vmin;
}

float lay2d_getvmin(lay2d l)
/*<TODO>*/
{return l->vmin;}

float lay2d_getvmax(lay2d l)
/*<TODO>*/
{return l->vmax;}

float lay2d_getvel(lay2d l)
/*< TODO >*/
{return l->v;}
