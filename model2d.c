/*
	 model2d.c (c)
	 
	 Purpose: Model 2D.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 10/02/2022

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#include <rsf.h>
#include "model2d.h"
#include "interface2d.h"
#include "layer2d.h"

#ifndef _model2d_h_

typedef struct Model2d *mod2d;
/* abstract data type */
/*^*/

#endif

struct Model2d{
	itf2d *itf;
	lay2d *lay;
	int nlay;
};
/* Concrete data type */

mod2d mod2d_init(int nlay, /* Number of layers */
		 float *sv,
		 float *minsv,
		 float *maxsv,
		 int nsz,
		 float osz,
		 float dsz,
		 float *sz
		 )
/*< TODO >*/
{
	int i, j;
	mod2d mod;
	float *szz;

	szz = sf_floatalloc(nsz/(nlay-1));

	mod = (mod2d) sf_alloc(1,sizeof(struct Model2d));

	mod->nlay = nlay;
	mod->lay = (lay2d*) sf_alloc(nlay,sizeof(lay2d));
	mod->itf = (itf2d*) sf_alloc(nlay-1,sizeof(itf2d));

	for(i=0;i<nlay;i++)
		mod->lay[i]=lay2d_init(sv[i],maxsv[i],minsv[i]);

	for(i=0;i<nlay-1;i++){
		for(j=0;j<nsz/(nlay-1);j++){
                	szz[j]=sz[j+((i-1)*nsz/(nlay-1))];
		}
		mod->itf[i]=itf2d_init(szz,nsz/(nlay-1),osz,dsz);
	}
	return mod;
}

float mod2d_getlayervmin(mod2d m, int n)
/*< TODO >*/
{return lay2d_getvmin(m->lay[n]);}

float mod2d_getlayervmax(mod2d m, int n)
/*< TODO >*/
{return lay2d_getvmax(m->lay[n]);}

float mod2d_getlayervel(mod2d m, int n)
/*< TODO >*/
{return lay2d_getvel(m->lay[n]);}

int mod2d_getnumlayers(mod2d m)
/*< TODO >*/
{return m->nlay;}

void mod2d_setlayervel(mod2d m, int n, float v)
/*< TODO >*/
{lay2d_setvel(m->lay[n],v);}
