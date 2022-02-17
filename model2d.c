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
		 float *sv, /* Layer's velocities */
		 float *minsv, /* Layer's minimum velocities */
		 float *maxsv, /* Layer's maximum velocities */
		 int nsz, /* Number of interfaces nodepoints */
		 float osz, /* Interfaces X coordinate origin */
		 float dsz, /* Interfaces sampling */
		 float *sz /* Interfaces nodepoints */
		 )
/*< Initialize model 2D struct >*/
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
                	szz[j]=sz[j+i*(nsz/(nlay-1))];
		}
		mod->itf[i]=itf2d_init(szz,nsz/(nlay-1),osz,dsz);
	}
	return mod;
}

float mod2d_getlayervmin(mod2d m, int n)

/*< Get a layer minimum velocity >*/
{return lay2d_getvmin(m->lay[n]);}

float mod2d_getlayervmax(mod2d m, int n)
/*< Get a layer maximum velocity >*/
{return lay2d_getvmax(m->lay[n]);}

float mod2d_getlayervel(mod2d m, int n)
/*< Get a layer velocity >*/
{return lay2d_getvel(m->lay[n]);}

int mod2d_getnumlayers(mod2d m)
/*< Get number of layers in model 2D >*/
{return m->nlay;}

void mod2d_setlayervel(mod2d m, int n, float v)
/*< Set a layer velocity >*/
{lay2d_setvel(m->lay[n],v);}

void mod2d_setlayervmin(mod2d m, int n, float vmin)
/*< Set a layer velocity minimum >*/
{lay2d_setvmin(m->lay[n],vmin);}

void mod2d_setlayervmax(mod2d m, int n, float vmax)
/*< Set a layer velocity maximum >*/
{lay2d_setvmax(m->lay[n],vmax);}

void mod2d_setinterfacesnodes(mod2d m, int n, float *z)
/*< Set a interface nodepoints >*/
{itf2d_setZNodepoints(m->itf[n],z);}

void mod2d_getinterfacesnodes(mod2d m, int n, float *z)
/*< Get all interfaces nodepoints >*/
{itf2d_getZNodepoints(m->itf[n],z);}

void mod2d_getallinterfacesnodes(mod2d m, float *z)
/*< Get all interfaces nodepoints >*/
{
	int i, j;
	int n;
	float *zz;

	n = itf2d_n(m->itf[0]);
	zz = sf_floatalloc(n);
	
	for(i=0;i<m->nlay-1;i++){
		itf2d_getZNodepoints(m->itf[i],zz);
		for(j=0;j<n;j++){
			z[i*n+j]=zz[j];
		}
	}
}

int mod2d_getnuminterfacesnodes(mod2d m)
/*< Get number of interfaces nodepoints >*/
{return itf2d_n(m->itf[0]);}