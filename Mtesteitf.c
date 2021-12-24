#include <rsf.h>
#include "interface2d.h"

int main(int argc, char* argv[]){
	
	itf2d itf;
	float *sz;
	int i;
	float f;

	sz = sf_floatalloc(11);

	for(i=0;i<11;i++)
		sz[i]=1.0;
	
	itf = itf2d_init(sz,11,-2.0,1.0);

	f = getZCoordinateOfInterface(itf,1.5);	
}
