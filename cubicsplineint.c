/*
	 cubicsplineint.c (c)
	 
	 Purpose: Cubic spline interpolation functions.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 19/09/2021

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/
#include <rsf.h>
#include "cubicsplineint.h"

void calculateSplineCoeficients(int n, /* Vectors (x,y) dimension */
				float* x, /* x coordinates */
				float* y, /* y coordinates */
				float** coef, /* Spline coefficients */
				int n_stripes /* Model x axis is divided in n stripes */)
/*< Function to calculate natural cubic spline coefficients

Note: It Receives n points and two vectors x and y with n dimension.
It returns a coefficients vector with 4 coefficients for each of the
n-1 natural cubic splines, coef[(n-1)*4].

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

		/* Calculate spline coefficients */
		for(i=0;i<n-1;i++){
			ha = x[i+1]-x[i];
			coef[k][0+i*4] = (s2[i+1]-s2[i])/(6*ha);
			coef[k][1+i*4] = s2[i]/2;
			coef[k][2+i*4] = (y[k*n+i+1]-y[k*n+i])/ha-(s2[i+1]+2*s2[i])*(ha/6);
			coef[k][3+i*4] = y[k*n+i];
		}
	}
}

void calcInterfacesZcoord(	float *zi, /* Interfaces depth coordinates */
				int nint, /* Number of interfaces */
				float xs, /* x coordinate */
				int si, /* Spline index */
				float **coef /* Cubic spline coefficients */)
/*< Calculate depth coordinates of the interfaces
 * Note: This function calculates interfaces depth coordinates and stores it
 * in the zi vector.
  >*/
{
	int i; // Loop counter

	for(i=0;i<nint;i++){
		zi[i] = coef[i][si*4+0]*xs*xs*xs+
			coef[i][si*4+1]*xs*xs+
			coef[i][si*4+2]*xs+
			coef[i][si*4+3];
	}
}

