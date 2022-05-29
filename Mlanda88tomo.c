/* Landa 1988 experiment: VFSA velocity inversion based on stereotomography and NIP tomography strategies

This program is a reproduction of the experiment in the article 'A method for determination of velocity and depth from seismic reflection data' avaliable in the doc directory of this repository

The initial velocity model and NIP sources position used in this program is set up using sfnipmodsetup, that program does the forward modeling by ray tracying, from NIP sources to acquisition surface and gets reflection traveltime for a set of reflection ray pairs.

The time misfit is calculated by the difference between the reflection traveltime obtained in the forward modeling and the traveltime calculated by CRE traveltime approximation formula, for RNIP and BETA parameters given. This time misfit is used as a convergence criteria for VFSA global optimization algorithm to obtain optimized velocity model.

*/

#include <math.h>
#include <rsf.h>
#include "tomography.h"
#include "vfsacrsnh_lib.h"
#include "velocity_lib.h"
#include "model2d.h"

#define RAD2DEG 180./SF_PI

int main(int argc, char* argv[])
{
	bool verb; // Verbose parameter
	int n[2]; // Velocity grid dimensions n[0]=n1, n[1]=n2
	float d[2]; // Velocity grid sampling d[0]=d1, d[1]=d2
	float o[2]; // Velocity grid origin o[0]=o1, o[1]=o2
	float** s; // NIP sources position (z,x)
	float* cnewv; // Temporary parameters vector used in VFSA
	float* cnewz; // Temporary parameters vector used in VFSA
	float* otsv; // Optimized parameters vector
	float* otsz; // Optimized parameters vector
	float tmis0; // Best time misfit
	float otmis=0; // Best time misfit
	float deltaE; // Delta (Metrópolis criteria in VFSA)
	float Em0=0; // Energy (VFSA algorithm)
	float PM; // Metrópolis criteria
	float temp=1; // Temperature for VFSA algorithm
	float u=0; // Random number between 0 and 1
	int nit; // Number of VFSA iterations
	float temp0; // Initial temperature for VFSA
	float c0; // Damping factor for VFSA
	int ndim; // n1 dimension in shotsfile, should be equal 2
	int nshot; // n2 dimensions in shotsfile, number of shots
	int nm; // Number of samples in velocity grid n1*n2
	float* a; // Normal Ray initial angle for each NIP source
	float* slow; // slowness model
	int im; // loop counter
	float v; // Velocity temporary variable
	float v0; // Near surface velocity
	int ns; // Number of NIP sources
	int q; // Loop counter for VFSA iteration
	float tmis; // data time misfit value
	float *m0; // CMP's for normal rays
	float *t0; // t0's for normal rays
	float *RNIP; // Rnip parameters vector
	float *BETA; // Beta parameters vector
	float* sz; // Interfaces depth coordinates (cubic spline function)
	int nsz; // Dimension of sz vector
	float dsz; // sz vector sampling
	float osz; // sz vector origin
	float* sv; // Layer's Velocity
	int nsv; // Dimension of sv vector
	float* mis; // Misfit of the current iteration
	int itf; // Interface to invert (index)
	int nx; // Nodepoints for each interface
	float ***data; // Prestack data A(m,h,t)
	int data_n[3]; // n1, n2, n3 dimension of data
	float data_o[3]; // o1, o2, o3 axis origins of data
	float data_d[3]; // d1, d2, d3 sampling of data
	float *otm0; // Normal rays new m0s vector
	float **otm02; // Normal rays new m0s 2D matrix
	float *ott0; // Normal rays new t0s
	float *otrnip; // Optimized RNIP parameter
	float *otbeta; // Optimized BETA parameter (radians)
	float *otangles; // Optimized Normal ray angles (degrees)
	mod2d mod;
        float minvel[5]={1.45,1.60,1.72,1.8,1.9};
        float maxvel[5]={1.55,1.75,1.80,1.9,2.2};
	float vf;
	float *t0p;
	float *m0p;
	float *RNIPp;
	float *BETAp;
	float *ota;
	sf_file shots; // NIP sources (z,x)
	sf_file vel; // background velocity model
	sf_file velinv; // Inverted velocity model
	sf_file angles; // Normal ray angles (degrees)
	sf_file m0s; // Central CMPs m0
	sf_file t0s; // Normal ray traveltimes
	sf_file rnips; // RNIP parameter for each m0
	sf_file betas; // BETA parameter for each m0
	sf_file sz_file; // interfaces z coordinates (cubic spline function)
	sf_file vz_file; // Layer's velocity (input)
	sf_file vspline; // Layers velocity (output)
	sf_file zspline; // Interfaces spline nodes 
	sf_file misfit; // Misfit of the previous iteration
	sf_file misinv; // Misfit result of this VFSA iteration
	sf_file datafile; // Prestack data A(m,h,t)
	sf_file rnips_out; // Optimized RNIP parameter
	sf_file betas_out; // Optimized BETA parameter (radians)
	sf_file angles_out; // Optmized normal ray angles (degrees)
	sf_file t0s_out; // Normal rays new t0s
	sf_file m0s_out; // Normal rays new m0s
	sf_file shots_out; // NIP sources new (z,x) coordinates

	sf_init(argc,argv);

	shots = sf_input("shotsfile");
	vel = sf_input("in");
	velinv = sf_output("out");
	vspline = sf_output("vspline");
	zspline = sf_output("zspline");
	angles = sf_input("anglefile");
	m0s = sf_input("m0s");
	t0s = sf_input("t0s");
	rnips = sf_input("rnips");
	betas = sf_input("betas");
	sz_file = sf_input("sz");
	vz_file = sf_input("vz");
	misfit = sf_input("misfit");
	misinv = sf_output("misinv");
	datafile = sf_input("data");
	rnips_out = sf_output("rnipsout");
	betas_out = sf_output("betasout");
	angles_out = sf_output("anglesout");
	t0s_out = sf_output("t0sout");
	m0s_out = sf_output("m0sout");
	shots_out = sf_output("shotsout");

	/* Velocity model: get 2D grid parameters */
	if(!sf_histint(vel,"n1",n)) sf_error("No n1= in input");
	if(!sf_histint(vel,"n2",n+1)) sf_error("No n2= in input");
	if(!sf_histfloat(vel,"d1",d)) sf_error("No d1= in input");
	if(!sf_histfloat(vel,"d2",d+1)) sf_error("No d2= in input");
	if(!sf_histfloat(vel,"o1",o)) o[0]=0.;
	if(!sf_histfloat(vel,"o2",o+1)) o[1]=0.;

	/* Prestack data cube */
	if(!sf_histint(datafile,"n1",data_n)) sf_error("No n1= in data");
	if(!sf_histint(datafile,"n2",data_n+1)) sf_error("No n2= in data");
	if(!sf_histint(datafile,"n3",data_n+2)) sf_error("No n3= in data");
	if(!sf_histfloat(datafile,"d1",data_d)) sf_error("No d1= in data");
	if(!sf_histfloat(datafile,"d2",data_d+1)) sf_error("No d2= in data");
	if(!sf_histfloat(datafile,"d3",data_d+2)) sf_error("No d3= in data");
	if(!sf_histfloat(datafile,"o1",data_o)) sf_error("No o1= in data");
	if(!sf_histfloat(datafile,"o2",data_o+1)) sf_error("No o2= in data");
	if(!sf_histfloat(datafile,"o3",data_o+2)) sf_error("No o3= in data");
	
	if(!sf_getbool("verb",&verb)) verb=true;
	/* verbose parameter (y/n) */

	if(!sf_getfloat("v0",&v0)) v0=1.5;
	/* Near surface velocity (Km/s) */

	if(!sf_getint("nit",&nit)) nit=1;
	/* Number of VFSA iterations */

	if(!sf_getfloat("temp0",&temp0)) temp0=5;
	/* Initial temperature for VFSA algorithm */

	if(!sf_getfloat("c0",&c0)) c0=0.1;
	/* Damping factor for VFSA algorithm */

	if(!sf_getint("itf",&itf)) itf=0;
	/* Interface to invert, others will be ignored */

	/* Shotsfile: get shot points */
	if(!sf_histint(shots,"n1",&ndim) || 2 != ndim)
		sf_error("Must have n1=2 in shotsfile");
	if(!sf_histint(shots,"n2",&nshot)) sf_error("No n2= in shotfile");
	s = sf_floatalloc2(ndim,nshot);
	sf_floatread(s[0],ndim*nshot,shots);
	sf_fileclose(shots);

	/* Cubic spline vector */
	if(!sf_histint(sz_file,"n1",&nsz)) sf_error("No n1= in sz file");
	if(!sf_histfloat(sz_file,"d1",&dsz)) sf_error("No d1= in sz file");
	if(!sf_histfloat(sz_file,"o1",&osz)) sf_error("No o1= in sz file");
	if(!sf_histint(vz_file,"n1",&nsv)) sf_error("No n1= in sv file");

	/* Build cubic spline velocity matrix */
	sv = sf_floatalloc(nsv);
	sz = sf_floatalloc(nsz);
	sf_floatread(sz,nsz,sz_file);
	sf_floatread(sv,nsv,vz_file);
	vf = sv[nsv-1];

	interfaceInterpolationFromNipSources(s,nshot,sz,nsz,osz,dsz,nsv);
	mod = mod2d_init(nsv,sv,minvel,maxvel,nsz,osz,dsz,sz);

	/* VFSA parameters vectors */
	cnewv = sf_floatalloc(nsv);
	cnewz = sf_floatalloc(nsz); 
	otsv = sf_floatalloc(nsv);
	otsz = sf_floatalloc(nsz);

	/* Read prestack data cube */
	data = sf_floatalloc3(data_n[0],data_n[1],data_n[2]);
	sf_floatread(data[0][0],data_n[0]*data_n[1]*data_n[2],datafile);

	/* Anglefile: get initial emergence angle */
	if(!sf_histint(angles,"n1",&ns)) sf_error("No n1= in anglefile");
	a = sf_floatalloc(ns);
	sf_floatread(a,ns,angles);
	if(ns!=nshot) sf_error("n1 in anglefile should be equal to n2 in shotsfile!");

	/* allocate parameters vectors */
	m0 = sf_floatalloc(ns);
	otm0 = sf_floatalloc(ns);
	otm02 = sf_floatalloc2(2,ns);
	sf_floatread(m0,ns,m0s);
	t0 = sf_floatalloc(ns);
	ott0 = sf_floatalloc(ns);
	sf_floatread(t0,ns,t0s);
	RNIP = sf_floatalloc(ns);
	otrnip = sf_floatalloc(ns);
	sf_floatread(RNIP,ns,rnips);
	BETA = sf_floatalloc(ns);
	otbeta = sf_floatalloc(ns);
	otangles = sf_floatalloc(ns);
	sf_floatread(BETA,ns,betas);
	t0p = sf_floatalloc(ns);
	m0p = sf_floatalloc(ns);
	RNIPp = sf_floatalloc(ns);
	BETAp = sf_floatalloc(ns);
	ota = sf_floatalloc(ns);

	for(im=0;im<ns;im++){
		otm0[im]=m0[im];
		otm02[im][0]=0.;
		otm02[im][1]=m0[im];
		ott0[im]=t0[im];
		otrnip[im]=RNIP[im];
		otbeta[im]=BETA[im];
		ota[im]=a[im];
		t0p[im]=t0[im];
		m0p[im]=m0[im];
		RNIPp[im]=RNIP[im];
		BETAp[im]=BETA[im];
	}

	/* get slowness squared (Background model) */
	nm = n[0]*n[1];
	slow =  sf_floatalloc(nm);
	sf_floatread(slow,nm,vel);

	for(im=0;im<nm;im++){
		v = slow[im];
		slow[im] = 1./(v*v);
	}

	if(verb){
		sf_warning("Command line Parameters");
		sf_warning("v0=%f nit=%d temp0=%f c0=%f itf=%d",v0,nit,temp0,c0,itf);
		sf_warning("Input file (Velocity model)");
		sf_warning("n1=%d d1=%f o1=%f",*n,*d,*o);
		sf_warning("n2=%d d2=%f o2=%f",*(n+1),*(d+1),*(o+1));
		sf_warning("Input file (Prestack data)");
		sf_warning("n1=%d d1=%f o1=%f",*data_n,*data_d,*data_o);
		sf_warning("n2=%d d2=%f o2=%f",*(data_n+1),*(data_d+1),*(data_o+1));
		sf_warning("n3=%d d3=%f o3=%f",*(data_n+2),*(data_d+2),*(data_o+2));
		sf_warning("Input file (shotsfile)");
		sf_warning("n1=%d",ndim);
		sf_warning("n2=%d",nshot);
		sf_warning("Input file (anglefile, t0s, m0s, rnips, betas)");
		sf_warning("n1=%d",ns);
		sf_warning("Input file (sz)");
		sf_warning("nz=%d",nsz);
		sf_warning("Input file (vz)");
		sf_warning("nv=%d",nsv);
	}

	/* Use previous misfit as the initial misfit value */
	mis=sf_floatalloc(1);
	sf_floatread(mis,1,misfit);
	tmis0=mis[0];
	otmis=tmis0;
	free(mis);

	/* Velocity model from inversion */
	sf_putint(velinv,"n1",n[0]);
	sf_putint(velinv,"n2",n[1]);
	sf_putint(velinv,"n3",1);
	sf_putfloat(velinv,"d1",d[0]);
	sf_putfloat(velinv,"d2",d[1]);
	sf_putfloat(velinv,"o1",o[0]);
	sf_putfloat(velinv,"o2",o[1]);
	sf_putfloat(velinv,"d3",1);
	sf_putfloat(velinv,"o3",0);

	/* velocity and interfaces (output) */
	sf_putint(vspline,"n1",nsv);
	sf_putint(vspline,"n2",1);
	sf_putint(zspline,"n1",nsz);
	sf_putfloat(zspline,"o1",osz);
	sf_putfloat(zspline,"d1",dsz);
	sf_putint(zspline,"n2",1);

	/* Misfit value */
	sf_putint(misinv,"n1",1);
	sf_putint(misinv,"n2",1);
	sf_putint(misinv,"n3",1);

	/* t0s values */
	sf_putint(t0s_out,"n1",ns);
	sf_putint(t0s_out,"n2",1);
	sf_putint(t0s_out,"n3",1);

	/* m0s values */
	sf_putint(m0s_out,"n1",ns);
	sf_putint(m0s_out,"n2",1);
	sf_putint(m0s_out,"n3",1);

	/* shotsfiles out */
	sf_putint(shots_out,"n1",2);
	sf_putint(shots_out,"n2",ns);
	sf_putint(shots_out,"n3",1);

	/* RNIPs parameters out */
	sf_putint(rnips_out,"n1",ns);
	sf_putint(rnips_out,"n2",1);
	sf_putint(rnips_out,"n3",1);

	/* BETAs parameters out (Radians) */
	sf_putint(betas_out,"n1",ns);
	sf_putint(betas_out,"n2",1);
	sf_putint(betas_out,"n3",1);

	/* BETAS parameters out (Degrees) */
	sf_putint(angles_out,"n1",ns);
	sf_putint(angles_out,"n2",1);
	sf_putint(angles_out,"n3",1);
	
	/* Intiate optimal parameters vectors */
	for(im=0;im<nsz;im++)
		otsz[im]=sz[im];
	for(im=0;im<nsv;im++)
		otsv[im]=sv[im];

	// TODO Next layer's velocity will be the same
	// in order to avoid interference during inversion
	sv[itf+1]=sv[itf];

	/* Number of NIP sources for each interface */
	nx = ns/(nsv-1);

	//disturbParameters(temp0,otsv,otsz,mod,0.1,itf);
	/* Setup NIP sources */
	//modelSetup(s, nx,  m0, t0, a,  itf,  n,  d,  o,  slow);
	//tmis0=calculateTimeMisfit(s,otsv[0],t0,m0,RNIP,BETA,n,o,d,slow,a,nx,itf,data,data_n,data_o,data_d,otsz,nsz,osz,dsz,otsv,nsv);

	//for(im=0;im<ns;im++){
	//	tmis0 += fabs(t0[im]-t0p[im]);
	//	tmis0 += fabs(m0[im]-m0p[im]);
		//sf_warning("tmis0=%f t0=%f t0p=%f",tmis0,t0[im],t0p[im]);
	//}
	//tmis0=10;
	//sf_error("oi");
	//otmis=tmis0;

	/* Very Fast Simulated Annealing (VFSA) algorithm */
	for (q=0; q<nit; q++){
	
		/* calculate VFSA temperature for this iteration */
		temp=getVfsaIterationTemperature(q,c0,temp0);
						
		/* parameter disturbance */
		disturbParameters(temp,cnewv,cnewz,mod,1,itf);

		/* Function to update velocity model */
		buildSlownessModelFromVelocityModel(n,o,d,cnewv,nsv,cnewz,nsz,osz,dsz,slow,nm);

		/* Setup NIP sources */
		modelSetup(s, nx,  m0, t0, a,  itf,  n,  d,  o,  slow);

		/**/
		interfaceInterpolationFromNipSources(s,nshot,cnewz,nsz,osz,dsz,nsv);

		/* Function to update velocity model */
		buildSlownessModelFromVelocityModel(n,o,d,cnewv,nsv,cnewz,nsz,osz,dsz,slow,nm);

		tmis=0;
	
		/* Calculate time misfit through forward modeling */		
		tmis=calculateTimeMisfit(s,cnewv[0],t0,m0,RNIP,BETA,n,o,d,slow,a,nx,itf,data,data_n,data_o,data_d,cnewz,nsz,osz,dsz,cnewv,nsv);

		//for(im=0;im<ns;im++){
		//	tmis += fabs(t0[im]-t0p[im])*fabs(t0[im]-t0p[im]);
		//	tmis += fabs(m0[im]-m0p[im])*fabs(m0[im]-m0p[im]);
			//tmis += fabs(RNIP[im]-RNIPp[im])*fabs(RNIP[im]-RNIPp[im]);
			//tmis += fabs(BETA[im]-BETAp[im])*fabs(BETA[im]-BETAp[im]);
		//}
		if(fabs(tmis) > fabs(tmis0) ){
			//sf_warning("tmis=%f v=%f",tmis,cnewv[itf]);
			otmis = fabs(tmis);
			/* optimized parameters */
			for(im=0;im<nsz;im++)
				otsz[im]=cnewz[im];
			for(im=0;im<itf+1;im++)
				otsv[im]=cnewv[im];
			#define UPDATE_CRS_PARAMETERS
			#ifdef UPDATE_CRS_PARAMETERS
			for(im=0;im<ns;im++){
				otm0[im]=m0[im];
				otm02[im][1]=m0[im];
				ott0[im]=t0[im];
				otrnip[im]=RNIP[im];
				otbeta[im]=BETA[im];
				ota[im]=a[im];
			}
			#endif
			tmis0 = fabs(tmis);
		}else{
			for(im=0;im<ns;im++){
				m0[im]=otm0[im];
				t0[im]=ott0[im];
				RNIP[im]=otrnip[im];
				BETA[im]=otbeta[im];
				a[im]=ota[im];
			}
		}

		/* VFSA parameters update condition */
		deltaE = fabs(tmis) - Em0;
		
		/* Metrópolis criteria */
		PM = expf(-deltaE/temp);
		
		if (deltaE<=0){
			for(im=0;im<nsz;im++)
				sz[im]=cnewz[im];
			for(im=0;im<itf+1;im++)
				mod2d_setlayervel(mod,im,cnewv[im]);
			Em0 = fabs(tmis);
		} else {
			u=getRandomNumberBetween0and1();
			if (PM > u){
				for(im=0;im<nsz;im++)
					sz[im]=cnewz[im];
				for(im=0;im<itf+1;im++)
					mod2d_setlayervel(mod,im,cnewv[im]);
				Em0 = fabs(tmis);
			}/*else{
				for(im=0;im<ns;im++){
					m0[im]=otm0[im];
					t0[im]=ott0[im];
					RNIP[im]=otrnip[im];
					BETA[im]=otbeta[im];
					a[im]=ota[im];
				}
			}*/
		}	
			
		sf_warning("%d/%d interface=%d => Semblance(%f) v=%f v=%f %f;",q+1,nit,itf,otmis,otsv[itf],cnewv[itf],tmis);

	} /* loop over VFSA iterations */

	/* Generate optimal velocity model */
	otsv[nsv-1]=vf;
	updateVelocityModel(n,o,d,otsv,nsv,otsz,nsz,osz,dsz,slow,nm);

	/* Write velocity model file */
	sf_floatwrite(slow,nm,velinv);
	
	for(im=0;im<ns;im++)
		otangles[im]=otbeta[im]*RAD2DEG+180.;

	/* Write velocity cubic spline function */
	sf_floatwrite(otsv,nsv,vspline);
	sf_floatwrite(otsz,nsz,zspline);
	sf_floatwrite(&otmis,1,misinv);
	sf_floatwrite(ott0,ns,t0s_out);
	sf_floatwrite(otm0,ns,m0s_out);
	sf_floatwrite(otm02[0],2*ns,shots_out);
	sf_floatwrite(otrnip,ns,rnips_out);
	sf_floatwrite(otbeta,ns,betas_out);
	sf_floatwrite(otangles,ns,angles_out);

}
