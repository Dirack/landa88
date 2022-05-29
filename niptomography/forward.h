#include "grid2.h"
#include "atela.h"
#include "dynamic.h"
#include "interface2d.h"
#define DT 0.001

#ifndef __FORWARD_H__
#define __FORWARD_H__

typedef struct RayTrace* raytrace;
/* abstract data type */
/*^*/

struct RayTrace {
    bool sym;
    int dim, nt;
    float dt, z0;
    grid2 grd2;
};
/* concrete data type */

void rayEndpointWarning(float *x,float *p,float **traj,float t)
/*< Output error message if ray get to the model side or bottom >*/
{
        /* TODO to correct the way you treat side rays */
        sf_warning("Ray endpoint => x=%f y=%f p[0]=%f p[1]=%f",x[1],x[0],p[0],p[1]);
        sf_warning("Ray starting point=> x=%f y=%f",traj[0][1],traj[0][0]);
        sf_warning("Ray traveltime => t=%f",t);
        sf_warning("Bad ray angle, ray get to the model side/bottom");
}

void setInitialRayPointAndRayVector(float **s, /* NIP sources Matrix */
                                float *x, /* Initial ray point (z,x) */
                                float *p, /* Initial ray vector */
                                int is, /* Source index */
                                float a /* Initial ray angle (radians) */)
/*< Set the initial ray point and ray vector >*/
{
                x[0]=s[is][0];
                x[1]=s[is][1];
                p[0] = -cosf(a);
                p[1] = sinf(a);
}

void calculateEscapeVector(
                           float *x, /* Ray endpoint */
                           float **traj, /* Ray trajectory */
                           int it /* Endpoint index */)
/*< Calculate Escape vector from ray trajectory, x is changed inside the function >*/
{
        int i; 
        i = it >= 2 ? it - 2 : it - 1;
        x[0]=traj[it][0];
        x[1]=traj[it][1];
        x[0]-=traj[i][0];
        x[1]-=traj[i][1];
}

float calculateBetaWithRayTrajectory(
                                     float *x, /* Ray endpoint */
                                     float **traj /* Ray trajectory */,
                                     int it /* Endpoint index */)
/*< Calculate BETA parameter using dot product with unit vector pointing upward
Note: x is changed inside the function
>*/
{
        float xx;

        calculateEscapeVector(x,traj,it);
        xx=sqrt(x[0]*x[0]+x[1]*x[1]);
        xx=acos(-x[0]/xx);
        if(x[1]<0) xx = -xx;
        return xx;
}

float getTransmissionAngleSnellsLaw(
		float ei, /* Incident angle in radians */
		float vi, /* Velocity in incident ray layer */
		float vt /* Velocity in transmited ray layer */)
/*< Snells Law to get transmission angle >*/
{
	return asinf((vt/vi)*sinf(ei));
}

float calculateIncidentAngle(
			     void *par, /* Interface struct */
			     float **traj, /* Ray trajectory */
			     int ir /* Ray sample index */)
/*< Calculate incident angle 
Note: This function calculates a vector using the diffence between
traj[ir] and traj[ir-1]. The angle between this vector and the normal
to the interface will be the incident angle returned
>*/
{
	float x[2];
	float n[2];
	float ei;
	itf2d it2;

	it2 = (itf2d) par;

	n[0] = getZCoordinateOfInterface(it2,traj[ir][1]+0.01)
	-getZCoordinateOfInterface(it2,traj[ir][1]);
	n[1] = 0.01;

	ei = sqrtf(n[0]*n[0]+n[1]*n[1]);
	n[0]/=ei;
	n[1]/=ei;

	ei = n[0];
	n[0]=-n[1];
	n[1]=ei;

	x[0] = traj[ir][0]-traj[ir-1][0];
	x[1] = traj[ir][1]-traj[ir-1][1];

	ei = sqrtf(x[0]*x[0]+x[1]*x[1]);
	ei = acos((x[0]*n[0]+x[1]*n[1])/ei);

	return ei;
}


void first_deriv( float h /* step */,
		float *zx /* f(x) */,
		float *der /* Derivative */,
		int n /* vectors dim */)
/*< Calculate first derivative numerically >*/
{
	int i;

	sf_deriv_init(n, 6, 0.);
        sf_deriv(zx,der);

        for (i=0; i < n; i++) {
                der[i] /= h;
        }
}


void second_deriv( float h /* passo h */,
		    float *zx /* f(x) */,
		    float *der /* Derivative */,
		    int n /* vectors dim */)
/*< Calculate second derivative numerically >*/
{
	float *firstder;
	firstder = sf_floatalloc(n);
	first_deriv(h,zx,firstder,n);
	first_deriv(h,firstder,der,n);
	free(firstder);
}

float calculateInterfaceCurvature(
				  void *par, /* Interface struct */
				  float x /* Coordinate (z,x) */)
/*< Calculate interface curvature using curvature formula for 2D function >*/
{

	float kf; // Interface curvature
	float f1, f2; // tmp variables
	int is; // Spline index
	int i;
	itf2d it2; // Interface struct
	float coef[4]; // Cubic spline coefficients matrix
	float a, b, c, d; // Spline coefficients
	//float xp, xm;
	//float xm7h,xm6h,xm5h,xm4h,xm3h, xm2h, xmh;
	float zx[5],dz1dx[5],dz2dx[5];

	it2 = (itf2d) par;

	is = (x-itf2d_o(it2))/itf2d_d(it2);
	itf2d_getSplineCoefficients(it2,coef,is);

	/*a = coef[0];
	b = coef[1];
	c = coef[2];
	d = coef[3];*/

	//f1 = 3*a*x*x+2*b*x+c;
	//xp = x+0.001;
	//xm = x-0.001;
	//printf("x=%f xm=%f xp=%f\n",x,xm,xp);
	for(i=0;i<5;i++)
		zx[i]=getZCoordinateOfInterface(it2,x-(i+2)*0.001);
	first_deriv(0.001,zx,dz1dx,5);
	first_deriv(0.001,dz1dx,dz2dx,5);
	//f2 = fabs(6*a*x+2*b);
	/*xm7h = getZCoordinateOfInterface(it2,x-7*0.01);
	xm6h = getZCoordinateOfInterface(it2,x-6*0.01);
	xm5h = getZCoordinateOfInterface(it2,x-5*0.01);
	xm4h = getZCoordinateOfInterface(it2,x-4*0.01);
	xm3h = getZCoordinateOfInterface(it2,x-3*0.01);
	xm2h = getZCoordinateOfInterface(it2,x-2*0.01);
	xmh = getZCoordinateOfInterface(it2,x-0.01);*/
	//x = getZCoordinateOfInterface(it2,x);

	f1 = dz1dx[2];
	f2 = fabs(dz2dx[2]);
	//printf("is=%d f1=%f f2=%f ",is,f1,f2);
	//f2 = second_deriv(0.01,a*xp*xp*xp+b*xp*xp+c*xp+d,a*xm*xm*xm+b*xm*xm+c*xm+d,a*x*x*x+b*x*x+c*x+d);
	//f2 = 6.;
	//f2 = 6*a*x+2*b;
	f1 = 1+f1*f1;
	f1 = f1*f1*f1;
	kf = f2/(sqrtf(f1));
	printf("kf=%f\n",kf);

	/*zx = a*x*x*x+b*x*x+c*x+d;
	x = x+0.01;
	zxp = a*x*x*x+b*x*x+c*x+d;

	if((zx-zxp)<0) kf*=-1.;*/
	return kf;
}

void getTransmitedRNIPHubralTransmissionLaw(
			   float *rnip, /* RNIP in incident ray layer */
			   float vt, /* Velocity in transmited ray layer */
			   float vi, /* Velocity in incident ray layer */
			   float ei, /* Ray incident angle*/
			   float kf /* Interface curvature */)
/*< Calculate transmited RNIP through interface using Hubral's transmission law
Note: RNIP parameter is modified inside the function to transmited RNIP value
>*/
{
	float et;
	float ri;

	et = getTransmissionAngleSnellsLaw(ei,vi,vt);

	ri = *rnip;
	ri = (1./ri);
	ri *= (cosf(ei)*cosf(ei))/(cosf(et)*cosf(et));
	ri = (vt/vi)*ri;

	ri = ri+ (1./(cosf(et)*cosf(et))*((vt/vi)*cosf(ei)-cosf(et)))*kf;

	*rnip = 1./ri;
}

void transmitedRNIPThroughInterface(
					void *par, /* Raytrace struct */
					void *interface, /* Interface struct */
					int *ir, /* ray sample index */
					int length,
					float vi,
					float vt,
					float **traj, /* ray trajectory */
					float *rnip /* RNIP parameter */)
/*< Calculate transmited RNIP parameter through interface using Hubral laws
Note: Velocity model interpolation makes a transition velocity zone close to
the interface, in this zone is applied transmission law where the ray passes
through interface
>*/
{
	//float vi; // Velocity - incident ray layer
	//float vt=0.; // Velocity - transmited ray layer
	raytrace rt; // Raytrace struct
	itf2d it2; // Interface struct
	float zi; // Interface z coordinate
	float ei; // Incident angle (radians)
	int pass=false; // Ray passed through interface?
	float kf; // Interface curvature
        
	rt = (raytrace) par;
	it2 = (itf2d) interface;

	//vi = getVelocityForRaySampleLocation(rt,traj,*ir);

	//vt = getVelocityForRaySampleLocation(rt,traj,++(*ir));
	
	//while(vi!=vt){
	//	*rnip+=2*vt*rt->dt;

	//	vi = vt;

	//	vt = getVelocityForRaySampleLocation(rt,traj,++(*ir));
		zi = getZCoordinateOfInterface(it2,traj[*ir+length/2][1]);
	//	if(zi>traj[*ir][0] && pass==false){
			ei = calculateIncidentAngle(it2,traj,*ir-1+length/2);
			kf = calculateInterfaceCurvature(it2,traj[*ir+length/2][1]);
			getTransmitedRNIPHubralTransmissionLaw(rnip,vt,vi,ei,kf);
			pass = true;
			//vt = getVelocityForRaySampleLocation(rt,traj,++(*ir));
	//	}
	//}

}

float calculateRNIPWithHubralLaws(
				  void *par, /* Raytrace struct */
				  float** traj, /* Normal ray trajectory */
				  int nt, /* Normal ray times samples */
				  float *v, /* layers velocities */
				  int nv, /* Number of layers */
				  float t0, /* Normal ray traveltime */
				  int itf, /* Interface index */
				  float *sz, /* Splines nodepoints */
				  int nsz, /* Number of nodepoints */
				  float osz, /* Nodepoints origin */
				  float dsz /* Nodepoints sampling */)
/*< Calculate RNIP parameter using Hubral's transmission and propagation laws >*/
{

	int i, j; // loop counter
	float rnip=0.; // RNIP parameter
	raytrace rt; // raytrace struct
	float vt, vi; // velocities
	float *x; // (z,x) position vector
	float* szz; // Z coordinates of interface being inverted
	itf2d interface; // Interface struct
	float zr, zi, vel;
	int length=0, it;
	int pass=0;

	// XXX if nv=1, zero division error!
	if(nv<=1) sf_error("%s: %d: nv can't be 1 or less! It will cause a zero division error!",__FILE__,__LINE__);
	szz = sf_floatalloc(nsz/(nv-1));

	rt = (raytrace) par;
	x = sf_floatalloc(2);

	x[0] = traj[0][0];
	x[1] = traj[0][1];
	vi = sqrtf(1./grid2_vel(rt->grd2,x));
	vt = vi;
	rnip+=2*vt*rt->dt;

	for(i=1;i<nt;i++){
		x[0] = traj[i][0];
		x[1] = traj[i][1];
		vt = sqrtf(1./grid2_vel(rt->grd2,x));

		/* If the ray reaches interface use transmission law */
		if(vt!=vi){
			// TODO: This loop could be outside this if?
			for(j=0;j<nsz/(nv-1);j++){
				szz[j]=sz[j+((itf-1)*nsz/(nv-1))];
			}
			interface = itf2d_init(szz,nsz/(nv-1),osz,dsz);
			//zr=x[0];
			//zi=getZCoordinateOfInterface(interface,x[1]);
			vi=v[itf];
			vt=v[itf-1];

			length=0;
			// Get transition zone length
			for(it=i;it<i+20;it++){
				x[0]=traj[it][0];
				x[1]=traj[it][1];
				vel=sqrt(1./grid2_vel(rt->grd2,x));
				if(vel!=v[itf-1]){
					length++;
				}else{
					break;
				}
			}

			rnip+=2*vt*rt->dt;

			//transmitedRNIPThroughInterface(rt,interface,&i,length,vi,vt,traj,&rnip);
			//sf_warning("pass vi=%f vt=%f length=%d",vi,vt,length); pass++;
			rnip = (vi/vt)*(rnip);
			//sf_warning("RNIP=%f %d vi=%f vt=%f\n",rnip,length,vi,vt);
			vi=vt;
			i+=length;
		}else{

			/* Propagation law */
			rnip+=2*vt*rt->dt;
			//if(itf==1) sf_warning("vt=%f vi=%f",vt,vi);
		}
	}

	//if(itf==1) sf_error("capa");

	return rnip;
}



void forwardModeling(float **s, int ns, float *m0, float *t0, float *a, int *n, float *d, float *o, float *slow, float *BETA, float *RNIP,float *vv, int nv, int itf)
{

	float x[2];
	float p[2];
	raytrace rt;
	float **traj;
	float normalRayAngleRad;
	int nt=10000;
	float dt=0.001;
	int it;
	float t;
	int is;
	float sz[5]={3.,3.,3.,3.,3.};
	int nsz=5;
	float osz=-2;
	float dsz=0.5;

        /* initialize ray tracing object */
        rt = raytrace_init(2,true,nt,dt,n,o,d,slow,ORDER);
        traj = sf_floatalloc2(2,nt+1);

	for(is=0; is<ns; is++){

                normalRayAngleRad = a[is]*DEG2RAD;

                /* Set initial ray point and ray vector */
                setInitialRayPointAndRayVector(s,x,p,is,normalRayAngleRad);
		//printf("z0=%f\n",normalRayAngleRad);

                /* Ray tracing */
                it = trace_ray (rt, x, p, traj);

                if(it>0){ // Ray endpoint at acquisition surface

                        m0[is]= x[1];
                        t0[is] = 2*it*dt;
			//printf("m0=%f\n",x[1]);
			//printf("t0=%f\n",t0[is]);
			//t0[is]=2.1;

                        /* Calculate RNIP */
                        RNIP[is] = calculateRNIPWithHubralLaws(rt,traj,it,vv,nv,2*it*dt,itf,sz,nsz,osz,dsz);
                        //RNIP[is] = calculateRNIPWithDynamicRayTracing(rt,dt,nt,traj,v0);

                        //if(RNIP[is]<0.0)
                          //      sf_warning("ERROR: RNIP=%f",RNIP[is]);

                        /* Calculate BETA */
                        BETA[is] = calculateBetaWithRayTrajectory(x,traj,it);
                }else if(it == 0){ // Ray endpoint inside model
                        t = abs(nt)*dt;
                        rayEndpointWarning(x,p,traj,t);
                }else{ // Side or bottom ray
                        rayEndpointWarning(x,p,traj,t);
                }

        } /* Loop over NIP sources */

	raytrace_close(rt);
        free(traj);
}

static void iso_rhs(void* par, float* y, float* f){}

static int term(void* par, float* y)
/* grid termination */
{
    raytrace rt;
    
    rt = (raytrace) par;
	
    switch (rt->dim) {
		case 2:
			return grid2_term(rt->grd2,y);
		default:
			sf_error("%s: Cannot raytrace with dim=%d",__FILE__,rt->dim);
			return 0;
    }
}

raytrace raytrace_init(int dim            /* dimensionality (2 or 3) */, 
					   bool sym,          /* if symplectic */
					   int nt             /* number of ray tracing steps */, 
					   float dt           /* ray tracing step (in time) */,
					   int* n             /* slowness dimensions [dim] */, 
					   float* o, float* d /* slowness grid [dim] */,
					   float* slow2       /* slowness squared [n3*n2*n1] */, 
					   int order          /* interpolation order */)
/*< Initialize ray tracing object. 
 * Increasing order increases accuracy but
 decreases efficiency. Recommended values: 3 or 4.
 * slow2 can be changed or deallocated after
 raytrace_init.
 >*/
{
    raytrace rt;
    
    rt = (raytrace) sf_alloc (1,sizeof(*rt));
    
    rt->dim = dim;
    rt->sym = sym;
    rt->nt = nt;
    rt->dt = dt;
    rt->z0 = o[0];
    
    switch (dim) {
		case 2:
			rt->grd2 = grid2_init (n[0], o[0], d[0], 
								   n[1], o[1], d[1],
								   slow2, order);
			break;
		default:
			sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
    }
	
    return rt;
}

void raytrace_close (raytrace rt)
/*< Free internal storage >*/
{
    switch (rt->dim) {
		case 2:
			grid2_close (rt->grd2);
			break;
    }
    free (rt);
}

int trace_ray (raytrace rt  /* ray tracing object */, 
			   float* x     /* point location {z,y,x} [dim] */, 
			   float* p     /* ray parameter vector [dim] */, 
			   float** traj /* output ray trajectory [nt+1,dim] */)
/*< Trace a ray.
 * Values of x and p are changed inside the function.
 * The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D
 * Vector p points in the direction of the ray. 
 The length of the vector is not important.
 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*cos(a); p[2] = sin(b)*sin(a) in 3-D
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians
 * The output code for it = trace_ray(...)
 it=0 - ray traced to the end without leaving the grid
 it>0 - ray exited at the top of the grid
 it<0 - ray exited at the side or bottom of the grid
 * The total traveltime along the ray is 
 nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int i, dim, it=0, nt;
    float y[6], s2;
	
    dim = rt->dim;
    nt = rt->nt;
	
    if (!rt->sym) {
		switch (dim) {
			case 2:
				s2 = grid2_vel(rt->grd2,x);
				break;
			default:
				s2 = 0.;
				sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
		}
		
		for (i=0; i < dim; i++) {
			y[i] = x[i];
			y[i+dim] = p[i]*sqrtf(s2);
		}
		
		//sf_runge_init(2*dim, nt, rt->dt);
		//it = sf_ode23_step (y, rt,iso_rhs,term,traj);
		//sf_runge_close();
		
		for (i=0; i < dim; i++) {
			x[i] = y[i];
			p[i] = y[i+dim];
		}
    } else {
		switch (dim) {
			case 2:
				it = atela_step (dim, nt, rt->dt, true, x, p, 
								 rt->grd2, 
								 grid2_vgrad, grid2_vel, grid2_term, traj);
				break;
			default:
				sf_error("%s: cannot handle %d dimensions",__FILE__,rt->dim);
				break;
		}
    }
    
    if (it > 0 && x[0] > rt->z0) {
		return (-it); /* exit through the side or bottom */
    } else {
		return it;
    }
}

#endif
