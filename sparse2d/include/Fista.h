/******************************************************************************
**                   Copyright (C) 2010 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  11 Jan. 2010
**    
**    File:  Fista.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**    Fista Algorithm
**    		
******************************************************************************/


#ifndef _FISTA_H
#define _FISTA_H

#include <fstream>
#include "GlobalInc.h"
#include "FloatTrans.h"

extern Bool Verbose;

class Fista_params
{
public:
	void reset();			// Initialise the class with default parameters
public:
	char NameOut[256];
	int MaxNiter;			// Maximum number of iterations
	float Threshold;		// Thresholding level
	float TolVar;			// Stopping criterium : ||xt+1 - xt||_infty < TolVar
	float Mu;				// Gradient regularization parameter
	bool Fast;				// (F)ISTA
	bool Decreasing;		// if true : linearily decreasing threshold
	bool Positivity;		// The reconstruction must be >0
	bool No_coarse;			// if true : kills the coarse scale before apaplying fista, then thresholds it
};
void Fista_params::reset()
{
	MaxNiter=10;
	Threshold = 1;
	TolVar = 1e-6;
	Mu=1;
	Fast = false;
	Decreasing = false;
	Positivity = false;
	No_coarse = false;
}


class Fista
{
public:
	Fista(FloatTrans *domain);
	~Fista(){};
	
	Fista_params P;
	FloatTrans* Domain;
	
	void run(cfarray &b, cfarray &z, void (*_degrade)(cfarray &,bool,bool));
	// b : observed signal, in degraded domain
	// z : recovered signal, in degraded domain
	// _degrade : degradation operator. 1bool to apply the degradation, 2bool to transpose
	// _degrade(z,*,true) is assumed real (imaginary part ignored)
};


Fista::Fista(FloatTrans *domain)
{
	P.reset();
	Domain = domain;
}



void Fista::run(cfarray &b, cfarray &z, void (*_degrade)(cfarray&,bool,bool))
{
	char filename[64];
	float *x, *xold, *xtmp; // current estimate, previous estimate, and temp variable
	fltarray zreal;
	
// Initialization
	float mu = P.Mu;
	float lvl = P.Threshold;
	float TolVar = P.TolVar;

// Other variables
	float tk=1,told;
	int iter; bool done;
	float speed, old_speed;
	fltarray yreal;// previous reconstruction
	iter=0; done = false; speed=0; old_speed=0;
	bool allocD=true; // Domain allocation
	
// Initial point
	float* CS;
	if(P.No_coarse)
	{
		_degrade(b, false, true);// degraded to direct space, no degradation made
		zreal.alloc(b.nx(),b.ny(),b.nz());
		if(b.nz()) for(int k=0;k<b.nz();k++) for(int j=0;j<b.ny();j++) for(int i=0;i<b.nx();i++) zreal(i,j,k) = b(i,j,k).real();
		else for(int j=0;j<b.ny();j++) for(int i=0;i<b.nx();i++) zreal(j,i) = b(j,i).real();
		Domain->transform(zreal, xtmp, true); allocD=false; // allocate xtmp and localTB
		Domain->substract_coarse_scale(xtmp,CS,true);
		//Domain->add_coarse_scale(xtmp,CS);
		Domain->recons(xtmp, zreal);// allocate xtmp and localTB
		if(b.nz()) for(int k=0;k<b.nz();k++) for(int j=0;j<b.ny();j++) for(int i=0;i<b.nx();i++) b(i,j,k) = complex_f(zreal(i,j,k),0.);
		else for(int j=0;j<b.ny();j++) for(int i=0;i<b.nx();i++) b(j,i) = complex_f(zreal(j,i),0.);
		_degrade(b, false, false);// direct to degraded space, no degradation made
	}
	z = b; // degraded space
	_degrade(z, false, true); // degraded to direct space, no degradation made
	zreal.alloc(z.nx(),z.ny(),z.nz());
	if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j,k) = z(i,j,k).real();
	else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(j,i) = z(j,i).real();
	Domain->transform(zreal, xtmp, allocD);
	int n=Domain->size_transform();
	x = new float[n]; for(int k=0;k<n;k++) x[k] = xtmp[k];
	xold = new float[n];

	cerr<<"##########\nBegin ";
	if(P.Decreasing) cerr<<"MCA";
	else if(P.Fast) cerr<<"FISTA";
	else cerr<<"ISTA";
	cerr<<" with mu="<<mu<<", and "<<P.MaxNiter<<" iterations.\n##########"<<endl;

// Fista algorithm. 
	while( iter<P.MaxNiter && (!done || P.Decreasing) )
	{
	// Threshold update
		if(P.Decreasing) 
		{
			//lvl = P.Threshold * (P.Niter-i-1.)/(P.Niter-1.);
			lvl = 1-iter*1.0/float(P.MaxNiter-1);			// 1..0
			lvl = (pow((double)lvl,3)+lvl/25.)/1.04;
			lvl = P.Threshold*(1+ lvl*(10-1));
		}
		
	// Save current solution
		if(iter!=0) Domain->recons(xtmp,zreal);
		if(P.Positivity) for(int k=0;k<zreal.n_elem();k++) zreal(k) = zreal(k) *(zreal(k)>0);
		if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) z(i,j,k) = complex_f(zreal(i,j,k),0.);
		else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) z(j,i) = complex_f(zreal(j,i),0.);
		sprintf(filename,"%s_%05d.fits",P.NameOut,iter);fits_write_fltarr(filename, zreal);
		
	// Evaluate the evolution
		yreal = yreal - zreal; // error between the last two estimates
		speed = abs(yreal.maxfabs());
		done = (speed < TolVar) && (iter>0) && (speed<old_speed);
		old_speed=speed;
		if(Verbose) cerr<<" Step "<<iter<<", lvl="<<lvl<<", || z - zt ||_infty = "<<abs(speed)<<", TolVar="<<TolVar<<endl;
		yreal = zreal; // Save the new solution
		
	// Gradient step
		_degrade(z, true, false);// direct to degraded space
		z = b - z ;// degraded space
		_degrade(z,true,true); // degraded to direct space
		if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j,k) = z(i,j,k).real();
		else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(j,i) = z(j,i).real();
		Domain->transform(zreal, xold);
		for(int k=0;k<n;k++) xtmp[k] = xtmp[k] + mu*xold[k];
		
	// Save estimate
		if(P.Fast) for(int k=0;k<n;k++) xold[k] = x[k];
		
	// Proximal operator : soft thresholding
		for(int k=0;k<n;k++) x[k] = xtmp[k];
		Domain->soft_threshold(x,lvl,P.No_coarse);
		
	// New point
		if(P.Fast) 
		{
			told=tk;
			tk = (1.+sqrt(1.+4.*told*told))/2.;
			for(int k=0;k<n;k++) xtmp[k] = x[k] + (told-1)/tk * (x[k]-xold[k]);
		}
		else
			for(int k=0;k<n;k++) xtmp[k] = x[k];
		
		iter++;
	}
	if(P.No_coarse)
		Domain->add_coarse_scale(xtmp,CS);
	Domain->recons(x,zreal);
	if(P.Positivity) for(int k=0;k<zreal.n_elem();k++) zreal(k) = zreal(k) *(zreal(k)>0);
	sprintf(filename,"%s_recons.fits",P.NameOut);fits_write_fltarr(filename, zreal);
	cerr<<"##########\nEnd.\n##########"<<endl;
}


#endif


