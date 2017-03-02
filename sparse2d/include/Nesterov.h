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
**    File:  Nesterov.h
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**    Nesterov Algorithm
**    		
******************************************************************************/


#ifndef _NESTEROV_H
#define _NESTEROV_H

#include <fstream>
#include "GlobalInc.h"
#include "FloatTrans.h"

extern Bool Verbose;

class Nesterov_params
{
public:
	void reset();			// Initialize the class with default parameters
public:
	char NameOut[256];
	int MaxNiter,MinNiter;	// Maximum/Minimum number of iterations per continuation step
	float Threshold;		// Lambda = Threshold / (Mu/2)
	float TolVar;			// Stopping criterium : ||xt+1 - xt||_infty < TolVar
	float Mu;				// Gradient step, 0<Mu<2/(|degrade| |Domain|)
	int Continuation;		// Number of continuation steps
	bool Positivity;		// The reconstruction must be >0
};
void Nesterov_params::reset()
{
	MaxNiter=100;
	MinNiter=2;
	Threshold = 1;
	TolVar = 1e-6;
	Mu=1;	// Mu < 2 / ( ||degrade*Domain|| )
	Continuation = 1;
	Positivity = false;
}


class Nesterov
{
public:
	Nesterov(FloatTrans *domain);
	~Nesterov(){};
	
	Nesterov_params P;
	FloatTrans* Domain;
	
	void run(cfarray &b, cfarray &z, void (*_degrade)(cfarray &,bool,bool));
	// b : observed signal, in degraded domain
	// z : recovered signal, in degraded domain
	// _degrade : degradation operator. 1bool to apply the degradation, 2bool to transpose
	// _degrade(z,*,true) is assumed real (imaginary part ignored)
};


Nesterov::Nesterov(FloatTrans *domain)
{
	P.reset();
	Domain = domain;
}




void Nesterov::run(cfarray &b, cfarray &z, void (*_degrade)(cfarray &,bool,bool))
{
	char filename[64];
	float *xt, *xit, *wt; // current estimate, gradient, original estimate

// Get Parameters
	float tt = 0.;
	float at,att;
	float gamma, gammat;
	float mu=P.Mu, mu2=P.Mu/2.;
	float threshold = P.Threshold; float threshold0, lambda;
	float TolVar = P.TolVar; float TolVar0;
	fltarray y,zreal;
	
// Initial point
	z = b;
	_degrade(z, false, true); // degraded to direct space, no degradation made
	zreal.alloc(z.nx(),z.ny(),z.nz());
	if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j,k) = z(i,j,k).real();
	else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(j,i) = z(j,i).real();
	Domain->transform(zreal, xit, true);// allocate xit=x0 and localTB
	int n=Domain->size_transform();
	xt = new float[n];
	wt = new float[n];
	for(int k=0;k<n;k++) xt[k] = xit[k] ; // at t=0, xt = x0 = xit
	
// Continuation Initialization
	threshold0 = abs(zreal.maxfabs())*0.1; // threshold0 = 90% of the dynamic of the signal
	TolVar0 = threshold0/100.; // TolVar0 = 1% of threshold0;
	if(P.Continuation>1)
	{
		gamma = pow(threshold/threshold0,float(1./(P.Continuation-1)));
		gammat = pow(TolVar/TolVar0,float(1./(P.Continuation-1)));
		threshold = threshold0/gamma;
		TolVar = TolVar0/gammat;
	}
	else 
	{
		gamma = 1;
		gammat = 1;
	}
	
// Temp variables
	int i; float lvl; bool done;
	float speed, old_speed;
	
//	if(Verbose)
	cerr<<"##########\nBegin Nesterov with mu="<<mu<<", and "<<P.MinNiter<<":"<<P.MaxNiter<<" iterations, with "<<P.Continuation<<" continuations.\n##########"<<endl;
				
// Continuation loop
	for(int t=0;t<P.Continuation;t++)
	{
	// Initialization
		threshold *= gamma;
		TolVar *= gammat;
		tt = 0;
		i=0; done = false; speed=0; old_speed=0;
		if(Verbose) cerr<<"Continuation step "<<t<<", threshold="<<threshold<<", TolVar="<<TolVar<<", MaxNiter="<<P.MaxNiter<<endl;
		
	// Nesterov algorithm. Input : xit=x0, xt=x0
		while( (i<P.MaxNiter && !done) || i<P.MinNiter )
		{
			lambda = threshold/mu2;
			
		// First proximal computation
			lvl = tt*lambda;
			for(int k=0;k<n;k++) wt[k] = xit[k];
			Domain->soft_threshold(wt, lvl);
			
			at = (mu + sqrt(mu*mu+4.*mu*tt))/2. ;
			att = tt+at;
			for(int k=0;k<n;k++) wt[k] = (tt*xt[k] + at*wt[k])/att ;

		// Second proximal computation
			Domain->recons(wt,zreal);
			if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) z(i,j,k) = complex_f(zreal(i,j,k),0.);
			else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) z(j,i) = complex_f(zreal(j,i),0.);
			_degrade(z, true, false);
			z = b - z;
			_degrade(z, true, true);
			if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j,k) = z(i,j,k).real();
			else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(j,i) = z(j,i).real();
			Domain->transform(zreal,xt);
			for(int k=0;k<n;k++) xt[k] *= mu2;
			for(int k=0;k<n;k++) xt[k] += wt[k];
			lvl = mu2*lambda;
			Domain->soft_threshold(xt, lvl);
			
		// Save current solution
			Domain->recons(xt,zreal); 
			if(P.Positivity) for(int k=0;k<zreal.n_elem();k++) zreal(k) = zreal(k) *(zreal(k)>0);
			if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) z(i,j,k) = complex_f(zreal(i,j,k),0.);
			else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) z(j,i) = complex_f(zreal(j,i),0.);
			sprintf(filename,"%s_%02d_%05d.fits",P.NameOut,t,i);fits_write_fltarr(filename, zreal);
			
		// Evaluate the evolution
			y = y-zreal; // error between the last two estimates
			speed = abs(y.maxfabs());
			done = (speed < TolVar) && (i>0) && (speed<old_speed);
			old_speed=speed;
			if(Verbose) cerr<<" Step "<<i<<", || z - zt ||_infty = "<<abs(speed)<<", TolVar="<<TolVar<<endl;
			y=zreal; // Save the new solution
			
		// Accumulate gradient directions
			_degrade(z, true, false);
			z = b - z;
			_degrade(z, true, true);
			if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j,k) = z(i,j,k).real();
			else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(j,i) = z(j,i).real();
			Domain->transform(zreal,wt);
			for(int k=0;k<n;k++) wt[k] *= at;
			for(int k=0;k<n;k++) xit[k] += wt[k];

		// Update the step tt
			tt += at;
			
			i++;
		}
		// Save the current state as the initial step xit
		for(int k=0;k<n;k++) xit[k] = xt[k];
	}
	delete [] xt;
	delete [] xit;
	delete [] wt;
	
	if(Verbose) cerr<<"##########\nEnd Nesterov.\n##########"<<endl;
}


#endif


