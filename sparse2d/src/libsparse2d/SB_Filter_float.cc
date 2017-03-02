#include "SB_Filter_float.h"


// ***********************************************************************
// *   Orthogonal 2D wavelet transform
// ***********************************************************************
 
void Ortho_2D_WT_float::alloc(fltarray& in, int _NbrScale)
{
	DataNx = in.nx();
	DataNy = in.ny();
	NbrScale = _NbrScale;
}

int Ortho_2D_WT_float::size_transform()
{	
	return DataNx*DataNy;
}

void Ortho_2D_WT_float::transform(fltarray& in, float* &out, bool alloc)
{
	if(alloc) out = new float [size_transform()];
	Ifloat f_in; f_in.alloc(in.buffer(),DataNy,DataNx);
	Ifloat f_out; f_out.alloc(out,DataNy,DataNx);
	Ortho_2D_WT::transform (f_in, f_out, NbrScale);
//	io_write_ima_float("out.fits",f_out);
}
void Ortho_2D_WT_float::recons(float* in, fltarray& out)
{
	out.resize(DataNy,DataNx);
	Ifloat f_in; f_in.alloc(in,DataNy,DataNx);
	Ifloat f_out; f_out.alloc(out.buffer(),DataNy,DataNx);
	Ortho_2D_WT::recons(f_in, f_out, NbrScale);
}
void Ortho_2D_WT_float::soft_threshold(float* in, float lvl, bool threshold_coarse)
{
	// threshold_coarse=false not implemented yet
	for(int i=0;i<size_transform();i++)
		in[i] = ::soft_threshold(in[i], lvl);
}
void Ortho_2D_WT_float::substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse)
{
	
}
void Ortho_2D_WT_float::add_coarse_scale(float* in, float* coarse)
{
	
}



// ***********************************************************************
// *   Orthogonal Undecimated 2D wavelet transform
// ***********************************************************************

void PAVE_2D_WT_float::alloc(fltarray& in, int _NbrScale)
{
	DataNx = in.nx();
	DataNy = in.ny();
	NbrScale = _NbrScale;
}

int PAVE_2D_WT_float::size_transform()
{
	return DataNx*DataNy*(3*NbrScale-2);
}

void PAVE_2D_WT_float::transform(fltarray& in, float* &out, bool alloc)
{
	if(alloc) out = new float [size_transform()];
	int nplan = 3*NbrScale-2;
	Ifloat f_in; f_in.alloc(in.buffer(),DataNy,DataNx);
	Ifloat * tab_f_out = new Ifloat[nplan];
	for(int i=0;i<nplan;i++)
		tab_f_out[i].alloc(out+i*DataNy*DataNx,DataNy,DataNx);
	
	PAVE_2D_WT::transform(f_in, tab_f_out, NbrScale);
	delete [] tab_f_out;
}
void PAVE_2D_WT_float::recons(float* in, fltarray& out)
{
	int nplan = 3*NbrScale-2;
	out.resize(DataNy,DataNx);
	Ifloat * tab_f_in = new Ifloat[nplan];
	for(int i=0;i<nplan;i++)
		tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
	Ifloat f_out; f_out.alloc(out.buffer(),DataNy,DataNx);
	
	PAVE_2D_WT::recons(tab_f_in, f_out, NbrScale);
	delete [] tab_f_in;
}
void PAVE_2D_WT_float::soft_threshold(float* in, float lvl, bool threshold_coarse)
{
	int nplan = 3*NbrScale-2;
	for(int i=0;i<(nplan-int(!threshold_coarse))*DataNx*DataNy;i++)
		in[i] = ::soft_threshold(in[i], lvl);
}
void PAVE_2D_WT_float::substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse)
{
	int start = 3*(NbrScale-1)*DataNx*DataNy;
	if(alloc_coarse) coarse = new float[DataNx*DataNy];
	for(int i=start;i<size_transform();i++)
	{
		coarse[i-start] = in[i];
		in[i] = 0.;
	}
}
void PAVE_2D_WT_float::add_coarse_scale(float* in, float* coarse)
{
	int start = 3*(NbrScale-1)*DataNx*DataNy;
	for(int i=start;i<size_transform();i++)
		in[i] += coarse[i-start];
}

		


// ***********************************************************************
// *   Isotropic a trous wavelet transform
// ***********************************************************************

void ATROUS_2D_WT_float::alloc(fltarray& in, int _NbrScale)
{
	DataNx = in.nx();
	DataNy = in.ny();
	NbrScale = _NbrScale;
}

int ATROUS_2D_WT_float::size_transform()
{
	return DataNx*DataNy*NbrScale;
}

void ATROUS_2D_WT_float::transform(fltarray& in, float* &out, bool alloc)
{
	if(alloc) out = new float [size_transform()];
	Ifloat f_in; f_in.alloc(in.buffer(),DataNy,DataNx);
	Ifloat * tab_f_out = new Ifloat[NbrScale];
	for(int i=0;i<NbrScale;i++)
		tab_f_out[i].alloc(out+i*DataNy*DataNx,DataNy,DataNx);
	
	ATROUS_2D_WT::transform(f_in, tab_f_out, NbrScale);

// Normalization	
	for(int i=0;i<size_transform();i++)
	{
		float norm = float(norm_band(int(floor(i/(DataNx*DataNy)))));
		out[i] /= norm;
	}
//for(int s=0;s<NbrScale;s++)
//cerr<<tab_f_out[s].sigma()<<","<<endl;
	
	delete [] tab_f_out;
}
void ATROUS_2D_WT_float::recons(float* in, fltarray& out)
{
	out.resize(DataNy,DataNx);
	Ifloat * tab_f_in = new Ifloat[NbrScale];
	for(int i=0;i<NbrScale;i++)
		tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
	Ifloat f_out; f_out.alloc(out.buffer(),DataNy,DataNx);
	
//cerr<<"RECONS "<<endl;
//for(int s=0;s<NbrScale;s++)
//cerr<<tab_f_in[s].sigma()<<","<<endl;

// Normalization	
	for(int i=0;i<size_transform();i++)
	{
		float norm = float(norm_band(int(floor(i/(DataNx*DataNy)))));
		in[i] *= norm;
	}
//for(int s=0;s<NbrScale;s++)
//cerr<<tab_f_in[s].sigma()<<","<<endl;
	
	ATROUS_2D_WT::recons(tab_f_in, f_out, NbrScale);

// re-normalization	(because in may be used again ouside this function)
	for(int i=0;i<size_transform();i++)
	{
		float norm = float(norm_band(int(floor(i/(DataNx*DataNy)))));
		in[i] /= norm;
	}
	delete [] tab_f_in;
}
void ATROUS_2D_WT_float::soft_threshold(float* in, float lvl, bool threshold_coarse)
{
	// There seems to be problems when thresholding the coarse scale
	threshold_coarse=false;
	for(int i=0;i<(NbrScale-int(!threshold_coarse))*DataNx*DataNy;i++)
		in[i] = ::soft_threshold(in[i], lvl);
	if(positive_coef)
	for(int i=0;i<size_transform();i++)
		if(in[i]<0) in[i] = 0;
}
void ATROUS_2D_WT_float::substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse)
{
//	int start = (NbrScale-1)*DataNx*DataNy;
//	if(alloc_coarse) coarse = new float[DataNx*DataNy];
//	for(int i=start;i<size_transform();i++)
//	{
//		coarse[i-start] = in[i];
//		in[i] = 0.;
//	}
}
void ATROUS_2D_WT_float::add_coarse_scale(float* in, float* coarse)
{
//	int start = (NbrScale-1)*DataNx*DataNy;
//	for(int i=start;i<size_transform();i++)
//		in[i] += coarse[i-start];
}

		
