//
// Floating version of the FCur, adapted to FloatTrans class, for
// algorithms like Fista/Nesterov
// This transform is fully normalized
//

#ifndef _SB_FILTER_FLOAT_H
#define _SB_FILTER_FLOAT_H

#include "SB_Filter.h"
#include "FloatTrans.h"

// ***********************************************************************
// *   Orthogonal 2D wavelet transform
// ***********************************************************************
 
class Ortho_2D_WT_float: public Ortho_2D_WT, public FloatTrans
{
	int DataNx;
	int DataNy;
	int NbrScale;
public:
	Ortho_2D_WT_float(SubBand1D &SB1D):Ortho_2D_WT(SB1D),DataNx(0),DataNy(0),NbrScale(0){}
	~Ortho_2D_WT_float(){}
	
	void alloc(fltarray& in, int _NbrScale);
	
// FloatTrans methods
	int size_transform();
	void transform(fltarray& in, float* &out, bool alloc=false);
	void recons(float* in, fltarray& out);
	void soft_threshold(float* in, float lvl, bool threshold_coarse=false);
	void substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse=false);
	void add_coarse_scale(float* in, float* coarse);
};


// ***********************************************************************
// *   Orthogonal Undecimated 2D wavelet transform
// ***********************************************************************

class PAVE_2D_WT_float: public PAVE_2D_WT, public FloatTrans
{
	int DataNx;
	int DataNy;
	int NbrScale;
public:
	PAVE_2D_WT_float(SubBand1D &SB1D):PAVE_2D_WT(SB1D),DataNx(0),DataNy(0),NbrScale(0){}
	~PAVE_2D_WT_float(){}
	
	void alloc(fltarray& in, int _NbrScale);
	
// FloatTrans methods
	int size_transform();
	void transform(fltarray& in, float* &out, bool alloc=false);
	void recons(float* in, fltarray& out);
	void soft_threshold(float* in, float lvl, bool threshold_coarse=false);
	void substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse=false);
	void add_coarse_scale(float* in, float* coarse);
};


/************************************************************************/
//    Isotropic a trous wavelet transform
/************************************************************************/

class ATROUS_2D_WT_float: public ATROUS_2D_WT, public FloatTrans
{
	int DataNx;
	int DataNy;
	int NbrScale;
	bool positive_coef;
public:
	ATROUS_2D_WT_float():ATROUS_2D_WT(),DataNx(0),DataNy(0),NbrScale(0),positive_coef(false){}
	~ATROUS_2D_WT_float(){}
	
	void alloc(fltarray& in, int _NbrScale);
	void set_positive_coef(bool pc){positive_coef = pc;}
// FloatTrans methods
	int size_transform();
	void transform(fltarray& in, float* &out, bool alloc=false);
	void recons(float* in, fltarray& out);
	void soft_threshold(float* in, float lvl, bool threshold_coarse=false);
	void substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse=false);
	void add_coarse_scale(float* in, float* coarse);
};

#endif
