/******************************************************************************
**                   Copyright (C) 2009 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  26/10/2009
**    
**    File:  FloatTrans.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION:  FloatTrans class 
**	  -----------	A purely virtual class defining the basis methods needed
**					for a general itetrative thresholding algorithm to work.
**					Needs a class/transform dealing with float* for coefficients
**					and fltarray for the data.
**					BFCurvelet3D, UOWT_float, and FCur_float derive from it.
**					Nesterov and Fista use this class
**					
******************************************************************************/


#ifndef _FLOATTRANS_H
#define _FLOATTRANS_H

#include <fstream>
#include "GlobalInc.h"

class FloatTrans
{
public:
	virtual int size_transform() = 0;
	// returns the number of elements of the transform
	
	virtual void transform(fltarray& in, float* &out, bool alloc=false) = 0; 
	// transforms the data in into the vector out[size_transform]
	
	virtual void recons(float* in, fltarray& out) = 0;
	// reconstructs a signal out from the coefficients in[size_transform]
	
	virtual void soft_threshold(float* in, float lvl, bool threshold_coarse=false) = 0;
	// soft_threshold the coefficients in at level lvl
	

// The two following functions are not always necessary and need not always be defined in the derived class
	virtual void substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse=false){};
	// saves the values of in[size_transform] which correspond to the coarse scale into coarse[unknown]
	// in is set to 0 on the coarse coefficients
	// the length of coarse is unknown by this virtual class, 
	//  the elements of the coarse scale are not necessarily contiguous
	
	virtual void add_coarse_scale(float* in, float* coarse){};
	// adds the values of the coarse scale to those inside in
};

#endif


