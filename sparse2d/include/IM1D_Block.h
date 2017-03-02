/******************************************************************************
**                   Copyright (C) 2003 by CEA 
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe Querre
**
**    Date:  10/06/03 
**    
**    File:  IM1D_Block1D.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  : CLASS Definition for 1D block manadgement
**    ----------- 
******************************************************************************/

#ifndef _BLOCK1D_H_
#define _BLOCK1D_H_

//#include "IM1D_Obj.h"
  
/***********************************************************************/

class Block1D { 
  
private:
   int _Nx;
   int _BlockSigSize;      // Sig block size  without taking account
                          // the overlapping aera
   int _BlockSize;         // Sig block size 	taking account
                          // the overlapping aera		  
   int _Nxb;              // Number of blocks in the x directions
   int _NbrBlock;   
   void reset_param() {_Nx=0;_BlockSigSize=_BlockSize=_NbrBlock; 
                       _Nxb=0;_BlockOverlap=False;
		       _Verbose=False;_WeightFirst=False;}   

public:
   // acessors
   Bool _Verbose;
   Bool _BlockOverlap;
   Bool _WeightFirst;
   inline int nx() { return _Nx;}
   inline int nbr_block_nx() { return _Nxb;}
   inline int nbr_block()    { return _NbrBlock;}
   inline int block_size()   { return _BlockSize;}
   inline int l_nx() {return _BlockSize*_Nxb;}
   
   Block1D() {reset_param();} 
   void alloc (int NSig, int ParamBlockSize);
   void get_block_sig (int Bx, fltarray& Sig, fltarray& SigBlock);
   void put_block_sig (int Bx, fltarray& Sig, fltarray& SigBlock);
   float get_weight (int PosPix, int CurrentBlock, int MaxBlockNumber); 
   void add_block_sig (int Bx, fltarray& Sig, fltarray& SigBlock);
   ~Block1D() {reset_param();} 
};

/***********************************************************************/

#endif
