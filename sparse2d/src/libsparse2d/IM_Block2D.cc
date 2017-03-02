/******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  28/03/2003 
**    
**    File:  IM_Block2D.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION 2D block manadgement
**    -----------  
**                 
******************************************************************************/

#include "IM_Block2D.h"

/***********************************************************************/

void Block2D::alloc(int Nlc, int Ncc, int ParamBlockSize)
{
    int MaxXY;
    Nl=Nlc;
    Nc=Ncc;
    BlockSize = ParamBlockSize;
    
    MaxXY = MAX(Nl,Nc);
        
    if (BlockSize==0) 
    {
       Nlb = Ncb = 1;
       BlockSize = MIN(Nl,Nc);
       BlockOverlap=False;
    }
    else if (BlockSize==MaxXY)
    {
       Nlb = Ncb = 1;
       BlockOverlap=False;
    }
    else
    {
       // cout << "BlockSize = " << BlockSize << endl;
       if (BlockOverlap == True) BlockImaSize = (BlockSize+1)/2;
       else BlockImaSize = BlockSize;

       // cout << "BlockImaSize = " << BlockImaSize << endl;
       Nlb =  (Nl %  BlockImaSize == 0) ? Nl / BlockImaSize: Nl / BlockImaSize + 1;
       Ncb =  (Nc %  BlockImaSize == 0) ? Nc / BlockImaSize: Nc / BlockImaSize + 1;
    }
    
    NbrBlock = Nlb*Ncb;
    if (BlockOverlap == False) BlockImaSize = 0;
    if (Verbose == True)
    {
       printf( "\nima size = (%2d,%2d), BlockSize = %2d\n", Nl,Nc,BlockSize);
       printf( "Nbr Blocks = (%2d,%2d)\n", Nlb,Ncb);
       if (BlockOverlap == True) printf( "Block Overlap: size = %d \n", BlockImaSize);
       else printf( "No Block Overlap \n\n");
    }
}

/***********************************************************************/

void Block2D::get_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock)
// Extract the block (Bi,Bj) from Ima and put it in ImaBlock
{
   int k,l;
   int Depi = (BlockSize-BlockImaSize)*Bi;
   int Depj = (BlockSize-BlockImaSize)*Bj;
   type_border Bord = I_MIRROR;
   
   if ((BlockOverlap == False) || (WeightFirst == False))
   {
      // cout << "get_block_ima NO Weight " << endl;
      for (k = 0; k < BlockSize; k++)
      for (l = 0; l < BlockSize; l++) ImaBlock(k,l) = Ima(Depi+k,Depj+l,Bord);  
   }
   else
   {
      // cout << "get_block_ima  Weight " << endl;
      for (k = 0; k < BlockSize; k++)
      for (l = 0; l < BlockSize; l++)
           ImaBlock(k,l) = Ima(Depi+k,Depj+l,Bord) * get_weight(k, Bi, Nlb) 
	                     * get_weight(l, Bj, Ncb);  
   }
   // cout << "BLOCK: " << Bi <<","<<Bj<<","<<Bk<<": Dep = "<< Depi<<","<<Depj<<","<<Depk<<endl;
}

/*************************************************************************/

void Block2D::put_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock)
// Extract the block (Bi,Bj,Bk) from ima and put it in imaBlock
{
   int k,l;
   int Depi = (BlockSize-BlockImaSize)*Bi;
   int Depj = (BlockSize-BlockImaSize)*Bj;
 
   for (k = 0; k < BlockSize; k++)
   for (l = 0; l < BlockSize; l++)
          if ((Depi+k >= 0) && (Depi+k < Nl) 
             && (Depj+l >= 0) && (Depj+l < Nc))
                               Ima(Depi+k,Depj+l) = ImaBlock(k,l);   
   // cout << "BLOCK: " << Bi <<","<<Bj<<","<<Bk<<": Dep = "<< Depi<<","<<Depj<<","<<Depk<<endl;
}

/*************************************************************************/

float Block2D::get_weight(int PosPix, int CurrentBlock, int MaxBlockNumber) 
{
    float ValReturn=1.;
    if ((BlockSize % 2 == 0) || 
        ((PosPix != 0)  && (PosPix != BlockSize-1))) 
    {
       int Center = BlockSize / 2;
       if ((CurrentBlock != 0)  &&  (PosPix < Center)) 
           ValReturn = (float) pow(sin((double) (PosPix) / (double) Center*PI/2.), 2.);
       else 
	if ((CurrentBlock != MaxBlockNumber-1)  &&  (PosPix > Center)) 
           ValReturn = (float) pow(cos((double)(PosPix-Center)/ (double) Center *PI/2.), 2.);
    } 
    else 
    {
       if  ((CurrentBlock != 0) && (PosPix==0))  ValReturn=0.;
       else if ((CurrentBlock != MaxBlockNumber-1) && (PosPix== BlockSize-1)) 
	   ValReturn=0.; 
    }
    return ValReturn;
}

/*************************************************************************/

void Block2D::add_block_ima(int Bi, int Bj, Ifloat &Ima, Ifloat &ImaBlock)
// Add the block (Bi,Bj) ImaBlock in Ima  with weighted values
{
   int k,l;
   int Depi = (BlockSize-BlockImaSize)*Bi;
   int Depj = (BlockSize-BlockImaSize)*Bj;
    
   if ((BlockOverlap == True) && (WeightFirst == False))
   {
     // cout << "add_block_ima Weight " << endl;
     // if (BlockOverlap == True) cout << "BlockOverlap True" << endl;
     // if (WeightFirst == False) cout << "WeightFirst False" << endl;
      for (k = 0; k < BlockSize; k++)
      for (l = 0; l < BlockSize; l++)
          if ((Depi+k >= 0) && (Depi+k < Nl) 
             && (Depj+l >= 0) && (Depj+l < Nc))
                   Ima(Depi+k,Depj+l) += ImaBlock(k,l)*
	             get_weight(k, Bi, Nlb) *  get_weight(l, Bj, Ncb);  
   }
   else
   {
      // cout << "add_block_ima NO Weight " << endl;
      for (k = 0; k < BlockSize; k++)
      for (l = 0; l < BlockSize; l++)
          if ((Depi+k >= 0) && (Depi+k < Nl) 
             && (Depj+l >= 0) && (Depj+l < Nc))
                   Ima(Depi+k,Depj+l) += ImaBlock(k,l);  
   }
}

/****************************************************************************/
