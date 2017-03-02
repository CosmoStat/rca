/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.3
**
**    Author: 02/29/00
**
**    Date:  00/07/13
**    
**    File:  mc2d_com.cc
**
*******************************************************************************/
#include "mr_com.h"








/******************************************************************************
classe to_Iter
******************************************************************************/

/*
void to_Iter::iter_Begin () {

  //cout << " ==> BEGIN" << endl;

  o_Residu.init(iter_getpResult()->res_GetpParam()->o_Imag);
  
  o_Init.init(iter_getpResult()->res_GetpParam()->o_Imag); 
  
  NbIter = 0;

  o_Init = iter_getpResult()->res_GetpParam()->o_Imag; // Einit
  o_Residu.init(0.0);  
                                         // I(0) = 0
//char File[80];
//sprintf(File, "ImagIn_%d.fits", NbIter+1);
//io_write_ima_float (File, iter_getpResult()->res_GetpParam()->o_Imag);         
  
  NbIter++;
  //cout << "Iter number : " << NbIter << endl;
  iter_getpResult()->res_FirstCompute ();  
  iter_getpResult()->res_Recons();
  iter_Residu();
}
   
Bool to_Iter::iter_End () {

   //cout << " ==> END" << endl;
   
   NbIter++;
      
   if (NbIter > iter_getpResult()->res_GetpParam()->i_NbIter) { 
      iter_getpResult()->res_GetpParam()->o_Imag = o_Residu;             
      return False;
   } else return True;

}

void to_Iter::operator++ () {(*this)++;}
void to_Iter::operator++ (int) {

   //cout << " ==> OP ++ " << endl;
   //cout << "Iter number : " << NbIter << endl;
//char File[80];
//sprintf(File, "ImagIn_%d.fits", NbIter);
//io_write_ima_float (File, iter_getpResult()->res_GetpParam()->o_Imag);            
   iter_getpResult()->res_CurrentCompute ();
   iter_getpResult()->res_Recons ();
   iter_Residu ();
}

void to_Iter::iter_Residu () {

   // I(i) = I(i) + E(i)
   o_Residu += iter_getpResult()->res_GetpParam()->o_Imag;
	                
   // E(i) = Einit - I(i)
   iter_getpResult()->res_GetpParam()->o_Imag = o_Init - o_Residu;
}
 
 */
 
 
 /******************************************************************************
classe to_SoftIter
******************************************************************************/
/*
 
void to_SoftIter::iter_Begin () {

  //cout << " ==> BEGIN" << endl;
  NbIter = 0;                                        
   
  NbIter++;
  
  //cout << "Iter number : " << NbIter << endl;
  iter_getpResult()->res_FirstCompute ();  
  iter_getpResult()->res_Recons();
  iter_ActionOnSol();
}


Bool to_SoftIter::iter_End () {

   //cout << " ==> END" << endl;
   
   NbIter++;
      
   if (NbIter > iter_getpResult()->res_GetpParam()->i_NbIter) {        
      return False;
   } else return True;
}
 
void to_SoftIter::operator++ () {(*this)++;}
void to_SoftIter::operator++ (int) {

   //cout << " ==> OP ++ " << endl;
   //cout << "Iter number : " << NbIter << endl;
   iter_getpResult()->res_CurrentCompute ();
   iter_getpResult()->res_Recons ();
   iter_ActionOnSol ();
}

 
 
 */
 
 
