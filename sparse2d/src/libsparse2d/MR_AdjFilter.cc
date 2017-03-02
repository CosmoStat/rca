
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Abaque.h"
#include "MR_Noise.h"
#include "MR_NoiseModel.h"
#include "MR_Psupport.h"
#include "MR_Filter.h"

/*********************************************************************/

void mr_grad_adj_filter(Ifloat &Imag, Ifloat &Result, MRNoiseModel &NoiseModel,
                        int Max_Iter, type_border Border, Bool KillLastScale)
{
   int i,j,Iter = 0;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   Bool UseLastScale = (KillLastScale == True) ? False: True;
   MultiResol MR_Data(Nl,Nc, NoiseModel.nbr_scale(), 
                      NoiseModel.type_trans(), "MR_Data");

   Ifloat Resi (Nl, Nc, "Residual");
   Ifloat Sol0 (Nl, Nc, "Sol0");
   Ifloat temp (Nl, Nc, "temp");
   float Sigma, Noise_Ima;
   float Delta;
   float num,den,w;

   Delta = Sigma = 1.e9;
   Noise_Ima=NoiseModel.SigmaNoise;

   MR_Data.transform (Imag, Border);
   NoiseModel.threshold(MR_Data);
   MR_Data.rec_adjoint(Sol0, UseLastScale, Border);

   Result = Sol0;
   Resi = Result;
   Sigma = sigma(Resi);
   if (Noise_Ima < FLOAT_EPSILON)  Noise_Ima = detect_noise_from_med (Resi);
   for (i=0; i< Nl; i++)
   for (j=0; j< Nc; j++)
   {
      temp(i,j) = (Result(i,j) > 0.) ? Result(i,j) : 0.;
   }

   while (True)
   {
        MR_Data.transform (temp, Border);
        NoiseModel.threshold(MR_Data);
        MR_Data.rec_adjoint (temp);
        Resi = Sol0 - temp;

        MR_Data.transform (Resi, Border);
        NoiseModel.threshold(MR_Data);
        MR_Data.rec_adjoint (temp);
        num = 0.;
        den = 0.;
        for (i=0; i< Nl; i++)
        for (j=0; j< Nc; j++) 
        {
            num += Resi(i,j)*temp(i,j);
            den += temp(i,j)*temp(i,j);
        }
        if(!den) w = 1.0;
        else w = MAX (1.0, num/den); 

        for (i=0; i< Nl; i++)
        for (j=0; j< Nc; j++) 
        {
           Result(i,j) +=  w*Resi(i,j);
           temp(i,j) = (Result(i,j) > 0.) ? Result(i,j) : 0.;
        }

        Sigma = sigma(Resi);
        Iter ++;
        cout << "Iter "<< Iter  << ": sigma(resi) = "<< Sigma <<"   convergence coeff  = " <<w<< endl;
        if (Iter >= Max_Iter)  break;
   }
   threshold(Result);
}


/****************************************************************************/

void mr_nmfilter(Ifloat &Imag, Ifloat &Result, MRNoiseModel &NoiseModel, 
                 int Max_Iter, float Epsilon, Bool Sup_Set, type_border Border,
                 Bool PositivIma, Bool KillLastScale)
{
   int Iter = 0;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   MultiResol MR_Data(Nl,Nc, NoiseModel.nbr_scale(), 
                      NoiseModel.type_trans(), "MR_Data");

   Ifloat Resi (Nl, Nc, "Residual");
   float Sigma, Noise_Ima;
   float Old_Sigma, Delta;

   if (NoiseModel.TransImag == True) NoiseModel.im_transform(Imag);

   Delta = Sigma = 1.e9;
   Noise_Ima=NoiseModel.SigmaNoise;
   Result.init();
   Resi = Imag;

   Sigma = energy(Resi);
   Old_Sigma = Sigma; 
   if (Noise_Ima < FLOAT_EPSILON)  Noise_Ima = detect_noise_from_med (Resi);
   while (True)
   {
//cout << "begin loop: " << Noise_Ima << ", " << Sigma << ", " <<  Old_Sigma << endl;    
//cout << "            nrj:" << energy(Resi) << endl;  
//cout << "            nrj(Out):" << energy(Result) << endl;          
//char File[80];
//sprintf(File, "ImagInOld_%d.fits", Iter+1);
//io_write_ima_float (File, Resi);         

        MR_Data.transform (Resi, Border);

        NoiseModel.threshold(MR_Data, Sup_Set);
        // MR_Data.recons (Resi);
        if (KillLastScale == True)
                         MR_Data.scale(MR_Data.nbr_scale()-1).init();
        MR_Data.rec_adjoint (Resi);
        Result += Resi;
	
//sprintf(File, "ResultOld_%d.fits", Iter+1);
//io_write_ima_float (File, Result);     	
	
        if (PositivIma == True) threshold(Result);
	
//sprintf(File, "ThreResultOld_%d.fits", Iter+1);
//io_write_ima_float (File, Result);     		

        Old_Sigma = Sigma; 
        Resi = Imag - Result;
	
//sprintf(File, "ImagOutOld_%d.fits", Iter+1);
//io_write_ima_float (File, Resi);     	
	
        Sigma = energy(Resi);
        Delta = (Old_Sigma - Sigma)  / (Nl*Nc*Noise_Ima*Noise_Ima);
        Iter ++;
//io_write_ima_float ("ResiPfilt.fits", Resi);         
        float Xi = Sigma / (Nl*Nc*Noise_Ima*Noise_Ima);
        if (NoiseModel.which_noise() == NOISE_GAUSSIAN) 
             cout << "Iter "<< Iter << ": Xi = "<< Xi <<" ==> Delta  = " <<Delta<< endl;
        else cout << "Iter "<< Iter  <<" ==> Delta  = " <<Delta<< endl;
//cout << "end loop: " << Noise_Ima << ", " << Sigma << ", " <<  Old_Sigma << endl; 	
//cout << "          nrj(In) :" << energy(Resi) << endl;   
//cout << "          nrj(Out):" << energy(Result) << endl;          
	
	if ((Iter >= Max_Iter) || (Delta < Epsilon)) break;
   }

   if (NoiseModel.TransImag == True) NoiseModel.im_invtransform(Imag);
   if (NoiseModel.TransImag == True) NoiseModel.im_invtransform(Result);
}

/********************************************************************/ 

 
