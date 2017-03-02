/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/02
**    
**    File:  IM_OperInt.cc
**
*******************************************************************************
**
**    DESCRIPTION 
**    -----------  
**                 
*******************************************************************************
**
** double sigma(Iint & Image)
**
** returns the standard deviation
**
*******************************************************************************
** 
** double average(Iint & Image)
** 
** returns the average
**
*******************************************************************************
**
** int min(Iint & Image)
** 
** returns the min
**
*******************************************************************************
**
** int max(Iint & Image)
** 
** returns the max
**
********************************************************************************
**
** int flux(Iint & Image, int N)
** 
** returns the flux
**
********************************************************************************
**
** int energy(Iint & Image, int N)
** 
** returns the energy
**
*******************************************************************************
**
** const Iint  operator +(const Iint &Im1, const Iint &Im2)
**
** returns the addition of two images
**
*******************************************************************************
**
** const Iint  operator -(const Iint &Im1, const Iint &Im2)
**
** returns the soustraction of two images
**
*******************************************************************************
**
** const Iint operator *(const Iint &Im1, const Iint &Im2)
**
** returns the multiplication of two images
**
*******************************************************************************
**
** const Iint  operator /(const Iint &Im1, const Iint &Im2)
**
** returns the division of two images
** 
*******************************************************************************
**
** const Iint conv_direct(const Iint &Im1, const Iint &Im2, 
**                     type_border Type_Border)
**
** returns the convolution product  Im1 * Im2
**
*******************************************************************************
**
** float lib_mat_correl (const Iint & Tab_Im1, const Iint & Tab_Im2)
** 
** returns the corelation between two images
**
*******************************************************************************
**
** float detect_noise_sigma (const Iint &Image, Bool Average_Non_Null, int Nit)
**
** returns an estimation of the standard deviation of the noise in Image
**
******************************************************************************/

// static char sccsid[] = "@(#)IM_OperInt.cc 3.1 96/05/02 CEA 1994 @(#)";


#include "IM_Obj.h"

/******************************************************************************/

double sigma(const Iint &Image)
{
    int i,j;
    double Sigma, Moy;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Moy = Sigma = 0.;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        Moy += Image(i,j);
        Sigma += Image(i,j)*Image(i,j);
    }
    Moy /= (float) (Nl*Nc);
    Sigma /= (float) (Nl*Nc);
    Sigma = sqrt(Sigma - Moy * Moy);
    return Sigma;
};

/******************************************************************************/

double average(const Iint &Image)
{
    int i,j;
    double Moy;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Moy = 0.;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        Moy += Image(i,j);
    }
    Moy /= (float) (Nl*Nc);
    return Moy;
}

/******************************************************************************/

int min (const Iint &Image) 
{
    int i,j;
    int Val;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Val = Image(0,0);
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
       if (Val > Image(i,j)) Val = Image(i,j);
    return Val;
};

/******************************************************************************/

int max (const Iint &Image)
{
    int i,j;
    int Val;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Val = Image(0,0);
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
       if (Val < Image(i,j)) Val = Image(i,j);
    return Val;
};

/******************************************************************************/

int flux (const Iint &Image)
{
    int i,j;
    int Val;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Val = Image(0,0);
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
              Val += Image(i,j);
    return Val;
};

/******************************************************************************/

int energy (const Iint &Image)
{
    int i,j;
    int Val;
    int Nl = Image.nl();
    int Nc = Image.nc();

    Val = Image(0,0);
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
            Val += Image(i,j)*Image(i,j);
    return Val;
};

/******************************************************************************/


const Iint  operator +(const Iint &Im1, const Iint &Im2)
{
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     Iint *Result = new Iint(Nl,Nc, "op+");
     int i,j;

     for (i=0; i < Nl; i++) {
         for (j=0; j < Nc; j++) {
             (*Result)(i,j) = Im1(i,j) + Im2(i,j);
         }
     }
     return (*Result);
}

/****************************************************************************/

const Iint   operator -(const Iint &Im1, const Iint &Im2)
{
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     Iint *Result = new Iint(Nl,Nc, "op-");
     int i,j;

     for (i=0; i < Nl; i++) {
         for (j=0; j < Nc; j++) {
             (*Result)(i,j) = Im1(i,j) - Im2(i,j);
         }
     }
     return (*Result);
}

/****************************************************************************/

const Iint operator *(const Iint &Im1, const Iint &Im2)
{
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     Iint *Result = new Iint(Nl,Nc, "op*");
     int i,j;

     for (i=0; i < Nl; i++) {
         for (j=0; j < Nc; j++) {
             (*Result)(i,j) = Im1(i,j) * Im2(i,j);
         }
     }
     return (*Result);
}

/****************************************************************************/

const Iint  operator /(const Iint &Im1, const Iint &Im2)
{
     int Nl = Im1.nl();
     int Nc = Im1.nc();
     Iint *Result = new Iint(Nl,Nc, "op/");
     int i,j;

     for (i=0; i < Nl; i++) {
         for (j=0; j < Nc; j++) {
             if ( ABS(Im2(i,j)) > FLOAT_EPSILON ) 
                  (*Result)(i,j) = Im1(i,j) / Im2(i,j);
             else (*Result)(i,j) = 0;
         }
     }
     return (*Result);
}

/****************************************************************************/

const Iint conv_direct(const Iint &Im1, const Iint &Im2, type_border Type_Border)
{
    int Nl1 = Im1.nl();
    int Nc1 = Im1.nc();
    int Nl2 = Im2.nl();
    int Nc2 = Im2.nc();
    int Nl_2,Nc_2;
    int i,j,k,l;
    int ii,jj;

    Iint *Result = new Iint(Nl1,Nc1, "conv_direct");

    Nl_2 = Nl2 / 2;
    Nc_2 = Nc2 / 2;
    for (i = 0; i < Nl1; i++)
    {
        ii = i + Nl_2;
        for (j = 0; j < Nc1; j++)
        {
            (*Result)(i,j) = 0;
            jj = j + Nc_2;
            for (k = 0; k < Nl2; k++)
            for (l = 0; l < Nc2; l++)
                (*Result)(i,j) += Im1(ii - k, jj - l, Type_Border) * Im2(k, l);
       }
    }
    return (*Result);
}

/****************************************************************************/

float correlation (const Iint & Tab_Im1, const Iint & Tab_Im2)
{
    float Sum_X2,Sum_Y2,Sum_XY;
    int i,j;
    int Nl = Tab_Im1.nl();
    int Nc = Tab_Im1.nc();
    float Coef;

    Sum_X2 = Sum_Y2 = Sum_XY = 0.;
    for (i = 0; i < Nl; i++)
    for (j = 0; j < Nc; j++)
    {
        Sum_X2 += (float) (Tab_Im1 (i,j) * Tab_Im1 (i,j));
        Sum_Y2 += (float) (Tab_Im2 (i,j) * Tab_Im2 (i,j));
        Sum_XY += (float) (Tab_Im1 (i,j) * Tab_Im2 (i,j));
    }
    Coef = Sum_XY / sqrt (Sum_X2 * Sum_Y2);
    return (Coef);
}

/***************************************************************************/

float detect_noise_sigma (const Iint &Image, Bool Average_Non_Null, int Nit)
{
    int It, i, j;
    float S0,S1,S2,Sm=0,x;
    int Nl = Image.nl();
    int Nc = Image.nc();
    float Average=0., Sigma=0.;

    for (It = 0; It < Nit; It++)
    {
       S0 = S1 = S2 = 0.;
       for (i = 0; i < Nl; i++)
       for (j = 0; j < Nc; j++)
       {
           x = (float) Image(i,j);
	   if ((It == 0) || (ABS(x - Average) < Sm))
	   { 
	       S0 ++;
	       S1 += x;
	       S2 += x*x;
	   }
       }
       if (Average_Non_Null)
       {
       	   Average = S1 / S0;
	   Sigma = sqrt(S2/S0- Average * Average);
       }
       else  Sigma = sqrt(S2/S0);
       Sm = 3. * Sigma;
    }
    return (Sigma);
}

/***************************************************************************/
