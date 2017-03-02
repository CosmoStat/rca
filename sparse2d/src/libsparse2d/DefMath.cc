/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  22/12/98
**    
**    File:  DefMath.cc
**
*******************************************************************************
**
** double b3_spline (double x)
** 
** Computes the value on a b3-spline 
** we have b3_spline (0) = 2 / 3 
** and if |x| >= 2 => b3_spline (x) = 0
** 
*******************************************************************************
**
**  double xerf (double X)
**
**  compute the inverse function of the repartition function of a Gaussian
**  law. Ref ABRAMOVITZ et STEGUN p. 933
**
*******************************************************************************
**
**  init_random (unsigned int Val)
** 
**  random values generator initialisation from secondes with the
**  unix function drand48
**
*******************************************************************************
**
** float get_random (float Min, float Max)
** 
** return a random value between min and max
**
*****************************************************************************/
 
#include "DefMath.h"
#include "OptMedian.h"

// #ifndef VMS
// extern "C" {
// extern void srand48(long);
// extern double drand48 ();
// }
// #endif

void sort(unsigned long n, float *arr);
void sort(unsigned long n, double *arr);

/***************************************************************************/

double b3_spline (double x)
{
    double A1,A2,A3,A4,A5,Val;

    A1 = ABS ((x - 2) * (x - 2) * (x - 2));
    A2 = ABS ((x - 1) * (x - 1) * (x - 1));
    A3 = ABS (x * x * x);
    A4 = ABS ((x + 1) * (x + 1) * (x + 1));
    A5 = ABS ((x + 2) * (x + 2) * (x + 2));
    Val = 1./12. * (A1 - 4. * A2 + 6. * A3 - 4. * A4 + A5);
    return (Val);
}
    
/***************************************************************************/

double xerf (double X)
/* 
   compute the inverse function of the repartition function of a Gaussian
   law. Ref ABRAMOVITZ et STEGUN p. 933
*/
{
    double Z,A0,A1,B1,B2,Val_Return;
    
    Z = X;
    if (X > 0.5) Z = 1 - X;
    if (Z < DOUBLE_EPSILON)
    {
       Val_Return  = 0.;
    }
    else
    {
       Z = sqrt (-2. * log (Z));
       A0 = 2.30753;
       A1 = 0.27061;
       B1 = 0.99229;
       B2 = 0.04481;
       Val_Return = Z - (A0 + A1 * Z) / (1 + ((B1 + B2 * Z) * Z));
       if (X > 0.5) Val_Return = - Val_Return;
    }
    return (Val_Return);
}

/***************************************************************************/

double inverfc(double y)
{
    double s, t, u, w, x, z;

    z = y;
    if (y > 1) {
        z = 2 - y;
    }
    w = 0.916461398268964 - log(z);
    u = sqrt(w);
    s = (log(u) + 0.488826640273108) / w;
    t = 1 / (u + 0.231729200323405);
    x = u * (1 - s * (s * 0.124610454613712 + 0.5)) - 
        ((((-0.0728846765585675 * t + 0.269999308670029) * t + 
        0.150689047360223) * t + 0.116065025341614) * t + 
        0.499999303439796) * t;
    t = 3.97886080735226 / (x + 3.97886080735226);
    u = t - 0.5;
    s = (((((((((0.00112648096188977922 * u + 
        1.05739299623423047e-4) * u - 0.00351287146129100025) * u - 
        7.71708358954120939e-4) * u + 0.00685649426074558612) * u + 
        0.00339721910367775861) * u - 0.011274916933250487) * u - 
        0.0118598117047771104) * u + 0.0142961988697898018) * u + 
        0.0346494207789099922) * u + 0.00220995927012179067;
    s = ((((((((((((s * u - 0.0743424357241784861) * u - 
        0.105872177941595488) * u + 0.0147297938331485121) * u + 
        0.316847638520135944) * u + 0.713657635868730364) * u + 
        1.05375024970847138) * u + 1.21448730779995237) * u + 
        1.16374581931560831) * u + 0.956464974744799006) * u + 
        0.686265948274097816) * u + 0.434397492331430115) * u + 
        0.244044510593190935) * t - 
        z * exp(x * x - 0.120782237635245222);
    x += s * (x * s + 1);
    if (y > 1) {
        x = -x;
    }
    return x;
}

/***************************************************************************/

// double xerf(double y) {return 1. - xerfc(y);}
double xerfc (double F)
{
   double Nu;
   double P = F;
   if (P > 1.) P = 1.;
   else if (P < 0) P = 0.;
   if (P > 0.5) Nu = sqrt(2.)*inverfc((double) (1.-P)/0.5);
   else Nu = - sqrt(2.)*xerfc((double) P/0.5);
   return Nu;
}      

/***************************************************************************/

void init_random (unsigned int Init)
/* 
  random values generator initialisation from secondes with the
  unix function drand48
*/
{    
// #ifdef RDNT
//     FILE *Fp = popen("date | tr \":\" \"\\12\" | tail -1 | tr \" \" \"\\12\" | head -1","r");
//     fscanf (Fp,"%d", &Init);
// #endif
     if (Init == 100) srand48 (time (0));   // Init=100 is the default value
     else srand48 (Init);
     if (Init == 100) srandom (time (0));   // Init=100 is the default value
     else srandom (Init);
     // srand(Init);    
}

/***************************************************************************/

float get_random()
{
   return rand()/(RAND_MAX+1.0);
}

/***************************************************************************/

float get_random (float Min, float Max)
{
    float Val;
          
    // Val = (float)( get_random() * (Max-Min) + Min);
    Val = drand48 () * (Max-Min) + Min;
    return(Val);
}

/***************************************************************************/

double entropy (float *Pict, int Size, float StepHisto)
{
    int i,Nbr_Val,ind;
    int *Tab_Histo=NULL;
    double Prob_Ev;
    float Min,Max;
    double Entr;

    Min = Max = Pict [0];
    for (i = 1; i < Size; i++)
    {
        if      (Pict [i] > Max) Max = Pict [i];
        else if (Pict [i] < Min) Min = Pict [i];
    }

    /* Calcul de l'entropie */
    Nbr_Val = (int) ( (Max - Min + 1) / StepHisto);
    Tab_Histo = new int [Nbr_Val];
    
    for (i = 0; i < Nbr_Val; i++) Tab_Histo [i] = 0;
    for (i = 0; i < Size; i++)
    {
        ind =  (int)(((Pict[i] - Min) / StepHisto));
	if ((ind < 0) || (ind >= Nbr_Val))
	{
	    cout << "Error in entropy  function ... " << endl;
	    cout << "Nbr_Val = " << Nbr_Val << " Ind = " << ind << endl;
	    exit(-1);
	}
        Tab_Histo[ind] ++;
    }
    Entr = 0.;
    for (i = 0; i < Nbr_Val; i++)
    {
        if (Tab_Histo [i] > 0)
        {
            Prob_Ev = (double) Tab_Histo [i] / (double) Size;
            Entr += - Prob_Ev * log (Prob_Ev) / log(2.);
        }
    }
    if (Tab_Histo != NULL) delete [] Tab_Histo;
    return (Entr);
}

/***************************************************************************/

float get_sigma_mad(float *Image, int N)
{
    float Noise, Med;
    int i;
    fltarray Buff(N);

    for (i = 0; i < N; i++) Buff(i) = Image[i];
    Med = get_median(Buff.buffer(), N);

    for (i = 0; i < N; i++) Buff(i) = ABS(Image[i] - Med);
    Noise = get_median(Buff.buffer(), N);
    return (Noise/0.6745);
}

/***************************************************************************/

float get_sigma_clip (float *Data, int N, int Nit, Bool Average_Non_Null,
            Bool UseBadPixel, float BadPVal)
{
    int It, i;
    double S0,S1,S2,Sm=0,x;
    double Average=0., Sigma=0.;

    for (It = 0; It < Nit; It++)
    {
       S0 = S1 = S2 = 0.;
       for (i = 0; i < N; i++)
        {
           x = Data[i];
           if ((UseBadPixel == False) || (ABS(x-BadPVal) > FLOAT_EPSILON))
           {
	      if ((It == 0) || (ABS(x - Average) < Sm))
	      { 
	         S0 += 1.;
	         S1 += x;
	         S2 += x*x;
	      }
           }
       }
       if (S0 == 0) S0=1;
       if (Average_Non_Null==True)
       {
       	   Average = S1 / S0;
       	   Sigma = S2/S0- Average * Average;
	   // printf("Sigma = %f\n", (float) Sigma);
       	   if (Sigma > 0) Sigma = sqrt(S2/S0- Average * Average);
	   else Sigma = 0.;
       }
       else  Sigma = sqrt(S2/S0);
       Sm = 3. * Sigma;       
    }
    return ((float) Sigma);
}

/***************************************************************************/

float get_sigma_clip_robust (float *Data, int N, int Nit, Bool Average_Non_Null,
            Bool UseBadPixel, float BadPVal)
{
    int It, i;
    double S0,S1,S2,Sm=0,x;
    double Average=0., Sigma=0.;

    for (It = 0; It < Nit; It++)
    {
       S0 = 0.;
       S1 = 0.;
       S2 = 0.;
       for (i = 0; i < N; i++)
       {
           x = Data[i];
           if ((UseBadPixel == False) || (ABS(x-BadPVal) > FLOAT_EPSILON))
           {
	      if ((It == 0) || (ABS(x - Average) < Sm))
	      { 
	         S0 += 1.;
	         S1 += x;
 	      }
           }
       }
       if (S0 == 0) S0=1;
       S1 /= S0;
       for (i = 0; i < N; i++)
       {
           x = Data[i] - S1;
           if ((UseBadPixel == False) || (ABS(x-BadPVal) > FLOAT_EPSILON))
 	      if ((It == 0) || (ABS(x - Average) < Sm)) S2 += x*x;
       }
       Average = S1;
       S2 /= S0;
       // printf("A = %f, S2 = %f\n", Average, S2);
       if (S2 > 0) S2 = sqrt(S2);
       else S2 = 0.;
       Sigma = S2;
       // printf("Sigma(%f) = %f\n", (float) S0, (float) Sigma);
       Sm = 3. * Sigma;       
    }
    return ((float) Sigma);
}

/***************************************************************************/

double skewness(float *Dat, int N)
{
   double Skew;
   double x1,x2,x3,Sigma;
   int i;
   x1=x2=x3=0.;
   for (i=0; i < N; i++)
   {
      x1 += Dat[i];
      x2 +=  pow((double) Dat[i], (double) 2.);
      x3 +=  pow((double) Dat[i], (double) 3.);
   }
   x1 /= (double) N;
   x2 /= (double) N;
   x3 /= (double) N;
   Sigma = x2 - x1*x1;
   if (Sigma > 0.)
   {
      Sigma = sqrt(Sigma);
      Skew = 1. / pow(Sigma,(double) 3.) * (x3 - 3.*x1*x2+2.*x1*x1*x1);
   }
   else Skew = 0.;
   return Skew;
}

/****************************************************************************/

double curtosis(float *Dat, int N)
{
   double Curt;
   double x1,x2,x3,x4,Sigma;
   int i;
   x1=x2=x3=x4=0.;
   for (i=0; i < N; i++)
   {
      x1 += Dat[i];
      x2 +=  pow((double) Dat[i], (double) 2.);
      x3 +=  pow((double) Dat[i], (double) 3.);
      x4 +=  pow((double) Dat[i], (double) 4.);
   }
   x1 /= (double) N;
   x2 /= (double) N;
   x3 /= (double) N;
   x4 /= (double) N;
   Sigma = x2 - x1*x1;
   if (Sigma > 0.)
   { 
      double x1_2 = x1*x1;
      double x1_4 = x1_2*x1_2;
      Sigma = sqrt(Sigma);
      Curt = 1. / pow(Sigma,(double) 4.) * (x4 -4*x1*x3 + 6.*x2*x1_2 -3.*x1_4 ) - 3.;
   }
   else Curt = 0.;
   return Curt; 
}

/****************************************************************************/

void moment4(float *Dat, int N, double &Mean, double &Sigma, 
             double &Skew, double & Curt, float & Min, float & Max)
{
   double x1,x2,x3,x4;
   int i;
   x1=x2=x3=x4=0.;
   Curt=Skew=0;
   Min = Max = Dat[0];
   
   for (i=0; i < N; i++)
   {
      x1 += Dat[i];
      if (Min > Dat[i]) Min = Dat[i];
      if (Max < Dat[i]) Max = Dat[i];
   }
   x1 /= (double) N;
   
   for (i=0; i < N; i++)
   {
      // double Coef=Dat[i];
      x2 +=  pow((double) (Dat[i]-x1), (double) 2.);
      x3 +=  pow((double) (Dat[i]-x1), (double) 3.);
      x4 +=  pow((double) (Dat[i]-x1), (double) 4.);
   }
   
   x2 /= (double) N;
   x3 /= (double) N;
   x4 /= (double) N;
   Sigma = x2;
   if (Sigma > 0.)
   { 
      // double x1_2 = x1*x1;
      // double x1_4 = x1_2*x1_2;
      Sigma = sqrt(Sigma);
      Skew = x3 / pow(Sigma,(double) 3.) ;
      Curt = x4 / pow(Sigma,(double) 4.) - 3.;
   }
   else Sigma = 0.;
   Mean = x1;
}


/****************************************************************************/

void hc_test(float *Dat, int N, float & HC1, float & HC2)
{
   double Sig = 0. ;
   double  Mean = 0.;
   for (int i=0; i < N; i++) 
   {
   Mean += 1/(double )N*Dat[i];
   Sig += 1/(double)N*Dat[i]*Dat[i];
   }
   Sig -= Mean*Mean;
   if (Sig > 0) Sig = sqrt(Sig);
   else Sig=0.;
   
   hc_test( Dat,  N,  HC1, HC2, (float) Sig,  (float) Mean);
}

/****************************************************************************/

void hc_test(float *Dat, int N, float & HC1, float & HC2,  float Mean)
{
   double  Sig = 0;
   for (int i=0; i < N; i++) 
   {
   Sig += 1/(double)N*Dat[i]*Dat[i];
   }
   Sig -= Mean*Mean;
   if (Sig > 0) Sig = sqrt(Sig);
   else Sig=0.;
   hc_test( Dat,  N,  HC1, HC2, (float) Sig,  (float) Mean);
}

/****************************************************************************/

void hc_test(float *Dat, int N, float & HC1, float & HC2, float Sigma, float Mean)
{
   int i;
   double *p = new double [N];
   double IndGen;
   double SqN = sqrt((double) N);
   double Ninv = 1./(double) N;
   int CptPI = 0;   
	 int CptPF = 0;
   HC1 = HC2 = 0;
    
   for (i=0; i < N; i++) 
   {
      double Coef = (Dat[i] - Mean) / Sigma;
      p[i] = 1. - (1. + erf(Coef/sqrt(2.))) / 2.;
      if (p[i] < Ninv) CptPI++;      
			if (p[i] > 1-Ninv) CptPF++;
    }
   double *Ptr = p - 1;
   sort((unsigned long) N, Ptr);
   for (i=0; i <= N-1; i++)
   {
      double W,Den = sqrt(p[i] - p[i]*p[i]);
      IndGen = (i+1.) / (double) N;
      if (Den > FLOAT_EPSILON) W = SqN * (IndGen - p[i]) / Den;
      else W = 0;
      if (abs(W) > HC1) HC1 = abs(W);
      if ((i >= CptPI) && (abs(W) > HC2) && (i <= N-CptPF-1)) HC2 = abs(W);
   }
   delete [] p;
}

/****************************************************************************/

static float root(float A)
{
   int Dim=9;
   fltarray Pol(Dim+1);
   Pol(0) = 1. - A;
   Pol(1) = 0.;
   Pol(2) = 2.- 2.*A;
   Pol(3) = -1.;
   Pol(4) = 1.1;
   Pol(5) = -1.;
   Pol(6) = 0.35;
   Pol(7) = -1./20.;
   Pol(8) = 1./400.;
   float X;
   X = 1.;
   for (int n=0; n < 20; n++)
   {
      float Num = 0.;
      float Den = 0.;
      for (int i=0; i <= 8; i++)
             Num += (float) (Pol(i)*pow( (double) X, (double) i));
       for (int i=1; i <= 8; i++)
             Den += i*Pol(i)*pow( (double) X, (double) (i-1));
       X = X  - Num / Den;
   }
   return X;
}

/*********************************************************************/

void gausstest(float *Band, int N, float &T1, float &T2)
{
   int i;
   float *X = new float[N];
   float M = 0.;
   for (i=0; i < N; i++) M += Band[i];
   M /= (float) N;
   for (i=0; i < N; i++) X[i] = Band[i] - M;
   float C1=0.;
   for (i=0; i < N; i++) C1 += Band[i] * Band[i];
   C1 /= (float) N;
   
   float C2 = get_sigma_mad(X, N);
   // cout << " Std = " << X.sigma() << " MAD = " << C2 << endl;
   C2 *= C2;
   float A = C2/C1;
   float Sol = root(A);
   float Tau = 0.;
   float Sigma;
   if (Sol < 0.18)
   {
      Tau = 0.;
      Sigma = sqrt(C1);
   }
   else
   {
      Tau = Sol;
      Sigma = sqrt(C1/(1.+2*Tau*Tau));
   } 
   // cout << " A = " << A << " Tau = " << Tau << " Sigma = " << Sigma << endl;
   float Coeff;
   T1=T2=0.;
   for (i=0; i < N; i++) 
   {
      X[i] /= Sigma;
      Coeff = ABS(X[i]);
      if (Coeff > 2.)
      {
         T1 += Coeff - 2.;
         T2 += (Coeff - 2.)*Coeff;
      }
   }
   float phi2 =  (1. + erf(2./sqrt(2.))) / 2.;
   float Mu = 2* (exp(-2.) / sqrt( 2.*PI) - 2. * (1. - phi2));
   float V =  2* (5.*(1. - phi2) - 2.*exp(-2.) / sqrt( 2.*PI));
   float Sd = sqrt(V - Mu*Mu);
   T1 = (T1 - N*Mu) / (sqrt( (float) N) * Sd);
   
   float t=2.;
   float u = exp(-t*t/2.) / sqrt( 2.*PI);
   V = 1. - phi2;
   float m2 = V;
   float n2 = (3.+t*t)*V-t*u;
   m2 = 2*m2;
   n2 = 2*n2;
   Sd = sqrt(n2 - m2*m2);
   T2 = (T2 - N*m2) / (sqrt( (float) N) * Sd);
//    cout << " phi2 = " << phi2  << " Mu = " << Mu << " V = " << V << endl;

//   cout << " T1 = " << T1 << " T2 = " << T2 << endl;
//    cout << " Sol 0.5 = " << root(0.5) << endl;
//    cout << " Sol 1 = " << root(1.) << endl;
//    cout << " Sol 1.5 = " << root(1.5) << endl;
//    exit(0);
  delete [] X;
}

/*********************************************************************/

float fdr_pvalue(float *TabPVal, int N, float Alpha)
{
   int i;
   fltarray PVal(N);
   float p_cutoff=0;
   
   for (i = 0; i < N; i++) PVal(i) = TabPVal[i];
   float *P_val_l_tri = PVal.buffer()-1;
   sort(N,P_val_l_tri);
   
   for(i = 0; i < N; i++)
   {
      float j_alpha = Alpha * (float) i / (float) N;
      if (PVal(i) <= j_alpha) p_cutoff = PVal(i);
   }
   return p_cutoff;
}


/*********************************************************************/

double fdr_pvalue(double *TabPVal, int N, float Alpha)
{
   int i;
   dblarray PVal(N);
   double p_cutoff=0;
   
   for (i = 0; i < N; i++) PVal(i) = TabPVal[i];
   double *P_val_l_tri=PVal.buffer()-1;
   sort(N,P_val_l_tri);
   
   for(i = 0; i < N; i++)
   {
      double j_alpha = Alpha * (double ) i / (double) N;
      if (PVal(i) <= j_alpha) p_cutoff = PVal(i);
   }
   return p_cutoff;
}

/*********************************************************************/

float  fdr_gauss_threshold(float *TabVal, int N, float Alpha, float SigmaNoise)
{
   void indexx(int n, double arr[], int indx[]);
   int i, Ind_Cut=-1;
   double S2 = sqrt(2.)*SigmaNoise;
   dblarray PVal(N);
   intarray TInd(N);
   float ValRet, p_cutoff=0;
   
   for (i = 0; i < N; i++) PVal(i) = 1. - erf(abs(TabVal[i])/S2);
   double *P_val_l_tri = PVal.buffer()-1;
   int  *P_ind_tri = TInd.buffer()-1;
   
   indexx(N, P_val_l_tri, P_ind_tri);
    
   for(i = 0; i < N; i++)
   {
      double j_alpha = Alpha * (float) i / (float) N;
      int Ind = TInd(i) - 1;
      if (PVal(Ind) <= j_alpha) 
      {
         p_cutoff = PVal(Ind);
	 Ind_Cut = Ind;
      }
   }
   if (Ind_Cut < 0) ValRet = 10. * TabVal[TInd(N-1)];
   else ValRet = TabVal[Ind_Cut];
   return ValRet;
}

/*********************************************************************/


// void hc_ima_test(float *Dat, float *HCout, int N, float & HC1, float & HC2)
// {
//    void indexx(int n, double arr[], int indx[]);
//    int i;
//    double *p = new double [N];
//    double IndGen;
//    double SqN = sqrt((double) N);
//    double Ninv = 1./(double) N;
//    int CptP = 0;
//    HC1 = HC2 = 0;
//    intarray TInd(N);
//    //cout << "IN " << endl;
//    for (i=0; i < N; i++) 
//    {
//       double Coef =  ABS(xerfc(0.5+(1-Dat[i])/2.));
//       p[i] = 1. - (1. + erf(Coef/sqrt(2.))) / 2.;
//       if (p[i] < Ninv) CptP++;
//    }
//    double *Ptr = p - 1;
//    sort((unsigned long) N, Ptr);
//    int  *P_ind_tri = TInd.buffer()-1;
//    indexx((unsigned long) N, Ptr, P_ind_tri);
//    
//    for (i=0; i < N; i++)
//    {
//       double W,Den = sqrt(p[i] - p[i]*p[i]);
//       IndGen = (i+1.) / (double) N;
//       if (Den > FLOAT_EPSILON) W = SqN * (IndGen - p[i]) / Den;
//       else W = 0;
//       HCout[ TInd(i)-1]=W;
//       if (i <= N/2)
//       {
//          if (W > HC1) HC1 = W;
//          if ((i >= CptP) && (W > HC2)) HC2 = W;
//       }
//    }
//    delete [] p;
//    // cout << "OUT " << endl;
// }


/****************************************************************************/
