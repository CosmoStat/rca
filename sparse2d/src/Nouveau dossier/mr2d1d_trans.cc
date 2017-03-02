/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  04/09/99
**    
**    File:  mr3d_trans.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform of cube
**    ----------- 
**                 
**    Usage: mr3d_trans options cube output
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "SB_Filter.h"
#include "MR1D_Obj.h"
#include "MR_Obj.h"
#include "IM3D_IO.h"
#include "MR3D_Obj.h"
#include "DefFunc.h"

/****************************************************************************/
#define TESTIND 1

class MR2D1D {
  int Nx, Ny, Nz;
  MultiResol *Tab;
  int NbrScale2D;
  int NbrScale1D;
  int NbrBand2D;
  int NbrBand1D;
  MultiResol WT2D;
  MR_1D WT1D;
  fltarray *TabBand;
  intarray TabFirstPosBandNz;
  intarray TabSizeBandNx;
  intarray TabSizeBandNy;
  intarray TabSizeBandNz;
  
  // WT 2D
  sb_type_norm Norm;
  type_sb_filter SB_Filter;
  type_border Bord;
  type_undec_filter U_Filter;
  FilterAnaSynt FAS;
  
  int mr_io_fill_header(fitsfile *fptr);
  public:
       Bool Verbose;
       MR2D1D (){ NbrBand2D=NbrBand1D=0;Verbose=False;}
       void alloc(int iNx, int iNy, int iNz, type_transform Trans2D, int Ns2D, int Ns1D);
       void transform (fltarray &Data);
       void recons (fltarray &Data);
       int nbr_band_2d () const { return NbrBand2D;}
       int nbr_band_1d () const { return NbrBand1D;}
       int size_band_nx(int s2d, int s1d) const { return TabSizeBandNx(s2d, s1d);}
       int size_band_ny(int s2d, int s1d) const { return TabSizeBandNy(s2d, s1d);}
       int size_band_nz(int s2d, int s1d) const { return TabSizeBandNz(s2d, s1d);}
       float & operator() (int s2, int s1, int i, int j, int k) const;
       fltarray get_band(int s2, int s1);
       void put_band(fltarray Band, int s2, int s1);
       void read(char *Name);
       void write(char *Name);
       ~MR2D1D() { delete [] TabBand;}
};

/****************************************************************************/
/****************************************************************************/

static void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 3) || (File_Name_In[L-1] != 'r')
                || (File_Name_In[L-2] != 'm')
                || (File_Name_In[L-3] != '.'))
    {
        strcat (File_Name_Out, ".mr");
    }
}

/****************************************************************************/

/*--------------------------------------------------------------------------*/
static void PrintError( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        /* get the error status description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  /* get first message; null if stack is empty */
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}

/*--------------------------------------------------------------------------*/


int MR2D1D::mr_io_fill_header( fitsfile *fptr)
{
  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/

// cout << " DDDD " << (long) NbrBand2D << " " << (long)NbrBand1D << endl;

if ( ffpkyj(fptr, "Nx", (long)Nx,"Nx of the input cube",&status))
     PrintError( status );  
if ( ffpkyj(fptr,"Ny",(long)Ny,"Ny of the input cube",&status))
     PrintError( status );  
if ( ffpkyj(fptr,"Nz",(long)Nz,"Nz of the input cube",&status))
     PrintError( status );  
 if ( ffpkyj(fptr, "NSCALE2D", (long) NbrScale2D, "Number of bands 2D", &status))
     PrintError( status );  
 if ( ffpkyj(fptr, "NSCALE1D", (long)NbrScale1D, "Number of bands 1D", &status))
     PrintError( status );  
 if ( ffpkyj(fptr, "Type_Tra", (long) WT2D.Type_Transform, 
                         (char*)StringTransform(WT2D.Type_Transform), &status))
     PrintError( status );  
 return(status);
} 


/****************************************************************************/

void  MR2D1D::read(char *Name)
{
    char filename[256];
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0, hdutype ;
    long hdunum;
    char comment[FLEN_COMMENT]; 
    int naxis;
    long naxes[3];
    long mon_long;
    int anynul = 0;
    long nulval = 0;
    long inc[3];
    void PrintError( int status);
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
    //long fpixels[3];
    //long int lpixels[3];
    
     // for multiresol
    float *Ptr;
    //int my_logical; // sais pas...

     mr_io_name (Name, filename);

    inc[0]=1;  inc[1]=1; inc[2]=1;
 
#if DEBUG_IO
    cout << "Read in " << filename << endl;
#endif
   
    /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
         PrintError( status );
                                    
    hdunum = 1;  /*read  table */
    if ( ffmahd(fptr, hdunum, &hdutype, &status) ) /* move to the HDU */
           PrintError( status );

    int simple, bitpix, extend;
    long pcount, gcount;
    if ( ffghpr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           PrintError( status );

    nelements = naxes[0];
     
    if (ffgkyj(fptr,"Nx", &mon_long, comment, &status)) PrintError( status );
    int Nxi = (int) mon_long;  
    if (ffgkyj(fptr,"Ny", &mon_long, comment, &status)) PrintError( status );
    int Nyi = (int) mon_long; 
    if (ffgkyj(fptr,"Nz", &mon_long, comment, &status)) PrintError( status );
    int Nzi = (int) mon_long; 
    
    if (ffgkyj(fptr,"NSCALE2D", &mon_long, comment, &status)) PrintError( status );
    int iNbrScale2D  = (int) mon_long;
    if (ffgkyj(fptr,"NSCALE1D", &mon_long, comment, &status)) PrintError( status );
    int iNbrScale1D  = (int) mon_long;

    type_transform  TT;
    if (ffgkyj(fptr,"Type_Tra", &mon_long, comment, &status))   PrintError( status );
    else TT = (type_transform) mon_long;
    
    alloc(Nxi,Nyi,Nzi,TT, iNbrScale2D, iNbrScale1D);
     
    fltarray Tab(nelements);
   //cout << "        Nl = " << (TabCF_Band[0][1]).nl() << " Nc = " << TabCF_Band[0][1].nc() << endl;

    Ptr = Tab.buffer();
    if ( ffgpve(fptr, 1, 1, nelements, nulval, Ptr, &anynul, &status))
             PrintError( status );

    int ind=2;
    for (int s=0; s <NbrBand2D; s++) 
    for (int s1=0; s1 < NbrBand1D; s1++)
    {
    
       int Nxb = (int) Tab(ind++);
       int Nyb = (int) Tab(ind++);
       int Nzb = (int) Tab(ind++);
       //  cout << s << " " << s1 << " " << Nxb << " " << Nyb << " " << Nzb << endl;

        //cout << " s = " << s+1 << " b = " << b+1 << " Nl = " << Ny << " Nc = " << Nx << endl;
       //cout << "        Nl = " << (TabCF_Band[s][b]).nl() << " Nc = " << TabCF_Band[s][b].nc() << endl;
       for (int k=0;  k < Nzb; k++) 
       for (int j=0;  j < Nyb; j++)
       for (int i=0;  i < Nxb; i++)
          (*this)(s,s1,i,j,k) = Tab(ind++);
   }
      
   if ( ffclos(fptr, &status) ) PrintError( status );
// cout << " end of read fits file " << endl;
}

/****************************************************************************/

void  MR2D1D::write(char *Name)
{

 char filename[256];
 Ifloat Ima;
 fitsfile *fptr;    
 int status;
 int simple;
 int bitpix;
 long naxis=0;
 long naxes[4];
 long nelements;
 long group = 1;  /* group to write in the fits file, 1= first group */
 long firstpixel = 1;    /* first pixel to write (begin with 1) */
 Ifloat Aux;
 
/* we keep mr as extension even if its fits ! */
 mr_io_name (Name, filename);

#if DEBUG_IO  
    cout << "Write on " << filename << endl;
#endif

 FILE *FEXIST = fopen(filename, "rb");
 if (FEXIST)
 {
    fclose(FEXIST);
    remove(filename);               /* Delete old file if it already exists */
 }

 status = 0;         /* initialize status before calling fitsio routines */

    /* open the file */
 if ( ffinit(&fptr, filename, &status) )     /* create the new FITS file */
     PrintError( status );           /* call PrintError if error occurs */
                                                                              
/* write  the header */
 simple   = True;
 bitpix   =  -32;   /* 32-bit real pixel values      */
 long pcount   =   0;  /* no group parameters */
 long gcount   =   1;  /* only a single image/group */
 int  extend   =   False;
  naxis = 1;

  
  int Nelem=2;
  for (int s=0; s < NbrBand2D; s++)
  for (int s1=0; s1 < NbrBand1D; s1++) 
  {
     Nelem += 3 +  TabSizeBandNx(s,s1)*TabSizeBandNy(s,s1)*TabSizeBandNz(s,s1);
  } 
  // cout << " NELEM = " << Nelem << endl;
  fltarray Data(Nelem);
  int ind=2;
  Data(0) = NbrBand2D;
  Data(1) = NbrBand1D;
 
  // Ifloat Band;
  for (int s=0; s < NbrBand2D; s++)
  for (int s1=0; s1 < NbrBand1D; s1++) 
  {
     Data(ind++) = TabSizeBandNx(s,s1);
     Data(ind++) = TabSizeBandNy(s,s1);
     Data(ind++) = TabSizeBandNz(s,s1);
    
     for (int k=0;  k < TabSizeBandNz(s,s1); k++) 
     for (int j=0;  j < TabSizeBandNy(s,s1); j++) 
     for (int i=0;  i < TabSizeBandNx(s,s1); i++) Data(ind++) = (*this)(s,s1,i,j,k);
  }
  // cout << " DATAOK = " <<   endl;
  naxes[0] = ind;
  
 // write first header part (parameters)
 if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
     PrintError( status );          /* call PrintError if error occurs */
  
   // write the header of the multiresolution file
  mr_io_fill_header(fptr);
 
 nelements = ind;
 if ( ffppre(fptr, group, firstpixel, nelements, Data.buffer(), &status) )
              PrintError( status );  
     
  /* close the FITS file */
  if ( ffclos(fptr, &status) )  PrintError( status ); 
}

/****************************************************************************/

float & MR2D1D::operator() (int s2, int s1, int i, int j, int k) const
{    
#ifdef TESTIND
   if ( (i < 0) || (i >= size_band_nx(s2,s1)) ||
        (j < 0) || (j >= size_band_ny(s2,s1)) ||
	(k < 0) || (k >= size_band_nz(s2,s1)) ||
	(s2 < 0) || (s2 >= nbr_band_2d()) ||
	(s1 < 0) || (s1 >= nbr_band_1d()))
   {
      printf("Error: (s2,s1,i,j,k)=(%d,%d,%d,%d,%d), size band = (%d,%d,%d) , scale=(%d,%d) \n", s2,s1,i,j,k,size_band_nx(s2,s1), size_band_ny(s2,s1), size_band_nz(s2,s1), nbr_band_2d(), nbr_band_1d());
      exit(-1);
   }
#endif
  
   return TabBand[s2] (i,j,k+TabFirstPosBandNz(s1));
}

/****************************************************************************/

fltarray MR2D1D::get_band(int s2, int s1)
{
    int Nxb = size_band_nx(s2,s1);
    int Nyb = size_band_ny(s2,s1);
    int Nzb = size_band_nz(s2,s1);
    fltarray *Cube_Return = NULL;

    Cube_Return = new fltarray(Nxb, Nyb, Nzb);
    for (int i=0; i < Nxb; i++)
    for (int j=0; j < Nyb; j++)
    for (int k=0; k < Nzb; k++) (*Cube_Return)(i,j,k) = (*this)(s2,s1,i,j,k);
    
    return (*Cube_Return);
}

/****************************************************************************/

void  MR2D1D::put_band(fltarray Band, int s2, int s1)
{
    int Nxb = size_band_nx(s2,s1);
    int Nyb = size_band_ny(s2,s1);
    int Nzb = size_band_nz(s2,s1);
    
    for (int i=0; i < Nxb; i++)
    for (int j=0; i < Nyb; j++)
    for (int k=0; k < Nzb; k++) (*this)(s2,s1,i,j,k) = Band(i,j,k);
}

/****************************************************************************/

void MR2D1D::alloc (int iNx, int iNy, int iNz, type_transform Trans2D, int Ns2D, int Ns1D)
{
   Nx = iNx;
   Ny = iNy;
   Nz = iNz;
   NbrScale2D = Ns2D;
   NbrScale1D = Ns1D;
   
   Norm = NORM_L2;
   SB_Filter = F_MALLAT_7_9;
   Bord = I_CONT;
   U_Filter = DEF_UNDER_FILTER; 
   FilterAnaSynt *PtrFAS = NULL;
    if ((Trans2D == TO_MALLAT) || (Trans2D == TO_UNDECIMATED_MALLAT))
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    int NbrUndec = -1;                     /*number of undecimated scale */
    type_lift LiftingTrans = DEF_LIFT;
    if (Trans2D == TO_LIFTING) WT2D.LiftingTrans = LiftingTrans;
    WT2D.Border = Bord;
    WT2D.Verbose = Verbose;    
    WT2D.alloc (Ny, Nx, Ns2D, Trans2D, PtrFAS, Norm, NbrUndec, U_Filter);
    NbrBand2D = WT2D.nbr_band();
    WT2D.ModifiedATWT = True;

     	       
    Bool Rebin=False;
    WT1D.U_Filter = U_Filter;
    type_trans_1d Trans1D = TO1_MALLAT;
    WT1D.alloc (Nz, Trans1D, Ns1D, PtrFAS, Norm, Rebin);   
    NbrBand1D = WT1D.nbr_band();
    
   TabBand = new fltarray [NbrBand2D];
   TabSizeBandNx.resize(NbrBand2D, NbrBand1D);
   TabSizeBandNy.resize(NbrBand2D, NbrBand1D);
   TabSizeBandNz.resize(NbrBand2D, NbrBand1D);
   TabFirstPosBandNz.resize(NbrBand1D);
   TabFirstPosBandNz(0) =0;
   for (int b=0; b < NbrBand2D; b++) 
   {
      TabBand[b].alloc(WT2D.size_band_nc(b), WT2D.size_band_nl(b),  WT1D.size_ima_np ());
      for (int b1=0; b1 < NbrBand1D; b1++)
      {
         TabSizeBandNx(b,b1) = WT2D.size_band_nc(b);
	 TabSizeBandNy(b,b1) = WT2D.size_band_nl(b);
	 TabSizeBandNz(b,b1) = WT1D.size_scale_np(b1);
      }
   }
   for (int b1=1; b1 < NbrBand1D; b1++) TabFirstPosBandNz(b1) = TabFirstPosBandNz(b1-1) + TabSizeBandNz(0,b1-1);
}

/****************************************************************************/

void MR2D1D::transform (fltarray &Data)
{
   int i,j,b,z;
   Nx = Data.nx();
   Ny = Data.ny();
   Nz = Data.nz();
   Ifloat Frame(Ny, Nx);
   fltarray Vect(Nz);
   
   // 2D wt transform per frame
   for (z=0; z < Nz; z++)
   {
      for (i=0; i < Ny; i++)
      for (j=0; j < Nx; j++) Frame(i,j) = Data(j,i,z);
      WT2D.transform(Frame);
      for (b=0; b < NbrBand2D; b++)
      {
         for (i=0; i < WT2D.size_band_nl(b); i++)
	 for (j=0; j < WT2D.size_band_nc(b); j++) TabBand[b](j,i,z) = WT2D(b,i,j);
      }
   }
 
   // 1D wt 
   if (NbrBand1D >= 2)
   {
     for (b=0; b < NbrBand2D; b++)
     for (i=0; i < WT2D.size_band_nl(b); i++)
     for (j=0; j < WT2D.size_band_nc(b); j++) 
     {
        for (z=0; z < Nz; z++) Vect(z) = TabBand[b](j,i,z);
        WT1D.transform(Vect);
        z = 0;
        for (int b1=0; b1 < NbrBand1D; b1++)
        {
         for (int p=0; p < WT1D.size_scale_np (b1); p++) TabBand[b](j,i,z++) = WT1D(b1,p); 
        }
      }
   }
}

/****************************************************************************/


/****************************************************************************/

void MR2D1D::recons (fltarray &Data)
{
   if ((Data.nx() != Nx) || (Data.ny() != Ny) || (Data.nz() != Nz)) Data.resize(Nx, Ny, Nz); 
   int i,j,b,z;
   Nx = Data.nx();
   Ny = Data.ny();
   Nz = Data.nz();
   Ifloat Frame(Ny, Nx);
   fltarray Vect(Nz);
     
   // 1D wt 
   if (NbrBand1D >= 2)
   {
      for (b=0; b < NbrBand2D; b++)
      for (i=0; i < WT2D.size_band_nl(b); i++)
      for (j=0; j < WT2D.size_band_nc(b); j++) 
      {
         // for (z=0; z < Nz; z++) Vect(z) = TabBand[b](j,i,z);
	z = 0;
        for (int b1=0; b1 < NbrBand1D; b1++)
        for (int p=0; p < WT1D.size_scale_np (b1); p++) WT1D(b1,p) = TabBand[b](j,i,z++); 
         Vect.init();
	 WT1D.recons(Vect);
         for (z=0; z < Nz; z++) TabBand[b](j,i,z) = Vect(z);
     }
   }
   
   // 2D wt 
   for (z=0; z < Nz; z++)
   {
      for (b=0; b < NbrBand2D; b++)
      {
         for (i=0; i < WT2D.size_band_nl(b); i++)
	 for (j=0; j < WT2D.size_band_nc(b); j++) WT2D(b,i,j) = TabBand[b](j,i,z);
      }   
      WT2D.recons(Frame);
      for (i=0; i < Ny; i++)
      for (j=0; j < Nx; j++) Data(j,i,z) = Frame(i,j);
   }
}
 

/****************************************************************************/

char Name_Cube_In[256];
char Name_Out[256];
 
int NbrScale2d = 5;
int Nbr_Plan=4;
type_transform  Transform=TO_MALLAT;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

// sb_type_norm Norm = NORM_L2;
// type_sb_filter SB_Filter = F_HAAR; // F_MALLAT_7_9;
// type_lift LiftingTrans = DEF_LIFT;
Bool Verbose=False;
Bool GetMax = False;
float N_Sigma=3;                // number of sigma (for the noise)  
Bool Normalize=False;       // normalize data in
 
type_border Bord = I_MIRROR;
Bool Reverse = False;

/*********************************************************************/

/*static int max_scale_number (int Nc)
{
    int Nmin = Nc;
    int ScaleMax;

    ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
    return (ScaleMax);
}*/

/*********************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options cube output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    all_transform_usage(Transform);
    manline();
    fprintf(OUTMAN, "         [-n number_of_scales_2D]\n");
    fprintf(OUTMAN, "              number of scales used in the multiresolution transform\n");
    manline();
    fprintf(OUTMAN, "         [-N number_of_scales_2D]\n");
    fprintf(OUTMAN, "              number of scales used in the multiresolution transform\n");
    manline();
    fprintf(OUTMAN, "         [-M]\n");
    fprintf(OUTMAN, "             Normalize the data. Default is no. \n");    
    manline();
    manline();
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Reconstruction. Default is no. \n");    
    manline();       
    verbose_usage();    
    manline();
    vm_usage();
    manline();    
    manline();
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif   

    /* get options */
    while ((c = GetOpt(argc,argv,"rN:t:n:MvzZ:")) != -1) 
    {
	switch (c) 
        {
           case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM+1)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}                
 		break;
	   case 'r': Reverse = True; break;
           case 'M': Normalize = (Normalize == True) ? False: True; break;
 	   case 'v': Verbose = True; break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&NbrScale2d) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((NbrScale2d <= 1) || (NbrScale2d > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
  		    exit(-1);
		}
		break;
	case 'N':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 0) || (Nbr_Plan > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
  		    exit(-1);
		}
		break;
#ifdef LARGE_BUFF
	    case 'z':
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: Z option already set...\n");
                   exit(-1);
                }
	        OptZ = True;
	        break;
            case 'Z':
	        if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1)
		{
		   fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
		   exit(-1);
		}
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: z option already set...\n");
                   exit(-1);
                }
		OptZ = True;
                break;
#endif
 	   case '?':
			usage(argv);
		}
	}
      
	
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Out, argv[OptInd++]);
        else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/
 
int main(int argc, char *argv[])
{
    fltarray Dat;
    /* Get command line arguments, open input file(s) if necessary */
   fitsstruct Header;
   char Cmd[512];
   Cmd[0] = '\0';
   for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
   
    lm_check(LIC_MR3);
    transinit(argc, argv);
 

    
    if (Reverse == False)
    {
        if (Verbose == True)
        {
           cout << "Filename in = " << Name_Cube_In << endl;
           cout << "Filename out = " << Name_Out  << endl;
           cout << "Transform = " << StringTransform((type_transform) Transform) << endl;
           cout << "NbrScale2d = " << NbrScale2d<< endl;
           cout << "NbrScale1d = " << Nbr_Plan<< endl;
        }
    
       io_3d_read_data(Name_Cube_In, Dat, &Header);
 
       int Nx = Dat.nx();
       int Ny = Dat.ny();
       int Nz = Dat.nz();
       if (Verbose == True) cout << "Nx = " << Dat.nx() << " Ny = " << Dat.ny() << " Nz = " << Dat.nz() << endl;
     
       if (Normalize == True)
       {
         double Mean = Dat.mean();
         double Sigma = Dat.sigma();
         for (int i=0;i<Nx;i++)
         for (int j=0;j<Ny;j++)
         for (int k=0;k<Nz;k++) Dat(i,j,k) = (Dat(i,j,k)-Mean)/Sigma;
       }    
    
   
       MR2D1D WT;
       if (Verbose == True) cout << "Alloc ...  " << endl;
       WT.alloc(Nx, Ny, Nz, Transform, NbrScale2d, Nbr_Plan);

       if (Verbose == True) cout << "Transform ...  " << endl;
       WT.transform (Dat);

       if (Verbose == True)cout << "Write result ...  " << endl;
       WT.write(Name_Out);
              
       if (Verbose == True)
       {
          for (int s2 = 0; s2 < WT.nbr_band_2d (); s2++)
	  for (int s1 = 0; s1 < WT.nbr_band_1d (); s1++)
	  {
	    cout << "  Band " << s2 << ", " << s1 << ": " << " Nx = " << WT.size_band_nx(s2,s1) << ", Ny = " << WT.size_band_ny(s2,s1) <<  ", Nz = " << WT.size_band_nz(s2,s1) << endl;
	    fltarray Band;
	    Band = WT.get_band(s2, s1);
	    cout << "  Sigma = " << Band.sigma() << " Min = " << Band.min() << " Max = " << Band.max() << endl;
        char Name[512];
        sprintf(Name, "Band_%d_%d.fits", s2+1, s1+1);
         fits_write_fltarr(Name, Band);

	 }
       }
//        cout << endl << "READ " << endl;
//        MR2D1D WT1;
//        WT1.read(Name_Out);
//        fltarray Dat1;
//        WT1.recons (Dat1);
//        fits_write_fltarr (Name_Out, Dat1);
    }
    else
    {
    
       MR2D1D WT;
       WT.read(Name_Cube_In);
       WT.recons (Dat);
 
       if (Verbose == True) cout << "Write result ...  " << endl;
       fits_write_fltarr(Name_Out, Dat);
    }
    // Header.origin = Cmd;	 
    // fits_write_fltarr (Name_Out, StaInfo, &Header);
     
    exit(0);
}
