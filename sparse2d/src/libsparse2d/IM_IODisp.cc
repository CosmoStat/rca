/*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck
**
**    Date:  96/06/13 
**    
**    File:  IM_IODisp.cc
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
*******************************************************************************
**
** void io_disp::read_type(char *File_Name)
** 
** read the type of data in image File_Name
**
*******************************************************************************
**
** void io_disp::read_size  (char *File_Name) 
**
** read the size (number of lines and columns) of the image File_Name
**
*******************************************************************************
**
** void io_disp::read_ima  (char *File_Name) 
**
** read the image from file
**
*******************************************************************************
**
** void io_disp::write_ima  (char *File_Step_Name) 
**
** write the image to file
**
*******************************************************************************
**
** void io_disp::read_ima_int (char *File_Name, Iint & Image)
**
** read int image from file
**
*******************************************************************************
**
** void io_disp::read_ima_float (char *File_Name, Ifloat & Image)
**
** read float image from file
**
*******************************************************************************
**
** void io_disp::read_ima_complex_f (char *File_Name, Icomplex_f & Image)
**
** read complex float image from file
**
*******************************************************************************
**
** void io_disp::write_ima_int (char *File_Name, Iint & Image)
**
** write int image to file
**
*******************************************************************************
**
** void io_disp::write_ima_float (char *File_Name, Ifloat & Image)
**
** write float image to file
** 
*******************************************************************************
**
** void io_disp::write_ima_complex_f (char *File_Name, Icomplex_f & Image)
**
** write complex float image to file
**
******************************************************************************/

// static char sccsid[] = "@(#)IM_IODisp.cc 3.2 96/06/13 CEA 1994 @(#)";
 
#include"GlobalInc.h"
#include "IM_IO.h"
#include "IM_BIO.h"

/*****************************************************************/
 
FILE *disp_file_des_in(char *fname)
{
   FILE *fp;
   
   if (std_inout(fname) == True)  fp = stdin;
   else 
   {
      fp = fopen(fname, "rb");
      if (!fp) 
      {
        cerr << "Unable to open file: " <<  fname  << endl;
        exit(-1);
      }
   }
   return fp;
}

FILE *disp_file_des_out(char *fname)
{
   FILE *fp;
   
   if (std_inout(fname) == True) fp = stdout;
   else 
   {
      fp = fopen(fname, "wb");
      if (!fp) 
      {
        cerr << "Unable to open file: " <<  fname  << endl;
        exit(-1);
      }
   }
   return fp;
}

/*****************************/

static void read_size_f (char *File_Name, int *nl, int *nc)
{
    int L, Nbr,Ok = 0;
    char *Ptr;

    L = strlen (File_Name);
    Ptr = File_Name;
    while ((L > 3) && ((Ptr[0] != '_') || (Ptr[1] != 'f') || (Ptr[2] != '_')))
    {
        Ptr ++;
        L --;
    }
    if ((Ptr[0] == '_') && (Ptr[1] == 'f') && (Ptr[2] == '_'))
    {
        Nbr = sscanf(Ptr,"_f_%d_%d.d",nl,nc);
        // printf("File_Name = %s, Nl = %d, Nc = %d\n",File_Name,*nl,*nc);
        if ((Nbr == 2) && (*nl > 0) && (*nc > 0) ) Ok = 1;
    }
    if (Ok == 0)
    {
        fprintf(stderr, "Error f: reading string  %s\n", File_Name);
        exit (-1);
    }
}

/*****************************************************************/

static void read_size_i (char *File_Name, int *nl, int *nc)
{
    int L, Nbr,Ok = 0;
    char *Ptr;

    L = strlen (File_Name);
    Ptr = File_Name;
    while ((L > 3) && ((Ptr[0] != '_') || (Ptr[1] != 'i') || (Ptr[2] != '_')))
    {
        Ptr ++;
        L --;
    }
    if ((Ptr[0] == '_') && (Ptr[1] == 'i') && (Ptr[2] == '_'))
    {
        Nbr = sscanf(Ptr,"_i_%d_%d.d",nl,nc);
        // printf("File_Name = %s, Nl = %d, Nc = %d\n",File_Name,*nl,*nc);
        if ((Nbr == 2) && (*nl > 0) && (*nc > 0) ) Ok = 1;
    }
    if (Ok == 0)
    {
        fprintf(stderr, "Error i: reading string  %s\n", File_Name);
        exit (-1);
    }
}


/*****************************************************************/

static void read_size_b (char *File_Name, int *nl, int *nc)
{
    int L, Nbr,Ok = 0;
    char *Ptr;

    L = strlen (File_Name);
    Ptr = File_Name;
    while ((L > 3) && ((Ptr[0] != '_') || (Ptr[1] != 'b') || (Ptr[2] != '_')))
    {
        Ptr ++;
        L --;
    }
    if ((Ptr[0] == '_') && (Ptr[1] == 'b') && (Ptr[2] == '_'))
    {
        Nbr = sscanf(Ptr,"_b_%d_%d.d",nl,nc);
        // printf("File_Name = %s, Nl = %d, Nc = %d\n",File_Name,*nl,*nc);
        if ((Nbr == 2) && (*nl > 0) && (*nc > 0) ) Ok = 1;
    }
    if (Ok == 0)
    {
        fprintf(stderr, "Error b: reading string  %s\n", File_Name);
        exit (-1);
    }
}

/*****************************************************************/

static void read_size_s (char *File_Name, int *nl, int *nc)
{
    int L, Nbr,Ok = 0;
    char *Ptr;

    L = strlen (File_Name);
    Ptr = File_Name;
    while ((L > 3) && ((Ptr[0] != '_') || (Ptr[1] != 's') || (Ptr[2] != '_')))
    {
        Ptr ++;
        L --;
    }
    if ((Ptr[0] == '_') && (Ptr[1] == 's') && (Ptr[2] == '_'))
    {
        Nbr = sscanf(Ptr,"_s_%d_%d.d",nl,nc);
        // printf("File_Name = %s, Nl = %d, Nc = %d\n",File_Name,*nl,*nc);
        if ((Nbr == 2) && (*nl > 0) && (*nc > 0) ) Ok = 1;
    }
    if (Ok == 0)
    {
        fprintf(stderr, "Error s: reading string  %s\n", File_Name);
        exit (-1);
    }
}

/*****************************************************************/

static void read_size_d (char *File_Name, int *nl, int *nc)
{
    int L, Nbr,Ok = 0;
    char *Ptr;

    L = strlen (File_Name);
    Ptr = File_Name;
    while ((L > 3) && ((Ptr[0] != '_') || (Ptr[1] != 'd') || (Ptr[2] != '_')))
    {
        Ptr ++;
        L --;
    }
    if ((Ptr[0] == '_') && (Ptr[1] == 'd') && (Ptr[2] == '_'))
    {
        Nbr = sscanf(Ptr,"_d_%d_%d.d",nl,nc);
        // printf("File_Name = %s, Nl = %d, Nc = %d\n",File_Name,*nl,*nc);
        if ((Nbr == 2) && (*nl > 0) && (*nc > 0) ) Ok = 1;
    }
    if (Ok == 0)
    {
        fprintf(stderr, "Error d: reading string  %s\n", File_Name);
        exit (-1);
    }
}

/*****************************************************************/

static void read_size_cf (char *File_Name, int *nl, int *nc)
{
    int L, Nbr,Ok = 0;
    char *Ptr;

    L = strlen (File_Name);
    Ptr = File_Name;
    while ((L > 4) && ((Ptr[0] != '_') || (Ptr[1] != 'c') || (Ptr[2] != 'f') || (Ptr[3] != '_')))
    {
        Ptr ++;
        L --;
    }
    if ((Ptr[0] == '_') && (Ptr[1] == 'c') && (Ptr[2] == 'f')&& (Ptr[3] == '_'))
    {
        Nbr = sscanf(Ptr,"_cf_%d_%d.d",nl,nc);
        // printf("File_Name = %s, Nl = %d, Nc = %d\n",File_Name,*nl,*nc);
        if ((Nbr == 2) && (*nl > 0) && (*nc > 0) ) Ok = 1;
    }
    if (Ok == 0)
    {
        fprintf(stderr, "Error cf: reading string  %s\n", File_Name);
        exit (-1);
    }
}

/*****************************************************************/

static void read_size_cd (char *File_Name, int *nl, int *nc)
{
    int L, Nbr,Ok = 0;
    char *Ptr;

    L = strlen (File_Name);
    Ptr = File_Name;
    while ((L > 4) && ((Ptr[0] != '_') || (Ptr[1] != 'c') || (Ptr[2] != 'd') || (Ptr[3] != '_')))
    {
        Ptr ++;
        L --;
    }
    if ((Ptr[0] == '_') && (Ptr[1] == 'c') && (Ptr[2] == 'd')&& (Ptr[3] == '_'))
    {
        Nbr = sscanf(Ptr,"_cd_%d_%d.d",nl,nc);
        // printf("File_Name = %s, Nl = %d, Nc = %d\n",File_Name,*nl,*nc);
        if ((Nbr == 2) && (*nl > 0) && (*nc > 0) ) Ok = 1;
    }
    if (Ok == 0)
    {
        fprintf(stderr, "Error cd: reading string  %s\n", File_Name);
        exit (-1);
    }
}

/*****************************************************************/

void io_disp::read_type(char *File_Name)
{
    const char * S_int = "_i_";
    const char *  S_float = "_f_";
    const char *  S_double = "_d_";
    const char *  S_complex_f = "_cf_";
    const char *  S_complex_d = "_cd_";
    const char *  S_byte = "_b_";
    const char *  S_short = "_s_";

    Type = UNKNOWN;
    if (strstr(File_Name, S_int) != NULL) Type = T_INT;
    else if (strstr(File_Name, S_float) != NULL) Type = T_FLOAT;
    else if (strstr(File_Name, S_double) != NULL) Type = T_DOUBLE;
    else if (strstr(File_Name, S_complex_f) != NULL) Type = T_COMPLEX_F;
    else if (strstr(File_Name, S_complex_d) != NULL) Type = T_COMPLEX_D;
    else if (strstr(File_Name, S_byte) != NULL) Type = T_BYTE;
    else if (strstr(File_Name, S_short) != NULL) Type = T_SHORT;
}

/*****************************************************************/

void io_disp::read_size  (char *File_Name) 
{
    switch (Type)
    {
        case T_INT:read_size_i (File_Name, &Nl_io, &Nc_io);break;
        case T_FLOAT:read_size_f (File_Name, &Nl_io, &Nc_io);break;
        case T_DOUBLE:read_size_d (File_Name, &Nl_io, &Nc_io);break;
        case T_COMPLEX_F:read_size_cf (File_Name, &Nl_io, &Nc_io);break;
        case T_COMPLEX_D:read_size_cd (File_Name, &Nl_io, &Nc_io);break;
        case T_BYTE:read_size_b (File_Name, &Nl_io, &Nc_io);break;
        case T_SHORT:read_size_s (File_Name, &Nl_io, &Nc_io);break;
        default:break;
    }
}

/*****************************************************************/

void io_disp::read_ima  (char *File_Name) 
{
    FILE *File_Des;
    int Size = Nl_io*Nc_io;
    int Nbr=0;


    File_Des = disp_file_des_in(File_Name);
    switch (Type)
    {
        case T_INT: 
               Data_int = MemInt.alloc (Size);
               Nbr = fread ((void *) Data_int, 
                     sizeof (int), Size , File_Des);
               break;
        case T_FLOAT:
               Data_float = MemFloat.alloc (Size);
               Nbr = fread ((void *) Data_float, 
                             sizeof (float), Size , File_Des);
               break;
        case T_COMPLEX_F:
               Data_complex_f = MemCF.alloc (Size);
               Nbr = fread ((void *) Data_complex_f, 
                            sizeof (complex_f), Size , File_Des);
               break;
        case T_DOUBLE:
               // Data_double = new double[Size];
	       Data_double = (double *) alloc_buffer(Size*sizeof(double));
               Nbr = fread ((void *) Data_double, 
                             sizeof (double), Size , File_Des);
               break;
        case T_COMPLEX_D:
               // Data_complex_d = new complex_d[Size];
	       Data_complex_d = (complex_d *) alloc_buffer(Size*sizeof(complex_d));
               Nbr = fread ((void *) Data_complex_d, 
                             sizeof (complex_d), Size , File_Des);
               break;
        case T_BYTE: 
               // Data_byte = new unsigned char[Size];
	       Data_byte= (unsigned char *) alloc_buffer(Size*sizeof(char));
               Nbr = fread ((void *) Data_byte, 
                     sizeof (unsigned char), Size , File_Des);
               break;
        case T_SHORT: 
               // Data_short = new short int [Size];
	       Data_short = (short *) alloc_buffer(Size*sizeof(short));
               Nbr = fread ((void *) Data_short, 
                     sizeof (short int), Size , File_Des);
               break;

        default:break;
    }
    if (Nbr <= 0)
    {
        fprintf(stderr,"Error reading file %s\n",File_Name);
        exit (-1);
    }
    if (File_Des != stdin) fclose (File_Des);
}

/*****************************************************************/

void io_disp::read_info_ima(char *File_Name, int &Nl, int &Nc, 
                            type_data & TypeDat)
{
    FILE *File_Des = disp_file_des_in (File_Name);
 
    read_type(File_Name);
    read_size(File_Name);
    Nl = Nl_io;
    Nc = Nc_io;
    TypeDat = Type;
    if (File_Des != stdin) fclose(File_Des);
}

/*****************************************************************/

void io_disp::read_block_ima(char *File_Name,Ifloat & Image,int Indi, int Indj)
{
    FILE *File_Des;
    unsigned char *Buff;
    unsigned char *buffer;
    int SizeElem;
    
    File_Des = disp_file_des_in(File_Name);
    int i;
    read_type(File_Name);
    read_size(File_Name);
    switch (Type)
    {
        case T_INT: SizeElem = sizeof(int); break;  
        case T_FLOAT: SizeElem = sizeof(float); break;
        case T_COMPLEX_F: SizeElem = sizeof(complex_f); break;
        case T_DOUBLE: SizeElem = sizeof(double); break;
        case T_COMPLEX_D: SizeElem = sizeof(complex_d); break;
        case T_BYTE: SizeElem = sizeof(char); break;
        case T_SHORT: SizeElem = sizeof(short int); break; 
        default:
	  cerr << "Error: unknown type of data ... " << endl;
	  exit(-1);
	  break;
    }
    if (Type != T_FLOAT)
           Buff = (unsigned char *) 
	          alloc_buffer((size_t)(SizeElem*Image.nl()*Image.nc()*sizeof(char)));
           // Buff = new unsigned char [SizeElem*Image.nl()*Image.nc()];
    else Buff = (unsigned char *) Image.buffer();
    buffer = Buff;
    long fp = Indi * Nc_io + Indj;
    fp *= SizeElem;
    fseek(File_Des, fp, SEEK_SET);
    fp = Nc_io - Image.nc();
    fp *= SizeElem;
    for (i=0; i < Image.nl(); i++)
    {
       if (fread (buffer, SizeElem, Image.nc(), File_Des) != (size_t) Image.nc())
       {
          cerr << "Error: Unable to read in file: " << File_Name << endl;
	  cerr << "       line " << i << endl;
	  cerr << "       Nc_io " << Nc_io << " nc = " << Image.nc() << endl;
          exit(-1);
       }
       if (i != Image.nl()-1)
       {
           fseek(File_Des, fp, SEEK_CUR);
	   buffer += Image.nc()*SizeElem;
       }
    } 
    switch (Type)
    {
       case T_INT:
        {
	 int *Dat = (int *) Buff;
         for (i=0; i < Image.n_elem(); i++) Image(i) = (float) Dat[i];
        }
        break;
       case T_BYTE:
        {
	 unsigned char *Dat = (unsigned char *) Buff;
         for (i=0; i < Image.n_elem(); i++) Image(i) = (float) Dat[i];
        }
        break;
       case T_SHORT:
        {
	 unsigned short *Dat = (unsigned short *) Buff;
         for (i=0; i < Image.n_elem(); i++) Image(i) = (float) Dat[i];
        }
	break;
       case T_FLOAT: break;
       case T_DOUBLE:
        {
	 double *Dat = (double *) Buff;
         for (i=0; i < Image.n_elem(); i++) Image(i) = (float) Dat[i];
        }
	break;
       case T_COMPLEX_F:
        {
	 complex_f *Dat = (complex_f  *) Buff;
         for (i=0; i < Image.n_elem(); i++) Image(i) = (float) Dat[i].real();
        }
	break;
       case T_COMPLEX_D:
       {
	 complex_d *Dat = (complex_d  *) Buff;
         for (i=0; i < Image.n_elem(); i++) Image(i) = (float) Dat[i].real();
        }
       default :
         fprintf (stderr, "Error: type unknown reading %s\n", File_Name);
         exit (-1);
    }
    // if (Type != T_FLOAT) delete Buff;
    if (Type != T_FLOAT) free_buffer( (char *) Buff);
    if (File_Des != stdin) fclose(File_Des);
}

/*****************************************************************/

void io_disp::read_block_ima(char *File_Name,Iint & Image,int Indi, int Indj)
{
    FILE *File_Des;
    unsigned char *Buff;
    unsigned char *buffer;
    int SizeElem;
    
    File_Des = disp_file_des_in(File_Name);
    int i;
    read_type(File_Name);
    read_size(File_Name);
    switch (Type)
    {
        case T_INT: SizeElem = sizeof(int); break;  
        case T_FLOAT: SizeElem = sizeof(float); break;
        case T_COMPLEX_F: SizeElem = sizeof(complex_f); break;
        case T_DOUBLE: SizeElem = sizeof(double); break;
        case T_COMPLEX_D: SizeElem = sizeof(complex_d); break;
        case T_BYTE: SizeElem = sizeof(char); break;
        case T_SHORT: SizeElem = sizeof(short int); break; 
        default: 
          cerr << "Error: unknown type of data ... " << endl;
	  exit(-1);
	  break;
    }
    if (Type != T_INT)
           Buff = (unsigned char *) 
	          alloc_buffer((size_t)(SizeElem*Image.nl()*Image.nc()*sizeof(char)));
    // Buff = new unsigned char [SizeElem*Image.nl()*Image.nc()];
    else Buff = (unsigned char *) Image.buffer();
    buffer = Buff;
    long fp = Indi * Nc_io + Indj;
    fp *= SizeElem;
    fseek(File_Des, fp, SEEK_SET);
    fp = Nc_io - Image.nc();
    fp *= SizeElem;
    for (i=0; i < Image.nl(); i++)
    {
       if (fread (buffer, SizeElem, Image.nc(), File_Des) != (size_t)  Image.nc())
       {
          cerr << "Unable to read in file: " << File_Name << endl;
          exit(-1);
       }
       if (i != Image.nl()-1)
       {
           fseek(File_Des, fp, SEEK_CUR);
	   buffer += Image.nc()*SizeElem;
       }
    } 
    int Nelem = Image.nl()*Image.nc();
    switch (Type)
    {
       case T_INT: break;
       case T_BYTE:
        {
	 unsigned char *Dat = (unsigned char *) Buff;
         for (i=0; i < Nelem; i++) Image(i) = (int) Dat[i];
        }
        break;
       case T_SHORT:
        {
	 unsigned short *Dat = (unsigned short *) Buff;
         for (i=0; i < Nelem; i++) Image(i) = (int) Dat[i];
        }
	break;
       case T_FLOAT:
        {
	 float *Dat = (float *) Buff;
         for (i=0; i < Nelem; i++) Image(i) = (int) Dat[i];
        }
	break;
       case T_DOUBLE:
        {
	 double *Dat = (double *) Buff;
         for (i=0; i < Nelem; i++) Image(i) = (int) Dat[i];
        }
	break;
       case T_COMPLEX_F:
        {
	 complex_f *Dat = (complex_f  *) Buff;
         for (i=0; i < Nelem; i++) Image(i) = (int) Dat[i].real();
        }
	break;
       case T_COMPLEX_D:
       {
	 complex_d *Dat = (complex_d  *) Buff;
         for (i=0; i < Nelem; i++) Image(i) = (int) Dat[i].real();
        }
       default :
         fprintf (stderr, "Error: type unknown reading %s\n", File_Name);
         exit (-1);
    }
    // if (Type != T_INT) delete Buff;
    if (Type != T_INT) free_buffer( (char *) Buff);
    if (File_Des != stdin) fclose(File_Des);
}

/*****************************************************************/

void io_disp::write_block_ima(char *File_Name,Ifloat & Image,int Indi, int Indj)
{
    FILE *File_Des;
    unsigned char *Buff=NULL;
    unsigned char *buffer;
    int SizeElem;
    
    File_Des = fopen (File_Name, "rb+");
    if (File_Des  == NULL)
    {
     cerr << "Unable to open file: " << File_Name  << endl;
     exit(-1);
    }
    int i;
    read_type(File_Name);
    read_size(File_Name);
    switch (Type)
    {
        case T_INT: SizeElem = sizeof(int); break;  
        case T_FLOAT: SizeElem = sizeof(float); break;
        case T_COMPLEX_F: SizeElem = sizeof(complex_f); break;
        case T_DOUBLE: SizeElem = sizeof(double); break;
        case T_COMPLEX_D: SizeElem = sizeof(complex_d); break;
        case T_BYTE: SizeElem = sizeof(char); break;
        case T_SHORT: SizeElem = sizeof(short int); break; 
        default:
          cerr << "Error: unknown type of data ... " << endl;
	  exit(-1);
	  break; 
    }
    if (Type != T_FLOAT)
           Buff = (unsigned char *) 
	          alloc_buffer((size_t)(SizeElem*Image.nl()*Image.nc()*sizeof(char)));    
    // Buff = new unsigned char [SizeElem*Image.nl()*Image.nc()];
    switch (Type)
    {
       case T_INT:
        {
	 int *Dat = (int *) Buff;
         for (i=0; i < Image.n_elem(); i++) Dat[i] = (int) Image(i);
        }
        break;
       case T_BYTE:
        {
	 unsigned char *Dat = (unsigned char *) Buff;
         for (i=0; i < Image.n_elem(); i++) Dat[i] = (byte) Image(i);
         }
        break;
       case T_SHORT:
        {
	 unsigned short *Dat = (unsigned short *) Buff;
          for (i=0; i < Image.n_elem(); i++) Dat[i] = (short) Image(i);
        }
	break;
       case T_FLOAT: Buff = (unsigned char *) Image.buffer(); break;
       case T_DOUBLE:
        {
	 double *Dat = (double *) Buff;
	 for (i=0; i < Image.n_elem(); i++) Dat[i] = (double) Image(i);
         }
	break;
       case T_COMPLEX_F:
        {
	 complex_f *Dat = (complex_f  *) Buff;
	 for (i=0; i < Image.n_elem(); i++) Dat[i] = (complex_f)(Image(i),0);
         }
	break;
       case T_COMPLEX_D:
       {
	 complex_d *Dat = (complex_d  *) Buff;
         for (i=0; i < Image.n_elem(); i++)  Dat[i] = (complex_d) ((double)Image(i),0);
        }
       default :
         fprintf (stderr, "Error: type unknown reading %s\n", File_Name);
         exit (-1);
    }
    
    
    buffer = Buff;
    long fp =  (Indi * Nc_io + Indj)*SizeElem;
    fseek(File_Des, fp, SEEK_SET);
    fp = Nc_io - Image.nc();
    fp *= SizeElem;
    for (i=0; i < Image.nl(); i++)
    {
       if (fwrite (buffer, SizeElem, Image.nc(),  File_Des) != (size_t) Image.nc())
       {
          cerr << "Unable to read in file: " << File_Name << endl;
          exit(-1);
       }
       if (i != Image.nl()-1)
       {
           fseek(File_Des, fp, SEEK_CUR);
	   buffer += Image.nc()*SizeElem;
       }
    } 
    // if (Type != T_FLOAT) delete Buff;
    if (Type != T_FLOAT) free_buffer( (char *) Buff);

    fclose (File_Des);
}

/*****************************************************************/

void io_disp::write_block_ima(char *File_Name, Iint & Image,int Indi, int Indj)
{
    FILE *File_Des;
    unsigned char *Buff=NULL;
    unsigned char *buffer;
    int SizeElem;
    
    File_Des = fopen (File_Name, "rb+");
    if (File_Des  == NULL)
    {
     cerr << "Unable to open file: " << File_Name  << endl;
     exit(-1);
    }
    int i;
    read_type(File_Name);
    read_size(File_Name);
    switch (Type)
    {
        case T_INT: SizeElem = sizeof(int); break;  
        case T_FLOAT: SizeElem = sizeof(float); break;
        case T_COMPLEX_F: SizeElem = sizeof(complex_f); break;
        case T_DOUBLE: SizeElem = sizeof(double); break;
        case T_COMPLEX_D: SizeElem = sizeof(complex_d); break;
        case T_BYTE: SizeElem = sizeof(char); break;
        case T_SHORT: SizeElem = sizeof(short int); break; 
        default: 
          cerr << "Error: unknown type of data ... " << endl;
	  exit(-1);
	  break;    }
    if (Type != T_INT)
        Buff = (unsigned char *) 
	          alloc_buffer((size_t)(SizeElem*Image.nl()*Image.nc()*sizeof(char)));    
    // Buff = new unsigned char [SizeElem*Image.nl()*Image.nc()];
    switch (Type)
    {
       case T_INT: Buff = (unsigned char *) Image.buffer(); break;
       case T_BYTE:
        {
	 unsigned char *Dat = (unsigned char *) Buff;
         for (i=0; i < Image.n_elem(); i++) Dat[i] = (byte) Image(i);
         }
        break;
       case T_SHORT:
        {
	 unsigned short *Dat = (unsigned short *) Buff;
          for (i=0; i < Image.n_elem(); i++) Dat[i] = (short) Image(i);
        }
	break;
       case T_FLOAT:
        {
	 float *Dat = (float *) Buff;
	 for (i=0; i < Image.n_elem(); i++) Dat[i] = (float) Image(i);
         }
	break;
       case T_DOUBLE:
        {
	 double *Dat = (double *) Buff;
	 for (i=0; i < Image.n_elem(); i++) Dat[i] = (double) Image(i);
         }
	break;
       case T_COMPLEX_F:
        {
	 complex_f *Dat = (complex_f  *) Buff;
	 for (i=0; i < Image.n_elem(); i++) Dat[i] = (complex_f)(Image(i),0);
         }
	break;
       case T_COMPLEX_D:
       {
	 complex_d *Dat = (complex_d  *) Buff;
         for (i=0; i < Image.n_elem(); i++)  Dat[i] = (complex_d) ((double)Image(i),0);
        }
       default :
         fprintf (stderr, "Error: type unknown reading %s\n", File_Name);
         exit (-1);
    }
    
    
    buffer = Buff;
    long fp =  (Indi * Nc_io + Indj)*SizeElem;
    fseek(File_Des, fp, SEEK_SET);
    fp = Nc_io - Image.nc();
    fp *= SizeElem;
    for (i=0; i < Image.nl(); i++)
    {
       if (fwrite (buffer, SizeElem, Image.nc(),  File_Des) != (size_t) Image.nc())
       {
          cerr << "Unable to read in file: " << File_Name << endl;
          exit(-1);
       }
       if (i != Image.nl()-1)
       {
           fseek(File_Des, fp, SEEK_CUR);
	   buffer += Image.nc()*SizeElem;
       }
    } 
    // if (Type != T_INT) delete Buff;
    if (Type != T_INT) free_buffer( (char *) Buff);
    fclose (File_Des);
}

/*****************************************************************/

void io_disp::create_file  (char *File_Step_Name) 
{
     FILE *File_Des;
    char File_Name[MAXCHAR];
    int L;
    char *Ptr;

    L = strlen (File_Step_Name);
    Ptr = File_Step_Name;

    if ((Ptr[L-1] != 'd') || (Ptr[L-2] != '.'))
    {
         switch (Type)
         {
            case T_INT: 
                sprintf(File_Name,"%s_i_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_BYTE: 
                sprintf(File_Name,"%s_b_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_SHORT: 
                sprintf(File_Name,"%s_s_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_FLOAT:
                sprintf(File_Name,"%s_f_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_COMPLEX_F:
                sprintf(File_Name,"%s_cf_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_DOUBLE:
                sprintf(File_Name,"%s_d_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_COMPLEX_D:
                sprintf(File_Name,"%s_cd_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            default:break;
         }
    }
    else sprintf(File_Name,"%s", File_Step_Name);
 
    File_Des = fopen (File_Name, "wb");
    fclose (File_Des);
}

   
/*****************************************************************/

void io_disp::write_ima  (char *File_Step_Name) 
{
    FILE *File_Des;
    char File_Name[MAXCHAR];
    int L,Nbr=0,Size=Nl_io*Nc_io;
    char *Ptr;

    L = strlen (File_Step_Name);
    Ptr = File_Step_Name;

    if ((Ptr[L-1] != 'd') || (Ptr[L-2] != '.'))
    {
         switch (Type)
         {
            case T_INT: 
                sprintf(File_Name,"%s_i_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_BYTE: 
                sprintf(File_Name,"%s_b_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_SHORT: 
                sprintf(File_Name,"%s_s_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_FLOAT:
                sprintf(File_Name,"%s_f_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_COMPLEX_F:
                sprintf(File_Name,"%s_cf_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_DOUBLE:
                sprintf(File_Name,"%s_d_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            case T_COMPLEX_D:
                sprintf(File_Name,"%s_cd_%d_%d.d",File_Step_Name,Nl_io,Nc_io);
                break;
            default:break;
         }
    }
    else sprintf(File_Name,"%s", File_Step_Name);
 
    File_Des = disp_file_des_out(File_Name);
    
    switch (Type)
    {
        case T_INT: 
              Nbr = fwrite ((void *) Data_int, sizeof (int), Size, File_Des);
              break;
        case T_BYTE: 
              Nbr = fwrite ((void *) Data_byte, sizeof (unsigned char), Size, File_Des);
              break;
        case T_SHORT: 
              Nbr = fwrite ((void *) Data_short, sizeof (short int), Size, File_Des);
              break;
        case T_FLOAT:
               Nbr = fwrite ((void *) Data_float, 
                             sizeof (float), Size, File_Des);
               break;
        case T_COMPLEX_F:
               Nbr = fwrite ((void *) Data_complex_f, 
                            sizeof (complex_f), Size , File_Des);
               break;
        case T_DOUBLE:
               Nbr = fwrite ((void *) Data_double, 
                             sizeof (double), Size, File_Des);
               break;
        case T_COMPLEX_D:
               Nbr = fwrite ((void *) Data_complex_d, 
                            sizeof (complex_d), Size , File_Des);
               break;

            default:break;
     }
    
    if (Nbr <= 0)
    {
        fprintf(stderr,"Error writting in file %s\n",File_Name);
        exit (-1);
    }
    if (File_Des != stdout) fclose (File_Des);

}

/**************************************************************************/

void io_disp::read_ima_int (char *File_Name, Iint & Image)
{
    int i,j;

    read_type(File_Name);

    switch (Type)
    {
       case T_INT:
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Data_int, Nl_io, Nc_io, File_Name, True);
         break;
       case T_BYTE:
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = (int) Data_byte[i*Nc_io+j];
         // delete Data_byte;
	 free_buffer((char *) Data_byte);
         break;
       case T_SHORT:
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = (int) Data_short[i*Nc_io+j];
         // delete Data_short;
	 free_buffer((char *) Data_short);
         break;
       case T_FLOAT:
         fprintf (stderr, "WARNING: convert to INT reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = iround(Data_float[i*Nc_io+j]);
         // delete Data_float;
	 MemFloat.free(Data_float);
         break;
       case T_DOUBLE:
         fprintf (stderr, "WARNING: convert to INT reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = iround(Data_double[i*Nc_io+j]);
 	 free_buffer((char *) Data_double);
         break;
       case T_COMPLEX_F:
         fprintf (stderr, "WARNING: convert to INT reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = iround((Data_complex_f[i*Nc_io+j]).real());
         // delete Data_complex_f;
	 MemCF.free(Data_complex_f);
         break;
       case T_COMPLEX_D:
         fprintf (stderr, "WARNING: convert to INT reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = iround((Data_complex_d[i*Nc_io+j]).real());
         // delete Data_complex_d;
	 free_buffer((char *)  Data_complex_d);
         break;
       default :
         fprintf (stderr, "Error: type unknown reading %s\n", File_Name);
         exit (-1);
    }
}

/**************************************************************************/

void io_disp::read_ima_float (char *File_Name, Ifloat & Image)
{
    int i,j;
    read_type(File_Name);

    switch (Type)
    {
       case T_INT:
         fprintf (stderr, "WARNING: convert to FLOAT reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for (i=0; i < Nl_io; i++)
         for (j=0; j < Nc_io; j++)
                       Image(i,j) = (float) (Data_int[i*Nc_io+j]);
         // delete Data_int;
	 MemInt.free(Data_int);
         break;
       case T_BYTE:
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = (float) Data_byte[i*Nc_io+j];
         // delete Data_byte;
	 free_buffer((char *)  Data_byte);
         break;
       case T_SHORT:
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = (float) Data_short[i*Nc_io+j];
         // delete Data_short;
	 free_buffer((char *) Data_short);
         break;
       case T_FLOAT:
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Data_float, Nl_io, Nc_io, File_Name, True);
         break;
       case T_DOUBLE:
         fprintf (stderr, "WARNING: convert to FLOAT reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = (float) (Data_double[i*Nc_io+j]);
         // delete Data_double;
	 free_buffer((char *) Data_double);
         break;
       case T_COMPLEX_F:
         fprintf (stderr, "WARNING: convert to FLOAT reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = (float)((Data_complex_f[i*Nc_io+j]).real());
         // delete Data_complex_f;
	 MemCF.free(Data_complex_f);
         break;
       case T_COMPLEX_D:
         fprintf (stderr, "WARNING: convert to FLOAT reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = (float)((Data_complex_d[i*Nc_io+j]).real());
         // delete Data_complex_d;
	 free_buffer((char *) Data_complex_d);
         break;
       default :
         fprintf (stderr, "Error: type unknown reading %s\n", File_Name);
         exit (-1);
    }
}

/**************************************************************************/

void io_disp::read_ima_complex_f (char *File_Name, Icomplex_f & Image)
{
    int i,j;
    read_type(File_Name);

    switch (Type)
    {
       case T_INT:
         fprintf (stderr, "WARNING: convert to COMPLEX_F reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for (i=0; i < Nl_io; i++)
         for (j=0; j < Nc_io; j++)
                       Image(i,j) = complex_f (Data_int[i*Nc_io+j], 0);
         // delete Data_int;
	 MemInt.free (Data_int);
         break;
       case T_BYTE:
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = complex_f ((int)Data_byte[i*Nc_io+j], 0);
         // delete Data_byte;
	 free_buffer((char *) Data_byte);
         break;
       case T_SHORT:
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                   Image(i,j) = complex_f ((int)Data_short[i*Nc_io+j], 0);
         // delete Data_short;
	 free_buffer((char *) Data_short);
         break;
       case T_FLOAT:
         fprintf (stderr, "WARNING: convert to COMPLEX_F reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for (i=0; i < Nl_io; i++)
         for (j=0; j < Nc_io; j++)
                       Image(i,j) = complex_f (Data_float[i*Nc_io+j], 0);
         // delete Data_float;
	 MemFloat.free(Data_float);
         break;
       case T_DOUBLE:
         fprintf (stderr, "WARNING: convert to COMPLEX_F reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                       Image(i,j) = complex_f (Data_double[i*Nc_io+j], 0.);
         // delete Data_double;
	 free_buffer((char *) Data_double);
         break;
       case T_COMPLEX_F:
         read_size (File_Name);
         read_ima (File_Name);
	 Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                     Image(i,j) = complex_f(
                                     (Data_complex_f[i*Nc_io+j]).real(),
                                     (Data_complex_f[i*Nc_io+j]).imag());
         // delete Data_complex_d;
	 free_buffer((char *) Data_complex_d);	 	 
         //Image.Data = Data_complex_f;
         //Image.Nl = Nl_io;
         //Image.Nc = Nc_io;
         //strcpy(Image.Name_Obj,File_Name);
         //break;
       case T_COMPLEX_D:
         fprintf (stderr, "WARNING: convert to COMPLEX_F reading %s\n", File_Name);
         read_size (File_Name);
         read_ima (File_Name);
         Image.alloc(Nl_io, Nc_io, File_Name);
         for ( i=0; i < Nl_io; i++)
         for ( j=0; j < Nc_io; j++)
                     Image(i,j) = complex_f(
                                     (Data_complex_d[i*Nc_io+j]).real(),
                                     (Data_complex_d[i*Nc_io+j]).imag());
         // delete Data_complex_d;
	 free_buffer((char *) Data_complex_d);
         break;
       default :
         fprintf (stderr, "Error: type unknown reading %s\n", File_Name);
         exit (-1);
    }
}

/**************************************************************************/

void io_disp::write_ima_int (char *File_Name, Iint & Image)
{
    Type = T_INT;
    Data_int = Image.buffer();
    Nl_io = Image.nl();
    Nc_io = Image.nc();
    int Size = Nl_io*Nc_io;
    write_ima (File_Name);   
    int L,i,j;
    char *Ptr;

    L = strlen (File_Name);
    Ptr = File_Name;
    if ((Ptr[L-1] != 'd') || (Ptr[L-2] != '.'))
    {
      read_type(File_Name);
      switch (Type)
      {
        case T_INT: 
               Data_int = Image.buffer();  // MemInt.alloc (Size);
               break;
        case T_FLOAT:
               Data_float = MemFloat.alloc (Size);
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                   Data_float[i*Nc_io+j] = (float) Image(i,j);
               break;
        case T_COMPLEX_F:
               Data_complex_f = MemCF.alloc (Size);
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                Data_complex_f[i*Nc_io+j] = (complex_f) ((float)Image(i,j), 0.);
               break;
        case T_DOUBLE:
               // Data_double = new double[Size];
	       Data_double = (double *) alloc_buffer(Size*sizeof(double));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                   Data_double[i*Nc_io+j] = (double) Image(i,j);
               break;
        case T_COMPLEX_D:
               // Data_complex_d = new complex_d[Size];
	       Data_complex_d = (complex_d *) alloc_buffer(Size*sizeof(complex_d));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                Data_complex_d[i*Nc_io+j] = (complex_d)((double)Image(i,j), 0.);
               break;
        case T_BYTE: 
              { //float Min,Max;
               //Min = min (Image);
               //Max = max (Image);
               // Data_byte = new unsigned char[Size];
	       Data_byte  = (unsigned char *) alloc_buffer(Size*sizeof(char));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                 Data_byte [i*Nc_io+j] = (unsigned char) Image(i,j);
                
                  // Data_byte [i*Nc_io+j] = (unsigned char) 
                  //               ((Image(i,j) - Min) / (Max-Min)*255.+0.5);
               }
               break;
        case T_SHORT: 
              { //float Min,Max;
               //Min = min (Image);
               //Max = max (Image);
               // Data_short = new short int [Size];
	       Data_short   = (short *) alloc_buffer(Size*sizeof(short));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
               { 
                 //if (Max > 32767.)
                 //  Data_short [i*Nc_io+j] = (short int) 
                 //       ((Image(i,j) - Min) / (Max-Min)*32767.+0.5); 
                 //else 
                 Data_short [i*Nc_io+j] = (short int) Image(i,j);
               }
               }
               break;
        default:break;
      }
    }
    else Data_int = Image.buffer();  
    write_ima (File_Name);   
    switch (Type)
    {
        case T_FLOAT:  MemFloat.free (Data_float); break;
        case T_COMPLEX_F: MemCF.free(Data_complex_f); break;
        case T_DOUBLE: free_buffer((char *) Data_double); break;
        case T_COMPLEX_D: free_buffer((char *) Data_complex_d); break;
        case T_BYTE: free_buffer((char *) Data_byte); break;
        case T_SHORT: free_buffer((char *)  Data_short); break;
        default:break;
    }
}

/**************************************************************************/

void io_disp::write_ima_float (char *File_Name, Ifloat & Image)
{
    Type = T_FLOAT;
    Nl_io = Image.nl();
    Nc_io = Image.nc();    
    int L,i,j;
    char *Ptr;
    int Size = Nl_io*Nc_io;

    L = strlen (File_Name);
    Ptr = File_Name;
    if ((Ptr[L-1] == 'd') && (Ptr[L-2] == '.'))
    {
      read_type(File_Name);      
      switch (Type)
      {
        case T_INT: 
               Data_int =  MemInt.alloc (Size);
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                   Data_int[i*Nc_io+j] = (int) Image(i,j);
               break;
        case T_FLOAT:
               Data_float = Image.buffer();
               break;
        case T_COMPLEX_F:
               Data_complex_f = MemCF.alloc (Size);
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                Data_complex_f[i*Nc_io+j] = (complex_f) ((float)Image(i,j), 0.);
               break;
        case T_DOUBLE:
               // Data_double = new double[Size];
	       Data_double    = (double *) alloc_buffer(Size*sizeof( double));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                   Data_double[i*Nc_io+j] = (double) Image(i,j);
               break;
        case T_COMPLEX_D:
               // Data_complex_d = new complex_d[Size];
	        Data_complex_d    = (complex_d *) alloc_buffer(Size*sizeof( complex_d ));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                Data_complex_d[i*Nc_io+j] = (complex_d)((double)Image(i,j), 0.);
               break;
        case T_BYTE: 
              { //float Min,Max;
               //Min = min (Image);
               //Max = max (Image);
               // Data_byte = new unsigned char[Size];
	        Data_byte    = (unsigned char *) alloc_buffer(Size*sizeof(char));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                   Data_byte [i*Nc_io+j] = (unsigned char) Image(i,j);
                 //  Data_byte [i*Nc_io+j] = (unsigned char) 
                 //                ((Image(i,j) - Min) / (Max-Min)*255.+0.5);
               }
               break;
        case T_SHORT: 
              { //float Min,Max;
               // Data_short = new short int [Size];
	        Data_short    = (short *) alloc_buffer(Size*sizeof(short));
               //Min = min (Image);
               //Max = max (Image);
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
               { 
                 //if (Max > 32767.)
                 //  Data_short [i*Nc_io+j] = (short int) 
                 //       ((Image(i,j) - Min) / (Max-Min)*32767.+0.5); 
                 //else 
                 Data_short [i*Nc_io+j] = (short int)Image(i,j);
               }}
               break;
        default:break;
      }
    }
    else Data_float = Image.buffer(); 
    write_ima (File_Name);
    switch (Type)
    {
        case T_INT:  MemInt.free (Data_int); break;
        case T_COMPLEX_F: MemCF.free(Data_complex_f); break;
        case T_DOUBLE: free_buffer((char *) Data_double); break;
        case T_COMPLEX_D: free_buffer((char *) Data_complex_d); break;
        case T_BYTE: free_buffer((char *) Data_byte); break;
        case T_SHORT: free_buffer((char *)  Data_short); break;
        default:break;
    }
}

/**************************************************************************/

void io_disp::write_ima_complex_f (char *File_Name, Icomplex_f & Image)
{
    Type = T_COMPLEX_F;
    Nl_io = Image.nl();
    Nc_io = Image.nc();    
    int L,i,j;
    char *Ptr;
    int Size = Nl_io*Nc_io;

    L = strlen (File_Name);
    Ptr = File_Name;
    if ((Ptr[L-1] == 'd') && (Ptr[L-2] != '.'))
    {
      read_type(File_Name);
      switch (Type)
      {
        case T_INT: 
               Data_int =  MemInt.alloc (Size);
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                   Data_int[i*Nc_io+j] = (int) Image(i,j).real();
               break;
        case T_FLOAT:
               Data_float = MemFloat.alloc (Size);
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                   Data_float[i*Nc_io+j] = (float) Image(i,j).real();
               break;
        case T_COMPLEX_F:
	        Data_complex_f = (complex_f *)alloc_buffer(Size*sizeof(complex_f));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                Data_complex_f[i*Nc_io+j] = (complex_f)((double)Image(i,j).real(),  Image(i,j).imag());	               
	       //Data_complex_f = Image.Data;
               break;
        case T_DOUBLE:
               // Data_double = new double[Size];
	        Data_double   = (double *) alloc_buffer(Size*sizeof(double ));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                   Data_double[i*Nc_io+j] = (double) Image(i,j).real();
               break;
        case T_COMPLEX_D:
               // Data_complex_d = new complex_d[Size];
	        Data_complex_d = (complex_d *) alloc_buffer(Size*sizeof(complex_d));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                Data_complex_d[i*Nc_io+j] = (complex_d)((double)Image(i,j).real(),  Image(i,j).imag());
               break;
        case T_BYTE: 
               // Data_byte = new unsigned char[Size];
	       Data_byte  = (unsigned char *) alloc_buffer(Size*sizeof(char ));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                   Data_byte [i*Nc_io+j] = (unsigned char) Image(i,j).real();
               break;
        case T_SHORT: 
               // Data_short = new short int [Size];
	        Data_short  = (short *) alloc_buffer(Size*sizeof( short ));
               for ( i=0; i < Nl_io; i++)
               for ( j=0; j < Nc_io; j++)
                   Data_byte [i*Nc_io+j] = (short int) Image(i,j).real();
               break;
        default:break;
      }
    }
    else  Data_complex_f = Image.buffer(); 
    write_ima (File_Name);
    switch (Type)
    {
        case T_INT:  MemInt.free (Data_int); break;
        case T_FLOAT:  MemFloat.free (Data_float); break;
        case T_DOUBLE: free_buffer((char *) Data_double); break;
        case T_COMPLEX_D: free_buffer((char *) Data_complex_d); break;
        case T_BYTE: free_buffer((char *) Data_byte); break;
        case T_SHORT: free_buffer((char *)  Data_short); break;
        default:break;
    }    
}

/**************************************************************************/
