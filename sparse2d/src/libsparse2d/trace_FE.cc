

#include "Mr_FE.h"
#include "IM_IO.h"

extern char *OptArg;

enum TestType {
   AUTOCONV_HISTO, 
   REDUCED_HISTO,
   DISTRIB,
   HISTO_THRESHOLD
};
unsigned int xNbTestType = 4;

unsigned int xNbAutoConv = 25;
TestType xTestType = AUTOCONV_HISTO;
bool xHelp = false;
double xEpsilon = 1e-03;



//---------------------------------------------------------------------------------
char* 
int2string_filter_type( TestType type ) {
//---------------------------------------------------------------------------------
   switch (type) {
      case AUTOCONV_HISTO : 
         return( "auto-convolued histogram" );
      case REDUCED_HISTO :
         return( "reduced auto-convolued histogram" );
      case DISTRIB :
         return( "distribution function" );
      case HISTO_THRESHOLD :
         return( "search for histogram threshold" );
      default :
         std::cout << "unknown test type" << std::endl;
	 exit(-1);
   }
}

//---------------------------------------------------------------------------------
static void 
usage(char *argv[]) {
//---------------------------------------------------------------------------------
    fprintf( OUTMAN, "Usage: %s \n\n", argv[0] );
    fprintf( OUTMAN, "   where options =  \n" );
    manline(); 
    
    fprintf( OUTMAN, "         [-A autoconv number]\n" );
    fprintf( OUTMAN, "             default is %d\n", xNbAutoConv );
    manline(); 

    fprintf( OUTMAN, "         [t test type]\n" ); 
    for( unsigned int i=0; i<xNbTestType; ++i )
       fprintf( OUTMAN, "             %d : %s \n", i+1, 
                int2string_filter_type( (TestType)i ) );
    fprintf( OUTMAN, "             default is %d : %s \n", 
             1, int2string_filter_type( xTestType ) );
    manline(); 
    
    fprintf( OUTMAN, "         [-e epsilon]\n" );
    fprintf( OUTMAN, "             default is %d\n", xEpsilon );
    manline(); 

    fprintf( OUTMAN, "         [-h]\n" );
    fprintf( OUTMAN, "             help mode \n" );
    manline(); 
    
   
    manline();
    manline();
    manline();
    exit(-1);
}


//---------------------------------------------------------------------------------
static void 
filtinit(int argc, char *argv[]) {
//---------------------------------------------------------------------------------

    int c;
    while( ( c = GetOpt( argc, argv, "A:t:e:h" ) ) != -1 ) {
   	
	switch (c) {
	
	   case 'A' : 
	       if( sscanf( OptArg, "%d", &xNbAutoConv) != 1) {
                  fprintf( OUTMAN, "bad autoconv number: %s\n", OptArg);
                  exit(-1);
	       }
	       break;
	
	   case 't' : 
	       if( sscanf( OptArg, "%d", &c) != 1) {
                  fprintf( OUTMAN, "bad test type: %s\n", OptArg );
                  exit(-1);
	       }
	       if( c>0 && c<=(int)xNbTestType ) xTestType = TestType(c-1);
	       else {
                  fprintf( OUTMAN, "bad test type: %s !!! \n", OptArg );
		      manline(); usage( argv );
	       }
	       break;
	       
	   case 'e' : 
	       if( sscanf( OptArg, "%lf", &xEpsilon) != 1) {
                  fprintf( OUTMAN, "bad epsilon: %s\n", OptArg);
                  exit(-1);
	       }
	       break;
	       
	   case 'h' :
	      xHelp = true;
	      break;
 
       }
    }

}

//---------------------------------------------------------------------------------
void autoconv_histo_test( FewEvent* pFE ) {
//---------------------------------------------------------------------------------
   pFE->bspline_histo_2D( True );
   pFE->histo_convolution( True );
}


//---------------------------------------------------------------------------------
void reduced_histo_test( FewEvent* pFE ) {
//---------------------------------------------------------------------------------
   pFE->bspline_histo_2D();
   pFE->histo_convolution();
   pFE->histo_reduced_wavelet_coef( True );
}


//---------------------------------------------------------------------------------
void distrib_test( FewEvent* pFE ) {
//---------------------------------------------------------------------------------
   pFE->bspline_histo_2D();
   pFE->histo_convolution();
   pFE->histo_reduced_wavelet_coef();
   pFE->histo_distribution( True );
}


//---------------------------------------------------------------------------------
void histo_threshold_test( double Epsilon, FewEvent* pFE ) {
//---------------------------------------------------------------------------------
   pFE->bspline_histo_2D();
   pFE->histo_convolution();
   pFE->histo_reduced_wavelet_coef();
   pFE->histo_distribution();
   pFE->histo_threshold( Epsilon, True );
}


//---------------------------------------------------------------------------------
int 
main(int argc, char *argv[]) {
//---------------------------------------------------------------------------------
    
   filtinit(argc, argv);
   
   if( xHelp ) usage( argv );
  
   std::cout << "test of FewEvent class" << std::endl;
   std::cout << "  autoconv number : " << xNbAutoConv << std::endl;
   std::cout << "  test type       : " 
             << int2string_filter_type( xTestType ) << std::endl;
    
   FewEvent* pFE = new FewEvent( xNbAutoConv );
   
   switch( xTestType ) {
   
      case AUTOCONV_HISTO : 
         autoconv_histo_test( pFE );
         break; 
      case REDUCED_HISTO : 
         reduced_histo_test( pFE );
         break;  
      case DISTRIB : 
         distrib_test( pFE );
         break; 
      case HISTO_THRESHOLD : 
         std::cout << "  epsilon         : " << xEpsilon << std::endl;
         histo_threshold_test( xEpsilon, pFE );
         break;
      default :  
         std::cout << "unknown test type" << std::endl;
	 exit(-1);
   }
}



