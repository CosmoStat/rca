#include "NR.h"

int ndatat=0;	/* defining declaration */
float *xt=0,*yt=0,aa=0.0,abdevt=0.0;	/* defining declaration */

void medfit(float *x,float *y,int ndata,float sig[],int mwt,float *a,float *b,float *abdev)
{
	int j, nbr=0;
	float bb,b1,b2,del,f,f1,f2,sigb,temp;
	float sx=0.0,sy=0.0,sxy=0.0,sxx=0.0,chisq=0.0;

	ndatat=ndata;
	xt=x;
	yt=y;
	if (mwt) {
	  for (j=0;j<ndata;j++) {
		 sx += x[j]*sig[j];
		 sy += y[j]*sig[j];
		 sxy += x[j]*y[j]*sig[j];
		 sxx += x[j]*x[j]*sig[j];
		 if (sig[j]) nbr++;
	  }
	}
	else {
	  for (j=0;j<ndata;j++) {
		 sx += x[j];
		 sy += y[j];
		 sxy += x[j]*y[j];
		 sxx += x[j]*x[j];
	  }
	  nbr = ndata;
	}
	del=nbr*sxx-sx*sx;
	aa=(sxx*sy-sx*sxy)/del;
	bb=(nbr*sxy-sx*sy)/del;
	for (j=0;j<ndata;j++)
		chisq += (temp=sig[j]*(y[j]-(aa+bb*x[j])),temp*temp);
	sigb=sqrt(chisq/del);
	b1=bb;
	f1=rofunc(b1);
	b2=bb+((f1 > 0.0) ? fabs(3.0*sigb) : -fabs(3.0*sigb));
	f2=rofunc(b2);
	while (f1*f2 > 0.0) {
		bb=2.0*b2-b1;
		b1=b2;
		f1=f2;
		b2=bb;
		f2=rofunc(b2);
	}
	sigb=0.01*sigb;
	while (fabs(b2-b1) > sigb) {
		bb=0.5*(b1+b2);
		if (bb == b1 || bb == b2) break;
		f=rofunc(bb);
		if (f*f1 >= 0.0) {
			f1=f;
			b1=bb;
		} else {
			f2=f;
			b2=bb;
		}
	}
	*a=aa;
	*b=bb;
	*abdev=abdevt/nbr;
}
