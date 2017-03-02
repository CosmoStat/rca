#include "NR.h"


#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

#define GET_PSUM for (j=0;j<ndim;j++) { for (i=0,sum=0.0;i<mpts;i++)\
						sum += p[i][j]; psum[j]=sum;}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define MAX(a,b) ( ((a) >= (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) >= (b)) ? (b) : (a) )

void amoeba(float **p,float y[],int ndim,float ftol,float (*funk)(float [],float [],float [],float [],int ,char []),float time[], float data[], float mask[],int size, char *method,int *nfunk, int NMAX)
{
	int i,j,ilo,ihi,inhi,mpts=ndim+1;
	float ytry,ysave,sum,rtol,*psum, swap;
	int compt=0;

	psum=vector(1,ndim);


	*nfunk=0;
	GET_PSUM

	for (;;) {
		ilo=0;
		compt++;
		/*
		printf("+++ %d iteration ######\n",compt);
		printf("Valeurs au sommets: \n");
		for (i=0;i<mpts;i++) {
		  printf("%.8f  ",y[i]);
		}
		printf("\n");
		*/

		ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
		for (i=0;i<mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if ((y[i] > y[inhi]) && (i != ihi)) inhi=i;
		}
		
		//printf("   le + mauvais est   %d\n",ihi);
		//printf("le 2e + mauvais est   %d\n",inhi);
		//printf("    le meilleur est   %d\n",ilo);
		

		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
		//printf("tolerance relative est %f \n ",rtol);

		if ((rtol < ftol) || (compt >= NMAX)) {
		  SWAP(y[0],y[ilo])
		  for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i])
		  break;
		}		
		/* if (*nfunk >= NMAX) nrerror("Too many iterations in AMOEBA"); */
		*nfunk += 2;

		ytry=amotry(p,y,psum,ndim,funk,time,data,mask,size,method,ihi,-1.0);
		//printf("Reflex: %.8f ,%f\n", ytry,p[ihi][3]);
		if (ytry <= y[ilo]){
			ytry=amotry(p,y,psum,ndim,funk,time,data,mask,size,method,ihi,2.0);
			//printf("Extrapolate: %.8f,%f\n", ytry,p[ihi][3]);
		}
		else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			ytry=amotry(p,y,psum,ndim,funk,time,data,mask,size,method,ihi,0.5);
			//if (ytry < ysave) printf("Contraction:  %.8f,%f\n", ytry,p[ihi][3]);
			if (ytry >= ysave) {
			  
			  for (i=0;i<mpts;i++) {
				 if (i != ilo) {
					for (j=0;j<ndim;j++) {
					  psum[j]=0.5*(p[i][j]+p[ilo][j]);
					  p[i][j]=psum[j];
					}
					y[i]=(*funk)(psum,time,data,mask,size,method);
				 }
			  }
			  //printf("Shrinkage:  %.8f,%f\n", ytry,p[ihi][3]);
			  *nfunk += ndim;
			  GET_PSUM
			}
		} else --(*nfunk);
	}
	free_vector(psum,1,ndim);
}

float amotry(float **p,float *y,float *psum,int ndim,float ( *funk)(float [],float [],float [],float [],int ,char []),float time[], float data[], float mask[],int size, char *method,int ihi,float fac)
{
	int j;
	float fac1,fac2,ytry,*ptry;


	ptry=vector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=0;j<ndim;j++) ptry[j]=MAX(0.0,psum[j]*fac1-p[ihi][j]*fac2);
	ytry=(*funk)(ptry,time, data,mask,size,method);

	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=0;j<ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=MAX(ptry[j],0.0);
		}
	}
	free_vector(ptry,1,ndim);
	return ytry;
}



void amoeba_overshoot(float **p,float y[],int ndim,float ftol,float (*funk)(float [],float [],float [],float [],int ,char []),float time[], float data[], float mask[],int size, char *method,int *nfunk, int NMAX)
{
	int i,j,ilo,ihi,inhi,mpts=ndim+1;
	float ytry,ysave,sum,rtol,*psum, swap;
	int compt=0;

	psum=vector(1,ndim);


	*nfunk=0;
	GET_PSUM

	for (;;) {
		ilo=0;
		compt++;
		/*
		printf("+++ %d iteration ######\n",compt);
		printf("Valeurs au sommets: \n");
		for (i=0;i<mpts;i++) {
		  printf("%.8f  ",y[i]);
		}
		printf("\n");
		*/

		ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
		for (i=0;i<mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if ((y[i] > y[inhi]) && (i != ihi)) inhi=i;
		}
		
		//printf("   le + mauvais est   %d\n",ihi);
		//printf("le 2e + mauvais est   %d\n",inhi);
		//printf("    le meilleur est   %d\n",ilo);
		

		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
		//printf("tolerance relative est %f \n ",rtol);

		if ((rtol < ftol) || (compt >= NMAX)) {
		  SWAP(y[0],y[ilo])
		  for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i])
		  break;
		}		
		/* if (*nfunk >= NMAX) nrerror("Too many iterations in AMOEBA"); */
		*nfunk += 2;

		ytry=amotry_overshoot(p,y,psum,ndim,funk,time,data,mask,size,method,ihi,-1.0);
		//printf("Reflex: %.8f ,%f\n", ytry,p[ihi][3]);
		if (ytry <= y[ilo]){
			ytry=amotry_overshoot(p,y,psum,ndim,funk,time,data,mask,size,method,ihi,2.0);
			//printf("Extrapolate: %.8f,%f\n", ytry,p[ihi][3]);
		}
		else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			ytry=amotry_overshoot(p,y,psum,ndim,funk,time,data,mask,size,method,ihi,0.5);
			//if (ytry < ysave) printf("Contraction:  %.8f,%f\n", ytry,p[ihi][3]);
			if (ytry >= ysave) {
			  
			  for (i=0;i<mpts;i++) {
				 if (i != ilo) {
					for (j=0;j<ndim;j++) {
					  psum[j]=0.5*(p[i][j]+p[ilo][j]);
					  p[i][j]=psum[j];
					}
					y[i]=(*funk)(psum,time,data,mask,size,method);
				 }
			  }
			  //printf("Shrinkage:  %.8f,%f\n", ytry,p[ihi][3]);
			  *nfunk += ndim;
			  GET_PSUM
			}
		} else --(*nfunk);
	}
	free_vector(psum,1,ndim);
}

float amotry_overshoot(float **p,float *y,float *psum,int ndim,float ( *funk)(float [],float [],float [],float [],int ,char []),float time[], float data[], float mask[],int size, char *method,int ihi,float fac)
{
	int j;
	float fac1,fac2,ytry,*ptry;


	ptry=vector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=0;j<ndim;j++) ptry[j]=MAX(0.0,psum[j]*fac1-p[ihi][j]*fac2);
	ptry[2]=MIN(0.999,ptry[2]);
	ptry[3]=MIN(0.999,ptry[3]);

	ytry=(*funk)(ptry,time, data,mask,size,method);

	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=0;j<ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=MAX(ptry[j],0.0);
			if ((j==2) || (j==3))
			  p[ihi][j]= MIN(0.999,p[ihi][j]);
		}
	}
	free_vector(ptry,1,ndim);
	return ytry;
}

#undef ALPHA
#undef BETA
#undef GAMMA
