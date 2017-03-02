#include "NR.h"

void mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **covar, float **alpha, float *chisq,
	void (*funcs)(float, float [], float *, float [], int), float *alamda)
{
	void covsrt(float **covar, int ma, int ia[], int mfit);
	void gaussj(float **a, int n, float **b, int m);
	void mrqcof(float x[], float y[], float sig[], int ndata, float a[],
		int ia[], int ma, float **alpha, float beta[], float *chisq,
		void (*funcs)(float, float [], float *, float [], int));
	int j,k,l,m;
	static int mfit;
	static float ochisq,*atry,*beta,*da,**oneda;

	if (*alamda < 0.0) {
		atry=vector(1,ma);
		beta=vector(1,ma);
		da=vector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]) mfit++;
		oneda=matrix(1,mfit,1,1);
		*alamda=0.001;
		mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
			for (j++,k=0,m=1;m<=ma;m++) {
				if (ia[m]) {
					k++;
					covar[j][k]=alpha[j][k];
				}
			}
			covar[j][j]=alpha[j][j]*(1.0+(*alamda));
			oneda[j][1]=beta[j];
		}
	}
	gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		covsrt(covar,ma,ia,mfit);
		free_matrix(oneda,1,mfit,1,1);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		free_vector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
	mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				for (j++,k=0,m=1;m<=ma;m++) {
					if (ia[m]) {
						k++;
						alpha[j][k]=covar[j][k];
					}
				}
				beta[j]=da[j];
				a[l]=atry[l];
			}
		}
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
}
#undef NRANSI

void old_mrqmin(float x[],float y[],float sig[],int ndata,float a[],int ma,int lista[],int mfit,float **covar,float **alpha,float *chisq,void (*funcs)(float,float *,float *,float *,int),float *alamda)
{
	int k,kk,j,ihit;
	static float *da,*atry,**oneda,*beta,ochisq;

	if (*alamda < 0.0) {
		oneda=matrix(1,mfit,1,1);
		atry=vector(1,ma);
		da=vector(1,ma);
		beta=vector(1,ma);
		kk=mfit+1;
		for (j=1;j<=ma;j++) {
			ihit=0;
			for (k=1;k<=mfit;k++)
				if (lista[k] == j) ihit++;
			if (ihit == 0)
				lista[kk++]=j;
			else if (ihit > 1) nrerror("Bad LISTA permutation in MRQMIN-1");
		}
		if (kk != ma+1) nrerror("Bad LISTA permutation in MRQMIN-2");
		*alamda=0.001;
		old_mrqcof(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,chisq,funcs);
		ochisq=(*chisq);
	}
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(*alamda));
		oneda[j][1]=beta[j];
	}
	gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++)
		da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		covsrt(covar,ma,lista,mfit);
		free_vector(beta,1,ma);
		free_vector(da,1,ma);
		free_vector(atry,1,ma);
		free_matrix(oneda,1,mfit,1,1);
		return;
	}
	for (j=1;j<=ma;j++) atry[j]=a[j];
	for (j=1;j<=mfit;j++)
		atry[lista[j]] = a[lista[j]]+da[j];
	old_mrqcof(x,y,sig,ndata,atry,ma,lista,mfit,covar,da,chisq,funcs);
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=1;j<=mfit;j++) {
			for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
			beta[j]=da[j];
			a[lista[j]]=atry[lista[j]];
		}
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
	return;
}
