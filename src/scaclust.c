/*****************************************************************/
/*
 *  Copyright (C)2000 Evgenia Dimitriadou
 *               
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdlib.h>
#include <math.h>
#include "R_ext/Linpack.h"
#include "R.h"
#include "S.h"

int  subcommon(int *xrows, int *xcols, double *x, int *ncenters,
	       double *centers, int *itermax, int *iter,
	       int *verbose,  double *U, double *UANT, 
	       double *beta, double *taf, double *theta, double *ermin)
{
    
    typedef enum {FALSE,TRUE} bool;
    bool control;
    double sum1, sum2;
    double  serror, temp, conv, epsi1;

        
    int k, i, i2, col, col1, col2, info, job;
    double hta, thsigma, product, summea, summeb, summeab;
    double f;


    /*pointers*/
    double *diafmatrix, *gin, *scatter, *scattermatrix, *a, *b;
    double *determinant, *product1;
    double *det;
    int *t;
    
    diafmatrix = (double *) R_alloc((*xrows)*(*xcols), sizeof(double));
    gin = (double *) R_alloc((*xcols)*(*xcols)*(*ncenters), sizeof(double));
    scatter = (double *) R_alloc((*xcols)*(*xcols)*(*ncenters), sizeof(double));
    scattermatrix = (double *) R_alloc((*xcols)*(*xcols), sizeof(double));
    a = (double *) R_alloc((*ncenters), sizeof(double));
    b = (double *) R_alloc((*xrows)*(*ncenters), sizeof(double));
    product1 = (double *) R_alloc((*xcols), sizeof(double));
    t = (int *) R_alloc((*xrows)*(*ncenters), sizeof(int));
    determinant = (double *) R_alloc((*ncenters), sizeof(double));
    det = (double *) R_alloc((2), sizeof(double));

    
/*  *ermin=0.0;*/
    serror=0.0;
    job=11;    
    
    f=2.0;
    hta=0.0;
    thsigma=0.0;
    epsi1=0.002;
    conv=0.0;
    
    /*initialize T_i*/
    
    for(i=0;i<*ncenters;i++)
	for(k=0;k<*xrows;k++)   
	    t[k+(*xrows)*i]=0;

    if (*iter!=0) {/*not predict*/
    /*update UANT*/
    for(i=0;i<*ncenters;i++){
	for(k=0;k<*xrows;k++){
	    UANT[k+(*xrows)*i]=U[k+(*xrows)*i];
	    /*  Rprintf("UANT: k: %5.2d, i:%5.2d,U:%5.2f\n",i,k,U[k+(*xrows)*i]);*/
	}}
    

    /*NEW CENTERS*/
    
    for(i=0;i<*ncenters;i++)
    {
	sum2=0.0;
	
	for(col=0;col<*xcols;col++)
	    centers[i+(*ncenters)*col]=0.0;
	for(k=0;k<*xrows;k++)
	{
	    temp=pow(UANT[k+(*xrows)*i],f);
	    sum2=sum2+temp;
	    
	    for(col=0;col<*xcols;col++)
	    {
		centers[i+(*ncenters)*col]+= temp*x[k+(*xrows)*col];
		
	    }
	}
	for(col=0;col<*xcols;col++)
	    centers[i+(*ncenters)*col]/=sum2;
	
    }
    }/*not predict*/
    
    /*initialize*/
    for(i=0;i<*ncenters;i++){
	a[i]=0.0;
	for(col1=0;col1<*xcols;col1++){
	    product1[col1]=0.0;
	    for(col2=0;col2<*xcols;col2++){
		gin[col1+(*xcols)*col2+(*xcols)*(*ncenters)*i]=0.0;
		scattermatrix[col1+(*xcols)*col2]=0.0;}} 
	for(k=0;k<*xrows;k++){
	    for(col=0;col<*xcols;col++)
		diafmatrix[k+(*xrows)*col]=0.0;
	    b[k+(*xrows)*i]=0.0;}}
    /*end initialize*/
    
    
    /*SCATTER MATRIX*/
    for(i=0;i<*ncenters;i++){
	
	for(col1=0;col1<*xcols;col1++)
	    for(col2=0;col2<*xcols;col2++)
		scatter[col1+(*xcols)*col2+(*xcols)*(*ncenters)*i]=0.0;
	
	/*sum2=0.0;*/
	sum1=0.0;
	for(k=0;k<*xrows;k++)
	{
	    temp=pow(UANT[k+(*xrows)*i],f);
	    
	    sum1+=UANT[k+(*xrows)*i];
	    /*sum2=sum2+temp;*/
	    
	    
	    for(col=0;col<*xcols;col++)
	    {
		diafmatrix[k+(*xrows)*col] =
		  x[k+(*xrows)*col]-centers[i+(*ncenters)*col];
		
	    }
	    
	    for(col1=0;col1<*xcols;col1++){
		for(col2=0;col2<*xcols;col2++){		
		    gin[col1+(*xcols)*col2+(*xcols)*(*ncenters)*i]=diafmatrix[k+(*xrows)*col1]*diafmatrix[k+(*xrows)*col2];
		    scatter[col1+(*xcols)*col2+(*xcols)*(*ncenters)*i] += temp* gin[col1+(*xcols)*col2+(*xcols)*(*ncenters)*i];
		} /*col2*/
	    }
	    
	}/*xrows*/
	
	for(col1=0;col1<*xcols;col1++)
	    for(col2=0;col2<*xcols;col2++)
		scatter[col1+(*xcols)*col2+(*xcols)*(*ncenters)*i]/=sum1;
	
    }
    
    
    /*Determinant and inverse of scatter matrix*/
    for(i=0;i<*ncenters;i++){
	for(col1=0;col1<*xcols;col1++){
	    for(col2=0;col2<*xcols;col2++){	
		scattermatrix[col1+(*xcols)*col2] =
			scatter[col1+(*xcols)*col2+(*xcols)*(*ncenters)*i];
		
	    }
	}
	
	F77_NAME(dpofa)(scattermatrix, xcols, xcols, &info);      
	F77_NAME(dpodi)(scattermatrix, xcols, xcols, det, &job);	
	
	
	for(col1=1;col1<*xcols;col1++){
	    for(col2=0;col2<col1;col2++){
		scattermatrix[col1+(*xcols)*col2]=scattermatrix[col2+(*xcols)*col1];}}
	
	for (col1=0;col1<*xcols;col1++){
	    for (col2=0;col2<*xcols;col2++){
		scatter[col1+(*xcols)*col2+(*xcols)*(*ncenters)*i]=scattermatrix[col1+(*xcols)*col2];}}
	determinant[i]=det[0]*pow(10,det[1]);
	
      
      /*A_t expression*/
      /*sum2=0.0;*/
      sum1=0.0;
      for(k=0;k<*xrows;k++)
      {
	  /*temp=pow(UANT[k+(*xrows)*i],f);
	  sum2=sum2+temp;}*/
	  sum1+=UANT[k+(*xrows)*i];}
      if (*beta != 0.0){ 
	  hta=pow(sum1,(*xcols)*(*beta)-(*taf+1));
	  thsigma=pow(theta[i]*determinant[i],*beta);
	  a[i]= 0.5*((*taf)/(*beta))*hta*thsigma;}
      else
	  a[i]=-0.5*log(determinant[i]);
	  
    
      /*B_it expression*/
       for(k=0;k<*xrows;k++){
	   for (col2=0;col2<*xcols;col2++){
	       product1[col2]=0.0;
	       for (col1=0;col1<*xcols;col1++){
		   product1[col2]+=diafmatrix[k+(*xrows)*col1]*scatter[col1+(*xcols)*col2+(*xcols)*(*ncenters)*i];}}
	   product=0.0;
	   for (col=0;col<*xcols;col++){
	       product+=product1[col]*diafmatrix[k+(*xrows)*col];}
	   if (*beta != 0.0) 
	       b[k+(*xrows)*i]=hta*thsigma*product;
	 else 
	     b[k+(*xrows)*i]=product; 
       }
	 
    }

    /*Membership*/
    control=FALSE;
    while (!control){
	control=TRUE;
	for(i=0;i<*ncenters;i++){
	    
	    for(k=0;k<*xrows;k++){
		summeb=0.0;
		summea=0.0;
		summeab=0.0;
		for(i2=0;i2<*ncenters;i2++){
		    if (t[k+(*xrows)*i2]>=0){
			summeb+= 1/b[k+(*xrows)*i2];
			summea+= a[i2];
			summeab+=a[i2]/b[k+(*xrows)*i2];}}
		if (t[k+(*xrows)*i]>=0)
		    U[k+(*xrows)*i]=((1/b[k+(*xrows)*i])/summeb)-(1/b[k+(*xrows)*i])*((summeab/summeb)-a[i]);
/*		    U[k+(*xrows)*i]=1/(b[k+(*xrows)*i]*summeb);*/
		/*Rprintf("UANT2: k: %5.2d, i:%5.2d, U:%5.2f,t:%5d\n",i,k,U[k+(*xrows)*i],t[k+(*xrows)*i]);*/
		
		if ( U[k+(*xrows)*i] < 0.0 ){
		    U[k+(*xrows)*i]=0.0;
		    t[k+(*xrows)*i]=-1;
		    /*Rprintf("NEGATIV\n");*/
		    control=FALSE;}
	    }
	}
    }

/*    for(i=0;i<*ncenters;i++){
      for(col=0;col<*xcols;col++)
      Rprintf("ce: %5.2f\n",centers[i+(*ncenters)*col]);}*/
    

    /*ERROR MINIMIZATION*/

    *ermin=0.0;
    for (i=0;i<*ncenters;i++){
	sum1=0.0;
	for(k=0;k<*xrows;k++)
	{
	    /*temp=pow(U[k+(*xrows)*i],f);
	      sum2=sum2+temp;*/
	    sum1+=U[k+(*xrows)*i];
	    conv += fabs(U[k+(*xrows)*i]-UANT[k+(*xrows)*i]);}
	if (*beta != 0.0){
	    hta=pow(sum1,(*xcols)*(*beta)-(*taf));
	    thsigma=pow(theta[i]*determinant[i],*beta);
	    *ermin+=hta*thsigma;}
	else
	    *ermin+=sum1*log(determinant[i]);
	
    }

    
/*for (m=0;m<*ncenters;m++){
  for (k=0;k<*xrows;k++){
  serror = 0.0;
  for(n=0;n<*xcols;n++){
    if(*dist == 0){
    serror += (x[k+(*xrows)*n] - centers[m
    +(*ncenters)*n])*(x[k+(*xrows)*n] - centers[m +(*ncenters)*n]);                                        
    }
    else if(*dist ==1){
    serror += fabs(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]);
    }
    
    }
    *ermin+=pow(U[k+(*xrows)*m],f)*serror;
    }
    }
    *ermin=*ermin/(*xrows));*/
    if (conv<= ((*xrows)*(*xcols)*epsi1)){
	if (*verbose){
	Rprintf("Iteration: %3d    converged, Error:   %13.10f\n",*iter,*ermin/(*xrows));}
	*iter=*itermax;}
    else
	if (*verbose){
	    Rprintf("Iteration: %3d    Error:   %13.10f\n",*iter,*ermin/(*xrows));
	}
   
  return 0;
}
    
int common(int *xrows, int *xcols, double *x, int *ncenters,
	   double *centers, int *itermax, int *iter, 
	   int *verbose,  double *U, double *beta, double *taf,
	   double *theta, double *ermin)
    
{
    int k;
    
    
    int i,j,col;
    double suma,tempu,exponente,tempu1,tempu2,f;
    double *UANT;
    UANT = (double *) R_alloc((*xrows)*(*ncenters), sizeof(double));

    *iter=0;
    f=2.0;


    /*Initialize Membership Matrix */
    exponente=2.0/(f-1.0);
    
    for(i=0;i<*ncenters;i++)
    {
	
	for(k=0;k<*xrows;k++)
	{
	    suma=0;
	    for(j=0;j<*ncenters;j++)
	    {
		tempu=0;
		tempu1=0;
		tempu2=0;
		for (col=0;col<*xcols;col++)
		{
		    tempu1+=(x[k+(*xrows)*col]-centers[i+(*ncenters)*col])*(x[k+(*xrows)*col]-centers[i+(*ncenters)*col]);
		    tempu2+=(x[k+(*xrows)*col]-centers[j+(*ncenters)*col])*(x[k+(*xrows)*col]-centers[j+(*ncenters)*col]);
		}

		tempu=sqrt(tempu1)/sqrt(tempu2);
		suma=suma+pow(tempu,exponente);
	    }
	    UANT[k+(*xrows)*i]=1.0/suma;
	}
    }
    
    for(i=0;i<*ncenters;i++)
    {
	
	for(j=0;j<*xrows;j++)
	    
	    U[j+(*xrows)*i]=UANT[j+(*xrows)*i];
	
    }
    

    while(((*iter)++ < *itermax)) {

	*ermin=0.0;
	
	subcommon(xrows, xcols, x, ncenters, centers, itermax,
		  iter, verbose, U, UANT, beta, taf, theta, ermin);
    }
    
    return 0;
}

 















