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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int  subcmeans(int *xrows, int *xcols, double *x, int *ncenters,
	       double *centers, int *itermax, int *iter,
	       int *verbose, int *dist, double *U, double *UANT, double
	       *f, double *ermin)
{
  int k, col, i, m, n ;

  double serror;

  
  double sum1;
  double sum2;
  double temp,tempu, tempu1, tempu2;
  int j;
  double suma;
  double exponente,epsi1,conv;
  epsi1=0.002;

  conv=0.0;
  
/*  *ermin=0.0;*/
  serror=0.0;
  
  
  sum1=0;
  sum2=0;
    

  /*update UANT*/
  for(i=0;i<*ncenters;i++){
      for(k=0;k<*xrows;k++){
	  UANT[k+(*xrows)*i]=U[k+(*xrows)*i];
      }}
    
  
  /*NEW CENTERS*/

  for(i=0;i<*ncenters;i++)
  {
      sum1=0;
      sum2=0;
      for(col=0;col<*xcols;col++)
	  centers[i+(*ncenters)*col]=0.0;
      for(k=0;k<*xrows;k++)
      {
	  temp=pow(UANT[k+(*xrows)*i],*f);
	  sum2=sum2+temp;

	  for(col=0;col<*xcols;col++)
	  {
	      centers[i+(*ncenters)*col]+= temp*x[k+(*xrows)*col];
	  }
      }
      for(col=0;col<*xcols;col++)
	      centers[i+(*ncenters)*col]/=sum2;
  }
  
  
  /*Membership Matrix */
  
  exponente=2.0/(*f-1.0);
  
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
		  if (*dist==0){
		      tempu1+=(x[k+(*xrows)*col]-centers[i+(*ncenters)*col])*(x[k+(*xrows)*col]-centers[i+(*ncenters)*col]);
		      tempu2+=(x[k+(*xrows)*col]-centers[j+(*ncenters)*col])*(x[k+(*xrows)*col]-centers[j+(*ncenters)*col]);
		  }
		   else if(*dist ==1){
		       tempu1+=fabs(x[k+(*xrows)*col]-centers[i+(*ncenters)*col]);
		       tempu2+=fabs(x[k+(*xrows)*col]-centers[j+(*ncenters)*col]);
		   }					     
	      }
	      if (*dist==0){
		  tempu=sqrt(tempu1)/sqrt(tempu2);
	      }
	      else if(*dist ==1){
		  tempu=tempu1/tempu2;
	      }
	      suma=suma+pow(tempu,exponente);
	  }
	  U[k+(*xrows)*i]=1.0/suma;
	  /* UANT[k+(*xrows)*i]=U[k+(*xrows)*i];*/
      }
      
  }
   
  /*ERROR MINIMIZATION*/

  for (m=0;m<*ncenters;m++){
      for (k=0;k<*xrows;k++){
	  serror = 0.0;
	  conv += fabs(U[k+(*xrows)*m]-UANT[k+(*xrows)*m]);
	  for(n=0;n<*xcols;n++){
	      if(*dist == 0){
		  serror += (x[k+(*xrows)*n] - centers[m
						      +(*ncenters)*n])*(x[k+(*xrows)*n] - centers[m +(*ncenters)*n]);                                        
	      }
	      else if(*dist ==1){
		  serror += fabs(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]);
	      }
	      
	  }
	  *ermin+=pow(U[k+(*xrows)*m],*f)*serror;
      }
  }
  /* *ermin=*ermin/(*xrows));*/

  if (conv<= ((*xrows)*(*xcols)*epsi1)){
      printf("Iteration: %3d    converged, Error:   %13.10f\n",*iter,*ermin/(*xrows));
      *iter=*itermax;}
  else
      if (*verbose){
	  printf("Iteration: %3d    Error:   %13.10f\n",*iter,*ermin/(*xrows));
      }
  
  return 0;
}
    
int cmeans(int *xrows, int *xcols, double *x, int *ncenters,
	   double *centers, int *itermax, int *iter, 
	   int *verbose, int *dist, double *U, double *UANT, double
	   *f, double *ermin)
    
{
    int k;
    
    
    int i,j,col;
    double suma,tempu,exponente,tempu1,tempu2;
    
    *iter=0;
    
    /*Initialize Membership Matrix */
    exponente=2.0/(*f-1.0);
    
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
		    if (*dist==0){
			tempu1+=(x[k+(*xrows)*col]-centers[i+(*ncenters)*col])*(x[k+(*xrows)*col]-centers[i+(*ncenters)*col]);
			tempu2+=(x[k+(*xrows)*col]-centers[j+(*ncenters)*col])*(x[k+(*xrows)*col]-centers[j+(*ncenters)*col]);
		    }
		    else if(*dist ==1){
			tempu1+=fabs(x[k+(*xrows)*col]-centers[i+(*ncenters)*col]);
			tempu2+=fabs(x[k+(*xrows)*col]-centers[j+(*ncenters)*col]);
		   }					     
		}

		if (*dist==0){
		    tempu=sqrt(tempu1)/sqrt(tempu2);
		}
		else if(*dist ==1){
		    tempu=tempu1/tempu2;
		}
		
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
	
	subcmeans(xrows, xcols, x, ncenters, centers, itermax,
		  iter, verbose, dist, U, UANT, f, ermin);
    }
    
    return 0;
}

 
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int  subufcl(int *xrows, int *xcols, double *x, int *ncenters,
	     double *centers, int *itermax, int *iter,
	     int *verbose, int *dist, double *U, double
	     *f, double *par, double *ermin)
{
  int k, col, i, m, n ;

  double serror;

  
  double tempu, tempu1, tempu2;
  int j;
  double suma;
  double exponente, lrate;

/*  *ermin=0.0;*/
  serror=0.0;
  
  
  
  /*Membership Matrix Computation*/
  
  exponente=2.0/(*f-1.0);
  
  for(k=0;k<*xrows;k++)
  {
      for(i=0;i<*ncenters;i++)
      {
	  suma=0;
	  for(j=0;j<*ncenters;j++)
	  {
	      tempu=0;
	      tempu1=0;
	      tempu2=0;
	      for (col=0;col<*xcols;col++)
	      {
		  if (*dist==0){
		      tempu1+=(x[k+(*xrows)*col]-centers[i+(*ncenters)*col])*(x[k+(*xrows)*col]-centers[i+(*ncenters)*col]);
		      tempu2+=(x[k+(*xrows)*col]-centers[j+(*ncenters)*col])*(x[k+(*xrows)*col]-centers[j+(*ncenters)*col]);
		  }
		   else if(*dist ==1){
		       tempu1+=fabs(x[k+(*xrows)*col]-centers[i+(*ncenters)*col]);
		       tempu2+=fabs(x[k+(*xrows)*col]-centers[j+(*ncenters)*col]);
		   }					     
	      }
	      if (*dist==0){
		  tempu=sqrt(tempu1)/sqrt(tempu2);
	      }
	      else if(*dist ==1){
		  tempu=tempu1/tempu2;
	      }
	      suma=suma+pow(tempu,exponente);
	  }
	  U[k+(*xrows)*i]=1.0/suma;


	  /*UPDATE CENTERS*/
	  lrate = (*par)*(1-(double)*iter/(*itermax));
	      for (n=0;n<*xcols;n++){
		  centers[i+(*ncenters)*n]+=lrate*pow(U[k+(*xrows)*i],*f)*(x[k+(*xrows)*n]-centers[i+(*ncenters)*n]);
	      }
      }
      
  }      
  /*Stochastic Approximate ERROR MINIMIZATION*/
  
  for (m=0;m<*ncenters;m++){
      for (k=0;k<*xrows;k++){
	  serror=0.0;
	  for(n=0;n<*xcols;n++){
		  
	      if(*dist == 0){
		  serror += (x[k+(*xrows)*n] - centers[m
						      +(*ncenters)*n])*(x[k+(*xrows)*n] - centers[m +(*ncenters)*n]);                                        
	      }
	      else if(*dist ==1){
		  serror += fabs(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]);
	      }
	      
	      }
	  *ermin+=pow(U[k+(*xrows)*m],*f)*serror;
      }
  }
  /* *ermin=*ermin/(*xrows));*/
      
  if (*verbose){
      printf("Iteration: %3d    Error:   %13.10f\n",*iter,*ermin/(*xrows));
  }
  
  return 0;
}

int ufcl(int *xrows, int *xcols, double *x, int *ncenters,
	 double *centers, int *itermax, int *iter, 
	 int *verbose, int *dist, double *U, double
	 *f, double *par, double *ermin)
    
{
    *iter=0;
        
    
    while(((*iter)++ < *itermax)) {

	*ermin=0.0;
	
	subufcl(xrows, xcols, x, ncenters, centers, itermax,
		iter, verbose, dist, U, f, par, ermin);
    }
    
    return 0;
}

















