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
#include "R.h"

int  subcshell(int *xrows, int *xcols, double *x, int *ncenters,
	       double *centers, int *itermax, int *iter,
	       int *verbose, int *dist, double *U, double *UANT, double
	       *f, double *ermin, double *radius, int *flag)
{
  int k, col, i, m, n ;

  double serror;

  /*convergence parameters*/
  double epsi1, epsi2, conv;
  
  double sum2;
  double temp,tempu, tempu1, tempu2, distance;
  int j;
  double suma;
  double exponente;
  
/*  *ermin=0.0;*/
  serror=0.0;
  
  
  sum2=0;
    
  if ((*flag==0) || (*flag==5)){
      
      /* UPDATE CENTERS*/

      for(i=0;i<*ncenters;i++)
      {
	  sum2=0;
	  for(col=0;col<*xcols;col++)
	      centers[i+(*ncenters)*col]=0.0;
	  for(k=0;k<*xrows;k++)
	  {
	      temp=pow(U[k+(*xrows)*i],*f);
	      sum2=sum2+temp;
	      
	      for(col=0;col<*xcols;col++)
	      {
		  centers[i+(*ncenters)*col]+= temp*x[k+(*xrows)*col];
	      }
	      
	      
	  }
	  for(col=0;col<*xcols;col++)
	      centers[i+(*ncenters)*col]/=sum2;
      }


      /*UPDATE radius*/
      for(i=0;i<*ncenters;i++)
      {
	  sum2=0;
	  
	  radius[i]=0.0;
	  for(k=0;k<*xrows;k++)
	  {
	      distance=0.0;
	      temp=pow(U[k+(*xrows)*i],*f);
	      sum2=sum2+temp;
	      
	      for(col=0;col<*xcols;col++)
	      {
		  if (*dist==0){
		      distance+= (x[k+(*xrows)*col]-centers[i+(*ncenters)*col])*(x[k+(*xrows)*col]-centers[i+(*ncenters)*col]);
		  }
		  else if(*dist ==1){
		      distance+=fabs(x[k+(*xrows)*col]-centers[i+(*ncenters)*col]);
		  }
	      }
	      if (*dist==0){
		  radius[i]+= temp*sqrt(distance);}
	      else if(*dist ==1){
		  radius[i]+= temp*distance;}
	  }
	  radius[i]/=sum2;
      }
      
  }/*flag=0*/
  
/*update UANT*/
  for(i=0;i<*ncenters;i++){
      for(k=0;k<*xrows;k++){
	  UANT[k+(*xrows)*i]=U[k+(*xrows)*i];}}
	      
  
  /* UPDATE Membership Matrix */
  
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
		  tempu=fabs(sqrt(tempu1)-radius[i])/fabs(sqrt(tempu2)-radius[j]);
	      }
	      else if(*dist ==1){
		  tempu=fabs((tempu1-radius[i])/(tempu2-radius[j]));
	      }
	      suma=suma+pow(tempu,exponente);
	  }
	  U[k+(*xrows)*i]=1.0/suma;
	  
      }
      
  }


  
  /*ERROR MINIMIZATION*/
  
  epsi1=0.002;
 epsi2=0.2;
  conv=0.0;


  for (m=0;m<*ncenters;m++){
      for (k=0;k<*xrows;k++){
	  serror = 0.0;
	  for(n=0;n<*xcols;n++){
	      if(*dist == 0){ 
		  serror += (x[k+(*xrows)*n] - centers[m+(*ncenters)*n])*(x[k+(*xrows)*n] - centers[m +(*ncenters)*n]);
	      }
	      else if(*dist ==1){
		  serror += fabs(x[k+(*xrows)*n] - centers[m + (*ncenters)*n]);
	      }
	      
	  }
	  if (*dist == 0){
	      serror=fabs(sqrt(serror)-radius[m]);}
	  else if(*dist ==1){
	      serror=fabs(serror-radius[m]);}
	 
	  *ermin+=pow(U[k+(*xrows)*m],*f)*pow(serror,2);
        
/* *ermin=*ermin/(*xrows));*/
      
	
         /*Convergence check*/
	  conv += fabs(U[k+(*xrows)*m]-UANT[k+(*xrows)*m]);
      }
  }
  
  if (conv<= ((*xrows)*(*xcols)*epsi1)){
       *flag=2;
       if (*verbose){
       Rprintf("Iteration: %3d    converged, Error:   %13.10f\n",*iter,conv);
       }}
  else if (conv<= ((*xrows)*(*xcols)*epsi2)){
      if (*verbose){
	  Rprintf("Iteration: %3d    Epsi2:   %13.10f\n",*iter,conv);}
      
      if (*flag==3)
	  *flag=4;
      else
	  *flag=1;
  }
  else if(*flag==3)
      *flag=5;
  
  
  if (*verbose){
      Rprintf("Iteration: %3d    Error:   %13.10f\n",*iter,*ermin/(*xrows));
  }
  
  return 0;
}
  
int cshell(int *xrows, int *xcols, double *x, int *ncenters,
	   double *centers, int *itermax, int *iter, 
	   int *verbose, int *dist, double *U, double *UANT, double
	   *f, double *ermin, double *radius, int *flag)
    
{
    int k;
    
    
    int i,j,col;
    double suma,tempu,exponente,tempu1,tempu2;
    
    exponente=2.0/(*f-1.0);
    
/*    *flag=0;*/

    if (*flag==0){
	*iter=0;
	
	/*Initialize Membership Matrix */

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
			tempu=fabs(sqrt(tempu1)-radius[i])/fabs(sqrt(tempu2)-radius[j]);
		    }
		    else if(*dist ==1){
			tempu=fabs((tempu1-radius[i])/(tempu2-radius[j]));
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
	
    }

    while(((*iter)++ < *itermax) && ((*flag)!=1 && (*flag)!=2) && (*flag)!=4) {
	
	*ermin=0.0;
	
	subcshell(xrows, xcols, x, ncenters, centers, itermax,
		  iter, verbose, dist, U, UANT, f, ermin, radius, flag);
    }
    
    return 0;
}

/*****************************************************************/
/*******only for prediction***************************************/
/*****************************************************************/

int  cshell_assign(int *xrows, int *xcols, double *x, int *ncenters,
		   double *centers, int *dist, double *U, double *f,
		   double *radius)
{
  int k, col, i;

  double tempu, tempu1, tempu2;
  int j;
  double suma;
  double exponente;

  
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
		  tempu=fabs(sqrt(tempu1)-radius[i])/fabs(sqrt(tempu2)-radius[j]);
	      }
	      else if(*dist ==1){
		  tempu=fabs((tempu1-radius[i])/(tempu2-radius[j]));
	      }
	      suma=suma+pow(tempu,exponente);
	  }
	  U[k+(*xrows)*i]=1.0/suma;
	  
      }
      
  }
  
  return 0;
}
  






