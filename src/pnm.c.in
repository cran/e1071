/*
 *  Copyright (C) 1998 Friedrich Leisch
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

#include <pnm.h>



void readpnm(char **filename, int *red, int *green, int *blue)
{
  FILE *fp;
  xelval xmaxval;
  xel **ximage;
  int cols, rows, k, xformat;

  fp = fopen(*filename, "r");

  if(fp==NULL){
    printf("Can't open %s for reading\n", *filename);
    return;
  }

  ximage = pnm_readpnm(fp, &cols, &rows, &xmaxval, &xformat);  

  if(PNM_FORMAT_TYPE(xformat) == PBM_TYPE ||
     PNM_FORMAT_TYPE(xformat) == PGM_TYPE){
      for(k=0; k<rows*cols; k++){
	  red[k] = (int) PNM_GET1((*ximage)[k]); 
      }
  }
  else if (PNM_FORMAT_TYPE(xformat) == PPM_TYPE){
      for(k=0; k<rows*cols; k++){
	  red[k] = (int) PPM_GETR((*ximage)[k]);
	  green[k] = (int) PPM_GETG((*ximage)[k]);
	  blue[k] = (int) PPM_GETB((*ximage)[k]);
      }
  }
  else {
      printf("Unknown PNM format\n");
  }

  fclose(fp);
  return;
}

void readpnminit(char **filename, int *cols, int *rows, int *maxval,
		char **type)
{
  FILE *fp;
  xelval xmaxval;
  int xformat;

  fp = fopen(*filename, "r");

  if(fp==NULL){
    printf("Can't open %s for reading\n", *filename);
    return;
  }

  pnm_readpnminit(fp, cols, rows, &xmaxval, &xformat);  
  *maxval = (int)xmaxval;

  if(PNM_FORMAT_TYPE(xformat) == PBM_TYPE)
      strcpy(*type, "pbm");
  else if(PNM_FORMAT_TYPE(xformat) == PGM_TYPE)
      strcpy(*type, "pgm");
  else if(PNM_FORMAT_TYPE(xformat) == PPM_TYPE)
      strcpy(*type, "ppm");


  fclose(fp);
}



void writepgm(char **filename, int *image, int *cols, int *rows, int
	     *maxval, int *forceplain)
{
  FILE *fp;
  gray gmaxval;
  gray **gimage;
  int k;

  fp = fopen(*filename, "w");

  if(fp==NULL){
    printf("Can't open %s for writing\n", *filename);
    return;
  }

  gmaxval = (gray)(*maxval);
  gimage = pgm_allocarray(*cols, *rows);
  for(k=0; k<(*rows) * (*cols); k++){
      (*gimage)[k] = (gray) image[k];
  }

  pgm_writepgm(fp, gimage, *cols, *rows, gmaxval, *forceplain);  
  fclose(fp);
}





