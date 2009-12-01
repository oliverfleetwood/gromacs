/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "recomb.h"
#include "futil.h"
#include "wgms.h"
#include "smalloc.h"

real *read_gammaf(char *fn,int nframes)
{
  FILE   *in;
  real   *gf;
  double y;
  int    i;
  
  snew(gf,nframes);
  in=ffopen(fn,"r");
  for(i=0; (i<nframes); i++) {
    fscanf(in,"%lf",&y);
    gf[i]=y;
  }
  ffclose(in);
  fprintf(stderr,"Succesfully read gamma\n");
  return gf;
}

void recombine(char *base,char *gammaf,int nskip,
	       int nframes,int nev,int natoms,
	       rvec *ev[],real *evprj[],
	       rvec yav[],atom_id all_index[])
{
  static char *format=
    "Recombined projection of Gamma trj (EV %d) in Cartesian Space\n";
  FILE *out;
  rvec *xxx,*evptr;
  real *gamma;
  real prj;
  char buf[256];
  int  i,j,n;
  real gt;
  
  gamma=read_gammaf(gammaf,nframes);
  snew(xxx,natoms);
  for(n=0; (n<nev); n++) {
    sprintf(buf,"%s%d",base,n+1);
    out=ffopen(buf,"w");
    fprintf(out,format,n+1);
    fprintf(stderr,format,n+1);
    evptr=ev[n];
    
    for(j=0; (j<nframes); j++) {
      if ((j % 50) == 0)
	fprintf(stderr,"\r frame %d",j);
      if ((nskip == 0) || ((j % nskip) == 0)) {
	gt=1.0/gamma[j];
	prj=evprj[n][j];
	for(i=0; (i<natoms); i++) {
	  xxx[i][XX]=(yav[i][XX]+prj*evptr[i][XX])*gt;
	  xxx[i][YY]=(yav[i][YY]+prj*evptr[i][YY])*gt;
	  xxx[i][ZZ]=(yav[i][ZZ]+prj*evptr[i][ZZ])*gt;
	}
	write_gms_ndx(out,natoms,all_index,xxx,NULL);
      }
    }
    ffclose(out);
    fprintf(stderr,"\r");
  }
  fprintf(stderr,"\n");
  sfree(xxx);
  sfree(gamma);
}
