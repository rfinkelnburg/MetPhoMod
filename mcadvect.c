/*
   IMPLEMENTATION MODULE mcadvect.h
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mchemparse.h"
#include "mpdata.h"

void Advect(int tinc, int chemtinc)
{
  double *x, *top;
  double *dcm;
  int i, j, k, fa, loc;
  double wm[103], wx[102], zfact, tmp;
  Entity et;
/* Horizontal Advection */
  if (!(top = NEW(layer, double)) ||
      !(dcm = (double *)calloc(maxentity*layer, sizeof(double))))  {
    fprintf(stderr, "Fatal Error: Unable to allocate arrays in Advect\n");
    exit (1);
  }
/*  for (et = SUBS; et < maxentity; et++)
    if (!subst[et-SUBS].do_transport)
      printf("\nno transport for \"%s\"", subst[et-SUBS].name); */
  for (k = 0; k < nz; k++)  {
    zfact = leveli[k];
    if (k < nzm)  {
      if (windadvection)
	for (et = TEMP - (pressuretype != NONHYDROSTATIC); et--; )
	  MpAdvect(tinc, flux[UWIND]+k*layer, flux[VWIND]+k*layer, g[et]+k*layer, g[et]+(k+1)*layer, flux[WWIND]+k*layer, dcm+et*layer,
	     zfact, density[k], density[k+1], nx, ny, smiord, smnonos);
      if (advection)
	for (et = TEMP; et < (chemtinc ? maxentity : SUBS); et++)
	  if (g[et])
	    if (et == TKE || et == EPS)  {
	      if (k+1 < nzm)
	        MpAdvect(tinc, flux[UWIND]+k*layer, flux[VWIND]+k*layer,
	           g[et]+(k+1)*layer, g[et]+(k+2)*layer, flux[WWIND]+k*layer, dcm+et*layer,
		   zfact, density[k], density[k+1], nx, ny, smiord, smnonos);
	    }
	    else
	      if (et < SUBS || subst[et-SUBS+1].do_transport)
		MpAdvect((et < SUBS ? tinc : chemtinc),
		         flux[UWIND]+k*layer, flux[VWIND]+k*layer,
		         g[et]+k*layer, g[et]+(k+1)*layer, flux[WWIND]+k*layer, dcm+et*layer,
			 zfact, density[k], density[k+1], nx, ny, smiord, smnonos);
    }
    else  {
      if (windadvection)
	for (et = TEMP - (pressuretype != NONHYDROSTATIC); et--; )  {
	  if (coriolistype == GEOSTROPHICCORIOLIS)  {
	    switch (et)  {			/* Use geostrophic wind as boundary condition */
	      case UWIND : tmp = ugeos[nzm]; break;
	      case VWIND : tmp = vgeos[nzm]; break;
	      case WWIND : tmp = 0.; break;
	    }
	    for (i = layer, x = top; i--; x++)  *x = tmp;
	  }
	  else
	    memcpy(top, g[et]+k*layer, layer*sizeof(double));   /* Die aeltere Formulierung */
	  MpAdvect(tinc, flux[UWIND]+k*layer, flux[VWIND]+k*layer, g[et]+k*layer, top, flux[WWIND]+k*layer, dcm+et*layer,
	     zfact, density[k], density[k], nx, ny, smiord, smnonos);
	}
      if (advection)  {
        for (i = layer; i--; )
          top[i] = g[TEMP][nzm*layer+i] + toptemp[i];
	MpAdvect(tinc, flux[UWIND]+k*layer, flux[VWIND]+k*layer, g[TEMP]+k*layer, top, flux[WWIND]+k*layer, dcm+TEMP*layer,
	   zfact, density[k], density[k], nx, ny, smiord, smnonos);
	for (et = HUMIDITY; et < (chemtinc ? maxentity : SUBS); et++)
	  if (g[et] && et != TKE && et != EPS &&
	      (et < SUBS || subst[et-SUBS+1].do_transport))  {
	    memcpy(top, g[et]+k*layer, layer * sizeof(double));
	    MpAdvect((et < SUBS ? tinc : chemtinc),
	       flux[UWIND]+k*layer, flux[VWIND]+k*layer,
	       g[et]+k*layer, top, flux[WWIND]+k*layer, dcm+et*layer,
	       zfact, density[k], density[k], nx, ny, smiord, smnonos);
	  }
      }
    }
  }
  free(dcm); free(top);
}
