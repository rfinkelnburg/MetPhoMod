/* IMPLEMENTATION MODULE MCLoop
   In diesem Modul sind die eigentlichen Iterationsschritte zusammengefasst. */

#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcpress.h"
#include "mcdynamics.h"
#include "mcfilter.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

double DampFact[] = {0., 0.1, 0.4, 0.9};

void SetAvgMountains(void)
{
  int i, k;
  Entity et;
  for (i = layer-row; i > row; i--)  {
    for (k = ground[i].firstabove; k--; )  {
      for (et = TEMP; et < maxentity; et++)
        if (g[et])  g[et][k*layer+i] = avg[et*nz+k];
      if (g[UWIND][k*layer+i] || g[VWIND][k*layer+i])  {
        printf("Wind is unequal 0 in mountain at point %i/%i/%i\n",
           i / row, i % row, k);
        plausible = FALSE;
      }
    }
  }
}

void InterpolateToFaces(void)
/* Interpolates the Wind-Speeds at the cell-boundaries.
   If a neighbour of a cell is in the mountain, the flux is zero. */
{
  int i, j, k, loc;
  /* Calculate normal points */
  for (i = nx; i--; )
    for (j = ny; j--; )
      for (k = ground[i*row + j].firstabove; k < nzm; k++)  {
        loc = k*layer + i*row + j;
        flux[UWIND][loc] = (pstat[loc] || pstat[loc+row] ?
        		    0. : 0.5 * (g[UWIND][loc] + g[UWIND][loc+row]));
        flux[VWIND][loc] = (pstat[loc] || pstat[loc+1] ?
        		    0. : 0.5 * (g[VWIND][loc] + g[VWIND][loc+1]));
        flux[WWIND][loc] = (pstat[loc] ? 0. : 0.5 * (g[WWIND][loc] + g[WWIND][loc+layer]));
      }
  /* Calculate points at the top */
  for (i = nx; i--; )
    for (j = ny; j--; )  {
      loc = i*row + j + nzm*layer;
      flux[UWIND][loc] = (pstat[loc] || pstat[loc+row] ?
			  0. : 0.5 * (g[UWIND][loc] + g[UWIND][loc+row]));
      flux[VWIND][loc] = (pstat[loc] || pstat[loc+1] ?
			  0. : 0.5 * (g[VWIND][loc] + g[VWIND][loc+1]));
      flux[WWIND][loc] = g[WWIND][loc];
    }
  /* Calculate points at left border if necessary. */
  if (westbordertype & DATABORDER)
    for (j = ny; j--; )
      for (k = ground[j].firstabove; k < nzm; k++)  {
        loc = k*layer + j;
        flux[UWIND][loc] = g[UWIND][loc];
      }
  /* Calculate points at right border, if necessary. */
  if (eastbordertype & DATABORDER)
    for (j = ny; j--; )
      for (k = ground[j + nx*row].firstabove; k < nzm; k++)  {
        loc = k*layer + j + nx*row;
        flux[UWIND][loc-row] = g[UWIND][loc];
      }
  /* Calculate points at "south"-border, if necessary. */
  if (southbordertype & DATABORDER)
    for (i = nx; i--; )
      for (k = ground[i*row].firstabove; k < nzm; k++)  {
        loc = k*layer + i*row;
        flux[VWIND][loc] = g[VWIND][loc];
      }
  /* Calculate points at "north"-border, if necessary. */
  if (northbordertype & DATABORDER)
    for (i = nx; i--; )
      for (k = ground[i*row+ny].firstabove; k < nzm; k++)  {
        loc = k*layer + i*row + ny;
        flux[VWIND][loc-1] = g[VWIND][loc];
      }
}

void InterpolateToCenter(double fact, int owfact)
{
  int i, j, k, loc;
  fact *= 0.5;
  for (i = nx; --i; )
    for (j = ny; --j; )
      for (k = ground[i*row+j].firstabove; k < nz; k++)  {
        loc = i*row+j+k*layer;
        g[UWIND][loc] += fact * (flux[UWIND][loc-row] + flux[UWIND][loc]);
        g[VWIND][loc] += fact * (flux[VWIND][loc-1]   + flux[VWIND][loc]);
        g[WWIND][loc] = owfact * g[WWIND][loc] +
           fact * ((k ? flux[WWIND][loc-layer] : 0.) + flux[WWIND][loc]);
      }
}

void ApplyCoriolis(int tinc)
/* This procedure will apply the Coriolis-Forces to the Wind-Speeds in the center of the
   Grid-Cells. */
{
  int i, j, k, loc, et;
  double d[3], sintma, costma;
  if (coriolistype == NOCORIOLIS)  return;
  costma = cos(turnmapangle); sintma = sin(turnmapangle);
  for (i = nx; --i; )
    for (j = ny; --j; )
      for (k = ground[i*row+j].firstabove; k < nz; k++)  {
        loc = i*row + j + k*layer;
	switch (coriolistype)  {
	case FULLCORIOLIS : 
	  d[UWIND] =  Coriol3 * g[VWIND][loc] - costma * Coriol2 * g[WWIND][loc];
	  d[VWIND] = -Coriol3 * g[UWIND][loc] + sintma * Coriol2 * g[WWIND][loc];
	  /*             d[WWIND] = (pressuretype == NONHYDROSTATIC ? Coriol3 * g[UWIND][loc] : 0.); */
	  break;
	case DIFFERENTIALCORIOLIS :
	  d[UWIND] =  Coriol3 *
	    (g[VWIND][loc] -
	     0.5 * (2 - pstat[loc-1] - pstat[loc+1]) * avg[VWIND*nz+k]) -
	    costma * Coriol2 * g[WWIND][loc];
	  d[VWIND] = -Coriol3 *
	    (g[UWIND][loc] -
	     0.5 * (2 - pstat[loc-row] - pstat[loc+row]) * avg[UWIND*nz+k]) +
	    sintma * Coriol2 * g[WWIND][loc];
	  /*             d[WWIND] = (pressuretype == NONHYDROSTATIC ? Coriol3 * g[UWIND][loc] : 0.); */
	  break;
	case GEOSTROPHICCORIOLIS :
	  d[UWIND] =  Coriol3 *
	    (g[VWIND][loc] - 0.5 * (2 - pstat[loc-1] - pstat[loc+1]) * vgeos[k]) -
	    costma * Coriol2 * g[WWIND][loc];
	  d[VWIND] = -Coriol3 *
	    (g[UWIND][loc] - 0.5 * (2 - pstat[loc-row] - pstat[loc+row]) * ugeos[k]) +
	    sintma * Coriol2 * g[WWIND][loc];
	  /*             d[WWIND] = (pressuretype == NONHYDROSTATIC ? Coriol3 * g[UWIND][loc] : 0.); */
	  break;
        }
        for (et = WWIND; et--; )  g[et][loc] += tinc * d[et];
      }
}

void ApplyBuoyancy(int tinc)
/* This Procedure will add the Bouyancy-Term to the vertical-Wind-Speed. This is
   only necessary in non-hydrostatic mode.
   The equation used is:
      dw     |theta'   Cv p'|
      -- = g |------ - -----|
      dt     |theta0   Cp p0|
   Unfortunately, the buoyancy has to be interpolated to the faces.
*/
{
  int i, j, k, hloc, loc;
  double lastb, ab;
  for (i = nx; i--; )
    for (j = ny; j--; )  {
      k = ground[hloc = i*row + j].firstabove;
      loc = hloc + k*layer;
      lastb = Grav * ((g[TEMP][loc] - avg[TEMP*nz+k])/avg[TEMP*nz+k] - CvbCp * press[loc] / pp[k]);
      for (k++; k < nz; k++)  {
        loc += layer;
        ab = Grav * ((g[TEMP][loc] - avg[TEMP*nz+k])/avg[TEMP*nz+k] - CvbCp * press[loc] / pp[k]);
        flux[WWIND][loc-layer] += 0.5 * tinc * (lastb + ab);
        lastb = ab;
      }
      flux[WWIND][loc] += tinc * ab;
    }
/*  for (i = nx; i--; )
    for (j = ny; j--; )  {
      for (k = ground[loc = i*row + j].firstabove; k < nz; k++, loc += layer)  {
        ab = Grav * ((g[TEMP][loc] - avg[TEMP*nz+k])/avg[TEMP*nz+k] - CvbCp * press[loc] / pp[k]);
        flux[WWIND][loc] += tinc * ab;
      }
    } */
}

void ApplyPressure(int tinc)
/* Will apply the wind acceleration, caused by pressure gradient, to the winds at
   the faces. */
{
  int i, j, k, hloc, loc, xstart, ystart, xend, yend;
  xstart = !!(westbordertype & DATABORDER);
  xend = nx - !!(eastbordertype & DATABORDER);
  ystart = !!(southbordertype & DATABORDER);
  yend = ny - !!(northbordertype & DATABORDER);
  for (i = xend; --i >= xstart; )
    for (j = yend; --j >= ystart; )
      for (k = ground[i*row+j].firstabove; k < nz; k++)  {
        loc = k*layer+i*row+j;
	if (!pstat[loc+row])
	  flux[UWIND][loc] -= tinc * dxi * (press[loc+row] - press[loc]) / density[k];
	if (!pstat[loc+1])
	  flux[VWIND][loc] -= tinc * dyi * (press[loc+1] - press[loc]) / density[k];
      }
}

void Continuity(void)
{
  int i, j, k, loc;
  double wflux;
  memset(flux[WWIND], 0, mesh * sizeof(double));
  for (i = nx; --i; )
    for (j = ny; --j; )  {
      wflux = 0;;
      for (k = ground[i*row+j].firstabove, loc = k*layer+i*row+j; k < nz; k++, loc += layer)  {
        wflux -= ((flux[UWIND][loc] - flux[UWIND][loc-row]) * dxi +
        	  (flux[VWIND][loc] - flux[VWIND][loc-1]) * dyi) *
        	 density[k] * level[k];
        flux[WWIND][loc] = wflux / (wflux >= 0. ? density[k] : density[k+1]);
      }
    }
}

#define PP(y)   (y > 0. ? y : 0.)
#define PM(y)   (y < 0. ? y : 0.)

void CheckTimeStep(long *tinc, long *chemtinc, double cfact)
/*  This routine controls, if the actual time-step fullfills the Courant-Condition:
   
    V*dt <= dx (resp. dy, dz)
*/
{
  int i, j, k, loc, maxt, tinctemp;
  const double eps = 1.e-3;
#ifdef DEBUG
  double tmp;
#endif
#ifdef PARALLEL
  if (parallel && master)  {
    ChangeTai(&tai);
    if (tai >= nta)  {
      fprintf(stderr, "Courant-condition is not fullfilled with shortest time-step (%li)\n", tincarray[nta-1]);
      plausible = FALSE;
      tai = nta-1;
    }
  }
  else  {
#endif
    maxt = 2 * *tincarray;
    for (i = nx; --i; )
      for (j = ny; --j; )
	for (k = ground[i*row+j].firstabove, loc = k*layer+i*row+j; k < nz; k++, loc += layer)  {
/* Es scheint, dies waere nur fuer Operator-Splitting korrekt! */
	  if (advectiontype == PPM_A)  {
            if ((tinctemp = int(cfact * dx / (PP(flux[UWIND][loc]) - PM(flux[UWIND][loc-row]) + eps))) < maxt)  maxt = tinctemp;
            if ((tinctemp = int(cfact * dy / (PP(flux[VWIND][loc]) - PM(flux[VWIND][loc-1]) + eps))) < maxt)  maxt = tinctemp;
            if (k < nzm)  {
              if ((tinctemp = int(cfact * level[k] / (PP(flux[WWIND][loc]) - (k ? PM(flux[WWIND][loc-layer]): 0) + eps))) < maxt)  maxt = tinctemp;
            }
          }
          else  {
/* ansonsten gilt: */
            if ((tinctemp = int(cfact /
        	 ( (PP(flux[UWIND][loc]) - PM(flux[UWIND][loc-row])) / dx + 
        	   (PP(flux[VWIND][loc]) - PM(flux[VWIND][loc-1])) / dy + 
        	   (PP(flux[WWIND][loc]) - (k ? PM(flux[WWIND][loc-layer]): 0)) / level[k] + 
        	  eps))) < maxt)
              maxt = tinctemp;
            }
	}
  /*  printf("maxt = %d  ", maxt); */
    for (tai = 0 ; tai < nta && tincarray[tai] > maxt; tai++);
    if (tai >= nta)  {
      fprintf(stderr, "Courant-condition is not fullfilled with shortest time-step (%li)\n", tincarray[nta-1]);
      plausible = FALSE;
      tai = nta-1;
    }
#ifdef PARALLEL
    if (parallel)
      ChangeTai(&tai);
  }
#endif
  *tinc = *chemtinc = tincarray[tai];
  if (*tinc > tincmax)  *tinc = tincmax;
  *tinc = *tinc - actime % *tinc;
  if (nsubs)  {
    if (*chemtinc > tchem)  *chemtinc = tchem;
    *chemtinc = *chemtinc - chemtime % *chemtinc;
  }
}

#define CopyField(net,i1,j1,i2,j2)  for (et = net; et--; )  {if (g[et])  g[et][k*layer+(i2)*row+(j2)] = g[et][k*layer+(i1)*row+(j1)];}
#define SetNorth(net)  for (et = net; et--; )  {if (g[et])  g[et][k*layer+i*row+ny] = (k < northborderval[et] ? avg[et*nz+k] : northborder[(et*nz+k)*xrow+i]);}
#define SetSouth(net)  for (et = net; et--; )  {if (g[et])  g[et][k*layer+i*row] = (k < southborderval[et] ? avg[et*nz+k] : southborder[(et*nz+k)*xrow+i]);}
#define SetWest(net)   for (et = net; et--; )  {if (g[et])  g[et][k*layer+j] = (k < westborderval[et] ? avg[et*nz+k] : westborder[(et*nz+k)*row+j]);}
#define SetEast(net)   for (et = net; et--; )  {if (g[et])  g[et][k*layer+nx*row+j] = (k < eastborderval[et] ? avg[et*nz+k] : eastborder[(et*nz+k)*row+j]);}

void DampValuesX(int from, int to, int direction)
{
  int i, j, k, lfrom, lto;
  Entity et;
  double *wdv, *sdv;
  static double scalardampval[] = {0.06, 0.04, 0.02, 0.005};
  static double winddampval[] = {0.2, 0.04, 0.008, 0.0016};
  wdv = winddampval; sdv = scalardampval;
  for (j = 4; j--; to += direction, wdv++, sdv++)
    for (i = nx; --i; )
      for (k = ground[i*row+from].firstabove; k < nz; k++)  {
        lfrom = k*layer+i*row+from;
        lto   = k*layer+i*row+to;
        if (!pstat[lto])  {
          for (et = SUBS; et-- > TEMP; )
            if (g[et])  g[et][lto] += *sdv * (g[et][lfrom] - g[et][lto]);
          for (et = TEMP; et--; )
            g[et][lto] += *wdv * (g[et][lfrom] - g[et][lto]);
        }
      }
}

void DampValuesY(int from, int to, int direction)
{
  int i, j, k, lfrom, lto;
  Entity et;
  static double scalardampval[] = {0.06, 0.04, 0.02, 0.005};
  static double winddampval[] = {0.2, 0.04, 0.008, 0.0016};
  double *wdv, *sdv;
  wdv = winddampval; sdv = scalardampval;
  for (i = 4; i--; to += direction, wdv++, sdv++)
    for (j = ny; --j; )
      for (k = ground[from*row+j].firstabove; k < nz; k++)  {
        lfrom = k*layer+from*row+j;
        lto   = k*layer+to*row+j;
        if (!pstat[lto])  {
          for (et = SUBS; et-- > TEMP; )
            if (g[et])  g[et][lto] += *sdv * (g[et][lfrom] - g[et][lto]);
          for (et = TEMP; et--; )
            g[et][lto] += *wdv * (g[et][lfrom] - g[et][lto]);
        }
      }
}

void DampTop(void)
{
/*  const double dfact[4] = {0.000244140625, 0.001953125, 0.015625, 0.125}, *dfp;  */
  const double dfact[4] = {0.001953125, 0.0078125, 0.03125, 0.125}, *dfp;
  int i, j, k, et, loc;
  double *p, tavg;
  for (k = nz-4, dfp = dfact; k < nz; k++, dfp++)  {
    p = g[TEMP] + k*layer;
    tavg = avg[TEMP*nz + k];
    for (i = nx; --i; )
      for (j = ny; --j; )  {
        loc = i*row + j;
        p[loc] += 4. * *dfp * (tavg - p[loc]);
      }
  }
  for (et = HUMIDITY; et--; )
    for (k = nz-4, dfp = dfact; k < nz; k++, dfp++)  {
      p = g[et] + k*layer;
      for (i = nx; --i; )
        for (j = ny; --j; )  {
          loc = i*row + j;
          tmplayer[loc] = p[loc-row] + p[loc+row] + p[loc-1] + p[loc+1] -
          		  4. * p[loc];
        }
      for (i = nx; --i; )
        for (j = ny; --j; )  {
          loc = i*row + j;
          p[loc] += *dfp * tmplayer[loc];
        }
    }
}


void SetBoundary(int net)
{
  int i, j, k, im, ip, jm, jp;
  Entity et;
  BorderType bt;
  switch (northbordertype)  {
    case AUTOCONSTANTBORDER :
      for (i = nx; --i; )
        for (k = ground[i*row+ny].firstabove; k < nz; k++)
          if (g[VWIND][k*layer+i*row+nym] > 0)  {CopyField(net, i, nym, i, ny) SetNorth(3)}
          else	SetNorth(net)
      break;
    case CONSTANTBORDER : case SPONGEBORDER :
      for (i = nx; --i; )
        for (k = ground[i*row+ny].firstabove; k < nz; k++)
          SetNorth(net)
      break;
    case FREEBORDER :
      for (i = nx; --i; )
        for (k = ground[i*row+ny].firstabove; k < nz; k++)
          CopyField(net, i, nym, i, ny)
      break;
    case CYCLIC :
      for (i = nx; --i; )
        for (k = ground[i*row+ny].firstabove; k < nz; k++)
          CopyField(net, i, 1, i, ny)
      break;
    case MIRRORED :
      for (i = nx; --i; )
        for (k = ground[i*row+ny].firstabove; k < nz; k++)  {
          CopyField(net, i, nym-1, i, ny)
          g[VWIND][k*layer+i*row+ny] = -g[VWIND][k*layer+i*row+ny];
        }
      break;
  }
  switch (southbordertype)  {
    case AUTOCONSTANTBORDER :
      for (i = nx; --i; )
        for (k = ground[i*row].firstabove; k < nz; k++)
          if (g[VWIND][k*layer+i*row+1] < 0)	{CopyField(net, i, 1, i, 0) SetSouth(3)}
          else				SetSouth(net)
      break;
    case CONSTANTBORDER : case SPONGEBORDER :
      for (i = nx; --i; )
        for (k = ground[i*row].firstabove; k < nz; k++)
          SetSouth(net)
      break;
    case FREEBORDER :
      for (i = nx; --i; )
        for (k = ground[i*row].firstabove; k < nz; k++)
          CopyField(net, i, 1, i, 0)
      break;
    case CYCLIC :
      for (i = nx; --i; )
        for (k = ground[i*row].firstabove; k < nz; k++)
          CopyField(net, i, nym, i, 0)
      break;
    case MIRRORED :
      for (i = nx; --i; )
        for (k = ground[i*row].firstabove; k < nz; k++)  {
          CopyField(net, i, 2, i, 0)
          g[VWIND][k*layer+i*row] = -g[VWIND][k*layer+i*row];
        }
      break;
  }
  switch (westbordertype)  {
    case AUTOCONSTANTBORDER :
      for (j = ny; --j; )
        for (k = ground[j].firstabove; k < nz; k++)
          if (g[UWIND][k*layer+row+j] < 0)  {CopyField(net, 1, j, 0, j)  SetWest(3)}
          else				SetWest(net)
      break;
    case CONSTANTBORDER : case SPONGEBORDER :
      for (j = ny; --j; )
        for (k = ground[j].firstabove; k < nz; k++)
          SetWest(net)
      break;
    case FREEBORDER :
      for (j = ny; --j; )
        for (k = ground[j].firstabove; k < nz; k++)
          CopyField(net, 1, j, 0, j)
      break;
    case CYCLIC :
      for (j = ny; --j; )
        for (k = ground[j].firstabove; k < nz; k++)
          CopyField(net, nxm, j, 0, j)
      break;
    case MIRRORED :
      for (j = ny; --j; )
        for (k = ground[j].firstabove; k < nz; k++)  {
          CopyField(net, 2, j, 0, j)
          g[UWIND][k*layer+j] = -g[UWIND][k*layer+j];
        }
      break;
  }
  switch (eastbordertype)  {
    case AUTOCONSTANTBORDER :
      for (j = ny+1; j--; )
        for (k = ground[nx*row+j].firstabove; k < nz; k++)
          if (g[UWIND][k*layer+nxm*row+j] > 0)	{CopyField(net, nxm, j, nx, j) SetEast(3)}
          else				SetEast(net)
      break;
    case CONSTANTBORDER : case SPONGEBORDER :
      for (j = ny+1; j--; )
        for (k = ground[nx*row+j].firstabove; k < nz; k++)
          SetEast(net)
      break;
    case FREEBORDER :
      for (j = ny+1; j--; )
        for (k = ground[nx*row+j].firstabove; k < nz; k++)
          CopyField(net, nxm, j, nx, j)
      break;
    case CYCLIC :
      for (j = ny+1; j--; )
        for (k = ground[nx*row+j].firstabove; k < nz; k++)
          CopyField(net, 1, j, nx, j)
      break;
    case MIRRORED :
      for (j = ny+1; j--; )
        for (k = ground[nx*row+j].firstabove; k < nz; k++)  {
          CopyField(net, nxm-1, j, nx, j)
          g[UWIND][k*layer+nx*row+j] = -g[UWIND][k*layer+nx*row+j];
        }
      break;
  }
  if (net > WWIND)  {
    if (northbordertype == SPONGEBORDER)
      DampValuesX(ny, nym, -1);
    if (southbordertype == SPONGEBORDER)
      DampValuesX(0, 1, 1);
    if (westbordertype == SPONGEBORDER)
      DampValuesY(0, 1, 1);
    if (eastbordertype == SPONGEBORDER)
      DampValuesY(nx, nxm, -1);
  }
#ifdef PARALLEL
  if (westbordertype == NEIGHBOUR)
    SendWallToWorker(net, 1 + OVERLAP, leftbrother, TOTHELEFT);
  if (eastbordertype == NEIGHBOUR)
    GetWallFromWorker(net, nx - OVERLAP, rightbrother, TOTHELEFT);
  if (eastbordertype == NEIGHBOUR)
    SendWallToWorker(net, nxm - 2*OVERLAP, rightbrother, TOTHERIGHT);
  if (westbordertype == NEIGHBOUR)
    GetWallFromWorker(net, 0, leftbrother, TOTHERIGHT);
#endif
  for (et = maxentity; et--; )
    if (g[et])
      for (k = (nz-1)*layer; k >= 0; k -= layer)  {
        if (westbordertype != NEIGHBOUR)  {
          g[et][k] = g[et][k+row+1];
          g[et][k+ny] = g[et][k+row+nym];
        }
        if (eastbordertype != NEIGHBOUR)  {
          g[et][k+nx*row] = g[et][k+nxm*row+1];
          g[et][k+nx*row+ny] = g[et][k+nxm*row+nym];
        }
      }
}

