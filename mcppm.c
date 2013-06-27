/*
   MODULE MCPPM
   Implementiert den PPM-Advektionsalgorithmus
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mcglobal.h"
#include "mcppm.h"
#include "sunpos.h"
#include "mcground.h"
#include "mchemparse.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

static double *al, *ar, *af, *da, *a6, *dz, *f1, *f2, *f3, *f4, *f5, *f6, *b,
	      *top;

void InitializePPM(void)
{
  int nmax, i;
  nmax = (nx > ny ? (nx > nz ? nx : nz) : (ny > nz ? ny : nz));
  if (!(al = NEW(nmax+1, double)) ||
      !(ar = NEW(nmax+1, double)) ||
      !(af = NEW(nmax+1, double)) ||
      !(da = NEW(nmax+1, double)) ||
      !(a6 = NEW(nmax+1, double)) ||
      !(top = NEW(layer, double))  ||
      !(b = NEW(mesh, double)) ||
      !(dz = NEW(nz, double)) ||
      !(f1 = NEW(nz, double)) ||
      !(f3 = NEW(nz, double)) ||
      !(f4 = NEW(nz, double)) ||
      !(f5 = NEW(nz, double)) ||
      !(f6 = NEW(nz, double)))  {
    fprintf(stderr, "Fatal Error: Unable to allocate memory in \"InitializePPM\"\n");
    exit (5);
  }
/*  for (i = mesh; i--; )
    b[i] = 1.e10; */
}

static void CalcPPMGeoFact(void)
{
  int i;
  double tmp;
  for (i = nz; i--; )
    dz[i] = density[i] * level[i];
  for (i = nz-2; i > 0; i--)  {
    f1[i] = dz[i] / (dz[i] + dz[i+1]);
    tmp = 1. / (dz[i-1] + dz[i] + dz[i+1] + dz[i+2]);
    f3[i] = tmp * (dz[i-1] + dz[i]) / (2. * dz[i] + dz[i+1]);
    f4[i] = tmp * (dz[i+1] + dz[i+2]) / (dz[i] + 2. * dz[i+1]);
    f1[i] += 2. * dz[i+1] * dz[i] / (dz[i] + dz[i+1]) * (f3[i] - f4[i]);
    f3[i] *= dz[i]; f4[i] *= dz[i+1];
  }
  for (i = nzm; i > 0; i--)  {
    tmp = dz[i] / (dz[i-1] + dz[i] + dz[i+1]);
    f5[i] = tmp * (2. * dz[i-1] + dz[i]) / (dz[i+1] + dz[i]);
    f6[i] = tmp * (dz[i] + 2. * dz[i+1]) / (dz[i-1] + dz[i]);
  }
}

static double Smooth1(double da, double am, double a, double ap)
{
  if ((ap - a) * (a - am) > 0.)  {
    am = 2. * fabs(a - am);
    ap = 2. * fabs(a - ap);
    ap = (ap < am ? ap : am);
    a = fabs(da);
    a = copysign((a < ap ? a : ap), da);
    return (a);
  }
  else
    return (0.);
}

static void SetBoundX(int nx, double *x, double *vx)
{
  int j, nxm;
  nxm = nx-1;
  switch (westbordertype)  {
    case AUTOCONSTANTBORDER :
       for (j = ny; --j; )
         if (vx[j] < 0.)
           x[j] = x[row+j];
       break;
    case FREEBORDER :
       for (j = ny; --j; )
         x[j] = x[row+j];
       break;
    case CYCLIC :
       for (j = ny; --j; )
         x[j] = x[nxm*row+j];
       break;
    case MIRRORED :
       for (j = ny; --j; )
         x[j] = x[2*row+j];
       break;
  }
  switch (eastbordertype)  {
    case AUTOCONSTANTBORDER :
       for (j = ny; --j; )
         if (vx[nxm*row+j] > 0.)
           x[nx*row+j] = x[nxm*row+j];
       break;
    case FREEBORDER :
       for (j = ny; --j; )
         x[nx*row+j] = x[nxm*row+j];
       break;
    case CYCLIC :
       for (j = ny; --j; )
         x[nx*row+j] = x[row+j];
       break;
    case MIRRORED :
       for (j = ny; --j; )
         x[nx*row+j] = x[(nx-2)*row+j];
       break;
  }
}

static void SetBoundY(int ny, double *x, double *vy)
{
  int i, nym;
  nym = ny-1;
  switch (southbordertype)  {
    case AUTOCONSTANTBORDER :
       for (i = nx; --i; )
         if (vy[i*row] < 0.)
           x[i*row] = x[i*row+1];
    case FREEBORDER :
       for (i = nx; --i; )
         x[i*row] = x[i*row+1];
       break;
    case CYCLIC : 
       for (i = nx; --i; )
         x[i*row] = x[i*row+nym];
       break;
    case MIRRORED :
       for (i = nx; --i; )
         x[i*row] = x[i*row+2];
       break;
  }
  switch (northbordertype)  {
    case AUTOCONSTANTBORDER :
       for (i = nx; --i; )
         if (vy[i*row+nym] > 0.)
           x[i*row+ny] = x[i*row+nym];
    case FREEBORDER :
       for (i = nx; --i; )
         x[i*row+ny] = x[i*row+nym];
       break;
    case CYCLIC : 
       for (i = nx; --i; )
         x[i*row+ny] = x[i*row+1];
       break;
    case MIRRORED :
       for (i = nx; --i; )
         x[i*row+ny] = x[i*row+ny-2];
       break;
  }
}

static void HorizontalPPM(int n, int m, double *v, double *a, PointStatus *p,
		   double dx, long dt, int iord, double *af,
		   short redleft, short redright)
/*
   Variables:
     n  : The number of points. Transport is calculated from point 1 to nx-1
          (general Metphomod definition).
     m  : The array multiplier.
     v  : The velocity field. Velocities are defined at the faces: v[0]
          id between c[0] and c[1]. (corresponds to Colella and Woodward.)
     a  : The concentration field.
     p  : The point status of the points.
     dx : dx
     dt : dt
     iord : order to use.
     cdiff : the concentration difference.
*/
{
  int i, j, loc;
  const double ofac = 1./12.;
  double t1, t2, courant, deltaa, deltap, tmp;
  /* Stage 1: Interpolate values to the faces */
  i = n;
  do  {
    while (i && p[i*m])  i--;
    af[i] = a[loc = m*i];
    if (i--)  {
      loc -= m;
      if (p[loc])
        af[i] = a[loc+m];
      else  {
	af[i] = 0.5 * (a[loc] + a[loc+m]);
	if (i && !p[loc-m])  {
          deltap = Smooth1(0.5 * (a[loc+m] - a[loc-m]), a[loc-m], a[loc], a[loc+m]);
	  while (--i > 0 && !p[m*(i-1)])  {
	    loc = i*m;
	    deltaa = Smooth1(0.5 * (a[loc+m] - a[loc-m]), a[loc-m], a[loc], a[loc+m]);
	    af[i] = 0.5 * (a[loc]+a[loc+m]) + ofac * (deltaa - deltap);
	    deltap = deltaa;
	  }
	}
	if (i >= 0)  {
	  af[i] = 0.5 * (a[m*i] + a[m*(i+1)]);
	  if (i > 0)  {
	    af[i-1] = a[m*i];
	    i--;
	  }
	}
      }
    }
  }  while (i > 0);
  /* Stage 2: Determine ends of parabels */
  for (i = n; --i; )
    if (!p[m*i])  {
      al[i] = af[i-1];
      ar[i] = af[i];
      tmp = a[m*i];
      if (iord == 0 || (ar[i] - tmp) * (tmp - al[i]) <= 0.)
        ar[i] = al[i] = tmp;
      else {
        t2 = ar[i] - al[i];
        t1 = t2 * (tmp - 0.5 * (al[i] + ar[i]));
        t2 = t2*t2 / 6.;
        if (t1 > t2)  
          al[i] = 3. * tmp - 2. * ar[i];
        if (t1 < -t2)
          ar[i] = 3. * tmp - 2. * al[i];
      }
      da[i] = ar[i] - al[i];
      a6[i] = 6. * (tmp - 0.5 * (al[i] + ar[i]));
    }
  /* Pass to first order on leaving boundary */
  if (redleft && v[0] < 0.)  {
    ar[1] = al[1] = a[m];
    da[1] = a6[1] = 0.;
  }
  if (redright && v[j = (n-1)*m] > 0.)  {
    ar[n-1] = al[n-1] = a[j];
    da[n-1] = a6[n-1] = 0.;
  }
  ar[0] = al[0] = a[0];
  ar[n] = al[n] = a[n*m];
  da[0] = a6[0] = da[n] = a6[n] = 0.;
  /* Stage 3: Calculate zone averages */
  dx = (double)dt / dx;
  i = n; loc = (i-1)*m;
  if (!p[loc])
    if (v[loc] < 0.)  {
      courant = -v[loc] * dx;
      af[i-1] = al[i] + 0.5 * courant * (da[i] + (1. - 0.6666666 * courant) * a6[i]);
    }
  for (i = n-1; i; i--)  {
    loc = m * i;
    if (!p[loc])  {
      if (v[loc] > 0.)  {
	courant = v[loc] * dx;
	af[i] = ar[i] - 0.5 * courant * (da[i] - (1. - 0.6666666 * courant) * a6[i]);
      }
      loc -= m;
      if (v[loc] < 0.)  {
	courant = -v[loc] * dx;
	af[i-1] = al[i] + 0.5 * courant * (da[i] + (1. - 0.6666666 * courant) * a6[i]);
      }
    }
  }
  loc = m * i;
  if (!p[loc])
    if (v[loc] > 0.)  {
      courant = v[loc] * dx;
      af[i] = ar[i] - 0.5 * courant * (da[i] - (1. - 0.6666666 * courant) * a6[i]);
    }
}    /* That's it ! */


static void VerticalPPM(double *v, double *a, double tp, long tinc, int fa, int iord, double *af)
{
  int i, j, loc;
  const double ofac = 1./12.;
  double t1, t2, courant, deltaa, deltap, tmp;
  /* Stage 1: Interpolate values to the faces */
  if (fa < nzm)
    af[fa] = (level[fa+1] * a[fa*layer] + level[fa] * a[(fa+1)*layer]) /
             (level[fa] + level[fa+1]);
  af[nzm] = 0.5 * (a[nzm*layer] + tp);
  if (nzm-1 > fa)
    af[nzm-1] = (level[nzm] * a[(nzm-1)*layer] + level[nzm-1] * a[nzm*layer]) /
  	        (level[nzm] + level[nzm-1]);
  if (fa + 2 < nzm)  {
    i = fa+1; loc = i*layer;
    deltaa = f5[i] * (a[loc+layer] - a[loc]) + f6[i] * (a[loc] - a[loc-layer]);
    deltaa = Smooth1(deltaa, a[loc-layer], a[loc], a[loc+layer]);
    while (i < nzm-1)  {
      loc = i * layer;
      deltap = f5[i+1] * (a[loc+2*layer] - a[loc+layer]) + f6[i+1] * (a[loc+layer] - a[loc]);
      deltap = Smooth1(deltap, a[loc], a[loc+layer], a[loc+2*layer]);
      af[i] = a[loc] + f1[i] * (a[loc+layer] - a[loc]) +
	      f3[i] * deltap + f4[i] * deltaa;
      deltaa = deltap;
      i++;
    }
  }
  /* Stage 2: Determine ends of parabels */
  i = fa;
  al[i] = ar[i] = a[layer*i];
  da[i] = a6[i] = 0.;
  for (i++; i < nz; i++)  {
    al[i] = af[i-1];
    ar[i] = af[i];
    tmp = a[layer*i];
    if (iord == 0 || (ar[i] - tmp) * (tmp - al[i]) <= 0.)
      ar[i] = al[i] = tmp;
    else {
      t2 = ar[i] - al[i];
      t1 = t2 * (tmp - 0.5 * (al[i] + ar[i]));
      t2 = t2*t2 / 6.;
      if (t1 > t2)  al[i] = 3. * tmp - 2. * ar[i];
      if (t1 < -t2) ar[i] = 3. * tmp - 2. * al[i];
    }
    da[i] = ar[i] - al[i];
    a6[i] = 6. * (tmp - 0.5 * (al[i] + ar[i]));
  }
  al[nz] = ar[nz] = tp;
  da[nz] = a6[nz] = 0.;
  /* Stage 3: Calculate zone averages */
  for (i = fa; i < nzm; i++)  {
    loc = i*layer;
    if (v[i] >= 0.)  {
      courant = v[i] * tinc / dz[i];
      af[i] = ar[i] - 0.5 * courant * (da[i] - (1. - 0.6666666 * courant) * a6[i]);
    }
    else  {
      j = i+1;
      courant = -v[i] * tinc / dz[j];
      af[i] = al[j] + 0.5 * courant * (da[j] + (1. - 0.6666666 * courant) * a6[j]);
    }
  }
  loc = i*layer;
  if (v[i] >= 0.)  {
    courant = v[i] * tinc / dz[i];
    af[i] = ar[i] - 0.5 * courant * (da[i] - (1. - 0.6666666 * courant) * a6[i]);
  }
}    /* That's it ! */


static void PPM_3D(double *a, double *toplayer, long tinc, int iord)
/*
   Transport "a" in three dimensions
*/
{
  int i, j, k, loc, loc2, fa;
  double vv[MAXZGRID];
  /* Transport in x-Richtung */
  for (k = 0; k < nz; k++)  {
    for (j = ny; --j; )  {
      loc = k*layer + j;
      HorizontalPPM(nx, row, flux[UWIND]+loc, a+loc, pstat+loc, dx, tinc,
      		    iord, af,
#ifdef PARALLEL
		    leftest, rightest
#else
		    1, 1
#endif
      		    );
      for (i = nx; --i; )  {
        loc2 = loc + row*i;
        if (!pstat[loc2])
          b[loc2] = a[loc2] - tinc / dx *
                    (flux[UWIND][loc2] * af[i] - flux[UWIND][loc2-row] * af[i-1] -
                     a[loc2] * (flux[UWIND][loc2] - flux[UWIND][loc2-row]));
      }
      b[loc] = a[loc]; b[loc+nx*row] = a[loc+nx*row];
    }
    for (i = nx; i >= 0; i--)  {
      loc = k*layer + i*row;
      b[loc] = a[loc]; b[loc+ny] = a[loc+ny];
    }
    SetBoundX(nx, b+k*layer, flux[UWIND]+k*layer);
    SetBoundY(ny, b+k*layer, flux[VWIND]+k*layer);
  }
  /* Transport in y-Richtung */
  for (k = 0; k < nz; k++)  {
    for (i = nx; --i; )  {
      loc = k*layer + i*row;
      HorizontalPPM(ny, 1, flux[VWIND]+loc, b+loc, pstat+loc, dy, tinc,
      		    iord, af, 1, 1);
      for (j = ny; --j; )  {
        loc2 = loc + j;
        if (!pstat[loc2])
          b[loc2] -= tinc / dy *
          	     (flux[VWIND][loc2] * af[j] - flux[VWIND][loc2-1] * af[j-1] -
          	      a[loc2] * (flux[VWIND][loc2] - flux[VWIND][loc2-1]));
      }
    }
    SetBoundX(nx, b+k*layer, flux[UWIND]+k*layer);
    SetBoundY(ny, b+k*layer, flux[VWIND]+k*layer);
  }
  /* Transport in z-Richtung */
  for (i = nx; --i; )
    for (j = ny; --j; )  {
      loc = i*row + j;
      fa = ground[loc].firstabove;
      for (k = fa; k < nz; k++)
        vv[k] = flux[WWIND][loc+k*layer] * (flux[WWIND][loc+k*layer] > 0. ? density[k] : density[k+1]);
      VerticalPPM(vv, b+loc, toplayer[loc], tinc, fa, iord, af);
      k = fa; loc2 = loc + k*layer;
      a[loc2] = b[loc2] - tinc / dz[k] * vv[k] * (af[k] - a[loc2]);
      for (k = fa+1; k < nz; k++)  {
        loc2 = loc + k*layer;
        a[loc2] = b[loc2] - tinc / dz[k] *
        	  (vv[k] * af[k] - vv[k-1] * af[k-1] -
		   a[loc2] * (vv[k] - vv[k-1]));
      }
    }
  for (k = fa; k < nz; k++) {
    SetBoundX(nx, a+k*layer, flux[UWIND]+k*layer);
    SetBoundY(ny, a+k*layer, flux[VWIND]+k*layer);
  }
}  /* End of PPM_3D */


void PPM_Transport(long tinc, long chemtinc, int iord)
{
  Entity et;
  double *x, tmp;
  int i;
  CalcPPMGeoFact();
  if (windadvection)  {
    for (et = TEMP - (pressuretype != NONHYDROSTATIC); et--; )  {
      if (coriolistype == GEOSTROPHICCORIOLIS)  {
        switch (et)  {
	  case UWIND : tmp = ugeos[nzm]; break;
	  case VWIND : tmp = vgeos[nzm]; break;
	  case WWIND : tmp = 0.; break;
	}
	for (i = layer, x = top; i--; x++)  *x = tmp;
      }
      else
	memcpy(top, g[et]+nzm*layer, layer*sizeof(double));   /* Die aeltere Formulierung */
      PPM_3D(g[et], top, tinc, iord);
    }
  }
  if (advection)  {
    for (i = layer; i--; )
      top[i] = g[TEMP][nzm*layer+i] + toptemp[i];
    PPM_3D(g[TEMP], top, tinc, iord);
    for (et = HUMIDITY; et < (chemtinc ? maxentity : SUBS); et++)
      if (g[et] && (et < SUBS || subst[et-SUBS+1].do_transport))  {
	if (et == TKE || et == EPS)  {
	  memcpy(top, g[et]+nzm*layer, layer * sizeof(double));
	  PPM_3D(g[et]+layer, top, tinc, iord);
	}
	else  {
	  memcpy(top, g[et]+nzm*layer, layer * sizeof(double));
	  PPM_3D(g[et], top, (et < SUBS ? tinc : chemtinc), iord);
	}
      }
  }
}  /* End of PPM_Transport */
