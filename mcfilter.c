/*
   IMPLEMENTATION MODULE mcfilter.c
   Dieses Modul implementiert ein raeumliches Filter. Es handelt sich um das
   auf Filter nach Pepper et al. (1979) - gemaess Pielke p. 329
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcglobal.h"
#include "mcfilter.h"

#define F(a,b)  ((a)*row+(b))

void ApplyFilter(double *lay, PointStatus *pstat, double delta)
{
  int i, j, k, from, to, n;
  double b, a, *smooth, *diag;
  b = 2.*(1 + delta);
  a = 1. - delta;
  n = (nx > ny ? nx : ny);
  smooth = (double *)calloc(n, sizeof(double));
  diag = (double *)calloc(n, sizeof(double));
  for (j = ny; --j; )  {
    from = 2;
    while (from < nxm)  {
      while (from < nxm && pstat[F(from-1,j)])  from++;
      if (from < nxm && !pstat[F(from,j)] && !pstat[F(from+1,j)])  {
        for (to = from; to < nxm && !pstat[F(to+1,j)]; to++)
          smooth[to] = lay[F(to-1,j)] + 2. * lay[F(to,j)] + lay[F(to+1,j)];
        to--;
        smooth[from] -= a * lay[F(from-1,j)];
        diag[from] = b;
        for (k = from + 1; k <= to; k++)  {
          smooth[k] -= a / diag[k-1] * smooth[k-1];
          diag[k] = b - a * a / diag[k-1];
        }
        smooth[to] -= a * lay[F(to+1,j)];
        for (k = to - 1; k >= from; k--)  {
          smooth[k] -= a / diag[k+1] * smooth[k+1];
        }
        for (k = from; k <= to; k++)
          lay[F(k,j)] = smooth[k] / diag[k];
        from = to + 2;
      }
      else  from++;
    }
  }
  for (i = nx; --i; )  {
    from = 2;
    while (from < nym)  {
      while (from < nym && pstat[F(i,from-1)])  from++;
      if (from < nym && !pstat[F(i,from)] && !pstat[F(i,from+1)])  {
        for (to = from; to < nym && !pstat[F(i,to+1)]; to++)
          smooth[to] = lay[F(i,to-1)] + 2. * lay[F(i,to)] + lay[F(i,to+1)];
        to--;
        smooth[from] -= a * lay[F(i,from-1)];
        diag[from] = b;
        for (k = from + 1; k <= to; k++)  {
          smooth[k] -= a / diag[k-1] * smooth[k-1];
          diag[k] = b - a * a / diag[k-1];
        }
        smooth[to] -= a * lay[F(i,to+1)];
        for (k = to - 1; k >= from; k--)  {
          smooth[k] -= a / diag[k+1] * smooth[k+1];
        }
        for (k = from; k <= to; k++)
          lay[F(i,k)] = smooth[k] / diag[k];
        from = to + 2;          
      }
      else  from++;
    }
  }
  free(smooth); free(diag);
}

void ShapiroStep(double fact, double *x, double *flx, PointStatus *pstat)
{
  int i, j, loc;
  if (nx > 3)  {
    for (i = nx; i--; )
      for (j = ny+1; j--; )  {
	loc = F(i,j);
	flx[loc] = (!pstat[loc+row] && !pstat[loc] ? x[loc+row] - x[loc] : 0.);
      }
    for (j = ny+1; j--; )  {
      for (i = nx; --i; )  {
	loc = F(i,j);
	if (!pstat[loc])
          x[loc] += fact * (flx[loc] - flx[loc-row]);
      }
      if (!pstat[j])  x[j] += fact * flx[j];
      if (!pstat[F(nx,j)])  x[F(nx,j)] -= fact * flx[F(nxm,j)];
    }
  }
  if (ny > 3)  {
    for (i = nx+1; i--; )
      for (j = ny; j--; )  {
	loc = F(i,j);
	flx[loc] = (!pstat[loc] && !pstat[loc+1] ? x[loc+1] - x[loc] : 0.);
      }
    for (i = nx+1; i--; )  {
      for (j = ny; --j; )  {
	loc = F(i,j);
	if (!pstat[loc])
          x[loc] += fact * (flx[loc] - flx[loc-1]);
      }
      if (!pstat[i*row])  x[i*row] += fact * flx[i*row];
      if (!pstat[F(i,ny)])  x[F(i,ny)] -= fact * flx[F(i,nym)];
    }
  }
}

void ShapiroStepFour(int level, double *x, double *flx, PointStatus *pstat, double fval)
{
  int i, j, loc;
  double ifv, ofv;
  const double fltfct[3] =
     { 0.1875, 
       -0.1660533905932737622,
       0.5410533905932737622 };
  ifv = fltfct[level] * fval;
  ofv = 0.0625 * fval;
  if (nx > 3)  {
    for (i = nx; i--; )
      for (j = ny; --j; )  {
	loc = F(i,j);
	flx[loc] = (!pstat[loc+row] && !pstat[loc] ? x[loc+row] - x[loc] : 0.);
      }
    for (j = ny+1; j--; )  {
      for (i = nxm; --i > 1; )  {
	loc = F(i,j);
	if (!pstat[loc])  {
          x[loc] += ifv * (flx[loc-row] - flx[loc]) +
        	     ofv * (flx[loc+row] - flx[F(i-2,j)]);
          if (pstat[loc+row])  x[loc] -= ofv * flx[loc-row];
          if (pstat[loc-row])  x[loc] += ofv * flx[loc];
	}
      }
      if (!pstat[j])
	x[j] += ofv * flx[row+j] + (ofv - ifv) * flx[j];
      if (!pstat[row+j])  {
	x[row+j] += ifv * (flx[j] - flx[row+j]) + ofv * flx[F(2,j)];
	if (pstat[F(2,j)])  x[row+j] -= ofv * flx[j];
      }
      loc = F(nx,j);
      if (!pstat[loc])
	x[loc] += (ifv - ofv) * flx[loc-row] - ofv * flx[F(nxm-1,j)];
      loc -= row;
      if (!pstat[loc])  {
	x[loc] += ifv * (flx[loc-row] - flx[loc]) - ofv * flx[F(nxm-2,j)];
	if (pstat[loc-row])  x[loc] += ofv * flx[loc];
      }
    }
  }
  if (ny > 3)  {
    for (i = nx; --i; )
      for (j = ny; j--; )  {
	loc = F(i,j);
	flx[loc] = (!pstat[loc] && !pstat[loc+1] ? x[loc+1] - x[loc] : 0.);
      }
    for (i = nx+1; i--; )  {
      for (j = nym; --j > 1; )  {
	loc = F(i,j);
	if (!pstat[loc])  {
          x[loc] += ifv * (flx[loc-1] - flx[loc]) +
        	     ofv * (flx[loc+1] - flx[loc-2]);
          if (pstat[loc+1])  x[loc] -= ofv * flx[loc-1];
          if (pstat[loc-1])  x[loc] += ofv * flx[loc];
	}
      }
      loc = i*row;
      if (!pstat[loc])
	x[loc] += ofv * flx[loc+1] + (ofv - ifv) * flx[loc];
      if (!pstat[loc+1])  {
	x[loc+1] += ifv * (flx[loc] - flx[loc+1]) + ofv * flx[loc+2];
	if (pstat[loc+2])  x[loc+1] -= ofv * flx[loc];
      }
      if (!pstat[loc+ny])
	x[loc+ny] += (ifv - ofv) * flx[loc+nym] - ofv * flx[loc+nym-1];
      if (!pstat[loc+nym])  {
	x[loc+nym] += ifv * (flx[loc+nym-1] - flx[loc+nym]) - ofv * flx[loc+nym-2];
	if (pstat[loc+nym-1])  x[loc+nym] += ofv * flx[loc+nym];
      }
    }
  }
}

void ShapiroFilter(double *x, PointStatus *pstat, int ord, double filtval)
{
  int i, j, loc;
  double oneminf;
  memcpy(tmplayer2, x, layer*sizeof(double));
  if (ord--)  ShapiroStep(0.25, x, tmplayer, pstat);
  if (ord--)  ShapiroStep(-0.25, x, tmplayer, pstat);
  if (ord--)  ShapiroStepFour(0, x, tmplayer, pstat, 1.);
  if (ord--)  {
    ShapiroStepFour(1, x, tmplayer, pstat, 1.);
    ShapiroStepFour(2, x, tmplayer, pstat, 1.);
  }
  oneminf = 1. - filtval;
  if (filtval < 1.)
    for (i = nx+1; i--; )
      for (j = ny+1; j--; )  {
        loc = F(i,j);
        if (!pstat[loc])
          x[loc] = filtval * x[loc] + oneminf * tmplayer2[loc];
      }
}

