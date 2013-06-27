/*
   IMPLEMENTATION MODULE mpdata.c
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "mcglobal.h"
#include "mpdata.h"
#ifdef PARALLEL
#include <pvm3.h>
#include "mcparallel.h"
#endif

void MinMax(double *min, double *max, ...)
{
  va_list ap;
  double *v;
  va_start(ap, max);
  *min = *max = *(va_arg(ap, double *));
  while (v = va_arg(ap, double *))  {
    if (*v > *max)  *max = *v;
    if (*v < *min)  *min = *v;
  }
  va_end(ap);
}

void SetBound1(int n, double *x)
{
  x[0] = x[1];
  x[n] = x[n-1];
}

void SetBound2(int nx, int ny, double *x, double *vx, double *vy)
{
  int i, j, nxm, nym;
  nxm = nx-1; nym = ny-1;
  switch (southbordertype)  {
    case AUTOCONSTANTBORDER :
       for (i = nx; --i; )
         if (vy[i*(row+1)+1] < 0.)
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
         if (vy[i*(row+1)+ny] > 0.)
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
  switch (westbordertype)  {
    case AUTOCONSTANTBORDER :
       for (j = ny; --j; )
         if (vx[row+j] < 0.)
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
         if (vx[nx*row+j] > 0.)
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

#define PROCNAME1 MpData1
#define PROCNAME2 MpData2
#define VMDIV(a,b)  ((a) - (b)) / ((a) + (b) + eps)
#define DONOR(y1,y2,a) (((a) > 0. ? (y1) : (y2))* (a))
#define VDYF(x1,x2,a)  ((fabs(a) - a*a)*VMDIV(x2,x1))
#define VCORR(a,b,y1,y2)  (-0.125 * (a) * (b) * (y1) / (y2))
#define VDIV1(a1,a2,a3)  (0.25 * (a2) * ((a3) - (a1)))
#define PP(y)   (y > 0. ? y : 0.)
#define PN(y)   (y < 0. ? -y : 0.)
#define GTP(i)  (x[i] > 0. ? cp[i] : cn[i])
#define GTM(i)  (x[i] < 0. ? cp[i] : cn[i])
#define GTCP(i,j)  (x[(i)*row+(j)] > 0. ? cp[(i)*row+(j)] : cn[(i)*row+(j)])
#define GTCM(i,j)  (x[(i)*row+(j)] < 0. ? cp[(i)*row+(j)] : cn[(i)*row+(j)])

double Min3(double a, double b, double c)
{
  return (a < b ? (a < c ? a : c) : (b < c ? b : c));
}

#define VDIV2(a,b1,b2,b3,b4)  (0.25*a*(b1+b2-b3-b4))

void SwapPtr(void **a, void **b)
{
  void *tmp;
  tmp = *a; *a = *b; *b = tmp;
}

/* #include "mpdata.incl" */

#undef VMDIV
#define VMDIV(a,b)  (fabs(a) - fabs(b)) / (fabs(a) + fabs(b) + eps)
#undef VDYF
#define VDYF(x1,x2,a)  ((fabs(a) - a*a)*VMDIV(x2,x1))
#undef PROCNAME1
#define PROCNAME1 NpData1
#undef PROCNAME2
#define PROCNAME2 NpData2

/* #include "mpdata.incl" */

void MpAdvect(int tinc,
	      double *ux, double *uy, double *x, double *xp, double *wwind, double *dcm,
              double zfact, double dens, double denskp, 
              int nx, int ny, int iord, BOOL nonos)
/* The Routine MpAdvect will calculate the transport of a spatial distibuted variable
   simultaneously in three directions. It is based on Smolarkievicz' mpdata-Advection
   Scheme. mpdata has to be invoked once for each level. The Variables mean: 
   tinc         : deltaT, the time step of integration.
   ux, vy	: Wind speeds at the faces of the grid cells, multiplied with the
   		  time step.
   x		: The Variable to be transported, on the actual level.
   xp		: The Variable to be transported, on the next higher level.
   wwind	: vertical wind speed, un the face of the grid cell above the values "x"
   dcm		: A help variable of the size nx*ny, to read the fluxes from the
   		  neighbour cells underneath the actual layer, and to store the
   		  fluxes into the upper neighbours.
   etc.
*/
{
  int i, j, loc;
  const double eps = 1.e-10;
  double *vx, *fx, *cp, *cn, *mx, *mn;
  double *vy, *fy;
  double hlp;
  for (i = nx+1; i--; )
    for (j = ny+1; j--; )
       x[i*row+j] *= dens;
  if (!(vx = NEW((nx+2)*row, double)) ||
      !(fx = NEW((nx+2)*row, double)) ||
      !(vy = NEW((nx+1)*(row+1), double)) ||
      !(fy = NEW((nx+1)*(row+1), double)))  {
    fprintf(stderr, "Fatal Error: Unable to allocate arrays in npdata2\n");
    exit (1);
  }
  if (nonos)  {
    if (!(cp = NEW(layer, double)) ||
        !(cn = NEW(layer, double)) ||
        !(mx = NEW(layer, double)) ||
        !(mn = NEW(layer, double)))  {
      fprintf(stderr, "Fatal Error: Unable to allocate arrays in npdata2\n");
      exit (1);
    }
    for (i = nx; --i; )
      for (j = ny; --j; )  {
        loc = i*row+j;
        MinMax(&mn[loc], &mx[loc], &x[loc-row], &x[loc], &x[loc+row],
               &x[loc-1], &x[loc+1], NULL);
      }
  }
  /* Copy speeds to fx, fy for continuing with the algorithm. Take special care
     of different formats of the arrays. Multiply by the time step. */
  for (i = nx+1; i--; )
    for (j = ny+1; j--; )  {
      loc = i*row+j;
      vx[loc+row] = tinc * dxi * ux[loc];
      vy[loc+i+1] = tinc * dyi * uy[loc];
    }
  /* Set Flux-Values at the boundaries */
  for (i = nx+1; i--; )  {
    loc = i*(row+1);
    vy[loc] = vy[loc+1];
    vy[loc+ny+1] = vy[loc+ny];
    loc -= i;
    vx[loc] = vx[loc+1];
    vx[loc+ny] = vx[loc+nym];
  }
  for (j = ny+1; j--; )  {
    vx[j] = vx[j+row];
    vx[j+(nx+1)*row] = vx[j+nx*row];
    vy[j] = vy[j+row+1];
    vy[j+nx*(row+1)] = vy[j+nxm*(row+1)];
  }
  for (i = nx+1; --i; )
    for (j = ny; --j; )  {
      loc = i*row+j;
      fx[loc] = DONOR(x[loc-row], x[loc],  vx[loc]);
    }
  for (i = nx; --i; )
    for (j = ny+1; --j; )  {
      loc = i*row+j;
      fy[loc+i] = DONOR(x[loc-1], x[loc], vy[loc+i]);
    }
  for (i = nx; --i; )
    for (j = ny; --j; )  {
      loc = i*row+j;
      x[loc] += fx[loc] - fx[loc+row] + fy[loc+i] - fy[loc+i+1] +
                    (dcm[loc] -
                     (hlp = DONOR(x[loc], xp[loc] * denskp, wwind[loc]))) * tinc * zfact;
      dcm[loc] = hlp;
    }
  SetBound2(nx, ny, x, vx, vy);
  while (iord--)  {
    SwapPtr((void *)&fx, (void *)&vx);
    SwapPtr((void *)&fy, (void *)&vy);
    for (j = ny; --j; )
      for (i = nx+1; --i; )  {
        loc = i*row+j;
	vx[loc] = VDYF(x[loc-row], x[loc], fx[loc]) -
	   0.125 * fx[loc] * (fy[(i-1)*(row+1)+j] + fy[(i-1)*(row+1)+j+1] + fy[loc+i+1] + fy[loc+i]) *
	   VMDIV(x[loc-row+1] + x[loc+1], x[loc-row-1] + x[loc-1]);
      }
    for (i = nx; --i; )
      for (j = ny+1; --j; )  {
        loc = i*row+j;
        vy[loc+i] = VDYF(x[loc-1], x[loc], fy[loc+i]) -
           0.125 * fy[loc+i] * (fx[loc-1] + fx[loc] + fx[loc+row] + fx[loc+row-1]) *
           VMDIV(x[loc+row-1] + x[loc+row], x[loc-row-1] + x[loc-row]);
      }
    if (nonos)  {
      for (i = nx; --i; )
        for (j = ny; --j; )  {
          loc = i*row+j;
          MinMax(&mn[loc], &mx[loc], &x[loc-row], &x[loc], &x[loc+row],
                 &x[loc-1], &x[loc+1] , &mn[loc], &mx[loc], NULL);
        }
      for (i = nx+1; --i; )
	for (j = ny; --j; )  {
	  loc = i*row+j;
	  fx[loc] = DONOR(x[loc-row], x[loc], vx[loc]);
	}
      for (i = nx; --i; )
	for (j = ny+1; --j; )  {
	  loc = i*row+j;
	  fy[loc+i] = DONOR(x[loc-1], x[loc], vy[loc+i]);
	}
      for (i = nx; --i; )
        for (j = ny; --j; )  {
          loc = i*row+j;
          cp[loc] = (mx[loc] - x[loc]) /
             (PN(fx[loc+row]) + PP(fx[loc]) + PN(fy[loc+i+1]) + PP(fy[loc+i]) + eps);
          cn[loc] = (x[loc] - mn[loc]) /
             (PP(fx[loc+row]) + PN(fx[loc]) + PP(fy[loc+i+1]) + PN(fy[loc+i]) + eps);
        }
      for (i = 2; i < nx; i++)
        for (j = 2; j < ny; j++)  {
          loc = i*row+j;
          vx[loc] = vx[loc] * (vx[loc] > 0. ? 
          			 Min3(1., GTCP(i, j), GTCM(i-1, j)) :
                                 Min3(1., GTCP(i-1, j), GTCM(i, j)));
          vy[loc+i] = vy[loc+i] * (vy[loc+i] > 0. ?
           			 Min3(1., GTCP(i, j), GTCM(i, j-1)) :
                     		 Min3(1., GTCP(i, j-1), GTCM(i, j)));
        }
    }
    for (i = nx+1; --i; )
      for (j = ny; --j; )  {
        loc = i*row+j;
	fx[loc] = DONOR(x[loc-row], x[loc], vx[loc]);
      }
    for (i = nx; --i; )
      for (j = ny+1; --j; )  {
        loc = i*row+j;
	fy[loc+i] = DONOR(x[loc-1], x[loc], vy[loc+i]);
      }
    for (i = nx - !!(eastbordertype & DATABORDER); --i > !!(westbordertype & DATABORDER); )
      for (j = ny - !!(northbordertype & DATABORDER); --j > !!(southbordertype & DATABORDER); )  {
        loc = i*row+j;
	x[loc] += fx[loc] - fx[loc+row] + fy[loc+i] - fy[loc+i+1];
      }
    SetBound2(nx, ny, x, vx, vy);
  }
  free(vx); free(fx); free(vy); free(fy);
  if (nonos)  {
    free(cp); free(cn); free(mx); free(mn);
  }
  for (i = nx+1; i--; )
    for (j = ny+1; j--; )
       x[i*row+j] /= dens;
}
