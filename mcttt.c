/*
   IMPLEMENTATION MODULE MCTTT
   Dieses Modul implementiert die Transilient Turbulent Theory.
*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "mcglobal.h"
#include "mcpress.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcttt.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

typedef double Row[101];

void PrintTMatrix(int n, Row *m, char *title)
{
  int i, j;
  printf("%s\n\n       ", title);
  for (j = 0; j < n; j++)
    printf("%8d ", j);
  printf("\n");
  for (i = 0; i < n; i++)  {
    printf("%6d ", i);
    for (j = 0; j < n; j++)  {
      printf("%8.4g ", m[i][j]);
    }
    printf("\n");
  }
}

void CalcTransilientMatrix(int ig, int jg, double tinc, Row *ttm, double *deltaz)
/* Berechnet die Transilient Turbulence Matrix an der Stelle (i/j). */
{
  const double t0 = 1000., Dy = 1., Yref = 1000., Rc = 0.21;
  int kg, i, j, fa, nte, loc;
  double dUWIND, dVWIND, dTEMP, dzl, sum, Ymax;
  GroundParam *gr;
  Row relmass;
  gr = ground+ig*row+jg;
  fa = gr->firstabove;
  loc = fa*layer+ig*row+jg;
/* Bestimme die Y(ij)-Werte, i!=j. (Stull p. 235)  */
  sum = 0;
  for (i = nte = nz - fa; i--; )  {
    dzl = 0.;
    for (j = i; j--; )  {
      dzl += ldiff[j+fa];
      dUWIND = g[UWIND][i*layer+loc] - g[UWIND][j*layer+loc];
      dVWIND = g[VWIND][i*layer+loc] - g[VWIND][j*layer+loc];
      dTEMP  = g[TEMP][i*layer+loc] - g[TEMP][j*layer+loc];
      ttm[i][j] = ttm[j][i] = 
          tinc * t0 / (dzl * dzl) *
          (dUWIND*dUWIND + dVWIND*dVWIND - Grav*dzl*dTEMP/(g[TEMP][j*layer+loc]*Rc)) -
          Dy * tinc / t0;
    }
    sum += density[i+fa] * dzl;
  }
/* Die Werte werden monoton steigend gemacht, und die Y(ii) Werte gerechnet. */
  for (i = nte; i--; )  {
    ttm[i][i] = 0.;
    for (j = 0; j <= i; j++)  {
      if (ttm[i][j] < 0.)  ttm[i][j] = ttm[j][i] = 0.;
      if (j && ttm[i][j-1] > ttm[i][j])  ttm[i][j] = ttm[j][i] = ttm[i][j-1];
      if (i < nte-1 && ttm[i+1][j] > ttm[i][j])  ttm[i][j] = ttm[j][i] = ttm[i+1][j];
    }
    ttm[i][i] += Yref;
  }
/* Berechne relative Massen */
  for (i = 0; i < nte; i++)
    deltaz[i] = 2. * level[i];
  sum = 0.5 * nte * nte / sum;
  for (i = nte; i--; )
    relmass[i] = deltaz[i] * density[i+fa] * sum;
/* Bestimme die maximale Spaltensumme. */
  Ymax = 0.;
  for (i = nte; i--; )  {
    sum = 0;
    for (j = nte; j--; )
      sum += ttm[i][j];
    if (sum > Ymax)  Ymax = sum;
  }
/* Skaliere die Matrix */
  Ymax = 1. / Ymax;
  for (i = nte; i--; )
    for (j = nte; j--; )
      if (i != j)  ttm[i][j] *= Ymax * relmass[j];
/* Berechne abschliessend die Diagonalelemente c(ii) */
  for (i = nte; i--; )  {
    ttm[i][i] = 1.;
    for (j = nte; j--; )
      if (i != j)  ttm[i][i] -= ttm[i][j];
  }
/*  if (ig == 74 && jg == 37)
    PrintTMatrix(nte, ttm, "Diagonal-Elemente berechnet");  */
}

void CalcTransilientExchange(int ig, int jg, Entity et, double *v, Row *ttm,
			     double *deltaz)
{
  Row nv;
  double *p[100];
  int i, j, nte, loc;
  loc = ig*row+jg;
  nte = nz - ground[loc].firstabove;
  if (nte < 2)  return;
  for (i = nte; --i; )
    p[i] = v+i*layer+ig*row+jg;
  for (i = nte; i--; )  {
    nv[i] = 0.;
    for (j = nte; j--; )
      nv[i] += ttm[i][j] * *p[j];
  }
  for (i = nte; i--; )
    *p[i] = nv[i];
}


void CalcTransilientTurbulence(double tinc)
{
  Row *ttm;
  int i, j, k, fa, loc;
  Entity et;
  double K, h, d;
  Row deltaz;
  int xs, xe;
#ifdef PARALLEL
  if (parallel && master) return; // Nothing to do in this case
  if (parallel && !master)  {
    xs = mfirstx + !mfirstx;
    xe = mlastx - rightest;
  }
  else  {
#endif
    xs = 1; xe = nx;
#ifdef PARALLEL
  }
#endif
  ttm = (Row *)calloc(nz+1, sizeof(Row));
  for (i = xs; i < xe; i++)
    for (j = ny; --j; )  {
      fa = ground[i*row+j].firstabove;
      CalcTransilientMatrix(i, j, tinc, ttm, deltaz);
      for (et = maxentity; et--; )
        if (g[et])  CalcTransilientExchange(i, j, et, g[et] + (fa - 1)*layer, ttm, deltaz);
      for (k = fa; k < nz; k++)
        Km[k*layer+i*row+j] = 0.2875 * sqr(deltaz[k-fa+1]) * (1. - ttm[k-fa+1][k-fa+1]) / tinc;
    }
  free(ttm);
}

void CalcGroundTurbulence(double tinc)
{
  int i, j, k, loc;
  GroundParam *gr;
  double de;
  Entity et;
  double h2;
  int xs, xe;
#ifdef PARALLEL
  if (parallel && master) return; // Nothing to do in this case
  if (parallel && !master)  {
    xs = mfirstx + !mfirstx;
    xe = mlastx - rightest;
  }
  else  {
#endif
    xs = 1; xe = nx;
#ifdef PARALLEL
  }
#endif
  for (i = xs; i < xe; i++)
    for (j = ny; --j; )  {
      gr = &ground[i*row+j];
      k = gr->firstabove;
      loc = k*layer + i*row + j;
      h2 = 1. / level[k];
      g[UWIND][loc] *= gr->uw;
      g[VWIND][loc] *= gr->vw;
      g[TEMP][loc] += tinc * gr->wtheta * h2;
      if (calcsoiltemp)  {
	g[HUMIDITY][loc] += tinc * gr->wq * h2;
	gr->X += tinc * (1. - gr->X - gr->r * gr->wq * Lc * density[k]) / 1000.;
	if (gr->X > 1.)  gr->X = 1.;
	if (gr->X < 0.)  gr->X = 0.;
      }
    }
}
