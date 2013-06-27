/*
   IMPLEMENTATION MODULE McKEps
   Implementiert Boundary-Layer Turbulenz nach dem k-epsilon Schema.
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mchemparse.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

#define PP(x) ((x) > 0. ? x : 0.)
#define PM(x) ((x) < 0. ? x : 0.)

double ThRoot(double x)
{
  if (x > 0.)  return (pow(x, 1./3.));
  if (x < 0.)  return (-pow(-x, 1./3.));
  return (0.);
}

double SolveCubic(double a2, double a1, double a0)
{
  char line[256];
  const double tthpi = 2./3.*M_PI; /* 2/3*pi */
  double p, q, d, d2, cosphi, y;
  p = a1 - a2*a2/3.;
  q = a0 - a1*a2/3. + 2./27.*a2*a2*a2;
  d = q*q / 4. + (d2 = p*p*p / 27.);
  if (d < 0.)  {
    cosphi = -q / (2. * sqrt(-d2));
    cosphi = cos(acos(cosphi)/3. - tthpi);
    y = 2. * sqrt(-p/3.) * cosphi;
  }
  else if (d > 0.)  {
    d = sqrt(d);
    q /= -2.;
    y = ThRoot(q+d) + ThRoot(q-d);
  }
  else  y = 0.;
  y -= a2 / 3.;
  return (y);
}

void SolveKEps(short kimp, int tinc, double sh, double bu, double tr,
   double Km, double *tke, double *eps)
{
  int it, nit = 0;
  double pr, a, b, c, d, ets;
  const double c1 = 1.44, c2 = 1.92, tol = 0.00000001, mu = 0.09;
  const double tkemin = 0.001, epsmin = 0.000001;
  double c3, sb, Ea, epsa;
  double t1, t4, t5, t6, t7, t8, t9, t11, t20;
  c3 = c2 - 1.;
  if (!kimp)  {
/*   Analytical Version:  Explicit K */
/*    if (sh + bu > 0. && Km < 1.)  Km = 1.;
    sh *= Km; bu *= Km; */
    pr = sh + bu/* + PP(tr)*/;
    ets = *tke + tinc*(sh + bu);
    a = tinc*c3;
    b = ets + tinc*(*eps - c1*pr);
    c = -*eps * ets;
    d = b*b - 4.*a*c;
    if (d < 0.)  {
      *eps = epsmin;
      *tke = tkemin;
    }
    else  {
      *eps = (-b + sqrt(d)) / (2.*a);
      *tke = ets - tinc * *eps;
      if (*eps < epsmin)  *eps = epsmin;
      if (*tke < tkemin)  *tke = tkemin;
    }
  }
  else  {
/*
  Analytical Version: Implicit K:  
  The following formulae have been obtained with MapleV:
*/
    pr = sh + PP(bu);
    tr = PP(tr);
    sb = sh + bu;
    if (!*tke || !*eps)  {
      *tke = tkemin;
      *eps = epsmin;
    }
    Ea = *tke; epsa = *eps;
    t1 = c1*pr;
    t5 = c1*c1;
    t6 = pr*pr;
    t8 = t1*sb;
    t9 = sb*sb;
    t11 = mu*mu;
    t20 = tinc*tinc;
    a = ((sb-t1)*mu*c3+((t5*t6-2.*t8+t9)*t11+(2.*(t9-t8)*t11+
	t11*t9*c3)*c3)*t20)*tinc;
    b = -c3*epsa+(-2.*c2*mu*sb*epsa+((2.-c2)*tr*mu*sb+pr*mu*(2.*epsa-
	tr*c1))*c1)*t20+(-mu*sb*c2+(2.0*pr*c2-pr)*mu*c1)*tinc*Ea;
    t1 = epsa*epsa;
    t4 = c1*c1;
    t5 = tr*tr;
    t7 = sb*mu;
    t8 = tinc*tinc;
    c = (t1-epsa*tr*c1-t4*t5*t7*t8)*tinc+((2.*c2-1.)*epsa+
	(tr*c2*t7+pr*mu*tr*c1)*c1*t8-tinc*c1*pr*mu*c2*Ea)*Ea;
    d = (epsa*tinc*c1*tr-Ea*c2*epsa)*Ea;
    if (a == 0.)  {
      fprintf(stderr, "a is equal to zero in SolveKEps. Leaving\n");
      exit (0);
    }
    *tke = SolveCubic(b/a, c/a, d/a);
    *eps = *tke *(*tke *mu * tinc *(c1*pr - c2*sb) + epsa) /
        	 (-*tke * c3 + c2*Ea - tinc * c1 * tr);
    if (*tke <= tkemin)  *tke = tkemin;
    if (*eps < epsmin)  *eps = epsmin;
  }
}

void CalcKTurb(long tinc, long chemtinc)
{
  double tmp, tm[MAXZGRID][3], Kf[MAXZGRID];
  int i, j, k, l, loc, hloc;
  long ltinc;
  Entity et;
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
      hloc = i*row + j;
      /* Berechnung der K-Turbulenz-Matrix */
      k = ground[hloc].firstabove;
      if (k >= nzm)  continue;
      loc = hloc + k*layer;
      if (groundinterface)  {
	Kf[k] = (ground[hloc].zL > 0. ? 10. : Km[loc]) /
	  (0.5 * (level[k] + ground[hloc].z)) * 
	  density[k];
	/* Falls stabile Lage werden unterstes Niveau und Bodenniveau zusammengelegt.
	   Hier mit einem rel. hohen Austauschparameter realisiert */
      }
      else  Kf[k] = 0.;
      for (k++; k < nz; k++)  {
        loc = hloc + k*layer;
        Kf[k] = Km[loc] / ldiff[k] * 0.5 * (density[k] + density[k-1]);
      }
      for (l = 2; (l--) - !(nsubs && chemtinc); )  {
        ltinc = (l ? tinc : chemtinc);
        k = ground[hloc].firstabove;
        loc = hloc + k*layer;
        tmp = ltinc / (ground[hloc].z * density[k] / ground[hloc].slope.z);
        tm[k][1] = 1. + tmp * Kf[k];
        tm[k][2] = -tmp * Kf[k];
        for (; k < nzm; k++)  {
          loc = hloc + k*layer;
          tmp = ltinc / (level[k] * density[k]);
          tm[k+1][0] = -tmp * Kf[k];
          tm[k+1][1] = 1. + tmp*(Kf[k] + Kf[k+1]);
          tm[k+1][2] = -tmp * Kf[k+1];
        }
        loc = hloc + k*layer;
        tmp = ltinc / (level[k] * density[k]);
        tm[k+1][0] = -tmp * Kf[k];
        tm[k+1][1] = 1. + tmp * Kf[k];
        /* Konditioniere Matrix zum schnellen Lösen des Systems */
        for (k = ground[hloc].firstabove+1; k <= nz; k++)  {
          tm[k][0] = tmp = -tm[k][0] / tm[k-1][1];
          tm[k][1] += tmp * tm[k-1][2];
        }
        for (k = nzm; k >= ground[hloc].firstabove; k--)  {
          tm[k][2] = -tm[k][2] / tm[k+1][1];
        }
        for (k = ground[hloc].firstabove; k <= nz; k++)
          tm[k][1] = 1. / tm[k][1];
        /* Rechne jetzt die Turbulenz fuer alle Meteo-Groessen mit einer
           impliziten Integration */
        for (et = (l ? UWIND : SUBS); et < (l ? HUMIDITY+1 : maxentity); et++)
          if (et < SUBS || subst[et-SUBS+1].do_transport)  {
            k = ground[hloc].firstabove;
            loc = hloc + k * layer;
            if (l && groundinterface)
              g[et][loc] += tm[k+1][0] * ground[hloc].a[et];
            for (k++ ; k < nz; k++)  {
              loc = hloc + k*layer;
              g[et][loc] += tm[k+1][0] * g[et][loc-layer];
            }
            for (k-- ; --k >= ground[hloc].firstabove;)  {
              loc = hloc + k*layer;
              g[et][loc] += tm[k+1][2] * g[et][loc+layer];
            }
            k++;
            loc = hloc + k * layer;
            if (l && groundinterface)  {
              ground[hloc].a[et] += tm[k][2] * g[et][loc];
              ground[hloc].a[et] *= tm[k][1];
            }
            for ( ; k < nz; k++)  {
              loc = hloc + k*layer;
              g[et][loc] *= tm[k+1][1];
            }
          }  /* if (et */
        Kf[ground[hloc].firstabove] = 0.;
      }  /* for (l */
    }
}

/*
void PrintESum(char *cmmt, int hloc)
{
  double sum;
  int i;
  i = ground[hloc].firstabove;
  sum = g[EPS][hloc+i*layer] * density[i] * (ground[hloc].z + level[i]*0.5);
  for (i++; i < nz; i++)
    sum += g[EPS][hloc+i*layer] * density[i] * ldiff[i];
  printf("%s %12lg\n", cmmt, sum);
}
*/

void CalcFaceTurb(long tinc, int hloc, double *tr)
{
  double tmp, tm[MAXZGRID][3], Kf[MAXZGRID];
  const double sigma[] = {1., 0.9};
  int k, loc;
  Entity et;
  k = ground[hloc].firstabove;
  if (k >= nzm)  return;  /* No Turbulence to calculate */
  for (; k < nz; k++)
    tr[k] = g[TKE][hloc+k*layer];
  for (et = TKE; et <= EPS; et++)  {
    /* Berechnung der K-Turbulenz-Matrix */
    k = ground[hloc].firstabove;
    loc = hloc + k*layer;
/*    Kf[k] = Km[loc] * leveli[k] * density[k] * sigma[et-TKE]; */
    for (; k < nz; k++)  {
      Kf[k] = 0.5 * (Km[loc] + (k < nzm ? Km[loc+layer] : 0)) *
                leveli[k] * density[k] * sigma[et-TKE];
      loc += layer;
    }
    k = ground[hloc].firstabove;
    tmp = tinc / ((0.5 * level[k] + ground[hloc].z) * 0.5 * (density[k-1] + density[k]));
    tm[k][1] = 1. + tmp*Kf[k];
    tm[k][2] = -tmp * Kf[k];
    for (k++; k < nz; k++)  {
      tmp = tinc / (ldiff[k] * 0.5 * (density[k-1] + density[k]));
      tm[k][0] = -tmp * Kf[k-1];
      tm[k][1] = 1. + tmp*(Kf[k-1] + Kf[k]);
      tm[k][2] = -tmp * Kf[k];
    }
    /* Konditioniere Matrix zum schnellen Lösen des Systems */
    for (k = ground[hloc].firstabove+1; k < nz; k++)  {
      tm[k][0] = tmp = -tm[k][0] / tm[k-1][1];
      tm[k][1] += tmp * tm[k-1][2];
    }
    for (k = nzm-1; k >= ground[hloc].firstabove; k--)  {
      tm[k][2] = -tm[k][2] / tm[k+1][1];
    }
    for (k = ground[hloc].firstabove; k < nz; k++)
      tm[k][1] = 1. / tm[k][1];
    /* Rechne jetzt die Turbulenz fuer die entsprechende Turbulenz-
       Groesse mit einer impliziten Integration */
    k = ground[hloc].firstabove;
    loc = hloc + k*layer;
    for (k++; k < nz; k++)  {
      loc += layer;
      g[et][loc] += tm[k][0] * g[et][loc-layer];
     }
    for (k-- ; --k >= ground[hloc].firstabove; )  {
      loc = hloc + k*layer;
      g[et][loc] += tm[k][2] * g[et][loc+layer];
    }
    for (k++ ; k < nz; k++)  {
      loc = hloc + k*layer;
      g[et][loc] *= tm[k][1];
    }
  }
  k = ground[hloc].firstabove;
  /* tr[k] = tr[k+1] = 0.; */ 
  for (; k < nz; k++)  {
    tr[k] = (g[TKE][hloc+k*layer] - tr[k]) / (double)tinc;
  }
}

void CalcKm(int hloc)
{
  int k, loc;
  k = ground[hloc].firstabove;
  loc = hloc + k*layer;
  if (groundinterface)
    Km[loc] = karman * ground[hloc].ustar * ground[hloc].z / ground[hloc].phim;
  else
    Km[loc] = 0.;
  for (k++; k < nz; k++)  {
    loc += layer;
    if (g[EPS][loc] <= 0.)  g[EPS][loc] = 1.e-33;
    Km[loc] = 0.09 * sqr(g[TKE][loc]) / g[EPS][loc];
  }
}

void CalcAllKm(void)
{
  int i, j;
  for (i = nx; --i; )
    for (j = ny; --j; )
      CalcKm(i*row+j);
}

void CalcKEpsilon(long tinc)
/* Berechnet die neuen K-Parameter, sowie neue Werte fuer TKE, EPS
   und THV aufgrund der K-Eps Parametrisierung. */
{
  double uw[MAXZGRID], vw[MAXZGRID], wtheta[MAXZGRID],
         tmp, tmp2, shear, boyancy, transport[MAXZGRID];
  GroundParam *gr;
  int i, j, k, hloc, loc;
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
      hloc = i*row + j;
      gr = ground + hloc;
      k = ground[hloc].firstabove;
      loc = hloc + k*layer;
      /* Evaluate K,lambda */
      if (groundinterface)  {
	g[TKE][loc] = 3.3333 * sqr(gr->ustar);
	g[EPS][loc] = sqr(gr->ustar) * gr->ustar /
	  (gr->z * karman) * (gr->phim - gr->zL);
      }
      CalcFaceTurb(tinc, hloc, transport);
/* Vermutlich keine gute Idee!!      CalcKm(hloc); */
      k = ground[hloc].firstabove;
/*      loc = hloc + k*layer;
      tmp = -Km[loc] / (0.5 * level[k] + gr->z);
      uw[k] = g[UWIND][loc] * tmp;
      vw[k] = g[VWIND][loc] * tmp;
      wtheta[k] = (g[TEMP][loc] - gr->a[TEMP]) * tmp; */
      for (k++ ; k < nz; k++)  {
        loc = hloc + k*layer;
        tmp = -Km[loc] / ldiff[k];
        uw[k] = (g[UWIND][loc] - g[UWIND][loc-layer]) * tmp;
        vw[k] = (g[VWIND][loc] - g[VWIND][loc-layer]) * tmp;
        wtheta[k] = (g[TEMP][loc] - g[TEMP][loc-layer]) * tmp;
      }
      k = ground[hloc].firstabove;
/*      loc = hloc + k*layer;
      shear = (-uw[k] * g[UWIND][loc] - vw[k] * g[VWIND][loc]) / (0.5 * level[k] + gr->z);
      boyancy = Grav * 2.22 / ((g[TEMP][loc] + gr->a[TEMP]) *density[k]) * wtheta[k];
      SolveKEps(0, tinc, shear, boyancy, 0., Km[loc], &g[TKE][loc], &g[EPS][loc]); */
      for (k++; k < nz; k++)  {
        loc = hloc + k*layer;
        shear = (-uw[k] * (g[UWIND][loc] - g[UWIND][loc-layer])
                 -vw[k] * (g[VWIND][loc] - g[VWIND][loc-layer])) / ldiff[k];
        boyancy = Grav * 2.22 / ((g[TEMP][loc] + g[TEMP][loc-layer]) * density[k]) * wtheta[k];
        SolveKEps(0, tinc, shear, boyancy, 0., Km[loc], &g[TKE][loc], &g[EPS][loc]);
      }
      CalcKm(hloc);
    }
}
