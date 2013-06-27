/*
   IMPLEMENTATION MODULE mcclouds.c
   Implementiert Wolkenphysik
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sunpos.h"
#include "mcglobal.h"
#include "mcground.h"
#include "mcpress.h"


double SumClouds(int hloc, Entity et)
{
  double sum = 0.;
  int k;
  for (k = ground[hloc].firstabove; k < nz; k++)
    sum += g[et][k*layer + hloc] * density[k] * level[k];
  return (1000. * sum);
}

double TotalWaterClouds(int k, int i, int j, struct VarDesc *v)
{
  int hloc = i*row + j;
  return (SumClouds(hloc, CLOUDWATER) + SumClouds(hloc, RAINWATER));
}

double TotalIceClouds(int k, int i, int j, struct VarDesc *v)
{
  int hloc = i*row + j;
  return (SumClouds(hloc, CLOUDICE));
}

double TotalClouds(int k, int i, int j, struct VarDesc *v)
{
  return (TotalWaterClouds(k, i, j, v) + TotalIceClouds(k, i, j, v));
}

void Condensate(int loc, int k, double qsw, double qsi, double *T)
{
  const double T00 = 238.16, TD = 1. / 35.;
  double q, qs, qc, qi, C, D;
  q = g[HUMIDITY][loc];
  if (*T > T0 || !cloudice)  {
    C = 0.008 * tinc * (q - qsw) / (1. + 1.348e7 * qsw / (*T * *T));
    if (C + g[CLOUDWATER][loc] < 0.)  C = -g[CLOUDWATER][loc];
    g[CLOUDWATER][loc] += C;
    g[HUMIDITY][loc] -= C;
    *T += C * Lc / Cp;
  }
  else if (*T < T00)  {
    D = 0.008 * tinc * (q - qsi) / (1. + 1.727e7 * qsi / (*T * *T));
    if (D + g[CLOUDICE][loc] < 0.)  D = -g[CLOUDICE][loc];
    g[CLOUDICE][loc] += D;
    g[HUMIDITY][loc] -= D;
    *T += D * Ls / Cp;
  }
  else  {
    qc = g[CLOUDWATER][loc];
    qi = g[CLOUDICE][loc];
    qs = (qc || qi ? (qc * qsw + qi * qsi) / (qc + qi) : 0.5 * (qsw + qsi));
    C = 0.008 * tinc * (q - qs) * (*T - T00) * TD / (1. + 1.348e7 * qs / (*T * *T));
    D = 0.008 * tinc * (q - qs) * (T0 - *T) * TD / (1. + 1.727e7 * qs / (*T * *T));
    if (C + g[CLOUDWATER][loc] < 0.)  C = -g[CLOUDWATER][loc];
    if (D + g[CLOUDICE][loc] < 0.)  D = -g[CLOUDICE][loc];
    g[CLOUDWATER][loc] += C;
    g[CLOUDICE][loc] += D;
    g[HUMIDITY][loc] -= C + D;
    *T += C * Lc / Cp + D * Ls / Cp;
  }
}

void AutoConvert(int loc, int k)
{
  const double qc0 = 0.0005, alpha = 0.001;
  double C = 0., q, lambda, vr;
  if (cloudwater)  {
    q = g[CLOUDWATER][loc];
    if (q > qc0)
      C = alpha * (q - qc0);  /* Autoconversion */
    if (g[RAINWATER][loc] > 1.e-5)  {
      lambda = pow(2.5133e10 / (density[k] * g[RAINWATER][loc]), 0.25);
      C += 2.4834e10 * g[CLOUDWATER][loc] * pow(lambda, -3.8); /* Accretion */
    }
    g[RAINWATER][loc] += tinc * C;
    g[CLOUDWATER][loc] -= tinc * C;
  }
}

double TerminalVelocityOfRain(int k, int i, int j, struct VarDesc *v)
{
  int loc;
  double lambda;
  loc = k * layer + i * row + j;
  if (g[RAINWATER][loc] > 1.e-5)  {
    lambda = pow(2.5133e10 / (density[k] * g[RAINWATER][loc]), 0.25);
    return (2770. / (pow(lambda, 0.8) * sqrt(density[k])));
  }
  return (0);
}

void RainFall(int loc)
{
  int i, k, fa;
  double s, dist, qr[MAXZGRID], qro, lambda, grw;
  fa = ground[loc].firstabove;
  memset(qr, 0, nz * sizeof(double));
  for (k = nzm; k >= fa; k--)  {
    if ((qro = g[RAINWATER][loc + k*layer]) > 1.e-5)  {
      lambda = pow(2.5133e10 / (density[k] * g[RAINWATER][loc+k*layer]), 0.25);
      s = 2770. / (pow(lambda, 0.8) * sqrt(density[k])) * tinc;
      grw = 0.;
      for (dist = level[i = k]; i > fa && s >= dist; dist += level[--i]);
      if (s < dist)  {
        qr[i] += (dist - s) / level[i] * qro;
        s += level[k];
        if (s < dist)  qr[i] -= (dist - s) / level[i] * qro;
        else  {
          for (i--; i >= fa && (dist = dist + level[i]) < s; i--)
            qr[i] += qro;
          if (i >= fa)  qr[i] += (s - dist + level[i]) / level[i] * qro;
          else
            grw = (s - dist) * qro * density[k];
        }
      }
      else
        grw = level[k] * qro * density[k];
      g[RAINWATER][loc + k*layer] = qr[k];
      ground[loc].rain_intensity = grw / tinc;
      ground[loc].cum_rain += grw;
    }
  }
}

void SnowFall(int loc, double *vt)
{
  int i, k, fa;
  double s, dist, qi[MAXZGRID], qio, grw;
  fa = ground[loc].firstabove;
  memset(qi, 0, nz * sizeof(double));
  for (k = nzm; k >= fa; k--)  {
    if ((qio = g[CLOUDICE][loc + k*layer]) > 1.e-6)  {
      s = vt[k] * tinc;
      grw = 0.;
      for (dist = level[i = k]; i > fa && s >= dist; dist += level[--i]);
      if (s < dist)  {
        qi[i] += (dist - s) / level[i] * qio;
        s += level[k];
        if (s < dist)  qi[i] -= (dist - s) / level[i] * qio;
        else  {
          for (i--; i >= fa && (dist = dist + level[i]) < s; i--)
            qi[i] += qio;
          if (i >= fa)  qi[i] += (s - dist + level[i]) / level[i] * qio;
          else
            grw = (s - dist) * qio * density[k];
        }
      }
      else
        grw = level[k] * qio * density[k];
      g[CLOUDICE][loc + k*layer] = qi[k];
      ground[loc].snow_intensity = grw / tinc;
      ground[loc].cum_snow += grw;
    }
  }
}

void Evaporate(int loc, int k, double *T, double qsw)
{
  double E, lambda;
  if (g[HUMIDITY][loc] < qsw && g[RAINWATER][loc] > 1.e-5)  {
    lambda = pow(2.5133e10 / (density[k] * g[RAINWATER][loc]), 0.25);
    E = tinc * 5.02655e7 * (1. - g[HUMIDITY][loc] / qsw) *
        (0.78/sqr(lambda) + 4067.05 * pow(lambda, -2.9)) /
        (density[k] * 5.35288e11 / (*T * *T) +
         1. / (density[k] * qsw * 1.875e-5));
    if (g[RAINWATER][loc] < E)  E = g[RAINWATER][loc];
    g[RAINWATER][loc] -= E;
    g[HUMIDITY][loc] += E;
    *T -= E * Lc / Cp;
  }
}

void RimingMelting(int loc, int k, double *T, double *vt)
{
  double Ni, mi, sqmi, Di, vti, rhoi, k1, k2, R;
  if (cloudice && g[CLOUDICE][loc] > 1.e-6)  {
    Ni = 0.01 * exp(0.6 * (T0 - *T));
    sqmi = sqrt(mi = density[k] * g[CLOUDICE][loc] / Ni);
    if (mi < 1.e-8)  {
      if (mi < 1.7e-10)  {
        k1 = 16.3;
        k2 = 304.;
      }
      else  {
        k1 = 6.07;
        k2 = 1250.;
      }
      Di = k1 * sqmi;
      vti = k2 * Di * sqrt(1.225 / density[k]);
      rhoi = 3.39531 / (k1*k1*k1 * sqmi);
    }
    else  {
      Di = 1.58 * pow(mi, 0.417);
      vti = 4.48 * pow(Di, 0.25) * sqrt(1.225 / density[k]);
      rhoi = 0.860811 * pow(mi, -0.251);
    }
    vt[k] = vti;
    if (*T <= T0)  {   /* Riming occurs */
      R = 0.628319 * vti * Ni * g[CLOUDWATER][loc] * Di * Di * tinc;
      g[CLOUDICE][loc] += R;
      g[CLOUDWATER][loc] -= R;
      *T += R * Li / Cp;
    }
    else  {   /* Melting occurs */
      R = (*T - T0) * g[CLOUDICE][loc] * tinc * 0.04;
      if (R > g[CLOUDICE][loc])  R = g[CLOUDICE][loc];
      g[CLOUDICE][loc] -= R;
      g[RAINWATER][loc] += R;
      *T -= R * Li / Cp;
    }
  }
}

void CloudPhysics(void)
{
  int i, j, k, loc, hloc;
  /* Berechne qsat und lasse ueberschuessige Fluessigkeit kondensieren. */
  double T, Tv, qsw, qsi, vt[MAXZGRID];
  const double eps = 1.e-6;
  for (i = nx; i--; )
    for (j = ny; j--; )  {
      for (k = ground[hloc = i*row+j].firstabove; k < nz; k++)  {
        loc = k*layer + hloc;
        T = (Tv = g[TEMP][loc] * pow(pp[k] * 1.e-5, Kappa)) / (1 + 0.61 * g[HUMIDITY][loc]);
        qsw = 611.2 * exp(17.67*(T - T0) / (T - 29.66));
        qsw = qsw * 0.622 / (pp[k] - 0.378 * qsw);
        qsi = 611.2 * exp(21.9*(T - T0) / (T - 7.7));
        qsi = qsi * 0.622 / (pp[k] - 0.378 * qsi);
        Condensate(loc, k, qsw, qsi, &T);
        AutoConvert(loc, k);
        Evaporate(loc, k, &T, qsw);
        RimingMelting(loc, k, &T, vt);
        g[TEMP][loc] = VirtPot(T, g[HUMIDITY][loc], pp[k]);
      }
      RainFall(hloc);
      if (cloudice)  SnowFall(hloc, vt);
    }
}


