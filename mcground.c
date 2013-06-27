/*
   IMPLEMENTATION MODULE mcground
   In diesem Modul werden Bodenbilanzen gerechnet.
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mcglobal.h"
#include "mcpress.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcfloat.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif
#include "mctwostream/mctwostream.h"
#ifdef SVR4
#include "mc_group.hh"
#include "mc_group.t"  // Solaris scheint einen Fehler in der Behandlung

template class Group<PhotoDissReact>;
template class GroupEnumerator<PhotoDissReact>;

#endif                 // von Templates zu haben und braucht einen Workaround

static double visc = 1.5e-5,
	      grc = 0.0962;

typedef enum {DZL, DTG, DTF, MAXDERIV}  DerivType;
typedef double Gradient[MAXDERIV];


double St;           /* Korrigierte Sonnenstrahlung oberhalb
                        der Atmosphaere. */
double grlevel[NGROUNDLAYER] = {0., 0.04, 0.16, 0.36, 0.64, 1.},
/* Weite der Bodenlevel unter Grund  */
       dgrlev[NGROUNDLAYER] = {50., 8.333333, 5., 3.57142857, 2.7777778};
/* 1 / Maechtigkeit der Bodenlevel */
GroundParam *ground;
double Rd;

double VirtPot(double T, double q, double p)
{
  return (T * (1. + 0.61*q) * pow(1.e5 / p, Kappa));
}

void CalcSlope(Vector *v, int i, int j)
{
  double vxz, vyz;
  double r;
  vxz = topo[(i < nxm ? i+1 : i)*row+j] - topo[(i > 1 ? i-1 : i)*row+j];
  vyz = topo[i*row+(j < nym ? j+1 : j)] - topo[i*row+(j > 1 ? j-1 : j)];
  v->x = -vxz / dx; v->y = -vyz / dy; v->z = 1.;
  r = 1. / sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
  v->x *= r; v->y *= r; v->z *= r;
}

BOOL PrecipitableLength(double *P, int i, int j, Vector *vsun, double CosZa)
/* Berechnet Precipitable water in cm */
{
  BOOL shadow, inrange;
  int k, ip, jp, loc, hloc;
  double posx, posy, gposx, gposy, height;
  gposx = (double)i * dx; gposy = (double)j * dy;
  *P = 0.;
  shadow = 0;
  hloc = i*row+j;
  for (k = ground[hloc].firstabove, height = 0.5 * level[k]; k < nz; k++)  {
    posx = gposx + height * vsun->x;
    posy = gposy + height * vsun->y;
    ip = (int)(posx / dx + 0.5);
    jp = (int)(posy / dy + 0.5);
    loc = ip*row+jp;
    inrange = (ip > 0 && ip < nx && jp > 0 && jp < ny);
    shadow |= (inrange && height + topo[hloc] < topo[loc]);
    if (shadow || !inrange)
      *P += density[k] * avg[HUMIDITY*nz+k] * level[k] / CosZa;
    else
      *P += density[k] * g[HUMIDITY][k*layer+loc] * level[k] / CosZa;
    if (k < nzm)  height += ldiff[k+1];
  }
  *P = *P /10. + topprecipwater;
  return (shadow);
}

double EpsH2O(double dP)      /* Nach Pielke 8-35  */
{
  double logP;
  dP *= 0.1;
  logP = log10(dP+1.e-12);
  if (logP <= -4.)  return (0.113 * log10(1. + 12.6*dP));
  if (logP <= -3.)  return (0.104 * logP + 0.44);
  if (logP <= -1.5) return (0.121 * logP + 0.491);
  if (logP <= -1.)  return (0.146 * logP + 0.527);
  if (logP <= 0.)   return (0.161 * logP + 0.542);
  return (0.136 * logP + 0.542);
}

double EpsCO2(double pp)
{
  double Hc;
  Hc = 0.252*pp;
  if (Hc <= 0.)  Hc = 0.;
  return (0.185 * (1. - exp(-0.39 * pow(Hc, 0.4))));
}

double AbsoluteTemp(double theta, double press, double humid)
{
  return (theta * pow(press / 1.e5, Kappa) / (1. + 0.61 * humid));
}

double AbsolutePointTemp(int k, int i, int j, struct VarDesc *v)
{
  return AbsoluteTemp(g[TEMP][k*layer+i*row+j], pp[k], g[HUMIDITY][k*layer+i*row+j]);
}

double EquivalentPointTemp(int k, int i, int j, struct VarDesc *v)
{
  int loc;
  double quot;
  loc = k*layer + i * row + j;
  quot = Lc * g[HUMIDITY][loc] / (Cp * g[TEMP][loc] * pow(pp[k] / 1.e5, Kappa));
  return (g[TEMP][loc] * exp(quot));
}

double RelativeHumidity(int k, int i, int j, struct VarDesc *v)
{
  double qsat,T;
  int loc;
  loc = k*layer + i * row + j;
  T = AbsoluteTemp(g[TEMP][loc], pp[k], g[HUMIDITY][loc]);
  qsat =  611.2 * exp(17.67*(T - T0) / (T - 29.66));
  qsat = qsat * 0.622 / (pp[k] - 0.378 * qsat);
  return (100. * g[HUMIDITY][loc] / qsat);
}

double NewLongWaveRadiation(int i, int j, double tinc)
{
  int h, k, fa, hloc, loc;
  double T4g, T4f, T4m, Tt, tmp, tmp2,
	 u[MAXZGRID], dP[MAXZGRID], dPc[MAXZGRID], dH[MAXZGRID],
	 Tabs[MAXZGRID], dTdu[MAXZGRID], Rup[MAXZGRID], Rdown[MAXZGRID];
  /* Paremeters at ground */
  hloc = i*row + j;
  fa = ground[hloc].firstabove;
  T4g = sqr(sqr(ground[hloc].Tg[0]));
  T4f = sqr(sqr(ground[hloc].Tf));
  T4m = ground[hloc].sigf * T4f + (1. - ground[hloc].sigf) * T4g;
  Rup[fa] = Tabs[fa] = sigma * T4m;
  Tt = AbsoluteTemp(g[TEMP][nzm*layer+hloc], pp[nzm], avg[HUMIDITY*nz+nzm]);
  Rdown[nz] = sigma * sqr(sqr(Tt));
  u[fa] = dP[fa] = 0.;
  for (k = fa; k < nz; k++)  {
    loc = hloc + k*layer;
    tmp = density[k] * level[k];
    u[k+1] = u[k] + tmp;
    dP[k+1] = dP[k] + tmp * g[HUMIDITY][loc];
    dPc[k] = dP[k] + 0.5 * tmp * g[HUMIDITY][loc];
  }
  for (k = fa; k < nz; k++)
    dH[k] = 0.5 * (density[k-1] + density[k]);
  dH[nz] = density[k-1];
  loc = hloc + fa*layer;
  tmp = AbsoluteTemp(g[TEMP][loc], pp[fa], g[HUMIDITY][loc]);
  for (k = fa+1; k < nz; k++)  {
    loc = hloc + k*layer;
    Tabs[k] = sigma *
      sqr(sqr(0.5 * (tmp + (tmp2 = AbsoluteTemp(g[TEMP][loc], pp[k], g[HUMIDITY][loc])))));
    tmp = tmp2;
  }
  Tabs[nz] = sigma * sqr(sqr(tmp + 0.5 * toptemp[hloc]));
  for (k = fa; k < nz; k++)  {
    dTdu[k] = (Tabs[k+1] - Tabs[k]);
  }
  for (k = fa+1; k <= nz; k++)  {
    tmp = Rup[fa];
    for (h = fa; h < k; h++)
      tmp += (EpsH2O(dP[k] - dPc[h]) + EpsCO2(pp[h] - dH[k])) * dTdu[h];
    Rup[k] = tmp;
  }
  for (k = fa; k < nz; k++)  {
    tmp = Rdown[nz];
    for (h = k; h < nz; h++)
      tmp -= (EpsH2O(dPc[h] - dP[k]) + EpsCO2(dH[k] - pp[h])) * dTdu[h];
    Rdown[k] = tmp;
  }
  for (k = fa; k < nz; k++)  {
    loc = hloc + k*layer;
    irdiv[loc] = (Rup[k] - Rup[k+1] + Rdown[k+1] - Rdown[k]) /
      (density[k] * level[k] * Cp) *
      (1 + 0.61 * g[HUMIDITY][loc]) * pow(1.e5 / pp[k], Kappa);
  }
  return (Rdown[fa]);
}

void ApplyIRRadiation(int hloc, double tinc)
{
  int k, loc;
  for (k = ground[hloc].firstabove; k < nz; k++)  {
    loc = hloc + k * layer;
    g[TEMP][loc] += tinc * irdiv[loc];
  }
}

double LongWaveRadiation(int i, int j, double tinc)
{
  int k, fa, hloc, loc;
  double drdz, lower, higher, Rdown, Tt, T4g, T4t, T4f, T4m, eps, epslast, dP,
         Tabs[MAXZGRID];
  hloc = i*row+j;
  Tt = AbsoluteTemp(g[TEMP][nzm*layer+hloc], pp[nzm], avg[HUMIDITY*nz+nzm]);
  /* Aufwaertsgerichtete Strahlung */ 
  T4g = sqr(sqr(ground[hloc].Tg[0]));
  T4f = sqr(sqr(ground[hloc].Tf));
  T4m = ground[hloc].sigf * T4f + (1. - ground[hloc].sigf) * T4g;
  epslast = dP = 0.;
  fa = ground[hloc].firstabove;
  for (k = fa; k < nz; k++)  {
    loc = k*layer+hloc;
    Tabs[k] = sqr(sqr(AbsoluteTemp(g[TEMP][loc], pp[k], g[HUMIDITY][loc])));
  }
  for (k = fa; k < nz; k++)  {
    loc = k*layer+hloc;
    if (g[HUMIDITY][loc] < 0.)  {
      fprintf(stderr, "Humidity is negative at point %d/%d/%d.\n", i, j, k);
      plausible = FALSE;
      g[HUMIDITY][loc] = 0.;
    }
    dP += density[k] * g[HUMIDITY][loc] * level[k];
    eps = EpsH2O(dP) + EpsCO2(pp[k]);
    drdz = sigma * (eps - epslast) *
           (Tabs[k] - T4m) * leveli[k];
    g[TEMP][loc] -= tinc * drdz / (density[k] * Cp);
    epslast = eps;
  }
/*
   Funktioniert leider nicht!!!...trotzdem noch ein Versuch! */
  T4t = Tabs[nzm];
  dP = topprecipwater;
  epslast = eps = EpsH2O(dP) + EpsCO2(pp[nzm]);
  Rdown = sigma * T4t;
  for (k = nzm-1; k >= fa; k--)  {
    loc = k*layer+hloc;
    dP += 0.5 * (density[k+1] * g[HUMIDITY][loc+layer] * level[k+1] +
                 density[k] * g[HUMIDITY][loc] * level[k]);
    eps = EpsH2O(dP) + EpsCO2(pp[k]);
    drdz = sigma * (eps - epslast) *
           (T4t - Tabs[k]);
    g[TEMP][loc] += tinc * drdz * leveli[k]/ (density[k] * Cp);
    Rdown -= drdz;
    epslast = eps;
  }
/*
Stattdessen:
  Rdown = sqr(sqr(IRdown))*sigma;
*/
  return (Rdown);
}


void ShortWaveRadiation(int i, int j, int k, Vector vsun, GroundParam *gr)
{
  BOOL shadow;
  double asunslope, dP, m, awv, aoz, tt, CosZa, I, scatt, Dr, Dri, f, x, tmp;
  double dTg, gsig;
  Vector vslope;
  St = S0 * Rd;
  CosZa = vsun.z;
  k = gr->firstabove;
  gsig = 1. - gr->sigf;
  if (vsun.z > 0.01)  {     /* Tag: Es existiert kurzwellige Strahlung */
    /* Routinen zur Berechnung der direkten Einstrahlung */
    asunslope = AngleBetween(&vsun, &gr->slope);
    shadow = PrecipitableLength(&dP, i, j, &vsun, CosZa) || asunslope <= 0.;
    m = 1.003447 / sqrt(CosZa*CosZa + 0.00689);
    x = m * stratozone;
    aoz = 0.02118 * x / (1. + 0.042*x + 0.000323*x*x) +
	  1.082*x / pow(1. + 138.6*x, 0.805) + 0.0658*x / (1. + pow(103.6*x, 3.));
    m = 1. / CosZa;
    tt = 0.01936 + 0.583583/(m + 6.325633);
    awv = 2.9 * dP / (pow(1. + 141.5*dP, 0.635) + 5.925*dP);
    scatt = exp((-tt * pp[gr->firstabove] / 1.e5 - turbidity) * m);
    I = St * ((1. - aoz) * scatt - awv);
    if (I < 0.)  I = 0.;
    /* Berechnung der diffusen Strahlung */
    Dri = CosZa * I * ((0.4 + 0.6 * exp(-turbidity) - 0.28 / (1. + 6.43*CosZa)) / scatt - 1.);  /* Scatter von oben */
#ifdef RADFACT
    Dr = Dri + (Dri + Rdirfact[k*layer+i*row+j] * I * CosZa) * 0.0685 *
#else
    Dr = Dri + (Dri + Rdirfact[k] * I * CosZa) * 0.0685 *
#endif
      (gsig * gr->albedo - gr->sigf * gr->albedof);    /* Scatter von unten */
    f = sqr(cos(acos(gr->slope.z) / 2.));      /* Anteil des sichtbaren Himmels; Pielke p.406 */
#ifdef RADFACT
    gr->Rdiff = gr->Rg = Rdifffact[k*layer+i*row+j] * (Dr * f + (I + Dr) * (1. - f) * (1. - gr->absorb));
    /* Einfluss der schiefen Oberflaeche */
    if (!shadow)  {
      gr->Rg += (gr->Rdir = Rdirfact[k*layer+i*row+j] * I * asunslope);
#else
    gr->Rdiff = gr->Rg = Rdifffact[k] * (Dr * f + (I + Dr) * (1. - f) * (1. - gr->absorb));
    /* Einfluss der schiefen Oberflaeche */
    if (!shadow)  {
      gr->Rg += (gr->Rdir = Rdirfact[k] * I * asunslope);
#endif
      /* Direkte Strahlung */
    }
  }
  else  {
    gr->Rg = 0.;
  }
}

void CalcWTheta(int k, GroundParam *gr, double Vm,
		double *wthetag, double *wthetaf, double *newzL,
		double *dwthetag, double *dwthetaf,
		double *wqg, double *wqf,
		double *dwqg, double *dwqf, double *dnewzL)
{
  double ustar, zL, gsig, cd, cf, uaf, ra, rs, f, oldwtheta, Taf,
         oldzL, phim, phit, psit, psim, qsat, qf, qaf, chih, wq, tmp, tmp2, tmp3,
         thetastar, qstar, qg, qsf, dphim, dphit, dpsim, dpsit, dustar, grtheta,
         dgrtheta, dchih, dcd, duaf, dcf, logzz0, wtheta, dqgdg, dqsfdf, dra, df,
         h2, lpress;
  Gradient dthetastar, dqstar, dwtheta, dqsat, dqaf, dqf;
  int nit, loc;
  /* Die Variablennamen richten sich nach den entsprechenden Formeln in Pielke 84 */
  if (Vm < 0.1)  Vm = 0.1;
  loc = (gr - ground) + k * layer;
  logzz0 = log(0.5 * level[k] / gr->z0);
  if (calcsoiltemp) {
    grtheta = VirtPot(gr->Tg[0], g[HUMIDITY][loc], tmp = 0.5*(pp[k-1]+pp[k]));
    dgrtheta = (1. + 0.61*g[HUMIDITY][loc]) * pow(1.e5 / tmp, Kappa);
  }
  else {
    grtheta = gr->Tg[0];
    dgrtheta = 1.;
  }
  wtheta = gr->wtheta;
give_it_a_new_try :
  zL = gr->zL;
  gsig = 1. - gr->sigf;
  lpress = 0.5 * (pp[k] + pp[k-1]);
  if (zL <= 0.)  {
    phim = pow(1. - 15.*zL, -0.25);   /* Pielke: 7-51, 7-52 */
    dphim = 3.75 * pow(1. - 15.*zL, -1.25);
    phit = 0.74 / sqrt(1. - 9.*zL);
    dphit = 4.5 * phit / (1. - 9.*zL);
    psim = 2. * log((1.+1./phim)/ 2.) + log((1.+1./(phim*phim))/2.) -
	   2. * atan(1./phim) + 1.5707963;
    dpsim = -4. / (phim * (phim + 1.) * (phim*phim +1.)) * dphim;
    psit = 2. * log((1.+0.74/phit)/ 2.);
    dpsit = -0.74 / (sqr(phit) * ((1. + 0.74/phit)/2.)) * dphit;
  }
  else  {
    phim = 1. + 4.7*zL;
    dphim = dphit = 4.7;
    phit = 0.74 + 4.7*zL;
    psim = -4.7 * zL;
    dpsim = -4.7;
    psit = -6.35 * zL;
    dpsit = -6.35;
  }
  gr->ustar = ustar = karman * Vm / (logzz0 - psim);
  tmp = 0.74 * (logzz0 - psit);
  if (tmp < 0.01)  {
    gr->zL = 0.;
    goto give_it_a_new_try;
  }
  gr->ra = tmp / (karman * gr->ustar);
  gr->phim = phim; gr->phit = phit;
  dustar = ustar / (logzz0 - psim) * dpsim;
  chih = karman / tmp;
  dchih = chih / tmp * dpsit * 0.74;
  if (ustar <= 0.)  {
/*    printf("inplausible ustar : %le  (%i/%i/%i) zL = %lf\n", ustar, i, j, k, zL); */
    gr->zL /= 10.;
    goto give_it_a_new_try;
  }
  if (calcsoiltemp) {
    tmp3 = pow(ustar * gr->z0 / visc, 0.45);
    tmp = chih /         /* Pielke: 7-49, 7-50 */
          (1. + chih * grc / karman * tmp3);
    thetastar = tmp * (g[TEMP][loc] - grtheta);       /* Pielke: 7-49, 7-50 */
    tmp2 = dchih / (1. + chih * grc / karman * tmp3) -
           chih / sqr(1. + chih * grc / karman * tmp3) *
           grc / karman * (dchih * tmp3 +
                           0.45 * chih * tmp3 * dustar / ustar);
    dthetastar[DZL] = tmp2 * (g[TEMP][loc] - grtheta);
    dthetastar[DTG] =  - tmp * dgrtheta;
  }
  else {
    thetastar = (g[TEMP][loc] - grtheta) / tmp;
    dthetastar[DZL] = thetastar / tmp * dpsit * 0.74;
    dthetastar[DTG] = -dgrtheta / tmp;
  }
  dthetastar[DTF] = 0.;
  *wthetag = -gsig * ustar * thetastar;
  dwthetag[DZL] = -gsig * (dustar * thetastar + ustar * dthetastar[DZL]);
  dwthetag[DTG] = -gsig * ustar * dthetastar[DTG];
  dwthetag[DTF] = 0.;
  if (gr->sigf)  {
    cd = sqr(ustar / Vm);
    dcd = 2. * dustar * ustar / (Vm * Vm);
    uaf = 0.83 * ustar;
    duaf = 0.83 * dustar;
    cf = 0.01 * (1. + 0.3 / uaf);
    dcf = -0.003 / sqr(uaf) * duaf;
    Taf = 0.3 * AbsoluteTemp(g[TEMP][loc], lpress, g[HUMIDITY][loc]) +
	  0.6 * gr->Tf + 0.1 * gr->Tg[0];
    tmp = pow(1.e5 / lpress, Kappa);
    *wthetaf = -gr->sigf * 1.1 * gr->La * cf * uaf * tmp * (Taf - gr->Tf);
    dwthetaf[DZL] = -gr->sigf * 1.1 * gr->La * tmp *
                    (Taf - gr->Tf)*(cf * duaf + dcf * uaf);
    tmp2 = gr->sigf * gr->La * tmp * cf * uaf;
    dwthetaf[DTG] = -0.11 * tmp2;
    dwthetaf[DTF] = 0.44 * tmp2;
    tmp2 = gr->sigf * cd * uaf * tmp;
    *wthetag -= tmp2 * (Taf - gr->Tg[0]);
    dwthetag[DZL] -= (Taf - gr->Tg[0]) * tmp * gr->sigf * (cd * duaf + dcd * uaf);
    dwthetag[DTG] += tmp2 * 0.9;
    dwthetag[DTF] -= gr->sigf * tmp * cd * uaf * 0.6;
    wtheta = *wthetag + *wthetaf;
    dwtheta[DZL] = dwthetag[DZL] + dwthetaf[DZL];
    dwtheta[DTG] = dwthetag[DTG] + dwthetaf[DTG];
    dwtheta[DTF] = dwthetag[DTF] + dwthetaf[DTF];
  }
  else  {
    wtheta = *wthetag;
    dwtheta[DZL] = dwthetag[DZL];
    dwtheta[DTG] = dwthetag[DTG];
    *wthetaf = dwtheta[DTF] = 0.;
  }
  nit++;
  /* Pielke: 11-37 */
  /* Formeln zur Berechnung des latenten Waermeflusses */
  qsat =  611.2 * exp(17.67*(gr->Tg[0] - T0) / (gr->Tg[0] - 29.66));  /* Stull p.276 */
  dqsat[DTG] = qsat * 17.67 * (1. / (gr->Tg[0] - 29.66) - (gr->Tg[0] - T0) / sqr(gr->Tg[0] - 29.66));
  tmp = gr->X * 0.622 / (lpress - 0.378 * qsat);
  qg = tmp * qsat;
  dqgdg = (tmp + 0.378 * qg / (lpress - 0.378 * qsat)) * dqsat[DTG];
  tmp = 1. + chih * grc / karman * tmp3;
  qstar = chih *         /* Pielke: 7-49, 7-50 */
	  (g[HUMIDITY][loc] - qg) / tmp;
  dqstar[DZL] = (g[HUMIDITY][loc] - qg) * (dchih / tmp - chih / (tmp*tmp) *
                (grc / karman * (dchih * tmp3 +
                 chih * 0.45 * tmp3 / ustar * dustar)));
  dqstar[DTG] = -dqgdg * chih / tmp;
  dqstar[DTF] = 0.;
  *wqg = -gsig * qstar * ustar;
  dwqg[DZL] = -gsig * (dqstar[DZL] * ustar + qstar * dustar);
  dwqg[DTG] = -gsig * dqstar[DTG] * ustar;
  dwqg[DTF] = 0.;
  if (gr->sigf)  {
    qsat =  611.2 * exp(17.67*(gr->Tf - T0) / (gr->Tf - 29.66));
    dqsat[DTF] = qsat * 17.67 * (1. / (gr->Tf - 29.66) - (gr->Tf - T0) / sqr(gr->Tf - 29.66));
    tmp = 0.622 / (lpress - 0.378 * qsat);    
    qsf = tmp * qsat;
    dqsfdf = (tmp + 0.378 * qsf / (lpress - 0.378 * qsat)) * dqsat[DTF];
    if (gr->wq >= 0.)  {
      /*gr->ra = */ ra = 1. / (cf * uaf);    /* Pielke 11-61a */
      dra = -ra / (cf * uaf) * (dcf * uaf + cf * duaf);
      rs = gr->rcf * 900. / (27. + gr->Rg);
      f = 1. - rs / (rs + ra) * gr->fwatercontent;
      df = rs / sqr(rs + ra) * gr->fwatercontent * dra;
    }
    else  {
      f = 1.;
      df = 0.;
    }
    qf = (10.*f*qsf + (3.*g[HUMIDITY][loc] + qg)*(1. - f)) / (4. + 6.*f);
      /* Pielke: 11-60, 11-61 */
    dqf[DZL] = (40.*df*qsf - 30.*df*g[HUMIDITY][loc] - 10.*df*qg) /
               sqr(4. + 6.*f);
    dqf[DTG] = (1. - f) / (4. + 6.*f) * dqgdg;
    dqf[DTF] = 10.*f*dqsfdf / (4. + 6.*f);
    if (qf > qsf)  {qf = qsf; dqf[DZL] = dqf[DTG] = 0.; dqf[DTF] = dqsfdf;}
    qaf = 0.3 * g[HUMIDITY][loc] + 0.6 * qf + 0.1 * qg;
    dqaf[DZL] = 0.6 * dqf[DZL];
    dqaf[DTG] = 0.1 * dqgdg + 0.6 * dqf[DTG];
    dqaf[DTF] = 0.6 * dqf[DTF];
    tmp = gr->sigf * gr->La * f;
    *wqf = -tmp * cf * uaf * (qaf - qf);
    tmp2 = gr->sigf * gr->La * df;
    dwqf[DZL] = -cf * uaf * (tmp2 * (qaf - qf) + tmp * (dqaf[DZL] - dqf[DZL])) -
                 tmp * (qaf - qf) * (dcf * uaf + cf * duaf);
    dwqf[DTG] = -cf * uaf * tmp * (dqaf[DTG] - dqf[DTG]);
    dwqf[DTF] = -cf * uaf * tmp * (dqaf[DTF] - dqf[DTF]);
    *wqg -= gr->sigf * cd * uaf * (qaf - qg);
    dwqg[DZL] -= gr->sigf * (cd * uaf * dqaf[DZL] + (dcd*uaf + cd*duaf)*(qaf - qg));
    dwqg[DTG] -= gr->sigf * cd * uaf * (dqaf[DTG] - dqgdg);
    dwqg[DTF] -= gr->sigf * cd * uaf * dqaf[DTF];
  }
  else  {
    *wqf = dwqf[DZL] = dwqf[DTG] = dwqf[DTF] = 0.;
  }
  tmp = -2. * sqr(ustar / Vm) / level[k];
  gr->uw = gr->vw = 1. / (1. + ustar * ustar / (0.5 * level[k] * Vm) * tinc);
  if (k)
    h2 = (level[k-1] * avg[TEMP*nz+k] + level[k] * avg[TEMP*nz+k-1]) / (level[k] + level[k-1]);
  else
    h2 = avg[TEMP*nz+k];
  tmp = - 0.5 * level[k] * Grav / (h2 * sqr(ustar) * ustar);
  *newzL =  tmp * wtheta;
  dnewzL[DZL] = tmp * (dwtheta[DZL] - 3. / ustar * dustar * wtheta);
  dnewzL[DTG] = tmp * dwtheta[DTG];
  dnewzL[DTF] = tmp * dwtheta[DTF];
}

double GroundHeatFlux(int i, int j, GroundParam *gr, double tinc)
{
  int k;
  double Qg, dTg[NGROUNDLAYER], a[NGROUNDLAYER][3], f, tmp;
/* Neue Version: implizite Integration */
  tmp = 2. * gr->ks * tinc;
  for (k = 1; k < NGROUNDLAYER-1; k++)  {
    f = tmp / (grlevel[k+1] - grlevel[k-1]);
    a[k][0] = -f * dgrlev[k-1];
    a[k][1] = 1. + f * (dgrlev[k-1] + dgrlev[k]);
    a[k][2] = -f * dgrlev[k];
  }
  gr->Tg[1] -= a[1][0] * gr->Tg[0];
  for (k = 2; k < NGROUNDLAYER-1; k++)  {
    f = -a[k][0] / a[k-1][1];
    a[k][1] += f * a[k-1][2];
    gr->Tg[k] += f * gr->Tg[k-1];
  }
  for (k = NGROUNDLAYER-1; --k; )  {
    gr->Tg[k] = (gr->Tg[k] - a[k][2] * gr->Tg[k+1]) / a[k][1];
  }
/* Alte Version: explizie Integration
  dTg[0] = gr->ks * 2. * (gr->Tg[1] - gr->Tg[0]) / sqr(grlevel[0] - grlevel[1]);
  Qg = 0.;
  for (k = NGROUNDLAYER-1; --k; )
    dTg[k] = gr->ks * 2. * (grlevel[k] * (gr->Tg[k+1] - gr->Tg[k-1]) +
                            grlevel[k+1] * (gr->Tg[k-1] - gr->Tg[k]) +
                            grlevel[k-1] * (gr->Tg[k] - gr->Tg[k+1]))/
             ((grlevel[k+1] - grlevel[k])*(grlevel[k] - grlevel[k-1])*
              (grlevel[k+1] - grlevel[k-1]));
  for (k = NGROUNDLAYER-1; --k; )
    gr->Tg[k] += tinc * dTg[k];
*/
  Qg = gr->Cg * gr->ks * (gr->Tg[0] - gr->Tg[1]) / grlevel[1];
  return (Qg);
}

void TurnVect(Vector *v, double angle)
{
  double x, y, c, s;
  x = v->x; y = v->y;
  c = cos(angle); s = sin(angle);
  v->x = x * c + y * s;
  v->y = -x * s + y * c;
}

double CalcBoundaryLayerHeight(GroundParam *gr, int i, int j)
{
  int k, hloc, loc;
  double reftemp, blh;
  hloc = i*row+j;
  reftemp = g[TEMP][hloc+gr->firstabove*layer];
  for (k = gr->firstabove+1; k < nz && reftemp > g[TEMP][k*layer+hloc]; k++);
  loc = k*layer + hloc;
  if (k < nz)  {
    if (k > gr->firstabove+1)  {
      blh = zcenter[k] + (reftemp - g[TEMP][loc]) /
            (g[TEMP][loc-layer] - g[TEMP][loc]) * (zcenter[k-1] - zcenter[k]);
    }
    else  blh = zcenter[k-1];
  }
  else  blh = zcenter[nzm];
  if (gr->firstabove > 0)
    blh -= zfaces[gr->firstabove-1];
  return (blh);
}

double CalcMeanWind(GroundParam *gr)
{
  double Vm, wstar, grtheta, tmp;
  int k, loc;
  k = gr->firstabove;
  loc = (gr - ground) + k * layer;
  grtheta = VirtPot(gr->Tg[0], g[HUMIDITY][loc], 0.5*(pp[k-1]+pp[k])) - g[TEMP][loc];
  Vm = sqr(g[UWIND][loc]) + sqr(g[VWIND][loc]);
  tmp = grtheta * 0.01;
  if (grtheta > 0. && gr->wtheta > 0. && gr->blh > 0.)  {
    if (gr->wtheta < tmp)  gr->wtheta = tmp;
    wstar = pow(Grav / g[TEMP][loc] * gr->wtheta * gr->blh, 1./3.);
    Vm += wstar*wstar;
  }
  return (sqrt(Vm));
}

void GroundInterface(double timeofday, int dayofyear, double tinc, Vector *sunpos)
{
  BOOL verystable, novegetation;
  int i, j, k, nit;
  double K, L, R, Qg, Fg, Ff, Leg, Lef, dFgdg, dFgdf, dFfdg,
    dFfdf, emissityratio, T4g, T4f, qg, qsf, qsat,
    T3g, T3f, dTg, dTf, wthetag, wthetaf, wqg, wqf, tmp, tmp2, det,
    gsig, dGdg, zL, deltazL, dens, h, deviat, bestdeviat, Vm;
  Gradient dwthetag, dwthetaf, dwqg, dwqf, dzL, dFg, dFf;
  Vector vsun;
  GroundParam *gr, bestgr;
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
  dayofyear += (int)(timeofday / 24);
  dayofyear = (dayofyear-1) % 365 + 1;
  timeofday = fmod(timeofday, 24.);
  Rd = RatioDistMeanDist(dayofyear);
  SunVector(&vsun, dayofyear, timeofday, Xlong, Xlat, timezonediff);
  TurnVect(&vsun, turnmapangle);
  *sunpos = vsun;
  if (shortwaveradiation)
    if (radiationtype)
      twostreammodule.CalcTwoStream(vsun, Rd);
    else {
      for (i = nx; --i; )
	for (j = ny; --j; )  {
	  gr = ground+i*row+j;
	  ShortWaveRadiation(i, j, gr->firstabove, vsun, gr);
	}
    }
  if (groundinterface)
    for (i = xs; i < xe; i++)
      for (j = ny; --j; )  {
	gr = ground+i*row+j;
	gr->blh = CalcBoundaryLayerHeight(gr, i, j);
	novegetation = gr->sigf == 0.;
	k = gr->firstabove;
	dens = 0.5 * (density[k] + density[k-1]);
	Vm = CalcMeanWind(gr);
	K = (shortwaveradiation ? gr->Rg * gr->absorb : 0);
	if (timeofday > 23.5)  gr->Rgcum = 0.;
	else		     gr->Rgcum += K * tinc;
	if (calcsoiltemp) {
	  L = (irtinc ? (actime % irtinc == 0 ?
			 gr->IRdown = NewLongWaveRadiation(i, j, tinc) :
			 gr->IRdown) : 0.);
	  ApplyIRRadiation(i*row+j, tinc);
	  Qg = GroundHeatFlux(i, j, gr, tinc);
	  dGdg = gr->ks * gr->Cg / grlevel[1];
	  gsig = 1. - gr->sigf;
	  if (!novegetation)
	    emissityratio = sigma * gr->emissityf * gr->emissity /
	      (gr->emissity + gr->emissityf - gr->emissityf * gr->emissity);
	  else  emissityratio = 0.;
	}
	else {
	  gr->Tf = gr->Tg[0];
	  Fg = Ff = 0.;
	}
	nit = 0; bestdeviat = 1.e33;
	bestgr = *gr;
	verystable = FALSE;
	do  {
	  CalcWTheta(k, gr, Vm,
		     &wthetag, &wthetaf, &zL,
		     dwthetag, dwthetaf, &wqg, &wqf, dwqg, dwqf, dzL);
	  gr->wtheta = wthetag + wthetaf;
	  gr->wq = wqg + wqf;
	  if (calcsoiltemp) {
	    T4g = sqr(sqr(gr->Tg[0])); T4f = sqr(sqr(gr->Tf));
	    tmp = gr->sigf * emissityratio * (T4g - T4f);
	    Leg = gsig * sigma * gr->emissity * T4g + tmp;
	    Lef = gr->sigf * sigma * gr->emissityf * T4f - tmp;
	    T3g = 4. * sqr(gr->Tg[0]) * gr->Tg[0];
	    T3f = 4. * sqr(gr->Tf) * gr->Tf;
	    Fg = gsig * (K + L * gr->emissity) - Leg - Qg - dens* (Cp * wthetag + Lc * wqg);
	    dFg[DZL] = -dens * (Cp * dwthetag[DZL] + Lc * dwqg[DZL]); 
	    dFg[DTG] = -(gsig * sigma * gr->emissity + gr->sigf * emissityratio) * T3g -
	      dens * (Cp * dwthetag[DTG] + Lc * dwqg[DTG]) - dGdg;
	    dFg[DTF] = gr->sigf * emissityratio * T3f - dens * (Cp * dwthetag[DTF] + Lc * dwqg[DTF]);
	    if (!novegetation)  {
	      Ff = gr->sigf * (K + L * gr->emissityf) - Lef - dens* (Cp * wthetaf + Lc * wqf);
	      dFf[DZL] = -dens * (Cp * dwthetaf[DZL] + Lc * dwqf[DZL]);
	      dFf[DTG] = gr->sigf * emissityratio * T3g - dens * (Cp * dwthetaf[DTG] + Lc * dwqf[DTG]);
	      dFf[DTF] = -(gr->sigf * sigma * gr->emissityf + gr->sigf * emissityratio) * T3f -
		dens * (Cp * dwthetaf[DTF] + Lc * dwqf[DTF]);
	    }
	    else  Ff = 0.;
	    deviat = fabs(Fg) + fabs(Ff);
	  }
	  else {
	    deviat = fabs(zL - gr->zL);
	  }
	  zL -= gr->zL;
	  dzL[DZL] -= 1.;
	  if (deviat < bestdeviat)  {
	    bestdeviat = deviat;
	    bestgr = *gr;
	  }
	  if (deviat < 0.01)  break;
	  if (calcsoiltemp) {
	    verystable |= dzL[DZL] > 0. || gr->zL > 30.;
	    if (verystable)  {
	      if (novegetation)  {
		dTg = -Fg / dFg[DTG];
		dTf = deltazL = 0.;
	      }
	      else  {
		det = dFg[DTG] * dFf[DTF] - dFg[DTF] * dFf[DTG];
		dTg = (dFg[DTF] * Ff - dFf[DTF] * Fg) / det;
		dTf = (dFf[DTG] * Fg - dFg[DTG] * Ff) / det;
		deltazL = 0.;
	      }
	      gr->zL = 30.;
	    }
	    else  {
	      if (novegetation)  {
		det = dFg[DTG] * dzL[DZL] - dFg[DZL] * dzL[DTG];
		dTg = (dFg[DZL] * zL - dzL[DZL] * Fg) / det;
		deltazL = (dzL[DTG] * Fg - dFg[DTG] * zL) / det;
		dTf = 0.;
	      }
	      else  {
		det = dFg[DTG]*dFf[DTF]*dzL[DZL] + dFg[DTF]*dFf[DZL]*dzL[DTG] + dFg[DZL]*dFf[DTG]*dzL[DTF] -
		  dFg[DZL]*dFf[DTF]*dzL[DTG] - dFg[DTG]*dFf[DZL]*dzL[DTF] - dFg[DTF]*dFf[DTG]*dzL[DZL];
		dTg = -(Fg*dFf[DTF]*dzL[DZL] + dFg[DTF]*dFf[DZL]*zL + dFg[DZL]*Ff*dzL[DTF] -
			dFg[DZL]*dFf[DTF]*zL - Fg*dFf[DZL]*dzL[DTF] - dFg[DTF]*Ff*dzL[DZL]) / det;
		dTf = -(dFg[DTG]*Ff*dzL[DZL] + Fg*dFf[DZL]*dzL[DTG] + dFg[DZL]*dFf[DTG]*zL -
			dFg[DZL]*Ff*dzL[DTG] - dFg[DTG]*dFf[DZL]*zL - Fg*dFf[DTG]*dzL[DZL]) / det;
		deltazL = -(dFg[DTG]*dFf[DTF]*zL + dFg[DTF]*Ff*dzL[DTG] + Fg*dFf[DTG]*dzL[DTF] -
			    Fg*dFf[DTF]*dzL[DTG] - dFg[DTG]*Ff*dzL[DTF] - dFg[DTF]*dFf[DTG]*zL) / det;
	      }
	    }
	    Qg += dTg * dGdg;
	    gr->Tg[0] += dTg;
	    gr->Tf += dTf;
	    if (gr->Tg[0] < 230.)  {
	      Qg += (230. - gr->Tg[0]) * dGdg;
	      gr->Tg[0] = 230.;
 	    }
	    if (gr->Tg[0] > 340.)  {
	      Qg += (340. - gr->Tg[0]) * dGdg;
	      gr->Tg[0] = 340.;
	    }
	    if (gr->Tf < 250.)  gr->Tf = 250.;
	    if (gr->Tf > 340.)  gr->Tf = 340.;
	  }
	  else {     /* of "if (calcsoiltemp)" */
	    deltazL = -zL / dzL[DZL];
	  }
	  gr->zL += deltazL;
	}  while (++nit < 150);
	if (nit >= 150)  {
	  /*        printf("Problems in finding correct Tg[0]/Tf at %i/%i\n", i, j);  */
	  *gr = bestgr;
	}
	if (gr->wtheta > 0.4)  gr->wtheta = 0.4;
	if (gr->wtheta < -0.2)  gr->wtheta = -0.2;
	if (gr->wq > 2.4e-4)  gr->wq = 2.4e-4;
	if (gr->wq < -6.e-5)  gr->wq = -6.e-5;
	gr->Qg = Qg;
	gr->R = K + L * (gsig * gr->emissity + gr->sigf * gr->emissityf) - Leg - Lef;
      }
}
