/* IMPLEMENTATION MODULE MCServ
   Diverse Rechenroutinen. */

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sunpos.h"
#include "mcglobal.h"
#include "mcground.h"
#include "mcpress.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

/*
double CalcPileEnergy(int i, int j)
{
  int k;
  GroundParam *gr;
  double E;
  gr = &ground[i][j];
  E = density[-1] * dx * dy * (level[0] + gr->z) * gr->a[TEMP] * Cp * 0.5;
  for (k = nzm; k--; )
    E += density[k] * dx * dy * (level[k+1] - (k ? level[k-1] : 0.)) *
         g[TEMP][k][i][j] * Cp * 0.5;
  return (E);
}
*/

double SatHumidity(double T, double p)
{
  double qsat;
  qsat =  610.78 * exp(17.67*(T - 273.16) / (T - 35.86));
  return (0.622 / (p - 0.378 * qsat) * qsat);
}

double IntegrateVertHydroPress(double p, double T1, double T2, double q, double dz)
/* Integrates hydrostatic Pressure in downwards direction, along the height "dz".
   The virt pot. temperature is assumed to change linearly between the two levels.
   This routine is based on the exact analytical solution of the governing equations,
   with the only assumptions, that air humidity remains constant, and temperature
   is interpolated linearly. */
{
  double f, dt;
  dt = T2 - T1;
  f = R / ((1. + 0.61*q) * stdPKappa);
  p = pow(p, Kappa);
  if (fabs(dt) > 0.5)  {
    f *= dt;
    p = (log(T2/T1)*Grav*Kappa*dz+f*p) / f;
  }
  else  {
    f *= 0.5 * (T1 + T2);
    p = p + dz * Kappa * Grav / f;
  }
  p = pow(p, 1./Kappa);
  return (p);
}

void CalcMeanHydrostaticPressure(void)
{
  int i, j, k, etemp, ehumid, loc, hloc;
  double hydp;
  /* Berechne hydrostatische Druckniveaus. Starte oben. */
  pp[nzm] = toppp;
  etemp = TEMP*nz;
  ehumid = HUMIDITY*nz;
  /* density[nz] wird in Continuity und in Advect gebraucht. */
  density[nz] = density[nzm] =
     pp[nzm] / (R * AbsoluteTemp(avg[etemp+nzm], pp[nzm], avg[ehumid+nzm]));
  for (k = nzm; k--; )  {
    loc = etemp+k; hloc = ehumid+k;
    pp[k] = IntegrateVertHydroPress(pp[k+1], avg[loc], avg[loc+1],
    				    0.5*(avg[ehumid]+avg[ehumid+1]), ldiff[k+1]);
    density[k] = pp[k] / (R * AbsoluteTemp(avg[loc], pp[k], avg[hloc]));
  }
  pp[-1] = IntegrateVertHydroPress(pp[0], avg[loc], avg[loc], avg[ehumid], ldiff[0]);
  density[-1] = pp[-1] / (R * AbsoluteTemp(avg[loc], pp[-1], avg[hloc]));
/* Version abgeleitet von Beniston...
  for (k = nzm; k--; )  {
    loc = etemp+k;
    hydp = pow(pp[k+1], Kappa) + Cpress * (level[k+1] + level[k]) /
	   (avg[loc+1] + avg[loc]);
    pp[k] = pow(hydp, 1. / Kappa);
    density[k] = pp[k] / (R * avg[loc] * pow(pp[k] / P0, Kappa));
  }
  hydp = pow(pp[0], Kappa) + Cpress * level[0] / avg[etemp];
  pp[-1] = pow(hydp, 1. / Kappa);
  density[-1] = pp[-1] / (R * avg[etemp] * pow(pp[k] / P0, Kappa));
*/
  if (pp[0] >= 101400.)  {
    printf("middle pressure is too big at ground. (%lf)\n", pp[0]);
    plausible = FALSE;
  }
  if (pp[0] < 60000.)  {
    printf("middle pressure is too low at ground. (%lf)\n", pp[0]);
    plausible = FALSE;
  }
}

double CalcTopEnergy(void)
{
  int i, j;
  double sum = 0., *l;
#ifdef PARALLEL
  if (parallel)  {
    i = (mfirstx + !mfirstx) * row;
    l = flux[WWIND] + nzm*layer + (mfirstx + !mfirstx) * row;
    for (i = mnx-leftest-rightest; i--; )
      for (j = ny+1; j--; )
	sum += sqr(l[i*row + j]);
  }
  else
#endif
    for (i = nzm*layer; i < mesh; i++)
      sum += sqr(flux[WWIND][i]);
  return (sum);
}

void CalculateLayerAverage(double *l, PointStatus *stat, double *bar)
{
  int i, j, k;
  double sum;
  if (!l)  {
    for (k = nz; k--; )  *bar++ = 0.;
    return;
  }
#ifdef PARALLEL
  if (parallel)  {
    if (master)  {
      GetAveragesFromWorkers(bar);
      for (k = nz; k--; )
        if (nlevelpt[k])  bar[k] /= nlevelpt[k];
        else              bar[k] /= ((nxm-1)*(nym-1));
      SendAveragesToWorkers(bar);
    }
    else  {
      i = (mfirstx + !mfirstx) * row;
      l += i; stat += i;
      for (k = 0; k < nz; k++)  {
        bar[k] = 0.;
        for (i = mnx-leftest-rightest; i--; )
          for (j = ny; --j; )
            if (!stat[i*row+j])
                bar[k] += l[i*row+j];
        l += layer; stat += layer;
      }
      ChangeAveragesWithMaster(bar);
    }
  }
  else  {
#endif
    for (k = 0; k < nz; k++)  {
      sum = 0.;
      for (i = nx; --i; )
        for (j = ny; --j; )
          if (!stat[i*row+j])  
            sum += l[i*row+j];
      if (nlevelpt[k])  *bar++ = sum / nlevelpt[k];
      else              *bar++ = sum / ((nxm-1)*(nym-1));
      l += layer; stat += layer;
    }
#ifdef PARALLEL
  }
#endif
}

#ifdef PARALLEL

void ExchangeWalls(double *x)
{
  if (westbordertype == NEIGHBOUR)  {
    SendPressureWall(x, 1 + 2*OVERLAP, leftbrother, TOTHELEFT);
  }
  if (eastbordertype == NEIGHBOUR)  {
    GetPressureWall(x, nx, rightbrother, TOTHELEFT);
  }
  if (eastbordertype == NEIGHBOUR)  {
    SendPressureWall(x, nxm - 2*OVERLAP, rightbrother, TOTHERIGHT);
  }
  if (westbordertype == NEIGHBOUR)  {
    GetPressureWall(x, 0, leftbrother, TOTHERIGHT);
  }
}

#endif

void CalcFluxDivergency(int tinc, double *udiv)
/* This routine calculates the divergencies of the fluxes, in the form:

     integral(rho * U * dS) / dt

   where U and dS are considered as vectors. CalcFluxDivergency is narrowly related to
   "Convergency" (in mcdynamics). It serves to calculate non-hydrostatic pressure fields.
*/
{
  int i, j, k, loc;
  double dxdy;
  for (i = nx; --i; )
    for (j = ny; --j; )
      for (k = ground[i*row+j].firstabove; k < nz; k++)  {
        loc = i*row + j + k*layer;
        udiv[loc] = -( ( (flux[UWIND][loc] - flux[UWIND][loc-row]) * dy +
                         (flux[VWIND][loc] - flux[VWIND][loc-1]) * dx
                       ) * density[k] * level[k]
                       + ( flux[WWIND][loc] *
                           (flux[WWIND][loc] > 0. ? density[k] : density[k+1]) -
                           (k ? flux[WWIND][loc-layer] * 
                             (flux[WWIND][loc-layer] > 0. ? density[k-1] : density[k]) 
                              : 0.)
                         ) * dx * dy
                     ) / tinc;
      }
  /* Add contribution of top layer */
  dxdy = dx * dy / ldiff[nzm];
  for (i = nx; --i; )
    for (j = ny; --j; )  {
      loc = i*row + j;
      udiv[loc+nzm*layer] += toppress[loc] * dxdy;
    }
#ifdef PARALLEL
  if (parallel)
    ExchangeWalls(udiv);
#endif
}

void PressureMatMult(const double *x, double *y)
/* Multiplies the pressure vector "x" with the pressure matrix. This
   routine calculates matrix pressure multiplication, without actually
   storing the pressure matrix. */
{
  int i, j, k, loc, locxp, locxm, xs, xe;
  double dxdy, dxdyp, dxdz, dydz;
#ifdef PARALLEL
  if (parallel)  {
    xs = leftest; xe = nx - rightest;
  }
  else  {
#endif
    xs = 1 - (westbordertype == CYCLIC); 
    xe = nxm + (eastbordertype == CYCLIC);
#ifdef PARALLEL
  }
#endif
  for (k = nz; k--; )  {
    dxdz = level[k] * dx * dyi;
    dydz = level[k] * dy * dxi;
    dxdy = dx * dy / ldiff[k]; dxdyp = dx * dy / ldiff[k < nzm ? k+1 : k];
    for (i = nx; --i; )  {
      locxp = locxm = row;
#ifndef PARALLEL
      if (eastbordertype == CYCLIC && i == nxm)  locxp = -(nxm-1)*row;
      if (westbordertype == CYCLIC && i == 1)    locxm = -(nxm-1)*row;
#endif
      for (j = ny; --j; )  {
        loc = i*row + j + k*layer;
        if (!pstat[loc])  {
          y[loc] = dydz * ((i > xs && !pstat[loc-locxm]) * (x[loc] - x[loc-locxm]) +
                           (i < xe && !pstat[loc+locxp]) * (x[loc] - x[loc+locxp])) +
                   dxdz * ((j > 1 && !pstat[loc-1]) * (x[loc] - x[loc-1]) +
                           (j < nym && !pstat[loc+1]) * (x[loc] - x[loc+1]));
          if (k && !pstat[loc-layer])  y[loc] += dxdy * (x[loc] - x[loc-layer]);
          if (k < nzm)		y[loc] += dxdyp * (x[loc] - x[loc+layer]);
          else			y[loc] += dxdy * x[loc];
        }
      }
    }
  }
#ifdef PARALLEL
  if (parallel)
    ExchangeWalls(y);
#endif
}

double PressVectMult(double *a, double *b)
/* Multiplies two 3D meshes with vector-multiplication */
{
  int i, j, k, xs, xe, hloc, loc;
  double sum = 0;
#ifdef PARALLEL
  if (parallel)  {
    xs = mfirstx + !mfirstx;
    xe = mlastx - rightest;
    for (i = xs; i < xe; i++)
      for (j = ny; --j; )  {
        hloc = i*row + j;
        for (k = ground[hloc].firstabove; k < nz; k++)  {
          loc = hloc + k*layer;
          sum += a[loc] * b[loc];
        }
      }
    ChangeVectMult(&sum);
  }
  else
#endif
    for (i = mesh; i--; a++, b++)  {
      sum += *a * *b;
    }
  return (sum);
}

void PressVectAdd(double af, double *a, double bf, double *b)
/* Adds vector "b" multiplied by "f" to "a": b = b + f*a */
{
  int i;
  if (af == 1.)
    for (i = mesh; i--; a++, b++)
      *a = *a + bf * *b;
  else if (bf == 1.)
    for (i = mesh; i--; a++, b++)
      *a = af * *a + *b;
  else if (bf == -1.)
    for (i = mesh; i--; a++, b++)
      *a = af * *a - *b;
  else
    for (i = mesh; i--; a++, b++)
      *a = af * *a + bf * *b;
}

void SetPressureBoundary(void)
/* Sets the pressure at the boundary */
{
  int i, j, k, wloc, eloc, nloc, sloc, vloc;
  wloc = (westbordertype == CYCLIC ? nxm : 1) * row;
  eloc = (eastbordertype == CYCLIC ? 1 : nxm) * row;
  sloc = (southbordertype == CYCLIC ? nym : 1);
  nloc = (northbordertype == CYCLIC ? 1 : nym);
  for (k = nz; k--; )  {
    vloc = k*layer;
    for (i = nx; --i; )  {
      press[i*row+vloc] = press[i*row+sloc+vloc];
      press[i*row+ny+vloc] = press[i*row+nloc+vloc];
    }
#ifdef PARALLEL
    if (parallel)  {
      if (leftest)
	for (j = ny; j--; )  {
	  press[j+vloc] = press[j+wloc+vloc];
	}
      if (rightest)
	for (j = ny; j--; )  {
	  press[j+row*nx] = press[j+eloc+vloc];
	}
    }
    else
#endif
      for (j = ny+1; j--; )  {
	press[j+vloc] = press[j+wloc+vloc];
	press[j+vloc+row*nx] = press[j+eloc+vloc];
      }
  }
}

void SetPressBoundaryZero(double *m)
{
  int i, j, k, loc;
  for (k = nz; k--; )  {
    loc = k * layer;
    for (i = nx+1; i--; )  m[i*row+loc] = m[i*row+ny+loc] = 0.;
#ifdef PARALLEL
    if (parallel)  {
      if (leftest)  for (j = ny+1; j--; )  m[j+loc] = 0.;
      if (rightest) for (j = ny+1; j--; )  m[j+nx*row+loc] = 0.;
    }
    else  {
#endif
      if (westbordertype != CYCLIC)   for (j = ny+1; j--; )  m[j+loc] = 0.;
      if (eastbordertype != CYCLIC)  for (j = ny+1; j--; )  m[j+nx*row+loc] = 0.;
#ifdef PARALLEL
    }
#endif
  }
  for (i = mesh; i--; )
    if (pstat[i]) m[i] = 0.;
}

void TestWorse(double *a, double *b)
{
  int i;
  for (i = mesh; i--; )
    if (fabs(b[i]) / (fabs(a[i])+1.) > 3.)
      printf("worse at point %d/%d/%d, (a = %lg, b = %lg)\n",
         (i % layer) / row, i % row, i / layer, a[i], b[i]);
}


void TestEqual(double *x)
{
  int i;
  for (i = xrow*nz; i--; x+= row)
    if (x[1] != x[2])  {
      i = (nx+1)*nz - i;
      printf("Unequal points found at %d/%d: %lg vs. %lg\n",
        (i % xrow), i / xrow, x[1], x[2]);
      return;
    }
}

void SolveVerticalPressure(const double *b, double *x)
/* Solves for the pressure, using only the vertical interdepence of points, with
   a straightforward linear inversion routine. The original vector b is left
   unchanged */
{
  int i, j, k, fa, loc, hloc, xe;
  double dxdy;
  dxdy = 1. / (dx * dy);
#ifdef PARALLEL
  xe = (parallel ? nx - rightest + 1 : nx);
  for (i = leftest; i < xe; i++)
#else
  for (i = nx; --i; )
#endif
    for (j = ny; --j; )  {
      hloc = i*row + j;
      fa = ground[hloc].firstabove;
      loc = hloc+fa*layer;
      x[loc] = b[loc];
      for (k = fa+1; k < nz; k++)  {
	loc += layer;
	x[loc] = b[loc] + x[loc-layer];
      }
      x[loc] = x[loc] * ldiff[nzm] * dxdy;
      for (k = nzm; k > fa; k--)  {
	loc -= layer;
	x[loc] = x[loc] * dxdy * ldiff[k] + x[loc+layer];
      }
/*      for (k = fa; k < nz; k++)  {
        loc = hloc + k*layer;
        x[loc] = b[loc] * 2. * ldiff[k] * dxdy;
      } */
    }
}

void SolveForPressure(int tinc)
/* This routine applies the conjugated gradient method, to solve for the non-hydrostatic
   pressure. The conjugated gradient routine was implemented according to N.R. p. 84. */
{
  const double eps = 1.;
  double alpha, beta, rz, rpzp;
  static double *r = NULL, *p = NULL, *z = NULL, *swap;
  int it;
#ifdef PARALLEL
  if (parallel && master)  {
    ChangeVectMult(&alpha);
  }
  else  {
#endif
    if (!r)  {
      if (!(r = (double *)calloc(mesh, sizeof(double))) ||
          !(p = (double *)calloc(mesh, sizeof(double))) ||
          !(z = (double *)calloc(mesh, sizeof(double))))  {
	fprintf(stderr, "Unable to allocate %ld bytes in \"SolveForPressure\"..exiting.\n",
	   3 * mesh * sizeof(double));
	exit (0);
      }
    }
    SetPressBoundaryZero(press);
    CalcFluxDivergency(tinc, r);     /* r1 = b - A*x */
    PressureMatMult(press, z);
    PressVectAdd(1., r, -1., z);
    SolveVerticalPressure(r, z);      /* A'*z1 = r1 */
    memcpy(p, z, mesh * sizeof(double));   /* p1 = z1 */
    rz = PressVectMult(r, z);
    if (fabs(rz) < eps)  {
      it = 0;
    }
    else  {
      for (it = 1; it < mesh; it++)  {
	PressureMatMult(p, z);		/* alpha = (rk*zk) / (pk * A * pk) */
	alpha = rz / PressVectMult(p, z);
	PressVectAdd(1., r, -alpha, z);	/* rkp = rk - alpha * A * pk */
	PressVectAdd(1., press, alpha, p);  /* xkp = x + alpha * pk */
	SolveVerticalPressure(r, z);
	beta = (rpzp = PressVectMult(r, z)) / rz;  /* beta = (rkp * zkp) / (rk * zk) */
	if (printresid
#ifdef PARALLEL
	    && !parallel
#endif
	    )
	  printf("it %3d: alpha = %lg,  beta = %lg, r = %lg\n",
		 it, alpha, beta, (rpzp >= 0. ? sqrt(rpzp) : -sqrt(-rpzp)));
	if (fabs(rpzp) < eps)  break;
	PressVectAdd(beta, p, 1., z);
	rz = rpzp;
      }
    }
    SetPressureBoundary();
#ifdef PARALLEL
    if (parallel)  {
      alpha = -9999.;
      ChangeVectMult(&alpha);
    }
  }
  if (!parallel)
#endif
    printf("cj-gr it: %i  ", it);
}

double NextPress(double t, double t0, double d, double p, double p0, double dz)
{
  double hf, pn;
  t -= t0;
  hf = 0.5 * GCvbCp * d * dz / p0;
  pn = -(p - Grav*d*dz*t/t0 + hf*p) / (hf - 1.);
  return (pn);
}

void CalcHydrostaticPressure(void)
{
  int i, j, k, ik, jk, loc;
  double densmid, tempmid, pressmid;
  k = nzm;
  for (i = nx+1; i--; )
    for (j = ny+1; j--; )  {
      loc = k*layer+i*row+j;
      press[loc] = NextPress(g[TEMP][loc], avg[TEMP*nz+k], 
         density[k], toppress[i*row+j], pp[k], 100.);
      press[loc-layer] = NextPress((g[TEMP][loc] + g[TEMP][loc-layer]) * 0.5,
      				   (avg[TEMP*nz+k] + avg[TEMP*nz+k-1]) * 0.5,
      				   (density[k] + density[k-1]) * 0.5,
      				   press[loc],
      				   (pp[k] + pp[k-1]) * 0.5,
      				   ldiff[k]);
    }
  k--;
  while (--k >= 0)  {
    tempmid = (avg[TEMP*nz+k+1] + avg[TEMP*nz+k]) * 0.5;
    densmid = (density[k+1] + density[k]) * 0.5;
    pressmid = (pp[k+1] + pp[k]) * 0.5;
    for (i = nx+1; i--; )  {
      for (j = ny+1; j--; )  {
        loc = k*layer+i*row+j;
        ik = i; jk = j;
        if (!i && (westbordertype & DATABORDER))  ik = 1;
        if (i == nx && (eastbordertype & DATABORDER))  ik = nxm;
        if (!j && (southbordertype & DATABORDER))  jk = 1;
        if (j == ny && (northbordertype & DATABORDER))  jk = nym;
        if (!pstat[loc])  {
          press[loc] = NextPress((g[TEMP][(k+1)*layer+ik*row+jk] + g[TEMP][k*layer+ik*row+jk]) * 0.5,
            	    tempmid,
            	    densmid,
            	    press[loc+layer],
            	    pressmid,
            	    ldiff[k+1]);
	}
        else if (!pstat[loc+layer])
            press[loc] = press[loc] = NextPress((g[TEMP][(k+1)*layer+ik*row+jk] +
                    g[TEMP][k*layer+ik*row+jk]) * 0.5,
            	    tempmid,
            	    densmid,
            	    press[loc+layer],
            	    pressmid,
            	    0.5 * level[k]);
      }
    }
  }
}
