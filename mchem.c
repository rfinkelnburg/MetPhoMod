/*
   IMPLEMENTATION MODULE mchem.c
   Dieses Modul berechnet die chemischen Reaktionen
   des Modelles.
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
#include "mchem.h"
#include "matrix.h"
#include "mcparallel.h"

#define CONZ(s)  g[s][gli.loc]

struct  {
  int i, j, k, loc;
}  gli;

Matrix jcb = NULL;
int allfast;
ProductionDesc *prodd = NULL;

typedef enum {REACTION_RATES, REACTION_CONSTANTS}  PrintRateType;

void PrintRates(PrintRateType prt)
{
  int i, j;
  Reaction *r;
  Educt *e;
  Product *p;
  BOOL prpl;
  for (i = nreact, r = reaction; i--; r++)  {
    prpl = FALSE;
    for (j = r->neduct, e = r->educt; j--; e++)  {
      if (prpl)  printf(" + ");
      if (e->prefact > 1)  printf("%i ", e->prefact);
      printf("%s", e->comp->name);
      prpl = TRUE;
    }
    prpl = FALSE;
    printf(" --> ");
    for (j = r->nproduct, p = r->product; j--; p++)  {
      if (prpl)  printf(" + ");
      if (p->prefact != 1.)  printf("%lf ", p->prefact);
      printf("%s", p->comp->name);
      prpl = TRUE;
    }
    switch (prt)  {
      case REACTION_CONSTANTS :
/*         printf("\t\t\t%le\t\t%lf\n", r->K, r->EoR); */
	 printf("\t\t\t%le\n", r->Ke);
         break;
      case REACTION_RATES :
         printf(",\t\t\t%le\n", r->rate);
         break;
    }
  }
}

void PrintConcentrations(void)
{
  int i;
  Substance *s;
  for (i = nsubst, s = subst; i--; s++)  {
    printf("%s\t\t\t%le\n", s->name, CONZ(s->subs));
  }
  printf("\n");
}

double AirConz(int k, int i, int j, struct VarDesc *v)
{
  int loc;
  double pr;
  loc = k*layer + i*row + j;
  pr = pp[k] + press[loc];
  return (pr / (AbsoluteTemp(g[TEMP][loc], pr, g[HUMIDITY][loc]) * kBoltz) * 1.e-6);
}

void CalcRateConstants(double T, double dt, double density, double press, double cozen, double *photorate)
{
  double Ti, secans, conv, cloudfact, Rsum, M, k0t, kinft;
  Reaction *r;
  int i;
  BOOL daytime;
  Ti = -1./T;
  M = press / (T * kBoltz) * 1.e-6;
  conv = M * 1.e-9;
  if (daytime = (cozen > 0.05))
    secans = 1. / cozen;
  for (i = nreact, r = reaction; i--; r++)  {
    switch (r->type)  {
      case THERMAL :
         r->Ke = r->K * exp(r->EoR * Ti) * dt;
         break;
      case THERMAL2 :
         r->Ke = r->K * T * T * exp(r->EoR * Ti) * dt;
         break;
      case TROE :
         k0t = r->k0 * pow(T / 300., r->N);
         kinft = r->kinf * pow(T / 300., r->M);
         r->Ke = r->K * dt * k0t * M / (1. + k0t * M / kinft) *
            pow(0.6, 1. / (1. + spow(log10(k0t * M / kinft), 2)));
         break;
      case TROEQUIL :
         k0t = r->k0 * pow(T / 300., r->N);
         kinft = r->kinf * pow(T / 300., r->M);
         r->Ke = r->K * dt * exp(r->EoR * Ti) *
            k0t * M / (1. + k0t * M / kinft) *
            pow(0.6, 1. / (1. + spow(log10(k0t * M / kinft), 2)));
         break;
      case SPECIAL :
         r->Ke = r->K * dt * SpecialConst(r->specidx, T, M);
         break;
      case PHOTODISS :
         r->Ke = (daytime ? r->K * exp(r->EoR * secans) * dt : 0.);
         break; 
      case PHOTODISS3 :
         r->Ke = (daytime ? r->K * pow(cozen, r->N) * exp(-r->EoR * secans) * dt : 0.);
         break;
      case TWOSTREAM :
	r->Ke = photorate[r->pos] * dt;
    }
    r->Ke *= spow(conv, r->netot-1);
  }
}

void NoteTransport(BOOL transport)
{
  int i;
  for (i = nsubst; i--; )
    subst[i].do_transport = !subst[i].ndecay || transport;
}

void CalcRates(void)
{
  int i, j, k;
  Reaction *r, **rptr;
  Educt *e;
  Product *p;
  Substance *s;
  double h;
  for (i = nreact, r = reaction; i--; r++)  {
    r->rate = r->Ke;
    for (j = r->neduct, e = r->educt; j--; e++)
      r->rate *= spow(CONZ(e->comp->subs), e->prefact);
/*    printf("rate %2i : %13.5le\n", nreact - i, r->rate);  */
  }
/*  printf("\n"); */
  nslow = nfast = 0;
  for (i = nsubst, s = subst;  i--; s++)  {
    s->ht = 0.;
    for (j = s->ndecay, rptr = s->decay; j--; rptr++)  {
      h = (*rptr)->Ke;
      for (k = (*rptr)->neduct, e = (*rptr)->educt; k--; e++)
        h *= spow(CONZ(e->comp->subs), e->prefact - (e->comp == s));
      s->ht += h;
    }
    if (allfast || s->take_as_fast || s->ht > 0.1)  {
      s->fastidx = nfast;
      fast[nfast++] = s;
    }
    else  {
      s->fastidx = -1;
      slow[nslow++] = s;
    }
    if (CONZ(s->subs) > 0.5 || 
        (s->ht < 1. && (s->ht < 1.e-20 || CONZ(s->subs) / s->ht > 1.e-3)))
      s->do_transport = 1;
  }
}

Matrix AdjustJacob(Matrix m, int n)
/* Kreiert eine Matrix oder vergroessert sie falls noetig auf die Groesse n*/
{
  static int size;
  if (m)  {
    if (n > size)  FreeMatrix(m);
    else	   return (m);
  }
  size = n;
  return (AllocateMatrix(n, n));
}

Matrix CalculateJacobian(Matrix jcb)
{
  int i, j, k, ns, sidx[MAXSUBST/2];
  double pr;
  Reaction *r;
  Educt *e;
  Product *p;
  double f[MAXSUBST/2], tmp;
  jcb = AdjustJacob(jcb, nfast);
  UnityMatrix(nfast, jcb);
  for (i = nfast; i--; )
    fast[i]->dcdt = 0.;
  for (i = nreact, r = reaction; i--; r++)  {
    ns = 0;
    pr = r->Ke;
    for (j = r->neduct, e = r->educt; j--; e++)  {
      if (e->comp->fastidx >= 0)  {
        f[ns] = pr * e->prefact * spow(CONZ(e->comp->subs), e->prefact-1);
        sidx[ns] = e->comp->fastidx;
      }
      tmp = spow(CONZ(e->comp->subs), e->prefact);
      pr *= tmp;
      for (k = ns; k--; )
        f[k] *= tmp;
      ns += (e->comp->fastidx >= 0);
      if (ns > MAXSUBST/2)  {
        fprintf(stderr, "Fatal error ns = %i\n", ns);
        exit (3);
      }
    }
    r->rate = pr;
    for (j = r->neduct, e = r->educt; j--; e++)  {
      if (e->comp->fastidx >= 0)  {
        e->comp->dcdt -= (double)e->prefact * pr;
        for (k = ns; k--; )
          jcb[e->comp->fastidx][sidx[k]] += (double)e->prefact * f[k];
      }
    }
    for (j = r->nproduct, p = r->product; j--; p++)  {
      if (p->comp->fastidx >= 0)  {
        p->comp->dcdt += p->prefact * pr;
        for (k = ns; k--; )
          jcb[p->comp->fastidx][sidx[k]] -= p->prefact * f[k];
      }
    }
  }
  return (jcb);
}

void CalcSlowSubst()
{
  int i, j;
  Reaction *r;
  Educt *e;
  Product *p;
  for (i = nsubst; i--; )
    subst[i].dcdt = 0.;
  for (i = nreact, r = reaction; i--; r++)  {
    for (j = r->neduct, e = r->educt; j--; e++)
      if (e->comp->fastidx < 0)
        e->comp->dcdt -= (double)e->prefact * r->rate;
    for (j = r->nproduct, p = r->product; j--; p++)
      if (p->comp->fastidx < 0)
        p->comp->dcdt += p->prefact * r->rate;
  }
  for (i = nslow; i--; )
    CONZ(slow[i]->subs) = slow[i]->oldc + slow[i]->dcdt;
}

void PrintMatrix(int n, Matrix m, double *v)
{
  int i, j;
  for (i = 0; i < n; i++)  {
    for (j = 0; j < n; j++)
      printf("%12.7le ", m[i][j]);
    printf("   %12.7le\n", v[i]);
  }
}

void BoxChemStep(int i, int j, int k, long dt,
		 const Vector *sunpos, double *photorate)
{
  double f[MAXSUBST], vv[MAXSUBST], fact, err, origh2o;
  int indx[MAXSUBST], il, jl, nit;
  ProductionDesc *p;
  BOOL repeatit = FALSE;
  gli.i = i; gli.j = j; gli.k = k; gli.loc = k*layer+i*row+j;
  origh2o = g[HUMIDITY][gli.loc];
  g[HUMIDITY][gli.loc] *= 1.608e9;
  CalcRateConstants(AbsoluteTemp(g[TEMP][gli.loc], pp[k], origh2o),
		    dt, density[k], pp[k], sunpos->z, photorate);
  for (il = nsubst; il--; )  {
    /*	  if (CONZ(subst[il].subs) < 0.)
	  printf("Warning: Substance %s got negative: %le\n",
	  subst[il].name, CONZ(subst[il].subs)); */
    if (CONZ(subst[il].subs) < 1.e-30)  /* can't really handle such small concentrations */
      CONZ(subst[il].subs) = 0.;
  }
  CalcRates();
  for (il = nsubst; il--; )  {
    subst[il].oldc = CONZ(subst[il].subs);
    subst[il].take_as_fast = 0;
  }
  /*	CalcSlowSubst();  Aenderung am 17.9.97! */
  do  {
    nit = 0;
    while (1)  {
      jcb = CalculateJacobian(jcb);
      err = 0.;
      for (il = nfast; il--; )  {
	err += fabs((f[il] = CONZ(fast[il]->subs) - fast[il]->dcdt - fast[il]->oldc) /
		    (CONZ(fast[il]->subs) > 1.e-7 ? CONZ(fast[il]->subs) : 1.));
      }
      if (!repeatit && err < 1.e-5)  break;
      if (++nit > 500)  {
	printf("Chemistry is not converging at point %i/%i/%i\n", i, j, k);
	PrintConcentrations();
	PrintRates(REACTION_RATES);
	break;
      }
      if (LUdecomp(nfast, jcb, indx, &fact, vv))  {
	printf("Matrix is singular (in ChemicalTimeStep) at point %d/%d/%d!\n",
	       i, j, k);
	plausible = FALSE;
	return;
      }
      /*	    if (i == 10 && j == 4 && k == 1)
		    PrintRates(REACTION_RATES);
		    PrintMatrix(nfast, jcb, f);  */
      LUbackSub(nfast, jcb, indx, f);
      repeatit = 0;
      for (il = nfast; il--; )  {
	CONZ(fast[il]->subs) -= f[il];
	if (CONZ(fast[il]->subs) < 0.)  {
	  /*		printf("Warning: Substance %s got negative: %le\n",
			fast[il]->name, CONZ(fast[il]->subs));  */
	  CONZ(fast[il]->subs) = 0.;
	  repeatit = 1;
	}
      }
    }
    repeatit = 0;
    CalcSlowSubst();
    for (il = nfast; il--; )  {
      if (!fast[il]->oldc || 
	  fabs(CONZ(fast[il]->subs) - fast[il]->oldc) < 0.1 * fast[il]->oldc)  {
	fast[il]->take_as_fast = 0;
      }
    }
    for (il = nslow; il--; )  {
      if (CONZ(slow[il]->subs) < 0.)  {
	slow[il]->take_as_fast = repeatit = 1;
      }
    }
    if (repeatit)  {
      printf("repeating chemistry at point %i/%i/%i\n", i, j, k);
      CalcRates();
    }
  }  while (repeatit);
  for (p = prodd; p; p = p->next)
    p->pps[gli.loc] = (CONZ(p->subs) - subst[p->subs - SUBS + 1].oldc) /
      (double)dt;
  g[HUMIDITY][k*layer+i*row+j] = origh2o;
}

void ChemicalTimeStep(long dt, const Vector *sunpos)
{
  int i, j, k;
  int xs, xe;
  ProductionDesc *p;
  BOOL repeatit = FALSE;
  NoteTransport(FALSE);
#ifdef PARALLEL
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
    for (j = ny; --j; )
      for (k = ground[i*row+j].firstabove; k < nz; k++)
	BoxChemStep(i, j, k, dt, sunpos, NULL);
}
