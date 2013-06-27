/*
   IMPLEMENTATION MODULE mcdep
   In diesem Modul wird die Deposition von Schadstoffen gerechnet.
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mcdep.h"
#ifdef PARALLEL
#include <pvm3.h>
#include "mcparallel.h"
#endif

typedef struct DepositionDesc {
  VarDesc v;
  Entity et;
  double diffusivity, depofact, depowater;
  double *deposition;
  struct DepositionDesc *next;
}  DepositionDesc;

DepositionDesc *deposit = NULL;

double *PlaceDepositionVariable(VarDesc *v)
{
  static char depunit[] = "mol/(sec*m^2)";
  DepositionDesc *d;
  VarDesc *nv, *vi;
  for (d = deposit; d && strcmp(d->v.name+2, v->name); d = d->next);
  if (!d)  {
    if (!(d = (DepositionDesc *)malloc(sizeof(DepositionDesc))) ||
        !(d->deposition = (double *)calloc(layer, sizeof(double))))  {
      fprintf(stderr, "ERROR: Error allocating memory in PlaceDepositionVariable. Exiting!\n");
      exit (1);
    }
    nv = &d->v;
    nv->name = (char *)malloc(2 * strlen(v->name) + 16);
    sprintf((char *)nv->name, "D-%s\0Deposition of %s", (char *)v->name, (char *)v->name);
    nv->comment = nv->name + strlen(v->name) + 1;
    nv->unit = depunit;
    nv->storetype = DOUBLE_PTR;
    nv->v.d = d->deposition;
    nv->dims = (X_DIM | Y_DIM);
    nv->init = CALCULATED_VALUES;
    nv->inputtype = NORMAL_NUM;
    for (vi = v; vi->next; vi = vi->next);
    nv->next = NULL;
    vi->next = nv;
    d->et = v->v.et;
    d->next = deposit;
    deposit = d;
  }
  return (&d->diffusivity);
}

VarDesc *DepositionVariable(VarDesc *v, VarDesc *vh)
{
  *vh = *v;
  if (actualsection->id == DEPOSITION && v->storetype == GRID_VAL)  {
    vh->storetype = DOUBLE_PTR;
    vh->dims = 3;
    vh->v.d = PlaceDepositionVariable(v);
    vh->section = DEPOSITION;
    vh->inputtype = VECTOR_VAL;
  }
  return (vh);
}

#ifdef PARALLEL

void SendDepositionData(void)
{
  DepositionDesc *d;
  if (deposit)
    printf("  deposition data\n");
  for (d = deposit; d; d = d->next)  {
    pvm_pkstr((char *)d->v.name+2);
    pvm_pkdouble(&d->diffusivity, 3, 1);
  }
  pvm_pkstr("");
}

void GetDepositionData(void)
{
  char name[80];
  pvm_upkstr(name);
  while (*name)  {
    PlaceDepositionVariable(GetNamedVar(name));
    pvm_upkdouble(&deposit->diffusivity, 3, 1);
    pvm_upkstr(name);
  }
}
  
#endif

void InitDeposition(void)
{
  DepositionDesc *dd;
  for (dd = deposit; dd; dd = dd->next)
    dd->diffusivity = 2. * pow(dd->diffusivity, 0.66666666666) / karman;
}

void Deposit(long tinc)
{
  int i, j, fa, loc;
  double rtot, rs, rc, T, dep;
  DepositionDesc *dd;
  GroundParam *gr;
  if (deposit)  {
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
	fa = gr->firstabove;
	loc = fa*layer+i*row+j;
	T = AbsoluteTemp(g[TEMP][loc], pp[fa], g[HUMIDITY][loc]);
	for (dd = deposit; dd; dd = dd->next)  {
	  if (!gr->ustar || (rs = dd->diffusivity / gr->ustar) > 100.)
	     rs = 100.;  /* rs wird in Realitaet kaum hoeher als 100. */
	  if (gr->Drmin && (gr->X < 1. || gr->wq >= 0.))  {
	    /* keine Wasseroberflaeche, keine feuchte Oberflaeche */
	    if (gr->Rg == 0.)
	      rc = gr->Drnight;
	    else if (gr->Rg < 400.)
	      rc = gr->Drmin + (gr->Drmax - gr->Drmin) *
		   (1. - pow(gr->Rg * 0.0025, 0.3333333333333333));
	    else
	      rc = gr->Drmin;
	    rc *= dd->depofact;
	  }
	  else  {    /* Wasseroberflaeche */
	    rc = dd->depowater + gr->Drwet * dd->depofact;
	  }
	  rtot = gr->ra + rs + rc;
	  dep = g[dd->et][loc] / rtot;
	  dd->deposition[i*row+j] = dep * pp[fa] / (Ru * T);
	  g[dd->et][loc] -= dep * leveli[fa] * tinc;
	}
      }
  }
}
