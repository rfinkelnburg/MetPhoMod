/*
   IMPLEMENTATION MODULE mcemmiss
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mcemiss.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

EmissionDesc *firstemiss = NULL;

Coordinate *AllocateCoordVar(int size)
{
  Coordinate *coord;
  coord = NEW(size, Coordinate);
  if (!coord)  {
    InputError(ERROR, "Unable to allocate memory (in AllocateCoordVar).");
    inputerror = TRUE;
    return (NULL);
  }
  else  return (coord);
}

double *PlaceEmissionVariable(VarDesc *v, VarDesc *vh, int pointvalues, Coordinate *coord)
{
  EmissionDesc *e, *n;
  for (e = firstemiss; e && e->vo != v; e = e->next);
  if (!(n = NEW(1, EmissionDesc))  ||
      !(n->em = NEW((pointvalues ? pointvalues : layer), double)))  {
    fprintf(stderr, "ERROR: allocating memory in PlaceEmissionVariable. Exiting!\n");
    exit (1);
  }
  n->v = vh;
  n->vo = v;
  n->reduction = (e ? e->reduction : NEW(1, double));
  *n->reduction = 1.;
  n->next = firstemiss;
  firstemiss = n;
  n->cem.n = pointvalues;
  n->et = v->v.et;
  if (pointvalues)  {
    n->cem.c = coord;
  }
  return (n->em);
}

VarDesc *GetEmVarWithID(int id)
{
  EmissionDesc *e;
  for (e = firstemiss; e && e->v->id != id; e = e->next);
  return (e->v);
}

EmissionDesc *FindEmissionVariable(VarDesc *v)
{
  EmissionDesc *e;
  for (e = firstemiss; e && e->vo != v; e = e->next);
  return (e);
}

VarDesc *EmissionVariable(VarDesc *v, int pointvalues, int id, Coordinate *coord)
{
  static double bloed;
  static int newid = 333;
  VarDesc *vh;
  vh = NEW(1, VarDesc);
  *vh = *v;
  /*  vh->name = strdup(v->name); */
  vh->unit = "mol/m^2*sec";
  vh->id = (id ? id : newid++);
  if (actualsection->id == EMISSIONS)  {
    vh->comment = "dynamically created emission variable";
    if (pointvalues)  {
      vh->storetype = DOUBLE_PTR;
      vh->dims = COUNT_DIM ;
      vh->rmax = 5.;
    }
    else  {
      vh->storetype = DOUBLE_PTR;
      vh->dims = X_DIM | Y_DIM;
      vh->rmax = 1.e-5;
    }
    vh->ncoord = pointvalues;
    vh->v.d = PlaceEmissionVariable(v, vh, pointvalues, coord);
    vh->init = InitType(vh->option >> 8);
    vh->defval = 0.;
    vh->rmin = 0.;
    vh->inputtype = NORMAL_NUM;
  }
  if (actualsection->id == REDUCTION)  {
    vh->storetype = DOUBLE_PTR;
    vh->section = REDUCTION;
    vh->dims = 0;
    vh->comment = "dynamically created reduction variable";
    if (!(vh->v.d = FindEmissionVariable(v)->reduction))  {
      InputError(ERROR, "%s not found in Emission-Inventory. Can not be reduced!\n",
         v->name);
      vh->v.d = &bloed;
      inputerror = 1;
    }
    vh->init = SET_DEFAULT;
    vh->defval = 1.;
    vh->rmin = 0.;
    vh->rmax = 1.e6;
    vh->inputtype = NORMAL_NUM;
  }
  return (vh);
}

static void OutsideWarning(Coordinate *c)
{
#ifdef PARALLEL
  if (master)
#endif
  InputError(WARNING, "Point source with coords %g/%g/%g is outside of modelling domain.",
  	     c->x, c->y, c->z);
}

void ConvertEmissionCoords(int n, Coordinate *coord, int xoff)
{
  int hloc, i, j, k, h;
  double z;
  Coordinate *c;
  for (h = n, c = coord; h--; c++)  {
    i = (int)(c->x / dx)+1 - xoff;
    j = (int)(c->y / dy)+1;
    if (i < 1 || i > nxm || j < 1 || j > nym)  {
      OutsideWarning(c);
      c->loc = 0;
    }
    else  {
      hloc = i*row + j;
      z = c->z + topo[hloc] - reflevel;
      for (k = ground[hloc].firstabove; k < nz && z > zfaces[k]; k++);
      if (k > nzm)  {
        OutsideWarning(c);
        c->loc = 0;
      }
      else
        c->loc = hloc + k*layer;
    }
  }
}

void Emit(double tinc)
{
  EmissionDesc *e;
  Coordinate *c;
  int h, i, j, k, loc, vloc;
  double T;
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
  if (!firstemiss)  return;
  for (i = xs; i < xe; i++)
    for (j = ny; --j; )  {
      loc = i*row+j;
      k = ground[loc].firstabove;
      vloc = loc+k*layer;
      T = AbsoluteTemp(g[TEMP][vloc], pp[k], g[HUMIDITY][vloc]) *
          tinc * Ru / (pp[k] * level[k]);
      for (e = firstemiss; e; e = e->next)
        if (!e->cem.n)
          g[e->et][vloc] += T * e->em[loc] * *e->reduction;
    }
  for (e = firstemiss; e; e = e->next)
    for (h = 0, c = e->cem.c; h < e->cem.n; h++, c++)
      if (c->loc)  {
	k = c->loc / layer;
	T = AbsoluteTemp(g[TEMP][c->loc], pp[k], g[HUMIDITY][c->loc]) *
            tinc * Ru / (pp[k] * level[k] * dx * dy);
	g[e->et][c->loc] += T * e->em[h] * *e->reduction;
      }
}
