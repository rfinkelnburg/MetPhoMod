/*
   MODULE McGroundClasses
   Dieses Modul ermöglicht eine Behandlung der Bodenparameter mit Hilfe
   von Klassen.
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mcgroundclasses.h"

#define MAXCLASS 40
#define MAXVAR 20

typedef double classrow[MAXVAR];

int nclass = 0, nclvar = 0, valuepos = 0;

int classidx[MAXCLASS];
VarDesc *classvar[MAXVAR];
classrow *classval;


int InitializeGroundClasses(void)
{
  if (nclass)  {
    InputError(ERROR, "You can not define Groundclasses twice.");
    return (1);
  }
  classval = (classrow *)calloc(MAXCLASS, sizeof(classrow));
  nclass = nclvar = valuepos = 0;
  return (0);
}

int AddClassVar(char *name)
{
  VarDesc *v;
  if (nclvar >= MAXVAR)  {
    InputError(ERROR, "Too many groundclass variables. This version supports %d.", MAXVAR);
    return (1);
  }
  for (v = variable; v && strcmp(name, v->name); v = v->next);
  if (!v)  {
    InputError(ERROR, "%s is not a known variable or keyword.", name);
    return (1);
  }
  classvar[nclvar++] = v;
  return (0);
}

int AddClassIdx(int idx)
{
  if (nclass >= MAXCLASS)  {
    InputError(ERROR, "Too many ground-classes. This version supports only %d classes!", MAXCLASS);
    return (1);
  }
  if (!idx)  {
    InputError(ERROR, "Class-Indexes must be != 0.\n");
    return (1);
  }
  classidx[nclass] = idx;
  return (0);
}

int AddClassVal(double val)
{
  if (valuepos >= nclvar)  {
    InputError(ERROR, "Too many numbers on this line.");
    return (1);
  }
  classval[nclass][valuepos++] = val;
  return (0);
}

int EndOfClassLine(void)
{
  if (valuepos < nclvar)  {
    InputError(ERROR, "Too few numbers on this line.");
    return (1);
  }
  valuepos = 0;
  nclass++;
  return (0);
}

void FillGroundValues(void)
{
  int i, j, k, l, loc, gc;
  InitType initstate = WAS_SET_BY_USER;
  if (!nclass)  return;
  for (i = nx; --i; )
    for (j = ny; --j; )  {
      loc = i*row+j;
      gc = (int)(ground[loc].groundclass + 0.5);
      if (gc)  {
        for (k = nclass; --k >= 0 && classidx[k] != gc; );
        if (k < 0)  {
          InputError(WARNING, "Variable \"GroundClass\" includes class #%d, which is not part of the table.",
             gc);
          initstate = WAS_PARTITALLY_SET;
        }
        else
          for (l = nclvar; l--; )
            classvar[l]->v.g[loc].Tg[0] = classval[k][l];
      }
      else  initstate = WAS_PARTITALLY_SET;
    }
  for (l = nclvar; l--; )
    classvar[l]->init = initstate;
  free(classval);
}
