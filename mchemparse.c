/*
   IMPLEMENTATION MODULE mchemparse.c
   Dieses Modul interpretiert und berechnet die chemischen Reaktionen
   des Modelles.
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mchemparse.h"

Educt educt[MAXEDUCT];
Product product[MAXPRODUCT];
Substance subst[MAXSUBST] = {HUMIDITY, 0., 0., 0., "H2O", NULL, NULL, 0, 0, 0},
  *fast[MAXSUBST], *slow[MAXSUBST];
Reaction reaction[MAXREACT];
int neduct, netot, nproduct, oldneduct, oldnproduct, nsubst, nreact, lineno, chemerror,
    nfast, nslow, nspecial;
double konzunitfact, timeunitfact, convertorfact, fixedconc;
char lasttoken[256];
Reaction *reactionptr[MAXEDUCT+MAXPRODUCT];
KinReference reference;

int chiparse(void);

double spow(double x, int y)
{
  int i;
  double p;
  if (!y)  return (1.);
  p = x;
  while (--y)  p *= x;
  return (p);
}

void PrintChemStatus()
{
  int i;
  char txt[30];
  printf("Substances:\n\nname		     sources       sinks\n");
  for (i = 0; i < nsubst; i++)  {
    printf("%-16s%12i%12i\n", subst[i].name, subst[i].nproduction, subst[i].ndecay);
  }
}

void CalcConvertor(void)
{
  convertorfact = kBoltz * reference.T / (1.e-15 * reference.press);
}

void ConvertUnits()
{
  reaction[nreact].K *= spow(konzunitfact, reaction[nreact].netot-1) * timeunitfact * fixedconc;
}

Substance *PlaceSubstance(char *name, BOOL mustbenew)
{
  int i;
  Substance *s;
  for (i = nsubst, s = subst; i && strcmp(name, s->name); i--, s++);
  if (i > 0)  {
    if (mustbenew)  {
      fprintf(stderr,
	 "CHEM-ERROR: Line %i. Substance \"%s\" was declared twice.\n", name);
      return (NULL);
    }
    return (subst + nsubst - i);
  }
  if (nsubst >= MAXSUBST)  {
    fprintf(stderr,
       "FATAL CHEM-ERROR: Too much substances.\n"
       "This version of meteochem supports a maximum of %i.\n", MAXSUBST);
    return (NULL);
  }
  strcpy(subst[nsubst].name, name);
  subst[nsubst].subs = SUBS + nsubst - 1;
  subst[nsubst].inert = mustbenew;
  return (subst + nsubst++);
}

int PlaceEduct(char *name, int fact)
{
  Substance *newsubs;
  int i;
  if (!fact)  {
    fprintf(stderr,
       "CHEM-ERROR: Line %i.0 is not allowed as a prefactor. Set to 1.\n",
       lineno);
    chemerror = 1;
    fact = 1;
  }
  netot += fact;
  if (!strcmp(name, "M"))  {
    fixedconc *= spow(1.e9, fact);
    return (0);
  }
  if (!strcmp(name, "O2"))  {
    fixedconc *= spow(780.8e6, fact);  /* Nach Weischet p.37 */
    return (0);
  }
  if (!strcmp(name, "N2"))  {
    fixedconc *= spow(209.5e6, fact);
    return (0);
  }
  if (neduct >= MAXEDUCT)  {
    fprintf(stderr,
       "FATAL CHEM-ERROR: Too many reaction educts.\n"
       "This version of meteochem supports a maximum of %i.\n", MAXEDUCT);
    return (1);
  }
  if (!(educt[neduct].comp = newsubs = PlaceSubstance(name, FALSE)))  return (1);
  for (i = neduct; i-- > oldneduct; )
    if (newsubs == educt[i].comp)  {
      fprintf(stderr,
         "CHEM-ERROR: Line %i. Substance %s occurs twice in the educts of this reaction.\n",
         lineno, name);
      chemerror = 1;
    }
  educt[neduct].prefact = fact;
  neduct++;
  return (0);
}

int PlaceProduct(char *name, double fact)
{
  Substance *newprod;
  int i;
  if (!fact)  {
    fprintf(stderr,
       "CHEM-ERROR: Line %i. 0 is not allowed as a prefactor. Set to 1.\n",
       lineno);
    chemerror = 1;
    fact = 1;
  }
  if (!strcmp(name, "M") || !strcmp(name, "O2") || !strcmp(name, "N2"))
    return (0);
  if (nproduct >= MAXPRODUCT)  {
    fprintf(stderr,
       "FATAL CHEM-ERROR: Too many reaction products.\n"
       "This version of meteochem supports a maximum of %i.\n", MAXPRODUCT);
    return (1);
  }
  if (!(product[nproduct].comp = newprod = PlaceSubstance(name, FALSE)))  return (1);
  for (i = nproduct; i-- > oldnproduct; )
    if (newprod == product[i].comp)  {
      fprintf(stderr,
         "CHEM-ERROR: Line %i. Substance %s occurs twice in the products of this reaction.\n",
         lineno, name);
      chemerror = 1;
    }
  product[nproduct].prefact = fact;
  nproduct++;
  return (0);
}

void chierror(char *s)
{
  extern char *chitext;
  fprintf(stderr,
     "CHEM-ERROR: %s at line %i. Offending token: \"%s\"\n",
     s, lineno, lasttoken);
}

void AssignReactToSubst()
{
  int i, j, k;
  Reaction **rptr, *r;
  Substance *s;
  Educt *e;
  Product *p;
  rptr = reactionptr;
  for (i = nsubst, s = subst; i--; s++)  {
    s->ndecay = s->nproduction = 0;
    s->decay = rptr;
    for (j = nreact, r = reaction; j--; r++)  {
      for (k = r->neduct, e = r->educt; k && e->comp != s; k--, e++);
      if (k)  {*rptr++ = r; s->ndecay++;}
    }
    s->production = rptr;
    for (j = nreact, r = reaction; j--; r++)  {
      for (k = r->nproduct, p = r->product; k && p->comp != s; k--, p++);
      if (k)  {*rptr++ = r; s->nproduction++;}
    }
  }
}

int EndReaction()
{
  if (nreact >= MAXREACT)  {
    fprintf(stderr, 
       "Too many reactions. Maximum of this version is %i\n",
       MAXREACT);
       return (1);
  }
  reaction[nreact].educt = educt + oldneduct;
  reaction[nreact].product = product + oldnproduct;
  reaction[nreact].neduct = neduct - oldneduct;
  reaction[nreact].nproduct = nproduct - oldnproduct;
  reaction[nreact].netot = netot;
#ifdef FULLM
  reaction[nreact].nm = nm;
  reaction[nreact].nox = nox;
  reaction[nreact].nn2 = nn2;
  nm = nox = nn2 = 0;
#endif
  ConvertUnits();
  oldneduct = neduct; oldnproduct = nproduct;
  netot = 0; fixedconc = 1.;
  nreact++;
  return (0);
}

int CreateChemTable(void)
{
  int i;
  VarDesc *v, *vh;
  static char conzunit[] = "ppb";
  if (!(v = (VarDesc *)calloc(nsubs, sizeof(VarDesc))))  return (1);
  for (i = nsubs; i--; )  {
    v[i].name = subst[i+1].name;
    v[i].comment = NULL;
    v[i].unit = conzunit;
    v[i].storetype = GRID_VAL;
    v[i].v.et = subst[i+1].subs;
    v[i].dims = ALL_DIM;
    v[i].section = v[i].origsection = CHEMISTRY;
    v[i].init = SET_DEFAULT;
    v[i].defval = 0.;
    v[i].rmin = 0.;
    v[i].rmax = 100000.;
    v[i].inputtype = NORMAL_NUM;
    v[i].option = SET_DEFAULT << 8;
    v[i].next = v + i + 1;
  }
  v[nsubs-1].next = NULL;
  for (vh = variable; vh->next; vh = vh->next);
  vh->next = v;
  return (0);
}

extern "C" int chiwrap()
{
  return (1);
}

int ReadChemFile(char *fname)
{
  extern FILE *chiin;
  int i;
  short int analyse = 0;
  neduct = netot = oldneduct = nproduct = oldnproduct = nreact = 
  nfast = nslow = 0;
  lineno = nsubst = 1;
  timeunitfact = konzunitfact = 1.;
  reference.T = 298.15; reference.press = 101325.;
  fixedconc = 1.;
  CalcConvertor();
  if (!(chiin = fopen(fname, "r")))  {
    InputError(ERROR, "Unable to open \"%s\".", fname);
    return (1);
  }
  chemerror |= chiparse();
  if (chemerror)  return (1);
  if (nsubs != nsubst-1)  {
    InputError(ERROR, "The number of species found in chem-file (%i) is unequal to variable \"nsubs\" (%i)",
       nsubst-1, nsubs);
    return (1);
  }
  if (CreateChemTable())  return (1);
  PrintChemStatus();
  return (0);
}

extern double SpecialConst(int idx, double T, double M);
extern char *specialformula[];
extern const int nspecialformula;

int GenerateSpecial(FILE *f, char *chemname)
{
  int i, j = 0;
  extern FILE *chiin;
  short int analyse = 0;
  neduct = netot = oldneduct = nproduct = oldnproduct = nreact = 
  nfast = nslow = 0;
  lineno = nsubst = 1;
  timeunitfact = konzunitfact = 1.;
  reference.T = 298.15; reference.press = 101325.;
  fixedconc = 1.;
  CalcConvertor();
  if (!(chiin = fopen(chemname, "r")))  {
    fprintf(stderr, "Unable to open \"%s\".", chemname);
    return (1);
  }
  groundinterface = TRUE;
  chemerror |= chiparse();
  if (chemerror)  {
    fprintf(stderr, "ERROR: Reading chem-file\n");
    return (1);
  }
  PrintChemStatus();
  fprintf(f, "#include <stddef.h>\n"
             "#include <stdlib.h>\n"
             "#include <math.h>\n\n"
             "int nspecialformula = %i;\n\nchar *specialformula[] = {\n", nspecial);
  if (nspecial)  {
    for (i = 0; i < nreact; i++)
      if (reaction[i].type == SPECIAL)  {
        if (j++)  fprintf(f, ",\n");
        fprintf(f, "   \"%s\"", reaction[i].formula);
      }
    }
  else
    fprintf(f, "   \"dummy\"");
  fprintf(f, "};\n\ndouble SpecialConst(int idx, double T, double M)\n"
  	     "{\n  double K, k1, k2, k3;\n");
  if (nspecial)  {
    fprintf(f, "  switch (idx)  {\n");
    for (i = j = 0; i < nreact; i++)
      if (reaction[i].type == SPECIAL)  {
        fprintf(f, "#line %d \"%s\"\n", reaction[i].specialine, chemname);
        fprintf(f, "    case %d : %s;\n             break;\n", j++, reaction[i].formula);
      }
    fprintf(f, "  }\n");
  }
  else
    fprintf(f, "  K = 0.;\n");
  fprintf(f, "  return (K);\n"
             "}\n\n");
  return (0);
}

int TestSpecial()
{
  int i, j;
  if (nspecial && nspecialformula != nspecial)  goto error_end;
  for (i = j = 0; i < nreact; i++)
    if (reaction[i].type == SPECIAL)
      if (strcmp(reaction[i].formula, specialformula[j++]))  goto error_end;
  return (0);
error_end :
  fprintf(stderr, "ERROR in Chemfile: Special-Formulas are incompatible with\n                   the actual version of the program!\n");
  return (1);
}
