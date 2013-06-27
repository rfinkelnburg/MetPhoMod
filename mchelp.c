/*
   IMPLEMENTATION MODULE mchelp.c
   Gibt dem Benutzer Informationen zum Programm.
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mchelp.h"
#include "mc_module.hh"
#include "mctwostream/mctwostream.h"

char *InpTypeNames(int t)
{
  static char *names[] = {
    "integer", "double", "on-off option", "coriolis option",
    "pressure algorithm option",
    "bordertype option", "filtertype option",
    "time of day (hh:mm)", "date MM-DD-YYYY",
    "array", "\"AVG\"-Keyword", "Vector (xx.x, xx.x, xx.x, ...)",
    "int-vector (xx, xx, xx, ...)",
    "text \"abcd...\"", "turbulence type option", "advection type option",
    "Soil temperature option"
    },
    txt[256];
  int i;
  *txt = '\0';
  for (i = 0; t; i++, t >>= 1)
    if (t & 1)  {
      strcat(txt, names[i]);
      if (t > 1)  strcat(txt, "/");
    }
  return (txt);
}

char *DimNames(DimensionSet d)
{
  static char *names[] = {
    "X", "Y", "Z", "Time"
    },
    txt[32];
  int i;
  if (d)  {
    strcpy(txt, "[");
    for (i = 0; d; i++, d >>= 1)
      if (d & 1)  {
        strcat(txt, names[i]);
        if (d > 1)  strcat(txt, ",");
      }
    strcat(txt, "]");
    return (txt);
  }
  else  return ("-");
}

char *GetOptionFromVal(int val, InputType type)
{
  int i;
  for (i = MAXOPTION; i >=0 && (option[i].type != type || option[i].value != val) ; i--);
  return (i >= 0 ? option[i].name : (char *)"no valid default!");
}

void PrintHelp(int argc, char **argv)
{
  static char *helptxt[] =  {
    "Help is available for these topics:\n\n"
    "cmd          : Instructions how to invoke meteochem.\n"
    "syntax       : The syntax of the Input-File\n"
    "sections     : Comments about the sections in the Input-File\n"
    "about        : General info about meteochem\n"
    "photodiss    : A list o available photodissociation reactions\n"
    "               in the twostream module\n"
    "<section>    : Info about a specific section\n"
    "<variable>   : Info about a specific variable or option\n",
  
    "meteochem command-line syntax:\n\n"
    "meteochem -help [<keyword>] : Help about the indicated keyword\n"
    "                              will be shown\n\n"
    "meteochem -h [<keyword>]    : Same effect as command above\n\n"
    "meteochem [-RESTART restart-file [-REOPEN]] [-noexcpt] [-syntax] \\\n"
    "          [-nooutput] [-debug] [-sequential] inpfile     \n"
    "   starts the calculation defined in <inputfile>\n"
    "   <inputfile>              : The name of the input File. For the syntax\n"
    "                              of this file call \"-help syntax\".\n"
    "   -RESTART <dumpfile>      : Read <dumpfile> and start calculation\n"
    "                              at that point of a previous simulation.\n"
    "   -REOPEN                  : Tries to reopen and to append new results\n"
    "                              to the old cdf-Files after a program restart.\n"
    "   -noexcpt                 : Does not handle floating-point exceptions.\n"
    "   -sytnax                  : Does only syntax checking on the input\n"
    "                              File. No calculations will be performed.\n"
    "   -debug                   : Only parallel version: Will start workers\n"
    "                              in Debug-Mode (see pvm-Documentation).\n"
    "   -sequential              : Only parallel version: Use sequential mode.\n\n",
  
    "meteochem input-files expect the syntax (in yacc-notation)\n"
    "shown below:\n\n"
#include "mcsyntax.incl"
    ,

    "The input-File of meteochem can have these sections:\n\n"
    "TITLE           : The Title of the calculation\n"
    "TIME            : Values dealing with the time, such as time-step etc.\n"
    "GRID            : Values dealing with the grid-architecture\n"
    "OPTIONS         : Specify program options.\n"
    "ENVIRONMENT     : Topography, geographical coordinates, etc.\n"
    "INITIAL-DATA    : Initial meteorological values at the\n"
    "                  beginning of the calculation\n"
    "CHEMISTRY       : Chemical Mechanism, and initial concentrations of chemical species\n"
    "RADIATION       : Input parameters for the TwoStream radiation module.\n"
    "BORDER-DATA     : Values at the borders of the modelling domain.\n"
    "       In this Section there exist the subsections:\n"
    "       NORTH, SOUTH, WEST, EAST, GROUND, TOP\n"
    "EMISSIONS       : Emissions\n"
    "REDUCTION       : Emission reductions\n"
    "DEPOSITION      : Deposition-Parameters\n"
    "NESTING         : Create output-files for a subsequent nested runs\n"
    "GROUNDCLASSES   : Define groundclasses and assign soil parameters to each class\n"
    "OUTPUT          : Definition of the output files to be written\n"
    "                  by the program.\n\n"
    "In the input-File, the Sections should be called in the order\n"
    "shown above!\n\n",

    "meteochem is eulerian mesoscale(alpha) meteorological model.\n"
    "Sophisticated routines for the calculation of chemical kinetics are\n"
    "included in the model. For this reason meteochem is especially\n"
    "strong in modelling summer smog episodes in a prognostic manner.\n\n"
    "meteochem was developped independently by S.Perego, as a part of\n"
    "the Swiss POLLUMET investigation program.\n"};
  static char *helpkey[] =  {
    "cmd", "syntax", "sections", "about", "photodiss"};
  int i, j;
  VarDesc *v;
  if (argc == 0)  {
    printf("%s\n", helptxt[0]);
    return;
  }
  if (argc > 1)  {
    printf("%s\n", helptxt[1]);
    return;
  }
  for (i = 5; --i >= 0 && strcmp(*argv, helpkey[i]); );
  if (i >= 0)
    if (i < 4) {
      printf("%s\n", helptxt[i+1]);
      return;
    }
    else {
      twostreammodule.PrintPhotoDissReact();
      return;
    }
  SectionDesc *s;
  s = McInterface::mcsect.FindSection(*argv);
  if (s)  {
    printf("Section %s includes the variables:\n\n", *argv);
    for (v = variable; v; v = v->next)
      if (s->id == v->section)
        printf("%-2c%-22s: %s\n",
	       (v->init == MUST_BE_SET ? '*' : 
		(v->init == CALCULATED_VALUES ? ' ' : '+')),
	       v->name, v->comment);
    if (s->dependson)  printf("\nIt depends on Section(s):\n\n");
    for (const SectionDesc *d = s->dependson; d; d = d->dependson)
      printf("%s\n", d->name);
    printf("\n");
    return;
  }
  for (v = variable; v && strcmp(*argv, v->name); v = v->next);
  if (v)  {
    printf("NAME         : %s\n", *argv);
    printf("MEANING      : %s\n", v->comment);
    if (v->unit)
      printf("UNITS        : %s\n", v->unit);
    printf("SECTION      : %s\n", McInterface::mcsect.FindSection(v->section)->name);
    printf("TYPE         : %s\n", InpTypeNames(v->inputtype));
    if (v->inputtype & ARRAY_VAL)  {
      printf("DIMENSIONS   : %s\n", DimNames(v->dims));
    }
    printf("INIT-COND    : %s\n",
	   (v->init == MUST_BE_SET ? "must be set" : 
	    (v->init == CALCULATED_VALUES ? "Calculated value" : "can be set")));
    printf("CDF possible : %s\n", (v->dims && v->init != CALCULATED_VALUES ? "yes" : "no"));
    if (v->init != CALCULATED_VALUES)  {
      if (v->inputtype & ALL_OPTIONS)  {
        printf("\nOne of the following options is expected:\n\n");
        for (j = 0; j < MAXOPTION; j++)
          if (option[j].type & v->inputtype)
            printf("%-22s : %s\n", option[j].name, option[j].comment);
        if (!(v->init & (MUST_BE_SET | CALCULATED_VALUES)))
          printf("\nDEFAULT       : %s\n", GetOptionFromVal((int)v->defval, v->inputtype));
      }
      else  {
        printf("RANGE        : %lf..%lf\n", v->rmin, v->rmax);
        if (!(v->init & (MUST_BE_SET | CALCULATED_VALUES)))
          printf("DEFAULT      : %lf\n", v->defval);
      }
    }
    printf("\n");
    if (v->option & STARTING_WITH_ZERO) {
      printf("Data must contain values for boundaries! ((nx+2)*(ny+2) values in total)\n\n");
    }
    return;
  }
  printf("Keyword \"%s\" is not known to meteochem.\n\n", *argv);
}

