%{
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mcprint.h"
#include "mcsyntax.h"

Date ConvertToDate(char *t);

Time ConvertToTime(char *t);

void AssignItem(char *to, char *from);

%}

D [0-9]*
E [eE]-?[0-9]+

%%

\#.*$	;

TITLE  return (TITLE);

TIME   |
GRID   |
ENVIRONMENT |
OPTIONS   |
INITIAL-DATA |
BORDER-DATA |
GROUND  |
TOP     |
NORTH   |
SOUTH   |
WEST    |
EAST    |
CHEMISTRY   |
EMISSIONS   |
REDUCTION   |
DEPOSITION  |
NESTING     |
OUTPUT      {
	      AssignItem(yylval.item, yytext);
	      return (SECTION);
	    }

CDFIMPORT   return (CDFIMPORT);
AVG	    return (AVG);
CHEMFILE    return (CHEMFILE);
SUBDOMAIN   return (SUBDOMAIN);
ARRAY	    return (ARRAY);
GROUNDLEVEL	return (GROUNDLEVEL);
GROUNDCLASSES	return (GROUNDCLASSES);

ASCII	return (ASCII);

CDF	return (CDF);

AS	return (AS);

X	    {
	      yylval.ival = X_DIM;
	      return DIMENSION;
	    }

Y	    {
	      yylval.ival = Y_DIM;
	      return DIMENSION;
	    }

Z	    {
	      yylval.ival = Z_DIM;
	      return DIMENSION;
	    }

-?{D}\.[0-9]+ |
-?{D}\.{D}{E} |
-?{D}\./[ \n,\t);] {
	       yylval.dval = atof(yytext);
               return (FLOATNUM);
             }

-?[1-9]{D}   |
0            {
	       yylval.ival = atoi(yytext);
	       return (INTNUM);
	     }

[A-Za-z_][0-9A-Za-z_-]*  {
	       AssignItem(yylval.item, yytext);
               return (IDENTIFIER);
             }

:=	return (VARDEFINE);

[\;=\[\](),:] return (*yytext);

\.\.         return (RANGE);

[0-9]+-[0-9]+-[0-9]+   {
               yylval.date = ConvertToDate(yytext);
               return (IDATE);
             }

[0-9]+:[0-9]+ {
	       yylval.time = ConvertToTime(yytext);
	       return (DAYTIME);
	     }

[ \t]       ;

\n	    inplineno++;

\".*\"      {
              yytext[strlen(yytext)-1] = '\0';
	      AssignItem(yylval.item, yytext+1);
	      return (STRING);
	    }

.	    {
	       fprintf(stderr, "ERROR: Illegal character in Inputfile: %c\n",
	          *yytext);
	       return (NONSENSE);
	    }

%%

Date ConvertToDate(char *t)
{
  char *p;
  Date d;
  for (p = t; *p != '-'; p++);
  *p = '\0'; d.month = atoi(t); t = ++p;
  while (*p != '-')  p++;
  *p = '\0'; d.day = atoi(t);
  d.year = atoi(++p);
  if (d.day < 1 || d.day > 31 || d.month < 1 || d.month > 12 || d.year < 1900 ||
      d.year > 2100)  {
    fprintf(stderr, "Impossible date: %i-%i-%i!\n", d.month, d.day, d.year);
    d.year = -1;
  }
  return (d);
}

Time ConvertToTime(char *tname)
{
  char *p;
  Time t;
  for (p = tname; *p != ':'; p++);
  *p = '\0'; t.hour = atoi(tname);
  t.minute = atoi(++p);
  if (t.hour > 23 || t.minute > 59)  {
    fprintf(stderr, "Impossible Time: %i:%i\n", t.hour, t.minute);
    t.hour = -1;
  }
  return (t);
}

void AssignItem(char *to, char *from)
{
  if (strlen(from) >= (size_t)ITEMLEN)
    printf("Warning: Item *s is too long and will be cut to 63 Characters\n");
  strncpy(to, from, ITEMLEN);
}

