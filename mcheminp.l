%{
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mcglobal.h"
#include "mchemsyntax.h"
#include "mchemparse.h"

char yysubstname[TOKENLEN];
%}

D [0-9]*
DM [0-9]+
E [eE]-?[0-9]+

%%

\#.*$	;

thermal    {
	     strcpy(lasttoken, yytext);
	     return (THERMKIN);
	   }

thermal-T2 {
	     strcpy(lasttoken, yytext);
	     return (TERM2KIN);
	   }

troe	   {
	     strcpy(lasttoken, yytext);
	     return (TROEKIN);
	   }

troe-equil {
	     strcpy(lasttoken, yytext);
	     return (TROEQUILKIN);
	   }

photodiss  {
	     strcpy(lasttoken, yytext);
	     return (PHOTOKIN);
	   }

photodiss3 {
	     strcpy(lasttoken, yytext);
	     return (PHOTOKIN3);
	   }

twostream  {
             strcpy(lasttoken, yytext);
             return (TWOSTREAMARG);
           }

special    {
	     strcpy(lasttoken, yytext);
	     return (SPECIALKIN);
	   }

ppm	   {
	     strcpy(lasttoken, yytext);
	     konzunitfact = 0.001 * convertorfact;
	     return (A_UNIT);
	   }

ppb	   {
	     strcpy(lasttoken, yytext);
	     konzunitfact = convertorfact;
	     return (A_UNIT);
	   }

molec\/cm3  {
	     strcpy(lasttoken, yytext);
	     konzunitfact = 1.;
	     return (A_UNIT);
	   }

sec	   {
	     strcpy(lasttoken, yytext);
	     timeunitfact = 1.;
	     return (A_UNIT);
	   }

min	   {
	     strcpy(lasttoken, yytext);
	     timeunitfact = 1. / 60.;
	     return (A_UNIT);
	   }

units	   {
	     strcpy(lasttoken, yytext);
	     return (UNIT);
	   }

refs       {
	     strcpy(lasttoken, yytext);
	     return (REF);
	   }

press      {
	     strcpy(lasttoken, yytext);
	     yylval.dref = &reference.press;
	     return (A_REF);
	   }

temp      {
	     strcpy(lasttoken, yytext);
	     yylval.dref = &reference.T;
	     return (A_REF);
	   }

REACTIONS  {
	     strcpy(lasttoken, yytext);
	     return (REACTIONS);
           }

SUBSTANCES {
	     strcpy(lasttoken, yytext);
	     return (SUBSTANCES);
           }

END  	   {
	     strcpy(lasttoken, yytext);
	     return (END);
           }

-?{D}\.{D} |
-?{D}\.{D}{E} |
-?{DM}{E} |
-{DM}        {
	       yylval.dval = atof(yytext);
	       strcpy(lasttoken, yytext);
               return (FLOATNUM);
             }

[1-9]{D}   |
0            {
	       yylval.ival = atoi(yytext);
	       strcpy(lasttoken, yytext);
	       return (INTNUM);
	     }

[A-Za-z_][0-9A-Za-z_]*  {
			 yylval.item = yysubstname;
	       		 strcpy(lasttoken, yytext);
	       		 if ((int)strlen(yytext) >= TOKENLEN)  {
    			   fprintf(stderr, 
       			      "WARNING : %s is longer than %i characters and will be cut to ",
       			      yytext, TOKENLEN);
	       		   yytext[TOKENLEN-1] = '\0';
    			   fprintf(stderr, "%s\n", yytext);
	       		 }
	       		 strcpy(yysubstname, yytext);
			 return (SUBSTNAME);
		       }

\;		{
	          strcpy(lasttoken, yytext);
		  return (DELIMITER);
		}

[(,)+=]		{
	          strcpy(lasttoken, yytext);
 		  return (*yytext);
 		}

[ \t]+	;

\n	lineno++;

\".*\"	    {
	      yylval.item = (char *)malloc(strlen(yytext) - 1);
	      strncpy(yylval.item, yytext + 1, strlen(yytext) - 2);
	      yylval.item[strlen(yytext)-2] = 0;
	      return (STRING);
	    }

.	    {
	      fprintf(stderr, "Illegal characters in inputfile: %s\nUse # for comments\n",
	             yytext);
	      return (NONSENSE);
	    }

