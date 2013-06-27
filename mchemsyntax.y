%{

#include <stdio.h>
#include <string.h>
#include "mcglobal.h"
#include "mchemparse.h"

%}

%start inpfile

%union  {
  int ival;
  double dval;
  double *dref;
  char *item;
}

%token <ival> INTNUM DELIMITER NONSENSE PHOTOKIN THERMKIN TERM2KIN TROEKIN
	      PHOTOKIN3 TROEQUILKIN SPECIALKIN UNIT A_UNIT
%token <ival> REACTIONS SUBSTANCES END REF
%token <dval> FLOATNUM
%token <item> SUBSTNAME STRING
%token <dref> A_REF

%type <dval> floatnum
%type <ival> reactions substances substance fullreaction units inpfile

%%      /* beginning of rules */

inpfile : REACTIONS reactions END
        | SUBSTANCES substances END
        | SUBSTANCES substances REACTIONS reactions END

reactions : fullreaction
          | reactions fullreaction
            {
	      AssignReactToSubst();
            }
          ;

fullreaction : educts '=' products kindata DELIMITER
           {
             if (EndReaction())
               YYABORT;
           }
           | UNIT '(' units ')' DELIMITER
           | REF '(' refs ')' DELIMITER
           {
             CalcConvertor();
           }
           | error DELIMITER
           {
             chemerror = 1;
             neduct = oldneduct; nproduct = oldnproduct;
             yyerrok;
           }
           ;

educts : educt
       | educts '+' educt ;

products :
         | productitems ;

productitems : product
            | productitems '+' product ;

educt : SUBSTNAME
        {
          if (PlaceEduct($1, 1))
            YYABORT;
        }
      | INTNUM SUBSTNAME
        {
          if (PlaceEduct($2, $1))
            YYABORT;
        }
      ;

substances : substance
	   | substances substance

substance : SUBSTNAME DELIMITER
          {
            if (!PlaceSubstance($1, TRUE))  {
              chemerror = 1;
            }
          }

floatnum : FLOATNUM
         | INTNUM
         {
           $$ = (double)$1;
         }
         ;

product : SUBSTNAME
          {
            if (PlaceProduct($1, 1.))
              YYABORT;
          }
        | floatnum SUBSTNAME
          {
            if (PlaceProduct($2, $1))
              YYABORT;
          }
        ;

kindata : PHOTOKIN '(' floatnum ',' floatnum ')'
	  {
            if (!groundinterface)  {
              fprintf(stderr, "CHEM-ERROR: photodiss-reactions can only be calculated\n"
                              "            with groundinterface switched on.");
              chemerror = 1;
            }
            reaction[nreact].K = $3;
            reaction[nreact].EoR = $5;
            reaction[nreact].type = PHOTODISS;
          }
        | PHOTOKIN3 '(' floatnum ',' floatnum ',' floatnum ')'
	  {
            if (!groundinterface)  {
              fprintf(stderr, "CHEM-ERROR: photodiss-reactions can only be calculated\n"
                              "            with groundinterface switched on.");
              chemerror = 1;
            }
            reaction[nreact].K = $3;
            reaction[nreact].N = $5;
            reaction[nreact].EoR = $7;
            reaction[nreact].type = PHOTODISS3;
          }
        | THERMKIN '(' floatnum ',' floatnum ')'
          {
            reaction[nreact].K = $3;
            reaction[nreact].EoR = $5;
            reaction[nreact].type = THERMAL;
          }
        | TERM2KIN '(' floatnum ',' floatnum ')'
          {
            reaction[nreact].K = $3;
            reaction[nreact].EoR = $5;
            reaction[nreact].type = THERMAL2;
          }
        | TROEKIN '(' floatnum ',' floatnum ',' floatnum ',' floatnum ')'
          {
            reaction[nreact].K = 1.;
            reaction[nreact].k0 = $3;
            reaction[nreact].N = $5;
            reaction[nreact].kinf = $7;
            reaction[nreact].M = $9;
            reaction[nreact].type = TROE;
          }
        | TROEQUILKIN '(' floatnum ',' floatnum ',' floatnum ',' floatnum ',' floatnum ',' floatnum ')'
          {
            reaction[nreact].k0 = $3;
            reaction[nreact].N = $5;
            reaction[nreact].kinf = $7;
            reaction[nreact].M = $9;
            reaction[nreact].K = $11;
            reaction[nreact].EoR = $13;
            reaction[nreact].type = TROEQUIL;
          }
        | SPECIALKIN '(' STRING ')'
          {
            reaction[nreact].K = 1.;
            reaction[nreact].formula = $3;
            reaction[nreact].specidx = nspecial++;
            reaction[nreact].specialine = lineno+1;
            reaction[nreact].type = SPECIAL;
          }
        ;

units   : A_UNIT
        | units ',' A_UNIT
        ;

refs	: ref
	| refs ',' ref

ref	: A_REF '=' floatnum
	{
	  *$1 = $3;
	}

