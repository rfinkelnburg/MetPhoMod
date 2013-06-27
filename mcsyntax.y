%{

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcgroundclasses.h"
#include "mcparse.h"
#include "mccdfin.h"
#include "mcprint.h"
#include "mchemparse.h"
#include "mcnest.h"

#define SCALLD(proc)  if (proc)  {inputerror = 1; YYERROR;}
#define SCALL(proc)  if (proc)  inputerror = 1

%}

%start inpfile

%union  {
  int ival;
  double dval;
  char item[ITEMLEN];
  Date date;
  Time time;
  OutVar *outvar;
}

%token <item>  SECTION IDENTIFIER STRING
%token <dval>  FLOATNUM
%token <ival>  INTNUM TITLE CDFIMPORT AVG CHEMFILE RANGE NONSENSE ARRAY DIMENSION GROUNDLEVEL
%token <ival>  ASCII CDF VARDEFINE AS GROUNDCLASSES SUBDOMAIN
%token <date>  IDATE
%token <time>  DAYTIME

%type <dval> floatnum floatlist
%type <outvar> outitem outlist
%type <item> printident

%%

inpfile : title sections

title : TITLE STRING
      {
        ShowTitle($2);
      }

sections : section
	 | sections section

section : sectionkey statements
        {
          SCALL(EndOfSection());
        }

sectionkey : SECTION
	{
	  SCALL(ChangeSectionTo($1));
	}

statements :
	   | statements statement ';'
	   | statements error ';'
	   {
	     inputerror = 1;
	     yyerrok;
	   }
	   | statements groundclasstable;

statement :
          | variable '=' value
          {
            SCALL(AssignVariable());
            InitArrayVar();
          }
          | variable IDENTIFIER
          {
            SCALLD(FindOption($2));
            SCALL(AssignVariable());
            InitArrayVar();
          }
          | CDFIMPORT STRING
          {
            SCALL(ReadCDFFile($2, actualsection));
          }
          | CHEMFILE STRING
          {
            SCALL(ReadChemFile($2));
            SCALL(TestSpecial());
          }
          | SUBDOMAIN STRING
          {
            SCALL(DefineNestDomain($2, actualsection));
          }
          | outstatement ;
          | IDENTIFIER VARDEFINE value
          {
            SCALL(DefineAVariable($1));
          }
          
groundclasstable : groundclkey titleline datalines;

groundclkey : GROUNDCLASSES
	{
	  SCALL(InitializeGroundClasses());
	}

titleline : idlist ';';

idlist : IDENTIFIER
	{
	  SCALLD(AddClassVar($1));
	}
       | idlist IDENTIFIER
	{
	  SCALLD(AddClassVar($2));
	}

datalines : dataline
	  | datalines dataline;

dataline : clidx numlist ';'
	{
	  SCALL(EndOfClassLine());
	}

clidx : INTNUM
	{
	  SCALLD(AddClassIdx($1));
	}

numlist : FLOATNUM
	{
	  SCALL(AddClassVal($1));
	}
        | numlist FLOATNUM
        {
          SCALL(AddClassVal($2));
        }

variable : IDENTIFIER
	 {
	   SCALLD(SetVariable($1));
	 }
         | IDENTIFIER range
         {
           SCALLD(SetVariable($1));
           SCALLD(CopyRangeToVar());
         }

value : INTNUM
      {
        val.type = INT_NUM;
        val.v.ival = $1;
      }
      | FLOATNUM
      {
        val.type = DOUBLE_NUM;
        val.v.dval = $1;
      }
      | array
      | IDATE
      {
        val.type = DATE_VAL;
        val.v.date = $1;
      }
      | DAYTIME
      {
        val.type = TIME_VAL;
        val.v.time = $1;
      }
      | AVG
      {
        val.type = AVG_VAL;
        val.v.ival = nz;
      }
      | AVG '(' INTNUM ')'
      {
        val.type = AVG_VAL;
        val.v.ival = $3;
      }
      | IDENTIFIER
      {
        SCALLD(CopyVarToVal($1));
      }
      | numbers
      {
        val.type = VECTOR_VAL;
      }
      | intnumbers
      {
        val.type = INT_VECTOR_VAL;
      }
      | STRING
      {
        val.type = TEXT_VAL;
        strcpy(val.v.txt, $1);
      }

array : arraykey dims numbers
	{
	  if (val.v.a.nval < val.v.a.maxval)  {
	    InputError(ERROR, "There are only %i initialisers. You need %i.",
	       val.v.a.nval, val.v.a.maxval);
	    inputerror = 1;
	  }  /* of if */
	}

arraykey : ARRAY
         {
           InitArrayVar();
         }

range : '[' rangespecs ']'
      {
        AllocateArraySpace();
      }

dims : '[' dimspecs ']'
      {
        AllocateArraySpace();
      }


rangespecs : rangespec |
             rangespecs ',' rangespec

dimspecs : dimspec |
           dimspecs ',' dimspec

dimspec : DIMENSION
	{
	  SCALL(SetDimension($1));
	}
	| rangespec
	| DIMENSION '=' GROUNDLEVEL
	{
	  SCALL(SetDimension($1));
	  SCALL(SetRange($1, -1, -1));
	}

rangespec : DIMENSION '=' INTNUM RANGE INTNUM
	  {
	    SCALL(SetDimension($1));
	    SCALL(SetRange($1, $3, $5));
	  }
	  | DIMENSION '=' INTNUM
	  {
	    SCALL(SetDimension($1));
	    SCALL(SetRange($1, $3, $3));
	  }

numbers : '(' floatlist ')'

floatlist : floatnum
	  | floatlist ',' floatnum

floatnum : FLOATNUM
	 {
	   SCALLD(AddArrayVal($1));
	 }

intnumbers : '(' intlist ')'

intlist : INTNUM
	  {
	    AddIntArrayVal($1);
	  }
	  | intlist ',' INTNUM
	  {
	    AddIntArrayVal($3);
	  }

outstatement : cdfheader ':' outlist
	     {
	       AttachOutlist($3);
	     }
	     | asciiheader ':' outlist
	     {
	       AttachOutlist($3);
	     }

cdfheader : CDF '(' STRING ',' INTNUM ',' INTNUM ')'
	  {
	    SCALLD(PrepareOutput(CDF_PRINT, $3, $5, $7));
	  }
	    

asciiheader : ASCII '(' STRING ',' INTNUM ',' INTNUM ')'
	    {
	      SCALLD(PrepareOutput(ASCII_PRINT, $3, $5, $7));
	    }

outlist : outitem
	{
	  $$ = $1;
	}
	| outlist ',' outitem
	{
	  $$ = AppendOutVar($1, $3);
	  if (!$$)  {
	    inputerror = 1;
	    YYERROR;
	  }  /* of if */
	}

outitem : printident
	{
	  if (!($$ = OutItem($1, $1)))  {
	    inputerror = 1;
	    YYERROR;
	  }   /* of if */
	  InitArrayVar();
	}
	| printident AS IDENTIFIER
	{
	  if (!($$ = OutItem($1, $3)))  {
	    inputerror = 1;
	    YYERROR;
	  }   /* of if */
	  InitArrayVar();
	}

printident : IDENTIFIER
	{
	  strcpy($$, $1);
	}
	| IDENTIFIER dims
	{
	  strcpy($$, $1);
	}
