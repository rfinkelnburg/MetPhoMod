
# line 2 "mcsyntax.y"

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


# line 23 "mcsyntax.y"
typedef union   {
  int ival;
  double dval;
  char item[ITEMLEN];
  Date date;
  Time time;
  OutVar *outvar;
} YYSTYPE;
# define SECTION 257
# define IDENTIFIER 258
# define STRING 259
# define FLOATNUM 260
# define INTNUM 261
# define TITLE 262
# define CDFIMPORT 263
# define AVG 264
# define CHEMFILE 265
# define RANGE 266
# define NONSENSE 267
# define ARRAY 268
# define DIMENSION 269
# define GROUNDLEVEL 270
# define ASCII 271
# define CDF 272
# define VARDEFINE 273
# define AS 274
# define GROUNDCLASSES 275
# define SUBDOMAIN 276
# define IDATE 277
# define DAYTIME 278
#define mciclearin mcichar = -1
#define mcierrok mcierrflag = 0
extern int mcichar;
extern int mcierrflag;
#ifndef YYMAXDEPTH
#define YYMAXDEPTH 150
#endif
YYSTYPE mcilval, mcival;
typedef int mcitabelem;
#include <stdio.h>
# define YYERRCODE 256
mcitabelem mciexca[] ={
	-1, 1,
	0, -1,
	-2, 0,
-1, 10,
	59, 11,
	-2, 5,
	};
# define YYNPROD 73
# define YYLAST 239
mcitabelem mciact[]={

    55,    12,    89,    19,    35,   105,    98,    60,    15,   118,
    16,    81,    80,   125,   122,     3,    25,    24,   124,   123,
    23,    17,   117,   111,   110,   106,   105,    93,    64,    81,
    87,    73,    72,    65,    32,    31,    30,    28,     8,   109,
    70,    38,     7,    97,    79,   114,    83,    76,    68,    99,
    51,    75,    62,   115,    84,    27,    26,    40,    39,   102,
   100,   120,   103,   101,   119,    92,    91,    88,   127,   126,
   112,    67,    95,    74,    42,    41,    43,     5,    22,    21,
    78,    96,     9,    58,    54,    59,    52,    46,    34,    86,
    63,    37,    61,    36,   113,    82,    20,    18,    14,    13,
    11,    10,     6,     4,     2,    69,    77,     1,     0,     0,
    57,     0,    71,     0,    85,     0,     0,     0,     0,     0,
     0,     0,    90,     0,     0,     0,    94,     0,     0,     0,
     0,     0,     0,   104,     0,     0,     0,   108,     0,     0,
     0,     0,     0,     0,     0,     0,   116,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,   121,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,    33,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,    50,    53,
    45,    44,     0,     0,    49,     0,   107,     0,    56,     0,
     0,     0,    66,     0,    29,     0,     0,    47,    48 };
mcitabelem mcipact[]={

  -247,-10000000,  -215,  -221,  -215,-10000000,-10000000,-10000000,-10000000,-10000000,
  -255,    -3,    -4,-10000000,   -24,  -223,  -224,  -225,-10000000,   -87,
  -217,     0,    -1,-10000000,    35,    34,-10000000,-10000000,   -40,-10000000,
-10000000,-10000000,-10000000,   -40,-10000000,  -262,  -233,   -26,-10000000,  -218,
  -218,  -227,  -228,-10000000,-10000000,-10000000,-10000000,-10000000,-10000000,    33,
-10000000,-10000000,-10000000,-10000000,   -44,  -249,-10000000,-10000000,     2,-10000000,
    -7,  -233,-10000000,  -230,-10000000,-10000000,-10000000,    23,-10000000,  -272,
   -44,    23,    22,    21,  -234,    32,  -263,    19,    18,-10000000,
-10000000,-10000000,-10000000,  -262,  -235,-10000000,   -34,-10000000,  -218,  -219,
-10000000,  -237,  -238,    29,-10000000,  -231,     1,-10000000,    -8,-10000000,
-10000000,  -231,-10000000,  -239,-10000000,  -257,-10000000,-10000000,-10000000,-10000000,
    20,    17,-10000000,-10000000,  -263,  -256,-10000000,-10000000,  -242,  -243,
  -248,-10000000,-10000000,-10000000,    28,    27,-10000000,-10000000 };
mcitabelem mcipgo[]={

     0,   107,    44,   106,    48,    71,   105,   104,   103,    77,
   102,   101,   100,    99,    98,    76,    97,    96,    93,    92,
    91,    52,    90,    89,    88,    87,    50,    86,    84,    51,
    83,    81,    49,    43,    80,    79,    78 };
mcitabelem mcir1[]={

     0,     1,     7,     8,     8,     9,    10,    11,    11,    11,
    11,    12,    12,    12,    12,    12,    12,    12,    12,    13,
    17,    18,    20,    20,    19,    19,    21,    22,    23,    23,
    14,    14,    15,    15,    15,    15,    15,    15,    15,    15,
    15,    15,    15,    25,    28,    24,    29,    30,    30,    31,
    31,    33,    33,    33,    32,    32,    26,     3,     3,     2,
    27,    34,    34,    16,    16,    35,    36,     5,     5,     4,
     4,     6,     6 };
mcitabelem mcir2[]={

     0,     4,     5,     2,     4,     5,     3,     0,     6,     7,
     4,     0,     7,     5,     5,     5,     5,     2,     7,     6,
     3,     4,     3,     5,     2,     4,     7,     3,     3,     5,
     3,     5,     3,     3,     2,     3,     3,     3,     9,     3,
     3,     3,     3,     7,     3,     7,     7,     2,     6,     2,
     6,     3,     2,     7,    11,     7,     6,     2,     6,     3,
     6,     3,     7,     7,     7,    17,    17,     3,     7,     3,
     7,     3,     5 };
mcitabelem mcichk[]={

-10000000,    -1,    -7,   262,    -8,    -9,   -10,   257,   259,    -9,
   -11,   -12,   256,   -13,   -14,   263,   265,   276,   -16,   258,
   -17,   -35,   -36,   275,   272,   271,    59,    59,    61,   258,
   259,   259,   259,   273,   -24,    91,   -18,   -20,   258,    58,
    58,    40,    40,   -15,   261,   260,   -25,   277,   278,   264,
   258,   -26,   -27,   259,   -28,    40,   268,   -15,   -30,   -32,
   269,   -19,   -21,   -22,   261,    59,   258,    -5,    -4,    -6,
   258,    -5,   259,   259,    40,   -29,    91,    -3,   -34,    -2,
   261,   260,    93,    44,    61,   -21,   -23,   260,    44,   274,
   -29,    44,    44,   261,   -26,    40,   -31,   -33,   269,   -32,
    41,    44,    41,    44,   -32,   261,    59,   260,    -4,   258,
   261,   261,    41,    93,    44,    61,    -2,   261,   266,    44,
    44,   -33,   270,   261,   261,   261,    41,    41 };
mcitabelem mcidef[]={

     0,    -2,     0,     0,     1,     3,     7,     6,     2,     4,
    -2,     0,     0,    10,     0,     0,     0,     0,    17,    30,
     0,     0,     0,    20,     0,     0,     8,     9,     0,    13,
    14,    15,    16,     0,    31,     0,     0,     0,    22,     0,
     0,     0,     0,    12,    32,    33,    34,    35,    36,    37,
    39,    40,    41,    42,     0,     0,    44,    18,     0,    47,
     0,    19,    24,     0,    27,    21,    23,    63,    67,    69,
    71,    64,     0,     0,     0,     0,     0,     0,     0,    57,
    61,    59,    45,     0,     0,    25,     0,    28,     0,     0,
    72,     0,     0,     0,    43,     0,     0,    49,    51,    52,
    56,     0,    60,     0,    48,    55,    26,    29,    68,    70,
     0,     0,    38,    46,     0,     0,    58,    62,     0,     0,
     0,    50,    53,    54,     0,     0,    65,    66 };
typedef struct { char *t_name; int t_val; } mcitoktype;
#ifndef YYDEBUG
#	define YYDEBUG	0	/* don't allow debugging */
#endif

#if YYDEBUG

char * mcireds[] =
{
	"-no such reduction-",
      "inpfile : title sections",
      "title : TITLE STRING",
      "sections : section",
      "sections : sections section",
      "section : sectionkey statements",
      "sectionkey : SECTION",
      "statements : /* empty */",
      "statements : statements statement ';'",
      "statements : statements error ';'",
      "statements : statements groundclasstable",
      "statement : /* empty */",
      "statement : variable '=' value",
      "statement : variable IDENTIFIER",
      "statement : CDFIMPORT STRING",
      "statement : CHEMFILE STRING",
      "statement : SUBDOMAIN STRING",
      "statement : outstatement",
      "statement : IDENTIFIER VARDEFINE value",
      "groundclasstable : groundclkey titleline datalines",
      "groundclkey : GROUNDCLASSES",
      "titleline : idlist ';'",
      "idlist : IDENTIFIER",
      "idlist : idlist IDENTIFIER",
      "datalines : dataline",
      "datalines : datalines dataline",
      "dataline : clidx numlist ';'",
      "clidx : INTNUM",
      "numlist : FLOATNUM",
      "numlist : numlist FLOATNUM",
      "variable : IDENTIFIER",
      "variable : IDENTIFIER range",
      "value : INTNUM",
      "value : FLOATNUM",
      "value : array",
      "value : IDATE",
      "value : DAYTIME",
      "value : AVG",
      "value : AVG '(' INTNUM ')'",
      "value : IDENTIFIER",
      "value : numbers",
      "value : intnumbers",
      "value : STRING",
      "array : arraykey dims numbers",
      "arraykey : ARRAY",
      "range : '[' rangespecs ']'",
      "dims : '[' dimspecs ']'",
      "rangespecs : rangespec",
      "rangespecs : rangespecs ',' rangespec",
      "dimspecs : dimspec",
      "dimspecs : dimspecs ',' dimspec",
      "dimspec : DIMENSION",
      "dimspec : rangespec",
      "dimspec : DIMENSION '=' GROUNDLEVEL",
      "rangespec : DIMENSION '=' INTNUM RANGE INTNUM",
      "rangespec : DIMENSION '=' INTNUM",
      "numbers : '(' floatlist ')'",
      "floatlist : floatnum",
      "floatlist : floatlist ',' floatnum",
      "floatnum : FLOATNUM",
      "intnumbers : '(' intlist ')'",
      "intlist : INTNUM",
      "intlist : intlist ',' INTNUM",
      "outstatement : cdfheader ':' outlist",
      "outstatement : asciiheader ':' outlist",
      "cdfheader : CDF '(' STRING ',' INTNUM ',' INTNUM ')'",
      "asciiheader : ASCII '(' STRING ',' INTNUM ',' INTNUM ')'",
      "outlist : outitem",
      "outlist : outlist ',' outitem",
      "outitem : printident",
      "outitem : printident AS IDENTIFIER",
      "printident : IDENTIFIER",
      "printident : IDENTIFIER dims",
};
mcitoktype mcitoks[] =
{
	"SECTION",	257,
	"IDENTIFIER",	258,
	"STRING",	259,
	"FLOATNUM",	260,
	"INTNUM",	261,
	"TITLE",	262,
	"CDFIMPORT",	263,
	"AVG",	264,
	"CHEMFILE",	265,
	"RANGE",	266,
	"NONSENSE",	267,
	"ARRAY",	268,
	"DIMENSION",	269,
	"GROUNDLEVEL",	270,
	"ASCII",	271,
	"CDF",	272,
	"VARDEFINE",	273,
	"AS",	274,
	"GROUNDCLASSES",	275,
	"SUBDOMAIN",	276,
	"IDATE",	277,
	"DAYTIME",	278,
	"';'",	59,
	"'='",	61,
	"'('",	40,
	"')'",	41,
	"'['",	91,
	"']'",	93,
	"','",	44,
	"':'",	58,
	"-unknown-",	-1	/* ends search */
};
#endif /* YYDEBUG */

/* @(#)27       1.7.1.4  src/bos/usr/ccs/bin/yacc/yaccpar, cmdlang, bos430, 9737A_430 11/28/95 13:48:59 */
/*
 * COMPONENT_NAME: (CMDLANG) Language Utilities
 *
 * FUNCTIONS: mciparse
 * ORIGINS: 3
 */
/*
** Skeleton parser driver for yacc output
*/

/*
** yacc user known macros and defines
*/
#ifdef YYSPLIT
#   define YYERROR      return(-2)
#else
#   define YYERROR      goto mcierrlab
#endif
#ifdef YACC_MSG
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif
#include <nl_types.h>
nl_catd mciusercatd;
#endif
#define YYACCEPT        return(0)
#define YYABORT         return(1)
#ifndef YACC_MSG
#define YYBACKUP( newtoken, newvalue )\
{\
        if ( mcichar >= 0 || ( mcir2[ mcitmp ] >> 1 ) != 1 )\
        {\
                mcierror( "syntax error - cannot backup" );\
                YYERROR;\
        }\
        mcichar = newtoken;\
        mcistate = *mcips;\
        mcilval = newvalue;\
        goto mcinewstate;\
}
#else
#define YYBACKUP( newtoken, newvalue )\
{\
        if ( mcichar >= 0 || ( mcir2[ mcitmp ] >> 1 ) != 1 )\
        {\
                mciusercatd=catopen("yacc_user.cat", NL_CAT_LOCALE);\
                mcierror(catgets(mciusercatd,1,1,"syntax error - cannot backup" ));\
                YYERROR;\
        }\
        mcichar = newtoken;\
        mcistate = *mcips;\
        mcilval = newvalue;\
        goto mcinewstate;\
}
#endif
#define YYRECOVERING()  (!!mcierrflag)
#ifndef YYDEBUG
#       define YYDEBUG  1       /* make debugging available */
#endif

/*
** user known globals
*/
int mcidebug;                    /* set to 1 to get debugging */

/*
** driver internal defines
*/
#define YYFLAG          (-10000000)

#ifdef YYSPLIT
#   define YYSCODE { \
                        extern int (*_mcif[])(); \
                        register int mciret; \
                        if (_mcif[mcitmp]) \
                            if ((mciret=(*_mcif[mcitmp])()) == -2) \
                                    goto mcierrlab; \
                                else if (mciret>=0) return(mciret); \
                   }
#endif

/*
** global variables used by the parser
*/
YYSTYPE mciv[ YYMAXDEPTH ];      /* value stack */
int mcis[ YYMAXDEPTH ];          /* state stack */

YYSTYPE *mcipv;                  /* top of value stack */
YYSTYPE *mcipvt;                 /* top of value stack for $vars */
int *mcips;                      /* top of state stack */

int mcistate;                    /* current state */
int mcitmp;                      /* extra var (lasts between blocks) */

int mcinerrs;                    /* number of errors */
int mcierrflag;                  /* error recovery flag */
int mcichar;                     /* current input token number */

#ifdef __cplusplus
 #ifdef _CPP_IOSTREAMS
  #include <iostream.h>
  extern void mcierror (char *); /* error message routine -- iostream version */
 #else
  #include <stdio.h>
  extern "C" void mcierror (char *); /* error message routine -- stdio version */
 #endif /* _CPP_IOSTREAMS */
 extern "C" int mcilex(void);        /* return the next token */
#endif /* __cplusplus */


/*
** mciparse - return 0 if worked, 1 if syntax error not recovered from
*/
#ifdef __cplusplus
extern "C"
#endif /* __cplusplus */
int
mciparse()
{
        /*
        ** Initialize externals - mciparse may be called more than once
        */
        mcipv = &mciv[-1];
        mcips = &mcis[-1];
        mcistate = 0;
        mcitmp = 0;
        mcinerrs = 0;
        mcierrflag = 0;
        mcichar = -1;
#ifdef YACC_MSG
        mciusercatd=catopen("yacc_user.cat", NL_CAT_LOCALE);
#endif
        goto mcistack;
        {
                register YYSTYPE *mci_pv;        /* top of value stack */
                register int *mci_ps;            /* top of state stack */
                register int mci_state;          /* current state */
                register int  mci_n;             /* internal state number info */

                /*
                ** get globals into registers.
                ** branch to here only if YYBACKUP was called.
                */
        mcinewstate:
                mci_pv = mcipv;
                mci_ps = mcips;
                mci_state = mcistate;
                goto mci_newstate;

                /*
                ** get globals into registers.
                ** either we just started, or we just finished a reduction
                */
        mcistack:
                mci_pv = mcipv;
                mci_ps = mcips;
                mci_state = mcistate;

                /*
                ** top of for (;;) loop while no reductions done
                */
        mci_stack:
                /*
                ** put a state and value onto the stacks
                */
#if YYDEBUG
                /*
                ** if debugging, look up token value in list of value vs.
                ** name pairs.  0 and negative (-1) are special values.
                ** Note: linear search is used since time is not a real
                ** consideration while debugging.
                */
                if ( mcidebug )
                {
                        register int mci_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                        cout << "State " << mci_state << " token ";
                        if ( mcichar == 0 )
                                cout << "end-of-file" << endl;
                        else if ( mcichar < 0 )
                                cout << "-none-" << endl;
#else
                        printf( "State %d, token ", mci_state );
                        if ( mcichar == 0 )
                                printf( "end-of-file\n" );
                        else if ( mcichar < 0 )
                                printf( "-none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        else
                        {
                                for ( mci_i = 0; mcitoks[mci_i].t_val >= 0;
                                        mci_i++ )
                                {
                                        if ( mcitoks[mci_i].t_val == mcichar )
                                                break;
                                }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                cout << mcitoks[mci_i].t_name << endl;
#else
                                printf( "%s\n", mcitoks[mci_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        }
                }
#endif /* YYDEBUG */
                if ( ++mci_ps >= &mcis[ YYMAXDEPTH ] )    /* room on stack? */
                {
#ifndef YACC_MSG
                        mcierror( "yacc stack overflow" );
#else
                        mcierror(catgets(mciusercatd,1,2,"yacc stack overflow" ));
#endif
                        YYABORT;
                }
                *mci_ps = mci_state;
                *++mci_pv = mcival;

                /*
                ** we have a new state - find out what to do
                */
        mci_newstate:
                if ( ( mci_n = mcipact[ mci_state ] ) <= YYFLAG )
                        goto mcidefault;         /* simple state */
#if YYDEBUG
                /*
                ** if debugging, need to mark whether new token grabbed
                */
                mcitmp = mcichar < 0;
#endif
                if ( ( mcichar < 0 ) && ( ( mcichar = mcilex() ) < 0 ) )
                        mcichar = 0;             /* reached EOF */
#if YYDEBUG
                if ( mcidebug && mcitmp )
                {
                        register int mci_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                        cout << "Received token " << endl;
                        if ( mcichar == 0 )
                                cout << "end-of-file" << endl;
                        else if ( mcichar < 0 )
                                cout << "-none-" << endl;
#else
                        printf( "Received token " );
                        if ( mcichar == 0 )
                                printf( "end-of-file\n" );
                        else if ( mcichar < 0 )
                                printf( "-none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        else
                        {
                                for ( mci_i = 0; mcitoks[mci_i].t_val >= 0;
                                        mci_i++ )
                                {
                                        if ( mcitoks[mci_i].t_val == mcichar )
                                                break;
                                }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                cout << mcitoks[mci_i].t_name << endl;
#else
                                printf( "%s\n", mcitoks[mci_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        }
                }
#endif /* YYDEBUG */
                if ( ( ( mci_n += mcichar ) < 0 ) || ( mci_n >= YYLAST ) )
                        goto mcidefault;
                if ( mcichk[ mci_n = mciact[ mci_n ] ] == mcichar )  /*valid shift*/
                {
                        mcichar = -1;
                        mcival = mcilval;
                        mci_state = mci_n;
                        if ( mcierrflag > 0 )
                                mcierrflag--;
                        goto mci_stack;
                }

        mcidefault:
                if ( ( mci_n = mcidef[ mci_state ] ) == -2 )
                {
#if YYDEBUG
                        mcitmp = mcichar < 0;
#endif
                        if ( ( mcichar < 0 ) && ( ( mcichar = mcilex() ) < 0 ) )
                                mcichar = 0;             /* reached EOF */
#if YYDEBUG
                        if ( mcidebug && mcitmp )
                        {
                                register int mci_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                cout << "Received token " << endl;
                                if ( mcichar == 0 )
                                        cout << "end-of-file" << endl;
                                else if ( mcichar < 0 )
                                        cout << "-none-" << endl;
#else
                                printf( "Received token " );
                                if ( mcichar == 0 )
                                        printf( "end-of-file\n" );
                                else if ( mcichar < 0 )
                                        printf( "-none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                else
                                {
                                        for ( mci_i = 0;
                                                mcitoks[mci_i].t_val >= 0;
                                                mci_i++ )
                                        {
                                                if ( mcitoks[mci_i].t_val
                                                        == mcichar )
                                                {
                                                        break;
                                                }
                                        }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                        cout << mcitoks[mci_i].t_name << endl;
#else
                                        printf( "%s\n", mcitoks[mci_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                }
                        }
#endif /* YYDEBUG */
                        /*
                        ** look through exception table
                        */
                        {
                                register int *mcixi = mciexca;

                                while ( ( *mcixi != -1 ) ||
                                        ( mcixi[1] != mci_state ) )
                                {
                                        mcixi += 2;
                                }
                                while ( ( *(mcixi += 2) >= 0 ) &&
                                        ( *mcixi != mcichar ) )
                                        ;
                                if ( ( mci_n = mcixi[1] ) < 0 )
                                        YYACCEPT;
                        }
                }

                /*
                ** check for syntax error
                */
                if ( mci_n == 0 )        /* have an error */
                {
                        /* no worry about speed here! */
                        switch ( mcierrflag )
                        {
                        case 0:         /* new error */
#ifndef YACC_MSG
                                mcierror( "syntax error" );
#else
                                mcierror(catgets(mciusercatd,1,3,"syntax error" ));
#endif
                                goto skip_init;
                        mcierrlab:
                                /*
                                ** get globals into registers.
                                ** we have a user generated syntax type error
                                */
                                mci_pv = mcipv;
                                mci_ps = mcips;
                                mci_state = mcistate;
                                mcinerrs++;
                        skip_init:
                        case 1:
                        case 2:         /* incompletely recovered error */
                                        /* try again... */
                                mcierrflag = 3;
                                /*
                                ** find state where "error" is a legal
                                ** shift action
                                */
                                while ( mci_ps >= mcis )
                                {
                                        mci_n = mcipact[ *mci_ps ] + YYERRCODE;
                                        if ( mci_n >= 0 && mci_n < YYLAST &&
                                                mcichk[mciact[mci_n]] == YYERRCODE)                                        {
                                                /*
                                                ** simulate shift of "error"
                                                */
                                                mci_state = mciact[ mci_n ];
                                                goto mci_stack;
                                        }
                                        /*
                                        ** current state has no shift on
                                        ** "error", pop stack
                                        */
#if YYDEBUG
                                        if ( mcidebug )
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                            cout << "Error recovery pops state "
                                                 << (*mci_ps)
                                                 << ", uncovers state "
                                                 << mci_ps[-1] << endl;
#else
#       define _POP_ "Error recovery pops state %d, uncovers state %d\n"
                                                printf( _POP_, *mci_ps,
                                                        mci_ps[-1] );
#       undef _POP_
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
#endif
                                        mci_ps--;
                                        mci_pv--;
                                }
                                /*
                                ** there is no state on stack with "error" as
                                ** a valid shift.  give up.
                                */
                                YYABORT;
                        case 3:         /* no shift yet; eat a token */
#if YYDEBUG
                                /*
                                ** if debugging, look up token in list of
                                ** pairs.  0 and negative shouldn't occur,
                                ** but since timing doesn't matter when
                                ** debugging, it doesn't hurt to leave the
                                ** tests here.
                                */
                                if ( mcidebug )
                                {
                                        register int mci_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                        cout << "Error recovery discards ";
                                        if ( mcichar == 0 )
                                            cout << "token end-of-file" << endl;
                                        else if ( mcichar < 0 )
                                            cout << "token -none-" << endl;
#else
                                        printf( "Error recovery discards " );
                                        if ( mcichar == 0 )
                                                printf( "token end-of-file\n" );
                                        else if ( mcichar < 0 )
                                                printf( "token -none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                        else
                                        {
                                                for ( mci_i = 0;
                                                        mcitoks[mci_i].t_val >= 0;
                                                        mci_i++ )
                                                {
                                                        if ( mcitoks[mci_i].t_val
                                                                == mcichar )
                                                        {
                                                                break;
                                                        }
                                                }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                                cout << "token " <<
                                                    mcitoks[mci_i].t_name <<
                                                    endl;
#else
                                                printf( "token %s\n",
                                                        mcitoks[mci_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                        }
                                }
#endif /* YYDEBUG */
                                if ( mcichar == 0 )      /* reached EOF. quit */
                                        YYABORT;
                                mcichar = -1;
                                goto mci_newstate;
                        }
                }/* end if ( mci_n == 0 ) */
                /*
                ** reduction by production mci_n
                ** put stack tops, etc. so things right after switch
                */
#if YYDEBUG
                /*
                ** if debugging, print the string that is the user's
                ** specification of the reduction which is just about
                ** to be done.
                */
                if ( mcidebug )
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                        cout << "Reduce by (" << mci_n << ") \"" <<
                            mcireds[ mci_n ] << "\"\n";
#else
                        printf( "Reduce by (%d) \"%s\"\n",
                                mci_n, mcireds[ mci_n ] );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
#endif
                mcitmp = mci_n;                   /* value to switch over */
                mcipvt = mci_pv;                  /* $vars top of value stack */
                /*
                ** Look in goto table for next state
                ** Sorry about using mci_state here as temporary
                ** register variable, but why not, if it works...
                ** If mcir2[ mci_n ] doesn't have the low order bit
                ** set, then there is no action to be done for
                ** this reduction.  So, no saving & unsaving of
                ** registers done.  The only difference between the
                ** code just after the if and the body of the if is
                ** the goto mci_stack in the body.  This way the test
                ** can be made before the choice of what to do is needed.
                */
                {
                        /* length of production doubled with extra bit */
                        register int mci_len = mcir2[ mci_n ];

                        if ( !( mci_len & 01 ) )
                        {
                                mci_len >>= 1;
                                mcival = ( mci_pv -= mci_len )[1]; /* $$ = $1 */
                                mci_state = mcipgo[ mci_n = mcir1[ mci_n ] ] +
                                        *( mci_ps -= mci_len ) + 1;
                                if ( mci_state >= YYLAST ||
                                        mcichk[ mci_state =
                                        mciact[ mci_state ] ] != -mci_n )
                                {
                                        mci_state = mciact[ mcipgo[ mci_n ] ];
                                }
                                goto mci_stack;
                        }
                        mci_len >>= 1;
                        mcival = ( mci_pv -= mci_len )[1]; /* $$ = $1 */
                        mci_state = mcipgo[ mci_n = mcir1[ mci_n ] ] +
                                *( mci_ps -= mci_len ) + 1;
                        if ( mci_state >= YYLAST ||
                                mcichk[ mci_state = mciact[ mci_state ] ] != -mci_n )
                        {
                                mci_state = mciact[ mcipgo[ mci_n ] ];
                        }
                }
                                        /* save until reenter driver code */
                mcistate = mci_state;
                mcips = mci_ps;
                mcipv = mci_pv;
        }
        /*
        ** code supplied by user is placed in this switch
        */

                switch(mcitmp){

case 2:
# line 48 "mcsyntax.y"
{
        ShowTitle(mcipvt[-0].item);
      } /*NOTREACHED*/ break;
case 5:
# line 56 "mcsyntax.y"
{
          SCALL(EndOfSection());
        } /*NOTREACHED*/ break;
case 6:
# line 61 "mcsyntax.y"
{
	  SCALL(ChangeSectionTo(mcipvt[-0].item));
	} /*NOTREACHED*/ break;
case 9:
# line 68 "mcsyntax.y"
{
	     inputerror = 1;
	     mcierrok;
	   } /*NOTREACHED*/ break;
case 12:
# line 76 "mcsyntax.y"
{
            SCALL(AssignVariable());
            InitArrayVar();
          } /*NOTREACHED*/ break;
case 13:
# line 81 "mcsyntax.y"
{
            SCALLD(FindOption(mcipvt[-0].item));
            SCALL(AssignVariable());
            InitArrayVar();
          } /*NOTREACHED*/ break;
case 14:
# line 87 "mcsyntax.y"
{
            SCALL(ReadCDFFile(mcipvt[-0].item, actualsection));
          } /*NOTREACHED*/ break;
case 15:
# line 91 "mcsyntax.y"
{
            SCALL(ReadChemFile(mcipvt[-0].item));
            SCALL(TestSpecial());
          } /*NOTREACHED*/ break;
case 16:
# line 96 "mcsyntax.y"
{
            SCALL(DefineNestDomain(mcipvt[-0].item, actualsection));
          } /*NOTREACHED*/ break;
case 18:
# line 101 "mcsyntax.y"
{
            SCALL(DefineAVariable(mcipvt[-2].item));
          } /*NOTREACHED*/ break;
case 20:
# line 108 "mcsyntax.y"
{
	  SCALL(InitializeGroundClasses());
	} /*NOTREACHED*/ break;
case 22:
# line 115 "mcsyntax.y"
{
	  SCALLD(AddClassVar(mcipvt[-0].item));
	} /*NOTREACHED*/ break;
case 23:
# line 119 "mcsyntax.y"
{
	  SCALLD(AddClassVar(mcipvt[-0].item));
	} /*NOTREACHED*/ break;
case 26:
# line 127 "mcsyntax.y"
{
	  SCALL(EndOfClassLine());
	} /*NOTREACHED*/ break;
case 27:
# line 132 "mcsyntax.y"
{
	  SCALLD(AddClassIdx(mcipvt[-0].ival));
	} /*NOTREACHED*/ break;
case 28:
# line 137 "mcsyntax.y"
{
	  SCALL(AddClassVal(mcipvt[-0].dval));
	} /*NOTREACHED*/ break;
case 29:
# line 141 "mcsyntax.y"
{
          SCALL(AddClassVal(mcipvt[-0].dval));
        } /*NOTREACHED*/ break;
case 30:
# line 146 "mcsyntax.y"
{
	   SCALLD(SetVariable(mcipvt[-0].item));
	 } /*NOTREACHED*/ break;
case 31:
# line 150 "mcsyntax.y"
{
           SCALLD(SetVariable(mcipvt[-1].item));
           SCALLD(CopyRangeToVar());
         } /*NOTREACHED*/ break;
case 32:
# line 156 "mcsyntax.y"
{
        val.type = INT_NUM;
        val.v.ival = mcipvt[-0].ival;
      } /*NOTREACHED*/ break;
case 33:
# line 161 "mcsyntax.y"
{
        val.type = DOUBLE_NUM;
        val.v.dval = mcipvt[-0].dval;
      } /*NOTREACHED*/ break;
case 35:
# line 167 "mcsyntax.y"
{
        val.type = DATE_VAL;
        val.v.date = mcipvt[-0].date;
      } /*NOTREACHED*/ break;
case 36:
# line 172 "mcsyntax.y"
{
        val.type = TIME_VAL;
        val.v.time = mcipvt[-0].time;
      } /*NOTREACHED*/ break;
case 37:
# line 177 "mcsyntax.y"
{
        val.type = AVG_VAL;
        val.v.ival = nz;
      } /*NOTREACHED*/ break;
case 38:
# line 182 "mcsyntax.y"
{
        val.type = AVG_VAL;
        val.v.ival = mcipvt[-1].ival;
      } /*NOTREACHED*/ break;
case 39:
# line 187 "mcsyntax.y"
{
        SCALLD(CopyVarToVal(mcipvt[-0].item));
      } /*NOTREACHED*/ break;
case 40:
# line 191 "mcsyntax.y"
{
        val.type = VECTOR_VAL;
      } /*NOTREACHED*/ break;
case 41:
# line 195 "mcsyntax.y"
{
        val.type = INT_VECTOR_VAL;
      } /*NOTREACHED*/ break;
case 42:
# line 199 "mcsyntax.y"
{
        val.type = TEXT_VAL;
        strcpy(val.v.txt, mcipvt[-0].item);
      } /*NOTREACHED*/ break;
case 43:
# line 205 "mcsyntax.y"
{
	  if (val.v.a.nval < val.v.a.maxval)  {
	    InputError(ERROR, "There are only %i initialisers. You need %i.",
	       val.v.a.nval, val.v.a.maxval);
	    inputerror = 1;
	  }  /* of if */
	} /*NOTREACHED*/ break;
case 44:
# line 214 "mcsyntax.y"
{
           InitArrayVar();
         } /*NOTREACHED*/ break;
case 45:
# line 219 "mcsyntax.y"
{
        AllocateArraySpace();
      } /*NOTREACHED*/ break;
case 46:
# line 224 "mcsyntax.y"
{
        AllocateArraySpace();
      } /*NOTREACHED*/ break;
case 51:
# line 236 "mcsyntax.y"
{
	  SCALL(SetDimension(mcipvt[-0].ival));
	} /*NOTREACHED*/ break;
case 53:
# line 241 "mcsyntax.y"
{
	  SCALL(SetDimension(mcipvt[-2].ival));
	  SCALL(SetRange(mcipvt[-2].ival, -1, -1));
	} /*NOTREACHED*/ break;
case 54:
# line 247 "mcsyntax.y"
{
	    SCALL(SetDimension(mcipvt[-4].ival));
	    SCALL(SetRange(mcipvt[-4].ival, mcipvt[-2].ival, mcipvt[-0].ival));
	  } /*NOTREACHED*/ break;
case 55:
# line 252 "mcsyntax.y"
{
	    SCALL(SetDimension(mcipvt[-2].ival));
	    SCALL(SetRange(mcipvt[-2].ival, mcipvt[-0].ival, mcipvt[-0].ival));
	  } /*NOTREACHED*/ break;
case 59:
# line 263 "mcsyntax.y"
{
	   SCALLD(AddArrayVal(mcipvt[-0].dval));
	 } /*NOTREACHED*/ break;
case 61:
# line 270 "mcsyntax.y"
{
	    AddIntArrayVal(mcipvt[-0].ival);
	  } /*NOTREACHED*/ break;
case 62:
# line 274 "mcsyntax.y"
{
	    AddIntArrayVal(mcipvt[-0].ival);
	  } /*NOTREACHED*/ break;
case 63:
# line 279 "mcsyntax.y"
{
	       AttachOutlist(mcipvt[-0].outvar);
	     } /*NOTREACHED*/ break;
case 64:
# line 283 "mcsyntax.y"
{
	       AttachOutlist(mcipvt[-0].outvar);
	     } /*NOTREACHED*/ break;
case 65:
# line 288 "mcsyntax.y"
{
	    SCALLD(PrepareOutput(CDF_PRINT, mcipvt[-5].item, mcipvt[-3].ival, mcipvt[-1].ival));
	  } /*NOTREACHED*/ break;
case 66:
# line 294 "mcsyntax.y"
{
	      SCALLD(PrepareOutput(ASCII_PRINT, mcipvt[-5].item, mcipvt[-3].ival, mcipvt[-1].ival));
	    } /*NOTREACHED*/ break;
case 67:
# line 299 "mcsyntax.y"
{
	  mcival.outvar = mcipvt[-0].outvar;
	} /*NOTREACHED*/ break;
case 68:
# line 303 "mcsyntax.y"
{
	  mcival.outvar = AppendOutVar(mcipvt[-2].outvar, mcipvt[-0].outvar);
	  if (!mcival.outvar)  {
	    inputerror = 1;
	    YYERROR;
	  }  /* of if */
	} /*NOTREACHED*/ break;
case 69:
# line 312 "mcsyntax.y"
{
	  if (!(mcival.outvar = OutItem(mcipvt[-0].item, mcipvt[-0].item)))  {
	    inputerror = 1;
	    YYERROR;
	  }   /* of if */
	  InitArrayVar();
	} /*NOTREACHED*/ break;
case 70:
# line 320 "mcsyntax.y"
{
	  if (!(mcival.outvar = OutItem(mcipvt[-2].item, mcipvt[-0].item)))  {
	    inputerror = 1;
	    YYERROR;
	  }   /* of if */
	  InitArrayVar();
	} /*NOTREACHED*/ break;
case 71:
# line 329 "mcsyntax.y"
{
	  strcpy(mcival.item, mcipvt[-0].item);
	} /*NOTREACHED*/ break;
case 72:
# line 333 "mcsyntax.y"
{
	  strcpy(mcival.item, mcipvt[-1].item);
	} /*NOTREACHED*/ break;
}


        goto mcistack;           /* reset registers in driver code */
}

# line 336 "mcsyntax.y"
