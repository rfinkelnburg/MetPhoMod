
# line 2 "mchemsyntax.y"

#include <stdio.h>
#include <string.h>
#include "mcglobal.h"
#include "mchemparse.h"


# line 12 "mchemsyntax.y"
typedef union   {
  int ival;
  double dval;
  double *dref;
  char *item;
} YYSTYPE;
# define INTNUM 257
# define DELIMITER 258
# define NONSENSE 259
# define PHOTOKIN 260
# define THERMKIN 261
# define TERM2KIN 262
# define TROEKIN 263
# define PHOTOKIN3 264
# define TROEQUILKIN 265
# define SPECIALKIN 266
# define UNIT 267
# define A_UNIT 268
# define REACTIONS 269
# define SUBSTANCES 270
# define END 271
# define REF 272
# define FLOATNUM 273
# define SUBSTNAME 274
# define STRING 275
# define A_REF 276
#define chiclearin chichar = -1
#define chierrok chierrflag = 0
extern int chichar;
extern int chierrflag;
#ifndef YYMAXDEPTH
#define YYMAXDEPTH 150
#endif
YYSTYPE chilval, chival;
typedef int chitabelem;
#include <stdio.h>
# define YYERRCODE 256
chitabelem chiexca[] ={
	-1, 1,
	0, -1,
	-2, 0,
	};
# define YYNPROD 37
# define YYLAST 159
chitabelem chiact[]={

     9,    12,    40,    25,    78,    24,     9,    12,    15,     9,
    12,     7,    34,    12,    51,    57,     8,     7,    11,    34,
     7,    16,     8,    23,    11,     8,    15,    11,    33,    31,
    11,    69,     2,     3,    68,    33,    43,    45,    46,    47,
    44,    48,    49,    37,    67,    58,    27,    22,    32,    19,
    56,    39,    30,    54,    52,   109,    55,    53,   107,   103,
   102,    97,    96,    93,    84,    83,    82,    18,    81,    80,
    79,    50,     5,   111,   106,   101,    95,    17,    94,    92,
    85,    65,    64,    63,    62,    61,    60,    59,    21,    20,
     4,    29,    14,    38,    42,    10,    28,     6,    36,    13,
     1,     0,     0,    66,     0,    71,    26,    70,    72,    73,
    74,    75,    76,    77,    17,    35,    41,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,    86,    87,
    88,    89,    90,    91,     0,     0,     0,     0,     0,     0,
     0,     0,    98,     0,     0,    99,   100,     0,     0,     0,
     0,   104,   105,     0,     0,     0,   108,     0,   110 };
chitabelem chipact[]={

  -237,-10000000,  -247,  -248,  -250,-10000000,     6,    49,    48,  -211,
-10000000,-10000000,  -251,  -266,-10000000,  -212,-10000000,-10000000,  -245,  -244,
  -225,  -274,-10000000,-10000000,-10000000,  -247,-10000000,-10000000,  -224,    28,
-10000000,-10000000,  -260,-10000000,-10000000,-10000000,    13,-10000000,    12,-10000000,
   -11,  -256,  -213,    47,    46,    45,    44,    43,    42,    41,
  -245,-10000000,  -214,  -234,  -227,  -274,  -238,-10000000,-10000000,  -238,
  -238,  -238,  -238,  -238,  -238,  -271,-10000000,-10000000,-10000000,-10000000,
-10000000,-10000000,    26,    25,    24,    22,    21,    20,    39,  -238,
  -238,  -238,  -238,  -238,  -238,-10000000,    38,    19,    37,    35,
    18,    17,-10000000,  -238,-10000000,-10000000,  -238,  -238,    34,    16,
    15,-10000000,  -238,  -238,    33,    14,-10000000,  -238,    11,  -238,
    32,-10000000 };
chitabelem chipgo[]={

     0,   100,    48,    90,    99,    92,    72,    98,    97,    96,
    94,    93,    95,    91,    52,    51 };
chitabelem chir1[]={

     0,     1,     1,     1,     3,     3,     6,     6,     6,     6,
     8,     8,     9,     9,    13,    13,    12,    12,     4,     4,
     5,     2,     2,    14,    14,    10,    10,    10,    10,    10,
    10,    10,     7,     7,    11,    11,    15 };
chitabelem chir2[]={

     0,     6,     6,    10,     2,     5,    11,    10,    11,     5,
     2,     6,     0,     2,     2,     6,     3,     5,     2,     4,
     5,     2,     3,     3,     5,    13,    17,    13,    13,    21,
    29,     9,     2,     6,     2,     6,     7 };
chitabelem chichk[]={

-10000000,    -1,   269,   270,    -3,    -6,    -8,   267,   272,   256,
   -12,   274,   257,    -4,    -5,   274,   271,    -6,    61,    43,
    40,    40,   258,   274,   271,   269,    -5,   258,    -9,   -13,
   -14,   274,    -2,   273,   257,   -12,    -7,   268,   -11,   -15,
   276,    -3,   -10,   260,   264,   261,   262,   263,   265,   266,
    43,   274,    41,    44,    41,    44,    61,   271,   258,    40,
    40,    40,    40,    40,    40,    40,   -14,   258,   268,   258,
   -15,    -2,    -2,    -2,    -2,    -2,    -2,    -2,   275,    44,
    44,    44,    44,    44,    44,    41,    -2,    -2,    -2,    -2,
    -2,    -2,    41,    44,    41,    41,    44,    44,    -2,    -2,
    -2,    41,    44,    44,    -2,    -2,    41,    44,    -2,    44,
    -2,    41 };
chitabelem chidef[]={

     0,    -2,     0,     0,     0,     4,     0,     0,     0,     0,
    10,    16,     0,     0,    18,     0,     1,     5,    12,     0,
     0,     0,     9,    17,     2,     0,    19,    20,     0,    13,
    14,    23,     0,    21,    22,    11,     0,    32,     0,    34,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,    24,     0,     0,     0,     0,     0,     3,     6,     0,
     0,     0,     0,     0,     0,     0,    15,     7,    33,     8,
    35,    36,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,    31,     0,     0,     0,     0,
     0,     0,    25,     0,    27,    28,     0,     0,     0,     0,
     0,    26,     0,     0,     0,     0,    29,     0,     0,     0,
     0,    30 };
typedef struct { char *t_name; int t_val; } chitoktype;
#ifndef YYDEBUG
#	define YYDEBUG	0	/* don't allow debugging */
#endif

#if YYDEBUG

char * chireds[] =
{
	"-no such reduction-",
      "inpfile : REACTIONS reactions END",
      "inpfile : SUBSTANCES substances END",
      "inpfile : SUBSTANCES substances REACTIONS reactions END",
      "reactions : fullreaction",
      "reactions : reactions fullreaction",
      "fullreaction : educts '=' products kindata DELIMITER",
      "fullreaction : UNIT '(' units ')' DELIMITER",
      "fullreaction : REF '(' refs ')' DELIMITER",
      "fullreaction : error DELIMITER",
      "educts : educt",
      "educts : educts '+' educt",
      "products : /* empty */",
      "products : productitems",
      "productitems : product",
      "productitems : productitems '+' product",
      "educt : SUBSTNAME",
      "educt : INTNUM SUBSTNAME",
      "substances : substance",
      "substances : substances substance",
      "substance : SUBSTNAME DELIMITER",
      "floatnum : FLOATNUM",
      "floatnum : INTNUM",
      "product : SUBSTNAME",
      "product : floatnum SUBSTNAME",
      "kindata : PHOTOKIN '(' floatnum ',' floatnum ')'",
      "kindata : PHOTOKIN3 '(' floatnum ',' floatnum ',' floatnum ')'",
      "kindata : THERMKIN '(' floatnum ',' floatnum ')'",
      "kindata : TERM2KIN '(' floatnum ',' floatnum ')'",
      "kindata : TROEKIN '(' floatnum ',' floatnum ',' floatnum ',' floatnum ')'",
      "kindata : TROEQUILKIN '(' floatnum ',' floatnum ',' floatnum ',' floatnum ',' floatnum ',' floatnum ')'",
      "kindata : SPECIALKIN '(' STRING ')'",
      "units : A_UNIT",
      "units : units ',' A_UNIT",
      "refs : ref",
      "refs : refs ',' ref",
      "ref : A_REF '=' floatnum",
};
chitoktype chitoks[] =
{
	"INTNUM",	257,
	"DELIMITER",	258,
	"NONSENSE",	259,
	"PHOTOKIN",	260,
	"THERMKIN",	261,
	"TERM2KIN",	262,
	"TROEKIN",	263,
	"PHOTOKIN3",	264,
	"TROEQUILKIN",	265,
	"SPECIALKIN",	266,
	"UNIT",	267,
	"A_UNIT",	268,
	"REACTIONS",	269,
	"SUBSTANCES",	270,
	"END",	271,
	"REF",	272,
	"FLOATNUM",	273,
	"SUBSTNAME",	274,
	"STRING",	275,
	"A_REF",	276,
	"'='",	61,
	"'('",	40,
	"')'",	41,
	"'+'",	43,
	"','",	44,
	"-unknown-",	-1	/* ends search */
};
#endif /* YYDEBUG */

/* @(#)27       1.7.1.4  src/bos/usr/ccs/bin/yacc/yaccpar, cmdlang, bos430, 9737A_430 11/28/95 13:48:59 */
/*
 * COMPONENT_NAME: (CMDLANG) Language Utilities
 *
 * FUNCTIONS: chiparse
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
#   define YYERROR      goto chierrlab
#endif
#ifdef YACC_MSG
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif
#include <nl_types.h>
nl_catd chiusercatd;
#endif
#define YYACCEPT        return(0)
#define YYABORT         return(1)
#ifndef YACC_MSG
#define YYBACKUP( newtoken, newvalue )\
{\
        if ( chichar >= 0 || ( chir2[ chitmp ] >> 1 ) != 1 )\
        {\
                chierror( "syntax error - cannot backup" );\
                YYERROR;\
        }\
        chichar = newtoken;\
        chistate = *chips;\
        chilval = newvalue;\
        goto chinewstate;\
}
#else
#define YYBACKUP( newtoken, newvalue )\
{\
        if ( chichar >= 0 || ( chir2[ chitmp ] >> 1 ) != 1 )\
        {\
                chiusercatd=catopen("yacc_user.cat", NL_CAT_LOCALE);\
                chierror(catgets(chiusercatd,1,1,"syntax error - cannot backup" ));\
                YYERROR;\
        }\
        chichar = newtoken;\
        chistate = *chips;\
        chilval = newvalue;\
        goto chinewstate;\
}
#endif
#define YYRECOVERING()  (!!chierrflag)
#ifndef YYDEBUG
#       define YYDEBUG  1       /* make debugging available */
#endif

/*
** user known globals
*/
int chidebug;                    /* set to 1 to get debugging */

/*
** driver internal defines
*/
#define YYFLAG          (-10000000)

#ifdef YYSPLIT
#   define YYSCODE { \
                        extern int (*_chif[])(); \
                        register int chiret; \
                        if (_chif[chitmp]) \
                            if ((chiret=(*_chif[chitmp])()) == -2) \
                                    goto chierrlab; \
                                else if (chiret>=0) return(chiret); \
                   }
#endif

/*
** global variables used by the parser
*/
YYSTYPE chiv[ YYMAXDEPTH ];      /* value stack */
int chis[ YYMAXDEPTH ];          /* state stack */

YYSTYPE *chipv;                  /* top of value stack */
YYSTYPE *chipvt;                 /* top of value stack for $vars */
int *chips;                      /* top of state stack */

int chistate;                    /* current state */
int chitmp;                      /* extra var (lasts between blocks) */

int chinerrs;                    /* number of errors */
int chierrflag;                  /* error recovery flag */
int chichar;                     /* current input token number */

#ifdef __cplusplus
 #ifdef _CPP_IOSTREAMS
  #include <iostream.h>
  extern void chierror (char *); /* error message routine -- iostream version */
 #else
  #include <stdio.h>
  extern "C" void chierror (char *); /* error message routine -- stdio version */
 #endif /* _CPP_IOSTREAMS */
 extern "C" int chilex(void);        /* return the next token */
#endif /* __cplusplus */


/*
** chiparse - return 0 if worked, 1 if syntax error not recovered from
*/
#ifdef __cplusplus
extern "C"
#endif /* __cplusplus */
int
chiparse()
{
        /*
        ** Initialize externals - chiparse may be called more than once
        */
        chipv = &chiv[-1];
        chips = &chis[-1];
        chistate = 0;
        chitmp = 0;
        chinerrs = 0;
        chierrflag = 0;
        chichar = -1;
#ifdef YACC_MSG
        chiusercatd=catopen("yacc_user.cat", NL_CAT_LOCALE);
#endif
        goto chistack;
        {
                register YYSTYPE *chi_pv;        /* top of value stack */
                register int *chi_ps;            /* top of state stack */
                register int chi_state;          /* current state */
                register int  chi_n;             /* internal state number info */

                /*
                ** get globals into registers.
                ** branch to here only if YYBACKUP was called.
                */
        chinewstate:
                chi_pv = chipv;
                chi_ps = chips;
                chi_state = chistate;
                goto chi_newstate;

                /*
                ** get globals into registers.
                ** either we just started, or we just finished a reduction
                */
        chistack:
                chi_pv = chipv;
                chi_ps = chips;
                chi_state = chistate;

                /*
                ** top of for (;;) loop while no reductions done
                */
        chi_stack:
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
                if ( chidebug )
                {
                        register int chi_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                        cout << "State " << chi_state << " token ";
                        if ( chichar == 0 )
                                cout << "end-of-file" << endl;
                        else if ( chichar < 0 )
                                cout << "-none-" << endl;
#else
                        printf( "State %d, token ", chi_state );
                        if ( chichar == 0 )
                                printf( "end-of-file\n" );
                        else if ( chichar < 0 )
                                printf( "-none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        else
                        {
                                for ( chi_i = 0; chitoks[chi_i].t_val >= 0;
                                        chi_i++ )
                                {
                                        if ( chitoks[chi_i].t_val == chichar )
                                                break;
                                }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                cout << chitoks[chi_i].t_name << endl;
#else
                                printf( "%s\n", chitoks[chi_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        }
                }
#endif /* YYDEBUG */
                if ( ++chi_ps >= &chis[ YYMAXDEPTH ] )    /* room on stack? */
                {
#ifndef YACC_MSG
                        chierror( "yacc stack overflow" );
#else
                        chierror(catgets(chiusercatd,1,2,"yacc stack overflow" ));
#endif
                        YYABORT;
                }
                *chi_ps = chi_state;
                *++chi_pv = chival;

                /*
                ** we have a new state - find out what to do
                */
        chi_newstate:
                if ( ( chi_n = chipact[ chi_state ] ) <= YYFLAG )
                        goto chidefault;         /* simple state */
#if YYDEBUG
                /*
                ** if debugging, need to mark whether new token grabbed
                */
                chitmp = chichar < 0;
#endif
                if ( ( chichar < 0 ) && ( ( chichar = chilex() ) < 0 ) )
                        chichar = 0;             /* reached EOF */
#if YYDEBUG
                if ( chidebug && chitmp )
                {
                        register int chi_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                        cout << "Received token " << endl;
                        if ( chichar == 0 )
                                cout << "end-of-file" << endl;
                        else if ( chichar < 0 )
                                cout << "-none-" << endl;
#else
                        printf( "Received token " );
                        if ( chichar == 0 )
                                printf( "end-of-file\n" );
                        else if ( chichar < 0 )
                                printf( "-none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        else
                        {
                                for ( chi_i = 0; chitoks[chi_i].t_val >= 0;
                                        chi_i++ )
                                {
                                        if ( chitoks[chi_i].t_val == chichar )
                                                break;
                                }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                cout << chitoks[chi_i].t_name << endl;
#else
                                printf( "%s\n", chitoks[chi_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        }
                }
#endif /* YYDEBUG */
                if ( ( ( chi_n += chichar ) < 0 ) || ( chi_n >= YYLAST ) )
                        goto chidefault;
                if ( chichk[ chi_n = chiact[ chi_n ] ] == chichar )  /*valid shift*/
                {
                        chichar = -1;
                        chival = chilval;
                        chi_state = chi_n;
                        if ( chierrflag > 0 )
                                chierrflag--;
                        goto chi_stack;
                }

        chidefault:
                if ( ( chi_n = chidef[ chi_state ] ) == -2 )
                {
#if YYDEBUG
                        chitmp = chichar < 0;
#endif
                        if ( ( chichar < 0 ) && ( ( chichar = chilex() ) < 0 ) )
                                chichar = 0;             /* reached EOF */
#if YYDEBUG
                        if ( chidebug && chitmp )
                        {
                                register int chi_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                cout << "Received token " << endl;
                                if ( chichar == 0 )
                                        cout << "end-of-file" << endl;
                                else if ( chichar < 0 )
                                        cout << "-none-" << endl;
#else
                                printf( "Received token " );
                                if ( chichar == 0 )
                                        printf( "end-of-file\n" );
                                else if ( chichar < 0 )
                                        printf( "-none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                else
                                {
                                        for ( chi_i = 0;
                                                chitoks[chi_i].t_val >= 0;
                                                chi_i++ )
                                        {
                                                if ( chitoks[chi_i].t_val
                                                        == chichar )
                                                {
                                                        break;
                                                }
                                        }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                        cout << chitoks[chi_i].t_name << endl;
#else
                                        printf( "%s\n", chitoks[chi_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                }
                        }
#endif /* YYDEBUG */
                        /*
                        ** look through exception table
                        */
                        {
                                register int *chixi = chiexca;

                                while ( ( *chixi != -1 ) ||
                                        ( chixi[1] != chi_state ) )
                                {
                                        chixi += 2;
                                }
                                while ( ( *(chixi += 2) >= 0 ) &&
                                        ( *chixi != chichar ) )
                                        ;
                                if ( ( chi_n = chixi[1] ) < 0 )
                                        YYACCEPT;
                        }
                }

                /*
                ** check for syntax error
                */
                if ( chi_n == 0 )        /* have an error */
                {
                        /* no worry about speed here! */
                        switch ( chierrflag )
                        {
                        case 0:         /* new error */
#ifndef YACC_MSG
                                chierror( "syntax error" );
#else
                                chierror(catgets(chiusercatd,1,3,"syntax error" ));
#endif
                                goto skip_init;
                        chierrlab:
                                /*
                                ** get globals into registers.
                                ** we have a user generated syntax type error
                                */
                                chi_pv = chipv;
                                chi_ps = chips;
                                chi_state = chistate;
                                chinerrs++;
                        skip_init:
                        case 1:
                        case 2:         /* incompletely recovered error */
                                        /* try again... */
                                chierrflag = 3;
                                /*
                                ** find state where "error" is a legal
                                ** shift action
                                */
                                while ( chi_ps >= chis )
                                {
                                        chi_n = chipact[ *chi_ps ] + YYERRCODE;
                                        if ( chi_n >= 0 && chi_n < YYLAST &&
                                                chichk[chiact[chi_n]] == YYERRCODE)                                        {
                                                /*
                                                ** simulate shift of "error"
                                                */
                                                chi_state = chiact[ chi_n ];
                                                goto chi_stack;
                                        }
                                        /*
                                        ** current state has no shift on
                                        ** "error", pop stack
                                        */
#if YYDEBUG
                                        if ( chidebug )
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                            cout << "Error recovery pops state "
                                                 << (*chi_ps)
                                                 << ", uncovers state "
                                                 << chi_ps[-1] << endl;
#else
#       define _POP_ "Error recovery pops state %d, uncovers state %d\n"
                                                printf( _POP_, *chi_ps,
                                                        chi_ps[-1] );
#       undef _POP_
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
#endif
                                        chi_ps--;
                                        chi_pv--;
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
                                if ( chidebug )
                                {
                                        register int chi_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                        cout << "Error recovery discards ";
                                        if ( chichar == 0 )
                                            cout << "token end-of-file" << endl;
                                        else if ( chichar < 0 )
                                            cout << "token -none-" << endl;
#else
                                        printf( "Error recovery discards " );
                                        if ( chichar == 0 )
                                                printf( "token end-of-file\n" );
                                        else if ( chichar < 0 )
                                                printf( "token -none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                        else
                                        {
                                                for ( chi_i = 0;
                                                        chitoks[chi_i].t_val >= 0;
                                                        chi_i++ )
                                                {
                                                        if ( chitoks[chi_i].t_val
                                                                == chichar )
                                                        {
                                                                break;
                                                        }
                                                }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                                cout << "token " <<
                                                    chitoks[chi_i].t_name <<
                                                    endl;
#else
                                                printf( "token %s\n",
                                                        chitoks[chi_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                        }
                                }
#endif /* YYDEBUG */
                                if ( chichar == 0 )      /* reached EOF. quit */
                                        YYABORT;
                                chichar = -1;
                                goto chi_newstate;
                        }
                }/* end if ( chi_n == 0 ) */
                /*
                ** reduction by production chi_n
                ** put stack tops, etc. so things right after switch
                */
#if YYDEBUG
                /*
                ** if debugging, print the string that is the user's
                ** specification of the reduction which is just about
                ** to be done.
                */
                if ( chidebug )
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                        cout << "Reduce by (" << chi_n << ") \"" <<
                            chireds[ chi_n ] << "\"\n";
#else
                        printf( "Reduce by (%d) \"%s\"\n",
                                chi_n, chireds[ chi_n ] );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
#endif
                chitmp = chi_n;                   /* value to switch over */
                chipvt = chi_pv;                  /* $vars top of value stack */
                /*
                ** Look in goto table for next state
                ** Sorry about using chi_state here as temporary
                ** register variable, but why not, if it works...
                ** If chir2[ chi_n ] doesn't have the low order bit
                ** set, then there is no action to be done for
                ** this reduction.  So, no saving & unsaving of
                ** registers done.  The only difference between the
                ** code just after the if and the body of the if is
                ** the goto chi_stack in the body.  This way the test
                ** can be made before the choice of what to do is needed.
                */
                {
                        /* length of production doubled with extra bit */
                        register int chi_len = chir2[ chi_n ];

                        if ( !( chi_len & 01 ) )
                        {
                                chi_len >>= 1;
                                chival = ( chi_pv -= chi_len )[1]; /* $$ = $1 */
                                chi_state = chipgo[ chi_n = chir1[ chi_n ] ] +
                                        *( chi_ps -= chi_len ) + 1;
                                if ( chi_state >= YYLAST ||
                                        chichk[ chi_state =
                                        chiact[ chi_state ] ] != -chi_n )
                                {
                                        chi_state = chiact[ chipgo[ chi_n ] ];
                                }
                                goto chi_stack;
                        }
                        chi_len >>= 1;
                        chival = ( chi_pv -= chi_len )[1]; /* $$ = $1 */
                        chi_state = chipgo[ chi_n = chir1[ chi_n ] ] +
                                *( chi_ps -= chi_len ) + 1;
                        if ( chi_state >= YYLAST ||
                                chichk[ chi_state = chiact[ chi_state ] ] != -chi_n )
                        {
                                chi_state = chiact[ chipgo[ chi_n ] ];
                        }
                }
                                        /* save until reenter driver code */
                chistate = chi_state;
                chips = chi_ps;
                chipv = chi_pv;
        }
        /*
        ** code supplied by user is placed in this switch
        */

                switch(chitmp){

case 5:
# line 37 "mchemsyntax.y"
{
	      AssignReactToSubst();
            } /*NOTREACHED*/ break;
case 6:
# line 43 "mchemsyntax.y"
{
             if (EndReaction())
               YYABORT;
           } /*NOTREACHED*/ break;
case 8:
# line 49 "mchemsyntax.y"
{
             CalcConvertor();
           } /*NOTREACHED*/ break;
case 9:
# line 53 "mchemsyntax.y"
{
             chemerror = 1;
             neduct = oldneduct; nproduct = oldnproduct;
             chierrok;
           } /*NOTREACHED*/ break;
case 16:
# line 70 "mchemsyntax.y"
{
          if (PlaceEduct(chipvt[-0].item, 1))
            YYABORT;
        } /*NOTREACHED*/ break;
case 17:
# line 75 "mchemsyntax.y"
{
          if (PlaceEduct(chipvt[-0].item, chipvt[-1].ival))
            YYABORT;
        } /*NOTREACHED*/ break;
case 20:
# line 85 "mchemsyntax.y"
{
            if (!PlaceSubstance(chipvt[-1].item, TRUE))  {
              chemerror = 1;
            }
          } /*NOTREACHED*/ break;
case 22:
# line 93 "mchemsyntax.y"
{
           chival.dval = (double)chipvt[-0].ival;
         } /*NOTREACHED*/ break;
case 23:
# line 99 "mchemsyntax.y"
{
            if (PlaceProduct(chipvt[-0].item, 1.))
              YYABORT;
          } /*NOTREACHED*/ break;
case 24:
# line 104 "mchemsyntax.y"
{
            if (PlaceProduct(chipvt[-0].item, chipvt[-1].dval))
              YYABORT;
          } /*NOTREACHED*/ break;
case 25:
# line 111 "mchemsyntax.y"
{
            if (!groundinterface)  {
              fprintf(stderr, "CHEM-ERROR: photodiss-reactions can only be calculated\n"
                              "            with groundinterface switched on.");
              chemerror = 1;
            }
            reaction[nreact].K = chipvt[-3].dval;
            reaction[nreact].EoR = chipvt[-1].dval;
            reaction[nreact].type = PHOTODISS;
          } /*NOTREACHED*/ break;
case 26:
# line 122 "mchemsyntax.y"
{
            if (!groundinterface)  {
              fprintf(stderr, "CHEM-ERROR: photodiss-reactions can only be calculated\n"
                              "            with groundinterface switched on.");
              chemerror = 1;
            }
            reaction[nreact].K = chipvt[-5].dval;
            reaction[nreact].N = chipvt[-3].dval;
            reaction[nreact].EoR = chipvt[-1].dval;
            reaction[nreact].type = PHOTODISS3;
          } /*NOTREACHED*/ break;
case 27:
# line 134 "mchemsyntax.y"
{
            reaction[nreact].K = chipvt[-3].dval;
            reaction[nreact].EoR = chipvt[-1].dval;
            reaction[nreact].type = THERMAL;
          } /*NOTREACHED*/ break;
case 28:
# line 140 "mchemsyntax.y"
{
            reaction[nreact].K = chipvt[-3].dval;
            reaction[nreact].EoR = chipvt[-1].dval;
            reaction[nreact].type = THERMAL2;
          } /*NOTREACHED*/ break;
case 29:
# line 146 "mchemsyntax.y"
{
            reaction[nreact].K = 1.;
            reaction[nreact].k0 = chipvt[-7].dval;
            reaction[nreact].N = chipvt[-5].dval;
            reaction[nreact].kinf = chipvt[-3].dval;
            reaction[nreact].M = chipvt[-1].dval;
            reaction[nreact].type = TROE;
          } /*NOTREACHED*/ break;
case 30:
# line 155 "mchemsyntax.y"
{
            reaction[nreact].k0 = chipvt[-11].dval;
            reaction[nreact].N = chipvt[-9].dval;
            reaction[nreact].kinf = chipvt[-7].dval;
            reaction[nreact].M = chipvt[-5].dval;
            reaction[nreact].K = chipvt[-3].dval;
            reaction[nreact].EoR = chipvt[-1].dval;
            reaction[nreact].type = TROEQUIL;
          } /*NOTREACHED*/ break;
case 31:
# line 165 "mchemsyntax.y"
{
            reaction[nreact].K = 1.;
            reaction[nreact].formula = chipvt[-1].item;
            reaction[nreact].specidx = nspecial++;
            reaction[nreact].specialine = lineno+1;
            reaction[nreact].type = SPECIAL;
          } /*NOTREACHED*/ break;
case 36:
# line 182 "mchemsyntax.y"
{
	  *chipvt[-2].dref = chipvt[-0].dval;
	} /*NOTREACHED*/ break;
}


        goto chistack;           /* reset registers in driver code */
}

# line 186 "mchemsyntax.y"
