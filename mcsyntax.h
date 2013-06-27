
typedef union   {
  int ival;
  double dval;
  char item[ITEMLEN];
  Date date;
  Time time;
  OutVar *outvar;
} YYSTYPE;
extern YYSTYPE mcilval;
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
