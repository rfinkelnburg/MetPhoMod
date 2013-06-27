/*
   DEFINITION MODULE mcparse.h
   Dieses Modul liest und interpretiert mit Hilfe des Modules mcsyntax.y
   das Input-File des Meteochem-Programmes.
*/

#define MCPARSE

#ifndef MCGROUND
#include "mcground.h"
#endif

#define NORMAL_NUM (INT_NUM | DOUBLE_NUM | ARRAY_VAL)
#define SINGLE_VAL (INT_NUM | DOUBLE_NUM)

typedef enum  {
   X_DIM = 1, Y_DIM = 2, Z_DIM = 4, TIME_DIM = 8, ALL_DIM = 7, COUNT_DIM = 16
}  Dimension;

typedef enum  {WARNING, ERROR}  ErrorLevel;

typedef enum  {TIME, GRID, OPTIONS, ENVIRONMENT, INITIAL_DATA,
   BORDER_DATA, GROUND_BORDER, TOP_BORDER, NORTH_BORDER, SOUTH_BORDER,
   WEST_BORDER, EAST_BORDER, CHEMISTRY, EMISSIONS, REDUCTION, DEPOSITION,
   NESTING, OUTPUT, NULL_SECTION
}  InputSection;

typedef enum  {
  INT_NUM = 1,
  DOUBLE_NUM = 2,
  ON_OFF = 4,
  CORIOLIS_OPTION = 8,
  PRESSUREALG_OPTION = 16,
  BORDERTYPE_OPTION = 32,
  FILTERTYPE_OPTION = 64,
  TURBTYPE_OPTION = 16384,
  ADVECTIONTYPE_OPTION = 32768,
  ALL_OPTIONS = CORIOLIS_OPTION | PRESSUREALG_OPTION | BORDERTYPE_OPTION |
                FILTERTYPE_OPTION | ON_OFF  | TURBTYPE_OPTION |
                ADVECTIONTYPE_OPTION,
  TIME_VAL = 128,
  DATE_VAL = 256,
  ARRAY_VAL = 512,
  AVG_VAL = 1024,
  VECTOR_VAL = 2048,
  INT_VECTOR_VAL = 4096,
  TEXT_VAL = 8192
}  InputType;

typedef enum  {GRID_VAL, LAYER_PTR, FILE_PTR, XFILE_PTR, DOUBLE_PTR, LONG_PTR,
               INT_PTR, CHAR_PTR, GROUND_PARAM, INT_VECT_PTR, PROC_VAL}  StoreType;

typedef enum {STARTING_WITH_ZERO = 1, NO_ERRORS_ALLOWED = 2}  VariableOption;

typedef struct  {
  char *name, *comment;
  InputType type;
  int value;
}  OptionDesc;

typedef struct  {
  int day, month, year;
}  Date;

typedef struct  {
  int hour, minute;
}  Time;

typedef struct  {
  int nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax;
  int has_dimension, has_range, nval, maxval;
  int deletable;
  double *a;
}  ArrayDesc;

typedef union  {
  int ival;
  double dval;
  char txt[ITEMLEN];
  Date date;
  Time time;
  ArrayDesc a;
}  ValueType;

typedef struct   {
  InputType type;
  ValueType v;
}  ValueDesc;

typedef enum  {CALCULATED_VALUES, SET_DEFAULT, MUST_BE_SET, WAS_SET_BY_USER, WAS_PARTITALLY_SET,
   MADE_AN_ERROR}  InitType;
typedef short unsigned DimensionSet;

typedef struct  {
  char *name;
  InputSection dependson;
  BOOL initialize;
  int (*startproc)(InputSection, InputSection);
  int (*endproc)();
  enum {MUST_BE_CALLED, CAN_BE_OMITTED, WAS_CALLED, CALL_WHEN_SUBSTANCES}  calltype;
}  SectionDesc;

typedef enum  {
  IS_RELEVANT, NAME_OF_RELVAR, VAL_OF_RELVAR
}  RelyCmmd;

typedef struct  {
  int *nval;
  long *v;
}  IntVectPtr;

typedef double (*PointProc)(int, int, int);
/* Eine PointProc berechnet aus den Koordinaten k, i, j einen Wert. */

typedef struct VarDesc  {
  char *name, *comment, *unit;
  StoreType storetype;
  union  {
    long *li;
    long et;
    double *d;
    int *i;
    GroundParam *g;
    IntVectPtr *ivp;
    char *txt;
    PointProc proc;
  }  v;
  DimensionSet dims;
  InputSection section, origsection;
  InitType init;
  double defval, rmin, rmax;   /* Range-min and Range-max */
  InputType inputtype;
  int option, ncoord, id;
  char *(*relieson)(RelyCmmd cmmd, struct VarDesc *v);
  struct VarDesc *next;
}  VarDesc;

extern VarDesc variable[];
extern SectionDesc section[];

VarDesc *BorderVariable(VarDesc *v, VarDesc *vh);

#define MAXOPTION 24

extern OptionDesc option[MAXOPTION];
extern ValueDesc val;
extern InputSection actualsection;
extern int inputerror, inplineno;
extern Date startdate;
extern Time starttime;

int EndOfTime(void);
int EndOfOptions(void);
int EndOfGrid(void);
int EndOfEnvironment(void);
int EndOfGround(void);

int ChangeSectionTo(char *sname);

int AssignVariable(void);

int SetVariable(char *name);

int CopyRangeToVar(void);

int FindOption(char *oname);

void InitArrayVar(void);

void AllocateArraySpace(void);

int SetDimension(Dimension d);

int SetRange(Dimension d, int min, int max);

void InitVarTable(void);

int ParseInput(char *fname);

void InputError(ErrorLevel level, const char *s, ...);

int DefineAVariable(char *name);

int CopyVarToVal(char *name);

void ShowTitle(char *title);

#ifdef PARALLEL

VarDesc *GetNamedVar(char *name);

#endif

void InitializeData(void);
