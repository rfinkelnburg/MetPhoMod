/*
   DEFINITION MODULE mcparse.h
   Dieses Modul liest und interpretiert mit Hilfe des Modules mcsyntax.y
   das Input-File des Meteochem-Programmes.
*/

#ifndef INCLUDE_MCPARSE
#define MCPARSE
#define INCLUDE_MCPARSE

#ifndef MCGROUND
#include "mcground.h"
#endif

#ifndef INCLUDE_MCGROUP
#include "mc_group.hh"
#endif

typedef enum  {
   X_DIM = 1, Y_DIM = 2, Z_DIM = 4, TIME_DIM = 8, ALL_DIM = 7, COUNT_DIM = 16
}  Dimension;

typedef enum  {WARNING, ERROR}  ErrorLevel;

typedef enum  {TIME, GRID, OPTIONS, ENVIRONMENT, INITIAL_DATA,
   BORDER_DATA, GROUND_BORDER, TOP_BORDER, NORTH_BORDER, SOUTH_BORDER,
   WEST_BORDER, EAST_BORDER, CHEMISTRY, EMISSIONS, REDUCTION, DEPOSITION,
   NESTING, OUTPUT, NULL_SECTION
}  InputSection;

typedef int (*StartprocTemplate)(InputSection, InputSection);
typedef int (*EndprocTemplate)(void);
typedef enum {MUST_BE_CALLED, CAN_BE_OMITTED, WAS_CALLED, CALL_WHEN_SUBSTANCES}  CallType;

typedef struct SectionDesc {
  const char *name;
  int id;
  const struct SectionDesc *dependson;
  BOOL initialize, found;
  StartprocTemplate startproc;
  EndprocTemplate endproc;
  CallType calltype;
  struct SectionDesc *next;
}  SectionDesc;

typedef enum  {
  INT_NUM = 1,
  DOUBLE_NUM = 2,
  ON_OFF = 4,
  CORIOLIS_OPTION = 8,
  PRESSUREALG_OPTION = 0x10,
  BORDERTYPE_OPTION = 0x20,
  FILTERTYPE_OPTION = 0x40,
  TURBTYPE_OPTION = 0x80,
  ADVECTIONTYPE_OPTION = 0x100,
  SOILTEMP_OPTION = 0x200,
  RADIATION_OPTION = 0x400,
  ALL_OPTIONS = CORIOLIS_OPTION | PRESSUREALG_OPTION | BORDERTYPE_OPTION |
                FILTERTYPE_OPTION | ON_OFF  | TURBTYPE_OPTION |
                ADVECTIONTYPE_OPTION | SOILTEMP_OPTION | RADIATION_OPTION,
  TIME_VAL = 0x800,
  DATE_VAL = 0x1000,
  ARRAY_VAL = 0x2000,
  AVG_VAL = 0x4000,
  VECTOR_VAL = 0x8000,
  INT_VECTOR_VAL = 0x10000,
  TEXT_VAL = 0x20000,
  VAR_GROUP = 0x40000,
  NORMAL_NUM = (INT_NUM | DOUBLE_NUM | ARRAY_VAL),
  SINGLE_VAL = (INT_NUM | DOUBLE_NUM)
}  InputType;

typedef enum  {GRID_VAL, DOUBLE_PTR, LONG_PTR, INT_PTR, CHAR_PTR,
               GROUND_PARAM, VECT_PTR, INT_VECT_PTR, PROC_VAL, GROUP_VAL}  StoreType;

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
  enum { IS_NORMAL, IS_INT_VECT, IS_VECT} vectype;
}  ValueDesc;

typedef enum  {CALCULATED_VALUES, SET_DEFAULT, MUST_BE_SET, WAS_SET_BY_USER, WAS_PARTITALLY_SET,
   MADE_AN_ERROR}  InitType;
typedef short unsigned DimensionSet;

typedef struct  {
  int *nval;
  long *v;
}  IntVectPtr;

typedef struct {
  int *nval;
  double *v;
}  VectPtr;

typedef double (*PointProc)(int, int, int, struct VarDesc *v);
/* Eine PointProc berechnet aus den Koordinaten k, i, j einen Wert. 
   Der vierte Parameter entspricht der id der Aufrufenden Variablen */

typedef enum  {
  IS_RELEVANT, NAME_OF_RELVAR, VAL_OF_RELVAR
}  RelyCmmd;

typedef char *(*RelyProc)(RelyCmmd cmmd, struct VarDesc *v);

struct VarDesc  {
  // Variables
  const char *name, *comment, *unit;
  StoreType storetype;
  union  {
    long *li;
    long et;
    double *d;
    int *i;
    GroundParam *g;
    IntVectPtr *ivp;
    VectPtr *vp;
    char *txt;
    PointProc proc;
    Group<VarDesc> *group;
  }  v;
  DimensionSet dims;
  InputSection section, origsection;
  InitType init;
  double defval, rmin, rmax;   /* Range-min and Range-max */
  InputType inputtype;
  int option, ncoord, id, userdata;
  RelyProc relieson;
  VarDesc *next;
  VarDesc(void) {}
  VarDesc(const char *name_in, const char *comment_in, const char *unit_in,
	  StoreType storetype_in,
	  void *data,
	  DimensionSet dims_in,
	  InputSection section_in,
	  InitType init_in,
	  double defval_in, double rmin_in, double rmax_in,
	  InputType it,
	  int option_in, RelyProc rely, 
	  int userdata_in = 0,
	  int ncoord_in = 0)
    {
      name = name_in;
      comment = comment_in;
      unit = unit_in;
      storetype = storetype_in;
      v.li = (long *)data;
      dims = dims_in;
      origsection = section = section_in;
      init = init_in;
      defval = defval_in;
      rmin = rmin_in; rmax = rmax_in;
      inputtype = it;
      option = option_in + (init << 8);
      relieson = rely;
      userdata = userdata_in;
      ncoord = ncoord_in;
    }
};

extern VarDesc variable[];

VarDesc *BorderVariable(VarDesc *v, VarDesc *vh);

#define MAXOPTION 28

extern OptionDesc option[MAXOPTION];
extern ValueDesc val;
extern SectionDesc *actualsection;
extern int inputerror, inplineno;
extern Date startdate;
extern Time starttime;

int EndOfTime(void);
int EndOfOptions(void);
int EndOfGrid(void);
int EndOfEnvironment(void);
int EndOfGround(void);
int EndOfSection(void);

int ChangeSectionTo(char *sname);

int AssignVariable(void);

int SetVariable(const char *name);
// Set variable to be set by name

int SetVariable(const char *groupname, const char *name);
// Set VariableGroup variable to be set

int SetVariable(VarDesc *v);
// Set variable to be set directly

int CopyRangeToVar(void);

int FindOption(char *oname);

void InitArrayVar(void);

void AllocateArraySpace(void);

int SetDimension(Dimension d);

int SetRange(Dimension d, int min, int max);

void InitVarTable(void);

int ParseInput(char *fname);

int AddArrayVal(double num);

int AddIntArrayVal(long num);

void InputError(ErrorLevel level, const char *s, ...);

int DefineAVariable(char *name);

int CopyVarToVal(char *name);

void ShowTitle(char *title);

#ifdef PARALLEL

VarDesc *GetNamedVar(char *name);

#endif

void InitializeData(void);

void mcierror(char *s);

extern "C" int mciwrap();

int mcilex(void);

double AbsolutePointTemp(int k, int i, int j, struct VarDesc *v);

double EquivalentPointTemp(int k, int i, int j, struct VarDesc *v);

double RelativeHumidity(int k, int i, int j, struct VarDesc *v);

#endif
