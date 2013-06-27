/*
   DEFINITION MODULE mcemmiss
*/

#ifdef MCGROUND

typedef struct  {
  float x, y, z;
  long loc;
}  Coordinate;

typedef struct CoordArray  {
  Coordinate *c;
  int n;
}  CoordArray;

typedef struct EmissionDesc {
  VarDesc *v, *vo;   /* Holds variable, refers to Grid-Variable */
  double *em;
  Entity et;    /* Number of species */
  CoordArray cem;
  double *reduction;
  int emvid;    /* Unique Emission-Variable ID */
  struct EmissionDesc *next;
}  EmissionDesc;

Coordinate *AllocateCoordVar(int size);

extern EmissionDesc *firstemiss;

VarDesc *EmissionVariable(VarDesc *v, int pointvalues, int id, Coordinate *coord);

VarDesc *GetEmVarWithID(int id);

EmissionDesc *FindEmissionVariable(VarDesc *v);

void ConvertEmissionCoords(int n, Coordinate *coord, int xoff);

#endif

void Emit(double tinc);
