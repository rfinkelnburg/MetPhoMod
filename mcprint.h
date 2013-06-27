/*
   DEFINITION MODULE mcprint
   Dieses Modul ermoeglicht die Erstellung behandelt die Ausgabe der Modell-
   resultate.
*/

#ifdef MCGROUND

typedef enum {ASCII_PRINT, CDF_PRINT}  PrintType;

typedef struct OutVar {
  DimensionSet has_range;
  int xmin, xmax, ymin, ymax, zmin, zmax, datid;
  VarDesc *v;
  struct OutVar *next;
  char printname[40];
}  OutVar;

typedef struct OutFile  {
  PrintType prtype;
  int fid, xid, yid, zid, timeid, tvarid;
  FILE *f;
  long nextprint, printinterval;
  OutVar *vptr;
  size_t ntimerec;
  char fname[80];
  struct OutFile *next;
}  OutFile;

extern OutFile *outfile;

OutVar *OutItem(char *name, char *prname);

OutVar *AppendOutVar(OutVar *first, OutVar *add);

int PrepareOutput(PrintType prtype, char *name, long firstprint, long dprint);

#endif

int OpenOutputFiles(char *inprow, BOOL reopen);

void CloseOutputFiles(void);

int WriteOutData(long actime, BOOL ignoretime);

void FlushOutputFiles(void);

int CreateNonExistingFile(char *fname);
