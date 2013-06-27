/*
   DEFINITION MODULE mchemparse.h
   Strukturdefinitionen fuer die Chemieberechnungen.
*/


#define MAXSUBST 120
#define MAXEDUCT 600
#define MAXPRODUCT 800
#define MAXREACT 300
#define TOKENLEN 16

typedef char Token[TOKENLEN];

typedef enum {THERMAL, THERMAL2, TROE, TROEQUIL, SPECIAL, PHOTODISS, PHOTODISS3}  KinType;

typedef struct  {
  int subs;
  double ht, dcdt, oldc;
  Token name;
  struct Reaction **decay, **production;
  int ndecay, nproduction, fastidx, take_as_fast, 
      inert, do_transport;
}  Substance;

typedef struct  {
  int prefact;
  Substance *comp;
}  Educt;

typedef struct  {
  double prefact;
  Substance *comp;
}  Product;

typedef struct Reaction  {
  Educt *educt;
  Product *product;
  double K, EoR, N, M, k0, kinf, rate, Ke;
  char *formula;
  int neduct, nproduct, netot, specidx, specialine;
#ifdef FULLM
  int nm, nox, nn2;
#endif
  KinType type;
}  Reaction;

typedef struct  {
  double press, T;
}  KinReference;

typedef struct ProductionDesc  {
  double *pps;
  int subs;
  struct ProductionDesc *next;
}  ProductionDesc;

extern Educt educt[MAXEDUCT];
extern Product product[MAXPRODUCT];
extern Substance subst[MAXSUBST], *fast[MAXSUBST], *slow[MAXSUBST];
extern Reaction reaction[MAXREACT];
extern int neduct, netot, nproduct, oldneduct, oldnproduct, 
   nreact, nsubst, lineno, chemerror, nfast, nslow, nspecial;
extern double konzunitfact, timeunitfact, convertorfact, fixedconc;
extern char lasttoken[];
extern Reaction *reactionptr[];
extern KinReference reference;
extern int allfast;
extern ProductionDesc *prodd;

double spow(double x, int y);

Substance *PlaceSubstance(char *name, BOOL mustbenew);

int PlaceEduct(char *name, int fact);

int PlaceProduct(char *name, double fact);

int EndReaction();

void SetSubstKonz(char *name, double conz);

void ConvertUnits(void);

int ReadChemFile(char *fname);

void CalcConvertor(void);

int CreateChemTable(void);

int GenerateSpecial(FILE *f, char *chemname);

int TestSpecial();

double SpecialConst(int idx, double T, double M);
/* Eine Prozedur in dieser Form wird durch "GenerateSpecial" erzeugt */

void chierror(char *s);
int chiwrap();

void AssignReactToSubst();
