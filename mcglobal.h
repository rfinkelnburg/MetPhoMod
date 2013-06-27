/* DEFINITION MODULE MCGlobal
   Globale Datendefinitionen des MeteoChem Programmes */

#ifndef INCLUDE_MCGLOBAL
#define INCLUDE_MCGLOBAL

#define UWIND 0
#define VWIND 1
#define WWIND 2
#define TEMP 3
#define HUMIDITY 4
#define CLOUDWATER 5
#define RAINWATER 6
#define CLOUDICE 7
/* Die folgenden Variablen dienen der k-eps Turbulenzparametrisierung */
#define TKE 8   /* Turbulent kinetic energy (e) */
#define EPS 9   /* Diffusion parameter of k-eps Param */
#define SUBS 10
#define NEW(n,t)  (t *)calloc((n), sizeof(t))

typedef int Entity;

#define MAXZGRID 200

#define BOOL int
#define TRUE 1
#define FALSE 0
#define MAXTINC 20
#define ITEMLEN 256

typedef double EntitySet[SUBS];
#define NORMALPOINT 0
#define UNDERGROUND 1
typedef long PointStatus;
typedef enum {
  CYCLIC = 1,
  MIRRORED = 2,
  FREEBORDER = 4,
  CONSTANTBORDER = 8,
  SPONGEBORDER = 16,
  AUTOCONSTANTBORDER = 32,
  NEIGHBOUR = 64,
  GEOMETRICBORDER = (FREEBORDER | CYCLIC | MIRRORED),
  DATABORDER = (CONSTANTBORDER | SPONGEBORDER | AUTOCONSTANTBORDER)
}  BorderType;
typedef enum  {
  FIXEDVALUE, AVERAGE}  BorderValue;
typedef enum  {
  NOPRESSURE, HYDROSTATIC, NONHYDROSTATIC}  PressureType;
typedef enum {
  NOCORIOLIS,
  FULLCORIOLIS,
  DIFFERENTIALCORIOLIS,
  GEOSTROPHICCORIOLIS
}  CoriolisType;
typedef enum  {ATPOINT, PILE, XAXIS, YAXIS, LAYER, TOTALGRID, MAXOUTTYPE}  OutType;
typedef enum  {NO_FILTER, PEPPER_FILTER, SHAPIRO_FILTER, DIFFUSION_FILTER}  FilterType;
typedef enum  {NO_TURB, TTTT, KEPS_T}  TurbType;
typedef enum  {MPDATA_A, PPM_A}  AdvectionType;

extern double dx, dy, dxi, dyi, dx2, dy2, dx2i, dy2i, ixx, iyy;
extern double **g;     /* Zeiger auf aktuelles Grid (g) */
extern BOOL *gactive;
extern double *Km, *flux[3], *oldu[3], *irdiv;
extern double *press;  /* Druck und Quellenterme */
extern double *pp, *density, *Rdirfact, *Rdifffact;
extern double *topo, *level, *leveli, *ldiff, *zcenter, *zfaces;
extern PointStatus *pstat;
extern double reflevel, blt, slt, sunelevation, photofact;
extern int nx, ny, nz, nxm, nym, nzm, nsubs, maxentity, avgentity, nta, tai,
	   row, xrow, layer, mesh;
extern int long actime, chemtime, tstart, tinc, tend, tincmax, dayofyear,
	   dumptime, chemtinc, irtinc, tchem, tincarray[MAXTINC], workerswanted;
extern double timeofday, topprecipwater, turnmapangle, toplevel;
extern PressureType pressuretype;
extern CoriolisType coriolistype;
extern FilterType filtertype;
extern TurbType turbtype;
extern AdvectionType advectiontype;
extern int *nlevelpt;
extern double toppp;
extern double *toppress, *toptemp, *ugeos, *vgeos;
extern double *avg;
extern double *tmplayer, *tmplayer2;
extern double Xlong, Xlat, Coriol3, Coriol2, timezonediff, topprecipwater,
	    stratozone, turbidity, spatialfilter, borderdamping;
extern BOOL groundinterface, printresid, plausible, withchem,
	    advection, windadvection, printwithborder, cloudwater, cloudice,
	    timeddump, dampinglayer, workeronthishost, shortwaveradiation,
            calcsoiltemp, radiationtype, horizontaldiffusion;
extern int smiord, smnonos;
extern double *westborder, *eastborder;
extern double *northborder, *southborder;
extern BorderType westbordertype, eastbordertype, northbordertype, southbordertype;
extern int *northborderval, *southborderval, *westborderval,
		   *eastborderval;
extern char dumppath[ITEMLEN], *worktitle;

typedef struct  {
  double *west, *east, *westdff, *eastdff;
  double *north, *south, *northdff, *southdff;
}  OldData;

extern OldData olddata;

extern double R, Ru, Cp, Kappa, Grav, P0, Cpress, Omega, OmegaM, Pi,
              CvbCp, GCvbCp, kBoltz, sigma, Lc, Ls, Li, T0, karman, w0,
              stdPKappa;
        /* Naturkonstanten. Bedeutung siehe MCGLOBAL.C */

#ifndef NOPROTO

double sqr(double x);

void MakeAFullDump(void);
/* Erstellt einen vollstaendigen binaeren Dump aller globalen Daten.*/

BOOL ReadFullDump(char *dumpname, BOOL reopen);
/* Liest einen vollstaendigen Dump wieder ein. */

#else

void InitData();
double sqr();
void MakeAFullDump();
BOOL ReadFullDump();

#endif

#endif

