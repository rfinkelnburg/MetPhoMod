/* IMPLEMENTATION MODULE MCGlobal
   Globale Datendefinitionen des MeteoChem Programmes */

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#ifdef PARALLEL
#include <pvm3.h>
#endif
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mcprint.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

double dx, dy, dxi, dyi, dx2, dy2, dx2i, dy2i, ixx, iyy;
double **g;           /* Zeiger auf aktuelles Grid (g) */
BOOL *gactive;
double *Km, *flux[3], *oldu[3], *irdiv;
double *pp, *density, *Rdirfact, *Rdifffact;
double *topo, *level, *leveli, *ldiff, *press, *tmplayer, *tmplayer2, *zcenter, *zfaces;
PointStatus *pstat;
double reflevel;
int nx, ny, nz, nxm, nym, nzm, nsubs, maxentity, avgentity,
	   row, xrow, layer, mesh;
int long actime, chemtime, tstart, tinc, tincmax, tend, dayofyear,
    dumptime, chemtinc, irtinc, tchem, tincarray[MAXTINC], workerswanted;
#ifdef PARALLEL
int worknum;
#endif
double timeofday, turnmapangle, sunelevation, photofact, toplevel;
BorderValue bordervalue;
PressureType pressuretype;
CoriolisType coriolistype;
FilterType filtertype;
int *nlevelpt, nta, tai;
double toppp;
double *toppress, *toptemp, *ugeos, *vgeos;
double *avg;
double Xlong, Xlat, Coriol3, Coriol2, timezonediff, topprecipwater, stratozone,
       turbidity, spatialfilter, borderdamping;
BOOL printresid, plausible, withchem, groundinterface,
     advection, windadvection, printwithborder, cloudwater, cloudice,
     timeddump, dampinglayer, workeronthishost, shortwaveradiation;
int smiord, smnonos;
double *westborder, *eastborder;
double *northborder, *southborder;
BorderType westbordertype, eastbordertype, northbordertype, southbordertype;
TurbType turbtype;
AdvectionType advectiontype;
int *northborderval, *southborderval, *westborderval,
		   *eastborderval;
char dumppath[ITEMLEN] = "", *worktitle;

double R  = 287.04,          /* Luft-Gaskonstante */
       Ru = 8.31441e9,	     /* Universelle Gaskonstante * 10^9 */
       Cp = 1005.,           /* isobare Waermekapazitaet der Luft */
       Kappa = 0.28561194,   /* hydrostatischer Exponent */
       Grav  = 9.80597,      /* Erdbeschleunigung */
       P0 = 100000.,         /* Standard Atmosphaerendruck */
       Cpress = 0.261447322, /* G*P0^Kappa/Cp */
       Lc = 2.5e6,           /* Verdampfungswaerme des Wassers */
       Ls = 2.83e6,	     /* Sublimationswaerme des Wassers */
       Li = 3.34e5,          /* Schmelzwaerme Eis/Wasser */
       T0 = 273.16,          /* Schmelztemparatur des Wassers */
       Omega  = 1.895,
       karman = 0.35,	     /* Karmans Konstante */
       OmegaM = -0.895,
       Pi = 3.14159265359,
       GCvbCp = 6.99590,     /* Cv / Cp * g */
       CvbCp = 0.714388,
       kBoltz = 1.380662e-23,/* Boltzmann Konstante */
       sigma = 5.67032e-8,   /* Stefan-Boltzmann Konstante */
       w0 = 0.017453293,       /* Umrechnung deg->rad */
       stdPKappa = 26.7954;  /* p0^Kappa */

double sqr(double x)
{
  return (x * x);
}

extern void InstallHandler(void);

#define dump(x)  fwrite(&x, sizeof(x), 1, df)


void MakeAFullDump(void)
{
  int i;
  FILE *df;
  char fname[128], buf[64];
  OutFile *of;
#ifdef PARALLEL
  if (parallel)
    if (!master)  {
      SendMeshToMaster(Km, 0, ALL_DIM);
      SendMeshToMaster(press, -1, ALL_DIM);
      SendGroundToMaster();
      for (i = 0; i < maxentity; i++)
	if (g[i])
	  SendMeshToMaster(g[i], 0, ALL_DIM);
      return;
    }
    else  {
      ReadMeshFromWorkers(Km, 0, ALL_DIM);
      ReadMeshFromWorkers(press, -1, ALL_DIM);
      ReadGroundFromWorker();
    }
#endif
  if (*dumppath)
    sprintf(fname, "%s/mcdump", dumppath);
  else
    strcpy(fname, "mcdump");
  if (timeddump)
    sprintf(buf, "%s.%li", fname, actime);
  else
    strcpy(buf, fname);
  strcpy(fname, buf);
  if (!(df = fopen(fname, "w")))  {
    fprintf(stderr, "Warning: Unable to open dump-file \"%s\"\n", fname);
    return;
  }
  dump(nx); dump(ny); dump(nz); dump(nsubs);
  fwrite(Km, sizeof(double), mesh, df);
  fwrite(press-layer, sizeof(double), (nz+1)*layer, df);
  fwrite(pp-1, sizeof(double), nz+1, df);
  fwrite(density-1, sizeof(double), nz+1, df);
  dump(actime); dump(chemtime); dump(tinc); dump(chemtinc);
  fwrite(ground, sizeof(GroundParam), layer, df);
  for (i = 0; i < maxentity; i++)
    if (g[i])  {
#ifdef PARALLEL
      if (parallel && master)
        ReadMeshFromWorkers(g[i], 0, ALL_DIM);
#endif
      fwrite(g[i], sizeof(double), mesh, df);
    }
  for (of = outfile; of; of = of->next)  {
    fwrite(of->fname, 1, 80, df);
  }
  fclose(df);
}

#define inr(x)  fread(&x, sizeof(x), 1, df)

/*
void ReadRange(void *one, int onelen, void *two, int twolen, FILE *f)
{
  char *from, *to;
  int len;
  if (two > one)  {from = (char *)one; to = (char *)two; len = twolen;}
  else            {from = (char *)two; to = (char *)one; len = onelen;}
  fread(from, to - from + len, 1, f);
}
*/
BOOL ReadFullDump(char *dumpname, BOOL reopen)
{
  int i, tnx, tny, tnz, tns;
  FILE *df;
  OutFile *of;
#ifdef PARALLEL
  if (parallel && !master)  {
    GetMeshFromMaster(Km, 0);
    GetMeshFromMaster(press, -1);
    pvm_recv(myparent, VALUES);
    pvm_upkdouble(pp-1, nz+1, 1);
    pvm_upkdouble(density-1, nz+1, 1);
    pvm_upklong(&actime, 1, 1);
    pvm_upklong(&chemtime, 1, 1);
    pvm_upklong(&tinc, 1, 1);
    pvm_upklong(&chemtinc, 1, 1);
    GetGroundFromMaster();
    for (i = 0; i < maxentity; i++)
      if(g[i])
        GetMeshFromMaster(g[i], 0);
    return (TRUE);
  }
#endif
  if (!(df = fopen(dumpname, "r")))  return (FALSE);
  printf("Reading dumpfile...");
  inr(tnx); inr(tny); inr(tnz); inr(tns);
  if (tnx != nx || tny != ny || tnz != nz || tns != nsubs)  {
    fprintf(stderr, "Size of dump-file is not correct!\n");
    exit (1);
  }
  fread(Km, sizeof(double), mesh, df);
  fread(press-layer, sizeof(double), (nz+1)*layer, df);
  fread(pp-1, sizeof(double), nz+1, df);
  fread(density-1, sizeof(double), nz+1, df);
  inr(actime); inr(chemtime); inr(tinc); inr(chemtinc);
  fread(ground, sizeof(GroundParam), layer, df);
#ifdef PARALLEL
  if (parallel)  {
    SendMeshToWorkers(Km, 0);
    SendMeshToWorkers(press, -1);
    pvm_initsend(PvmDataRaw);
    pvm_pkdouble(pp-1, nz+1, 1);
    pvm_pkdouble(density-1, nz+1, 1);
    pvm_pklong(&actime, 1, 1);
    pvm_pklong(&chemtime, 1, 1);
    pvm_pklong(&tinc, 1, 1);
    pvm_pklong(&chemtinc, 1, 1);
    SendToAll(VALUES);
    SendGroundToWorkers();
  }
#endif
  for (i = 0; i < maxentity; i++)
    if (g[i])  {
      fread(g[i], sizeof(double), mesh, df);
#ifdef PARALLEL
      SendMeshToWorkers(g[i], 0);
#endif
    }
  if (reopen)
    for (of = outfile; of; of = of->next)
      if (fread(of->fname, 1, 80, df) != 80)  {
	printf("error\n");
	return (FALSE);
      }
  printf("successful\n");
  fclose(df);
  return (TRUE);
}
