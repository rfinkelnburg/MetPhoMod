/*
   MODULE mcnest.c
   Dieses Modul unterstuetzt das Herstellen von Nest-Files. Berphomod in
   seiner aktuellen Version unterstuetzt zwar das Nesten nicht wirklich.
   Mit den Routinen in diesem Modul ist es aber einfach moeglich, Subdomains
   zu definieren, fuer die die entsprechenden Randwerte dann in entsprechende
   Dateien rausgeschrieben werden
   
   This modules supports writing of Nest-Files. While Berphomod in its actual
   version is not able to do a real time nesting, the routines defined here,
   support the creation of files, which describe the borders of a subdomain,
   and which can easily be used for a new run on this new domain.
   
   mcnest is a part of Berphomod, the comprehensive eulerian atmospheric
   air pollution simulation program.
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mcnest.h"
#include "mcprint.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

typedef struct CDFVarDescr  {
  int vid, et;
}  CDFVarDescr;

typedef struct CDFileDescr  {
  int fid, tid, nv, uvid, vvid;
  CDFVarDescr *v;
}  CDFileDescr;

typedef struct NestDescr {
  char name[80];
  long xorigin, yorigin,
       dx, dy,
       nnx, nny,
       ndt;
  double turnangle, sina, cosa;
  CDFileDescr west, east, south, north, top;
  struct NestDescr *next;
}  NestDescr;

char nestvarname[NNESTVAR][20] =
  {"Nest-xorigin", "Nest-yorigin", "Nest-nx", "Nest-ny", 
   "Nest-dx", "Nest-dy", "Nest-angle", "Nest-dt"};
double dummynestvar;

NestDescr *firstdomain = NULL;

void ResetVarReference(char *vname, long *var)
{
  VarDesc *v;
  for (v = variable; v && v->name != vname; v = v->next);
  if (!v)  {
    InputError(ERROR, "Internal implementation error in \"SearchVarReference\"");
    exit (1);
  }
  v->v.li = var;
  v->init = MUST_BE_SET;
}

int TestInside(NestDescr *domain, double x, double y)
{
  if (x < 0. || y < 0. || x > nx * dx || y > ny * dy)  {
    InputError(ERROR, "Nest-Domain \"%s\" is not inside modeling domain.", domain->name);
    return (1);
  }
  return (0);
}

int TestNestDomain(NestDescr *domain)
{
  VarDesc *v;
  int rc = 0, xmin, ymin, xmax, ymax, xe, ye;
  if (firstdomain)
    for (v = variable; v; v = v->next)
      if (v->section == NESTING && v->init == MUST_BE_SET)  {
	InputError(ERROR, "Variable \"%s\" was not set for Domain \"%s\"",
           v->name, domain->name);
	rc = 1;
      }
  domain->turnangle *= w0;  /* w0 ist in mcglobal.c definiert. */
  domain->sina = sin(domain->turnangle);
  domain->cosa = cos(domain->turnangle);
  rc |= TestInside(domain, domain->xorigin, domain->yorigin);
  rc |= TestInside(domain,
	  domain->xorigin + domain->cosa * (domain->nnx - 1) * domain->dx,
	  domain->yorigin + domain->sina * (domain->nnx - 1) * domain->dx);
  rc |= TestInside(domain,
          domain->xorigin - domain->sina * (domain->nny - 1) * domain->dy,
          domain->yorigin + domain->cosa * (domain->nny - 1) * domain->dy);
  rc |= TestInside(domain,
          domain->xorigin + domain->cosa * (domain->nnx - 1) * domain->dx - domain->sina * (domain->nny - 1) * domain->dy,
          domain->yorigin + domain->sina * (domain->nnx - 1) * domain->dx + domain->cosa * (domain->nny - 1) * domain->dy);
  return (rc);
}

int EndOfNesting(void)
{
  return (TestNestDomain(firstdomain));
}

int DefineNestDomain(char *basename, int section)
{
  NestDescr *newdomain, *hd;
  VarDesc *v;
  if (section != NESTING)  {
    InputError(ERROR, "Nesting-Subdomains have to be defined in section \"NESTING\"");
    return (1);
  }
  if (firstdomain)  TestNestDomain(firstdomain);
  newdomain = (NestDescr *)calloc(1, sizeof(NestDescr));
  if (!newdomain)  {
    InputError(ERROR, "Error allocating memory in \"NESTING\"\n");
    exit (1);
  }
  strncpy(newdomain->name, basename, 80);
  for (hd = firstdomain; hd; hd = hd->next)
    if (!strcmp(newdomain->name, hd->name))  {
      InputError(ERROR, "Subdomains must have differing names!\n"
      		"There are two domains called \"%s\"", hd->name);
      return (1);
    }
  newdomain->next = firstdomain;
  firstdomain = newdomain;
  ResetVarReference(nestvarname[0], &newdomain->xorigin);
  ResetVarReference(nestvarname[1], &newdomain->yorigin);
  ResetVarReference(nestvarname[2], &newdomain->nnx);
  ResetVarReference(nestvarname[3], &newdomain->nny);
  ResetVarReference(nestvarname[4], &newdomain->dx);
  ResetVarReference(nestvarname[5], &newdomain->dy);
  ResetVarReference(nestvarname[6], (long *)&newdomain->turnangle);
  ResetVarReference(nestvarname[7], &newdomain->ndt);
  return (0);
}

char *OnNestDomain(RelyCmmd cmd, VarDesc *v)
{
  switch (cmd)  {
    case IS_RELEVANT :
       return ((char *)!!firstdomain);
    case NAME_OF_RELVAR :
       if (!firstdomain)  {
         InputError(ERROR, "You have to establish a domain using the \"DOMAIN\"-satement before initializing variable \"%s\"!",
            v->name);
         inputerror = 1;
       }
       return ("Nesting-domain");
    case VAL_OF_RELVAR :
       return (char *)(firstdomain ? "defined" : "undefined");
  }
  return (NULL);
}

BOOL ReopenDomainFile(char *fname, char *dim1name, char *dim2name, CDFileDescr *f,
		      BOOL pressure_only)
{
  int net, i;
  VarDesc *v;
  if (nc_open(fname, NC_WRITE, &f->fid))  {
    fprintf(stderr, "FATAL ERROR: Unable to reopen file \"%s\" (a nesting domain file)\n", fname);
    return (TRUE);
  }
  printf("%s...\n", fname);
  if (nc_inq_varid(f->fid, "Time", &f->tid))
    goto not_all_vars;
  net = 0;
  if (pressure_only)  net = 1;
  else
    for (i = maxentity; i--; )
      net += (g[i] && i != WWIND);
  f->v = (CDFVarDescr *)calloc(net, sizeof(CDFVarDescr));
  if (!f->v)  {
    fprintf(stderr, "FATAL ERROR: Unable to allocate memory in \"ReopenDomainFile\"\n");
    exit (3);
  }
  if (pressure_only)  {
    if (nc_inq_varid(f->fid, "delta_TopPressure", &f->v[0].vid))
      goto not_all_vars;
    f->nv = 1;
  }
  else  {
    i = 0;
    for (v = variable; v; v = v->next)  {
      if (v->storetype == GRID_VAL && v->v.et != WWIND &&
	  g[v->v.et])  {
	if (nc_inq_varid(f->fid, v->name, &f->v[i].vid))  goto not_all_vars;
	f->v[i].et = v->v.et;
	if (v->v.et == UWIND)
	  f->uvid = f->v[i].vid;
	else if (v->v.et == VWIND)
	  f->vvid = f->v[i].vid;
	i++;
      }
    }
    f->nv = i;
  }
  return (FALSE);
not_all_vars :
  fprintf(stderr, "FATAL ERROR: Unable to find variables in file \"%s\"\n", fname);
  return (TRUE);
}

void CreateDomainFile(char *fname, char *dim1name, int dim1n, double *dim1k,
				   char *dim2name, int dim2n, double *dim2k,
		      CDFileDescr *f, BOOL pressure_only, BOOL reopen)
{
  int dim[3], i, xid, yid, net;
  char datename[128];
  VarDesc *v;
  if (reopen)  {
    if (ReopenDomainFile(fname, dim1name, dim2name, f, pressure_only))  {
      fprintf(stderr, "FATAL ERROR: Unable to reopen domain file \"%s\".\n", fname);
      inputerror = TRUE; plausible = FALSE;
    }
    return;
  }
  f->fid = CreateNonExistingFile(fname);
  nc_def_dim(f->fid, "Time", NC_UNLIMITED, dim);
  nc_def_dim(f->fid, dim2name, dim2n, dim+1);
  nc_def_dim(f->fid, dim1name, dim1n, dim+2);
  nc_def_var(f->fid, "Time", NC_LONG, 1, dim, &f->tid);
  sprintf(datename, "seconds since %d-%d-%d %d:%d:00 %+d:%02d",
          startdate.year, startdate.month, startdate.day,
          starttime.hour, starttime.minute,
          (int)(timezonediff + 0.5), (int)(60. * fmod(timezonediff, 1.) + 0.5));
  nc_put_att_text(f->fid, f->tid, "units", strlen(datename), datename);
  nc_def_var(f->fid, dim1name, NC_DOUBLE, 1, dim+2, &xid);
  nc_put_att_text(f->fid, xid, "units", 1, "m");
  nc_def_var(f->fid, dim2name, NC_DOUBLE, 1, dim+1, &yid);
  nc_put_att_text(f->fid, yid, "units", 1, "m");
  net = 0;
  if (pressure_only)  net = 1;
  else
    for (i = maxentity; i--; )
      net += (g[i] && i != WWIND);
  f->v = (CDFVarDescr *)calloc(net, sizeof(CDFVarDescr));
  if (!f->v)  {
    fprintf(stderr, "FATAL ERROR: Unable to allocate memory in \"CreateDomainFile\"\n");
    exit (3);
  }
  if (pressure_only)  {
    nc_def_var(f->fid, "delta_TopPressure", NC_DOUBLE, 3, dim, &f->v[0].vid);
    f->nv = 1;
  }
  else  {
    i = 0;
    for (v = variable; v; v = v->next)  {
      if (v->storetype == GRID_VAL && v->v.et != WWIND &&
	  g[v->v.et])  {
	nc_def_var(f->fid, v->name, NC_DOUBLE, 3, dim, &f->v[i].vid);
	if (v->unit)
	  nc_put_att_text(f->fid, f->v[i].vid, "units", strlen(v->unit), v->unit);
	f->v[i].et = v->v.et;
	if (v->v.et == UWIND)
	  f->uvid = f->v[i].vid;
	else if (v->v.et == VWIND)
	  f->vvid = f->v[i].vid;
	i++;
      }
    }
    f->nv = i;
  }
  nc_enddef(f->fid);
  nc_put_var_double(f->fid, xid, dim1k);
  nc_put_var_double(f->fid, yid, dim2k);
}

void CreateDomainFiles(BOOL reopen)
{
  NestDescr *d;
  char name[100];
  double xstep[1000], ystep[1000];
  int i;
  if (firstdomain)
    printf("\n%sing Domain files...\n", (reopen ? "Reopen" : "Creat"));
  for (d = firstdomain; d; d = d->next)  {
    for (i = d->nnx; i--; )  xstep[i] = i * d->dx;
    for (i = d->nny; i--; )  ystep[i] = i * d->dy;
    sprintf(name, "%s.west.cdf", d->name);
    CreateDomainFile(name, "Y", d->nny, ystep, "Z", nz, zcenter, &d->west, FALSE, reopen);
    sprintf(name, "%s.east.cdf", d->name);
    CreateDomainFile(name, "Y", d->nny, ystep, "Z", nz, zcenter, &d->east, FALSE, reopen);
    sprintf(name, "%s.south.cdf", d->name);
    CreateDomainFile(name, "X", d->nnx, xstep, "Z", nz, zcenter, &d->south, FALSE, reopen);
    sprintf(name, "%s.north.cdf", d->name);
    CreateDomainFile(name, "X", d->nnx, xstep, "Z", nz, zcenter, &d->north, FALSE, reopen);
    sprintf(name, "%s.top.cdf", d->name);
    CreateDomainFile(name, "X", d->nnx, xstep, "Y", d->nny, ystep, &d->top, TRUE, reopen);
  }
  printf("done.\n");
}

double *TestNestBuffer(int size)
{
  static double *buff = NULL;
  static int nbuff = 0;
  if (!nbuff)
    buff = (double *)calloc(nbuff = size, sizeof(double));
  else if (nbuff < size)
    buff = (double *)realloc(buff, size * sizeof(double));
  if (!buff)  {
    fprintf(stderr, "FATAL ERROR: Unable to allocate memory in \"TestNestBuff\"\n");
    exit (3);
  }
  return (buff);
}

double Interpolate(double *mesh, double x, double y, int k)
{
  double xf, yf, xfm, yfm, r;
  int xi, yi, loc;
  x /= dx; y /= dy;
  xi = (int)x; yi = (int)y;
  xf = x - xi; yf = y - yi;
  xfm = 1. - xf; yfm = 1. - yf;
  loc = yi + xi * row + k * layer;
  r = xfm * yfm * mesh[loc] +
      xf  * yfm * mesh[loc+row] +
      xfm * yf  * mesh[loc+1] +
      xf  * yf  * mesh[loc+row+1];
  return (r);
}

void TurnWinds(int n, double *u, double *v, double sina, double cosa)
{
  double tmp;
  for ( ; n--; u++, v++)  {
    tmp = cosa * *u + sina * *v;
    *v = -sina * *u + cosa * *v;
    *u = tmp;
  }
}

void WriteToDomainFile(long actime, long entry, CDFileDescr nc, NestDescr *d,
		       int n1,
		       long xorigin, long yorigin,
		       double dx1, double dy1)
{
  int i, k, h, et, nf;
  size_t start[3], count[3];
  double *buff, *u, *v, *bf, *mv;
  start[0] = entry; start[1] = start[2] = 0;
  count[0] = 1; count[1] = nz; count[2] = n1;
  nf = n1 * nz;
  buff = TestNestBuffer(3 * n1 * nz);
  u = buff + nf; v = u + nf;
  for (h = 0; h < nc.nv; h++)  {
    et = nc.v[h].et;
#ifdef PARALLEL
    if (parallel)
      mv = GetGridVarFromWorker(et);
    else
#endif
      mv = g[et];
    bf = (et > VWIND ? buff : buff + (et+1) * nf);
    for (i = 0; i < n1; i++)
      for (k = 0; k < nz; k++)
	bf[k*n1 + i] = Interpolate(g[et], xorigin + i*dx1, yorigin + i*dy1, k);
    if (et > VWIND)
      nc_put_vara_double(nc.fid, nc.v[h].vid, start, count, bf);
  }
  TurnWinds(nf, u, v, d->sina, d->cosa);
  nc_put_vara_double(nc.fid, nc.uvid, start, count, u);
  nc_put_vara_double(nc.fid, nc.vvid, start, count, v);
  nc_put_var1_long(nc.fid, nc.tid, start, &actime);
  nc_sync(nc.fid);
}

void WritePressureToFile(long actime, long entry, CDFileDescr nc,
			 int n1, int n2,
			 long xorigin, long yorigin,
			 double dx1, double dx2,
			 double dy1, double dy2)
{
  int i, j, h;
  size_t start[3], count[3];
  double *buff;
  start[0] = entry; start[1] = start[2] = 0;
  count[0] = 1; count[1] = n2; count[2] = n1;
  buff = TestNestBuffer(n1 * n2);
  for (i = 0; i < n1; i++)
    for (j = 0; j < n2; j++)
      buff[j*n1+i] = Interpolate(toppress,
			xorigin + i*dx1 + j*dy1, yorigin + i*dx2 + j*dy2, 0);
  nc_put_vara_double(nc.fid, nc.v->vid, start, count, buff);
  nc_put_var1_long(nc.fid, nc.tid, start, &actime);
  nc_sync(nc.fid);
}

void WriteToDomainFiles(long actime)
{
  NestDescr *d;
  long entry = 0;
  for (d = firstdomain; d; d = d->next)  {
    if (actime % d->ndt == 0)  {
      entry = actime / d->ndt;
      WriteToDomainFile(actime, entry, d->west, d, d->nny,
			d->xorigin, d->yorigin,
			- d->dy * d->sina, d->dy * d->cosa);
      WriteToDomainFile(actime, entry, d->east, d, d->nny,
			long(d->xorigin + (d->nnx-1) * d->dx * d->cosa),
			long(d->yorigin + (d->nnx-1) * d->dx * d->sina),
			- d->dy * d->sina, d->dy * d->cosa);
      WriteToDomainFile(actime, entry, d->south, d, d->nnx,
			d->xorigin, d->yorigin,
			d->dx * d->cosa, d->dy * d->sina);
      WriteToDomainFile(actime, entry, d->north, d, d->nnx,
			long(d->xorigin - (d->nny-1) * d->dy * d->sina),
			long(d->yorigin + (d->nny-1) * d->dy * d->cosa),
			d->dx * d->cosa, d->dy * d->sina);
      WritePressureToFile(actime, entry, d->top, d->nnx, d->nny,
			  d->xorigin, d->yorigin,
			  d->dx * d->cosa, d->dy * d->sina,
			  - d->dy * d->sina, d->dy * d->cosa);
    }
  }
}

void CloseDomainFiles(void)
{
  NestDescr *d;
  for (d = firstdomain; d; d = d->next)  {
    nc_close(d->west.fid);
    nc_close(d->east.fid);
    nc_close(d->south.fid);
    nc_close(d->north.fid);
    nc_close(d->top.fid);
  }
}
