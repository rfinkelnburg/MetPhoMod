/*
   IMPLEMENTATION MODULE mcprint
   Dieses Modul ermoeglicht die Erstellung behandelt die Ausgabe der Modell-
   resultate.
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
#include "mcprint.h"
#include "mchemparse.h"
#include "mchem.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

const static double FillValue = -999.;
OutFile *outfile;

VarDesc *AddProductionVariable(char *name)
{
  VarDesc *v, *vp;
  ProductionDesc *p;
  for (vp = variable; vp && strcmp(name+5, vp->name); vp = vp->next);
  if (vp)  {
    if (vp->storetype != GRID_VAL || vp->v.et < SUBS)  {
      InputError(ERROR, "\"%s\" is not a chemical substance.", name+5);
      return (NULL);
    }
    v = (VarDesc *)malloc(sizeof(VarDesc));
    v->name = strdup(name);
    v->unit = "ppb/sec";
    v->comment = "Production of Substance per second";
    v->storetype = DOUBLE_PTR;
    if (!(v->v.d = (double *)calloc(mesh, sizeof(double))))  {
      InputError(ERROR, "Error allocating memory for Variable \"%s\"",
         name);
      return (NULL);
    }
    v->dims = ALL_DIM;
    v->section = v->origsection = CHEMISTRY;
    v->init = CALCULATED_VALUES;
    v->inputtype = NORMAL_NUM;
    v->next = NULL;
    p = (ProductionDesc *)malloc(sizeof(ProductionDesc));
    p->pps = v->v.d;
    p->subs = vp->v.et;
    p->next = prodd;
    prodd = p;
  }
  for (vp = variable; vp->next; vp = vp->next);
  vp->next = v;
  return (v);
}

OutVar *OutItem(char *name, char *prname)
{
  int rdiff;
  OutVar *v;
  VarDesc *vp;
  for (vp = variable; vp && strcmp(name, vp->name); vp = vp->next);
  if (!vp)  {
    if (!strncmp(name, "prod-", 5))  vp = AddProductionVariable(name);
    if (!vp)  {
      InputError(ERROR, "There is no variable \"%s\".", name);
      return (NULL);
    }
  }
  if (!vp->dims && outfile->prtype == CDF_PRINT)  {
    InputError(ERROR, "The variable \"%s\" can not be written into a CDF-File", name);
    return (NULL);
  }
  if (vp->inputtype & ALL_OPTIONS)  {
    InputError(ERROR, "The Option-Variable \"%s\" is not printable.", name);
    return (NULL);
  }
  if (vp->storetype == GRID_VAL && !gactive[vp->v.et])  {
    InputError(ERROR, "Variable \"%s\" is not active in actual calculation.", name);
    return (NULL);
  }
  if (val.v.a.has_dimension & ~vp->dims)  {
    InputError(ERROR, "Inapproriate Dimensions for variable \"%s\".", name);
    return (NULL);
  }
  if (outfile->prtype == CDF_PRINT &&
      (((val.v.a.has_range & X_DIM) && (val.v.a.xmax != val.v.a.xmin)) ||
       ((val.v.a.has_range & Y_DIM) && (val.v.a.ymax != val.v.a.ymin)) ||
       ((val.v.a.has_range & Z_DIM) && (val.v.a.zmax != val.v.a.zmin))))  {
    InputError(ERROR, "Subranges are not allowed in CDF-Files.");
    return (NULL);
  }
  if ((val.v.a.has_range & X_DIM) && (val.v.a.xmin == -1) ||
      (val.v.a.has_range & Y_DIM) && (val.v.a.ymin == -1))  {
    InputError(ERROR, "The keyword \"GROUNDLEVEL\" is only allowed for the z-Dimension.");
    return (NULL);
  }
  if (v = (OutVar *)malloc(sizeof(OutVar)))  {
    v->has_range = val.v.a.has_range;
    rdiff = !printwithborder && !(vp->option & STARTING_WITH_ZERO);
    v->xmin = v->xmax = v->ymin = v->ymax = rdiff;
    v->zmin = v->zmax = 0;
    if (vp->dims & X_DIM)  v->xmax = nx - rdiff;
    if (vp->dims & Y_DIM)  v->ymax = ny - rdiff;
    if (vp->dims & Z_DIM)  v->zmax = nzm;
    if (v->has_range & X_DIM)  {
      v->xmin = val.v.a.xmin;
      v->xmax = val.v.a.xmax;
    }
    if (v->has_range & Y_DIM)  {
      v->ymin = val.v.a.ymin;
      v->ymax = val.v.a.ymax;
    }
    if (v->has_range & Z_DIM)  {
      v->zmin = val.v.a.zmin;
      v->zmax = val.v.a.zmax;
    }
    v->v = vp;
    strncpy(v->printname, prname, 39);
    v->next = NULL;
  }
  else  InputError(ERROR, "Error in allocating memory.");
  return (v);
}

#define ASSMIN(x, y)  if ((y) < (x))  x = y
#define ASSMAX(x, y)  if ((y) > (x))  x = y

OutVar *AppendOutVar(OutVar *first, OutVar *add)
{
  OutVar *o;
  if (first)  {
    if (outfile->prtype == ASCII_PRINT)  {
      if (first->v->dims != add->v->dims)  {
        InputError(ERROR, "Variables within the same ASCI-output File are not compatible.\n"
                   "             (Have different dimensions).");
        return (NULL);
      }
      if (first->has_range & add->has_range)  {
        InputError(ERROR, "The Range of a variable can only be defined once in one ASCII-File");
        return (NULL);
      }
      ASSMAX(first->xmin, add->xmin);
      ASSMIN(first->xmax, add->xmax);
      ASSMAX(first->ymin, add->ymin);
      ASSMIN(first->ymax, add->ymax);
      ASSMAX(first->zmin, add->zmin);
      ASSMIN(first->zmax, add->zmax);
      first->has_range |= add->has_range;
    }
    for (o = first; o->next; o = o->next)
      if (!strcmp(o->printname, add->printname))  {
        if (o->v == add->v)
          InputError(ERROR, "You can not print two variables with the same name \"%s\" within the same file.\n"
                            "(Use \"AS\" option if you want to have two ranges of the same variable!)",
                            add->printname);
        else
          InputError(ERROR, "You can not print two variables with the same name \"%s\" within the same file.\n"
                            "(Improper use of \"AS\" option!)",
                            add->printname);
        return (NULL);
      }
    o->next = add;
    return (first);
  }
  else  {
    return (add);
  }
}

int PrepareOutput(PrintType prtype, char *name, long firstprint, long dprint)
{
  OutFile *of;
  if (actualsection->id != OUTPUT)  {
    InputError(ERROR, "This command has to be called in the OUTPUT-section.");
    return (1);
  }
  of = (OutFile *)malloc(sizeof(OutFile));
  of->prtype = prtype;
  of->nextprint = firstprint;
  of->printinterval = dprint;
  of->vptr = NULL;
  of->ntimerec = 0;
  strcpy(of->fname, name);
  of->next = outfile;
  outfile = of;
  return (0);
}

void AttachOutlist(OutVar *v)
{
  outfile->vptr = v;
}

int CreateNonExistingFile(char *fname)
{
  char buff[255];
  int i;
  FILE *f;
  strcpy(buff, fname);
  i = 0;
  while (f = fopen(buff, "r"))  {
    fclose(f);
    sprintf(buff, "%s.%i", fname, ++i);
  }
  printf("%s ", buff);
  strcpy(fname, buff);
  nc_create(buff, NC_CLOBBER, &i);
  return (i);
}

void CopyWith_(char *a, const char *b)
{
  while (*a++ = (*b == '-' ? '_' : *b))  b++;
}

int OpenCDFFile(char *inpfile, OutFile *o)
{
  int i, j, topovar, boolvar, deltaid, xvar, yvar, zvar;
  int dims[4], deltadim[2];
  size_t count[4] = {3, 2}, startval[4] = {0, 0}, li;
  double numbuff;
  char namebuf[80], timestr[60];
  static float deltatempl[3][2] =  {
    0., 1000., 0., 1000., 0., 100.
  };
  OutVar *v;
  if ((o->fid = CreateNonExistingFile(o->fname)) < 0)  {
    printf("ERROR: Unable to open output-file %s (nccreate)\n", o->fname);
    return (1);
  }
  nc_def_dim(o->fid, "X", nx-1+2*printwithborder, &o->xid);
  nc_def_dim(o->fid, "Y", ny-1+2*printwithborder, &o->yid);
  nc_def_dim(o->fid, "Z", nz, &o->zid);
  nc_def_dim(o->fid, "Time", NC_UNLIMITED, &o->timeid);
  nc_def_dim(o->fid, "SpatialDims", 3, deltadim);
  nc_def_dim(o->fid, "AxIdx", 2, deltadim+1);
  nc_def_var(o->fid, "delta", NC_FLOAT, 2, deltadim, &deltaid);
  nc_put_att_text(o->fid, NC_GLOBAL, "Input", strlen(inpfile)+1, inpfile);
  nc_put_att_text(o->fid, NC_GLOBAL, "Title", strlen(worktitle)+1, worktitle);
  nc_put_att_double(o->fid, NC_GLOBAL, "dX", NC_FLOAT, 1, &dx);
  nc_put_att_double(o->fid, NC_GLOBAL, "dY", NC_FLOAT, 1, &dy);
  nc_put_att_text(o->fid, NC_GLOBAL, "X - Unit", 2, "m");
  nc_put_att_text(o->fid, NC_GLOBAL, "Y - Unit", 2, "m");
  nc_put_att_text(o->fid, NC_GLOBAL, "Z - Unit", 2, "m");
  nc_put_att_text(o->fid, NC_GLOBAL, "Time - Unit", 4, "sec");
  nc_def_var(o->fid, "Time", NC_LONG, 1, &o->timeid, &o->tvarid);
  nc_put_att_text(o->fid, o->tvarid, "Unit", 4, "sec");
  sprintf(timestr, "seconds since %d-%d-%d %d:%d:0 %+d:00",
     startdate.year, startdate.month, startdate.day,
     starttime.hour, starttime.minute, (int)(timezonediff+0.5));
  nc_put_att_text(o->fid, o->tvarid, "units", strlen(timestr), timestr);
  nc_put_att_text(o->fid, o->tvarid, "field", 20, "Time, scalar, series");
  nc_put_att_text(o->fid, o->tvarid, "positions", 14, "delta, compact");
#ifdef FERRET_FORMAT
  nc_def_var(o->fid, "X", NC_FLOAT, 1, &o->xid, &xvar);
  nc_put_att_text(o->fid, xvar, "units", 2, "km");
  nc_def_var(o->fid, "Y", NC_FLOAT, 1, &o->yid, &yvar);
  nc_put_att_text(o->fid, yvar, "units", 2, "km");
  nc_def_var(o->fid, "Z", NC_FLOAT, 1, &o->zid, &zvar);
  nc_put_att_text(o->fid, zvar, "units", 1, "m");
  dims[0] = o->zid; dims[1] = o->yid; dims[2] = o->xid;
  nc_def_var(o->fid, "Topography", NC_FLOAT, 2, dims+1, &topovar);
#else
  dims[0] = o->xid; dims[1] = o->yid; dims[2] = o->zid;
  nc_def_var(o->fid, "Topography", NC_FLOAT, 2, dims, topovar);
#endif
  nc_def_var(o->fid, "POINTSTATUS", NC_LONG, 3, dims, &boolvar);
  nc_put_att_text(o->fid, topovar, "Unit", 1, "m");
  nc_put_att_text(o->fid, topovar, "units",  1, "m");
  nc_put_att_text(o->fid, topovar, "field", 18, "Topography, scalar");
  nc_put_att_text(o->fid, topovar, "positions", 14, "delta, compact");
  nc_put_att_text(o->fid, boolvar, "field", 19, "POINTSTATUS, scalar");
  nc_put_att_text(o->fid, boolvar, "positions", 14, "delta, compact");
  for (v = o->vptr; v; v = v->next)  {
    i = 0; dims[i++] = o->timeid;
#ifdef FERRET_FORMAT
    if (v->v->dims & Z_DIM & ~v->has_range)  dims[i++] = o->zid;
    if (v->v->dims & Y_DIM & ~v->has_range)  dims[i++] = o->yid;
    if (v->v->dims & X_DIM & ~v->has_range)  dims[i++] = o->xid;
    CopyWith_(namebuf, v->printname);
    nc_def_var(o->fid, namebuf, NC_FLOAT, i, dims, &v->datid);
#else
    if (v->v->dims & X_DIM & ~v->has_range)  dims[i++] = o->xid;
    if (v->v->dims & Y_DIM & ~v->has_range)  dims[i++] = o->yid;
    if (v->v->dims & Z_DIM & ~v->has_range)  dims[i++] = o->zid;
    nc_def_var(o->fid, v->printname, NC_FLOAT, i, dims, &v->datid);
#endif
    if (v->v->unit)  {
      nc_put_att_text(o->fid, v->datid, "Unit", strlen(v->v->unit),
         v->v->unit);
      nc_put_att_text(o->fid, v->datid, "units", strlen(v->v->unit),
         v->v->unit);
      nc_put_att_double(o->fid, v->datid, "_FillValue", NC_FLOAT, 1, &FillValue);
      nc_put_att_double(o->fid, v->datid, "missing_value", NC_FLOAT, 1, &FillValue);
      sprintf(namebuf, "%s, scalar, series", v->v->name);
      nc_put_att_text(o->fid, v->datid, "field", strlen(namebuf), namebuf);
      nc_put_att_text(o->fid, v->datid, "positions", 14, "delta, compact");
    }
    o->ntimerec = 0;
  }
  nc_enddef(o->fid);
#ifdef FERRET_FORMAT
  deltatempl[0][1] = level[0]; deltatempl[1][1] = dy; deltatempl[2][1] = dx;
#else
  deltatempl[0][1] = dx; deltatempl[1][1] = dy; deltatempl[2][1] = level[0];
#endif
  nc_put_vara_float(o->fid, deltaid, startval, count, deltatempl[0]);
  startval[1] = 0;
#ifdef FERRET_FORMAT
  for (startval[1] = nxm+2*printwithborder; startval[1]--; )
    for (startval[0] = nym+2*printwithborder; startval[0]--; )
      nc_put_var1_double(o->fid, topovar, startval, topo+(startval[1]+!printwithborder)*row+startval[0]+!printwithborder);
  for (li = nxm+2*printwithborder; li--; )  {
    numbuff = li * dx * 0.001;
    nc_put_var1_double(o->fid, xvar, &li, &numbuff);
  }
  for (li = nym+2*printwithborder; li--; )  {
    numbuff = li * dy * 0.001;
    nc_put_var1_double(o->fid, yvar, &li, &numbuff);
  }
  numbuff = level[0] * 0.5 + reflevel;
  li = 0;
  nc_put_var1_double(o->fid, zvar, &li, &numbuff);
  for (li = 1; li < nz; li++)  {
    numbuff += 0.5 * (level[li-1] + level[li]);
    nc_put_var1_double(o->fid, zvar, &li, &numbuff);
  }
  for (startval[0] = nz; startval[0]--; )
    for (startval[1] = nym+2*printwithborder; startval[1]--; )
      for (startval[2] = nxm+2*printwithborder; startval[2]--; )
        nc_put_var1_long(o->fid, boolvar, startval,
           &pstat[startval[0]*layer + (startval[2]+!printwithborder)*row + startval[1]+!printwithborder]);
#else
  count[0] = 1; count[1] = nym + 2*printwithborder; count[2] = 1;
  for (*startval = nxm+2*printwithborder; (*startval)--; )
    nc_put_vara_double(o->fid, topovar, startval, count, topo+(startval[0]+!printwithborder)*row+!printwithborder);
  for (startval[0] = nxm+2*printwithborder; startval[0]--; )
    for (startval[1] = nym+2*printwithborder; startval[1]--; )
      for (startval[2] = nz; startval[2]--; )
        nc_put_var1_long(o->fid, boolvar, startval,
           &pstat[startval[2]*layer + (startval[0]+!printwithborder)*row + startval[1]+!printwithborder]);
#endif
  return (0);
}

int ReopenCDFFile(char *inpfile, OutFile *o)
{
  int i, j;
  int dims[4];
  char txt[MAX_NC_NAME], namebuf[80];
  long count[4], startval[4];
  OutVar *v;
  if (nc_open(o->fname, NC_WRITE, &o->fid))  {
    printf("Unable to append to file \"%s\". Creating a new one.\n", o->fname);
    return (OpenCDFFile(inpfile, o));
  }
  if (nc_inq_dimid(o->fid, "X", &o->xid) ||
      nc_inq_dimid(o->fid, "Y", &o->yid) ||
      nc_inq_dimid(o->fid, "Z", &o->zid) ||
      nc_inq_dimid(o->fid, "Time", &o->timeid))  {
    fprintf(stderr, "ERROR: %s has wrong dimensions\n", o->fname);
    return (1);
  }
  if (nc_inq_varid(o->fid, "Time", &o->tvarid))  {
    fprintf(stderr, "ERROR: I couldn't find Variable \"Time\" in File %s.\n", o->fname);
    return (1);
  }
  for (v = o->vptr; v; v = v->next)  {
    CopyWith_(namebuf, v->printname);
    if (nc_inq_varid(o->fid, namebuf, &v->datid))  {
      fprintf(stderr, "ERROR: I couldn't find Variable \"%s\" in File %s.\n", v->v->name, o->fname);
      return (1);
    }
  }
  nc_inq_dim(o->fid, o->timeid, NULL, &o->ntimerec);
  o->ntimerec--;
  nc_get_var1_long(o->fid, o->tvarid, &o->ntimerec, &o->nextprint);
  o->ntimerec++;
  o->nextprint += o->printinterval;
  printf("Appending data to %s.\n", o->fname);
  return (0);
}

int OpenOutputFiles(char *inpfile, BOOL reopen)
{
  OutFile *o;
  OutVar *v;
  int len, i;
  for (o = outfile; o; o = o->next)  {
    if (o->prtype == CDF_PRINT)  {
      if ((reopen && ReopenCDFFile(inpfile, o)) ||
          (!reopen && OpenCDFFile(inpfile, o)))  return (1);
    }
    else  {
      if (!(o->f = fopen(o->fname, "w")))  {
        fprintf(stderr, "\nERROR opening output-file \"%s\"\n", o->fname);
        return (1);
      }
      len = strlen(worktitle);
      for (i = len+4; i--; )
        fprintf(o->f, "*");
      fprintf(o->f, "\n");
      fprintf(o->f, "* %s *\n", worktitle);
      for (i = len+4; i--; )
        fprintf(o->f, "*");
      fprintf(o->f, "\n\n");
      printf("%s ", o->fname);
      fprintf(o->f, "%s : %s\n\n", inpfile, o->fname);
      for (v = o->vptr; v; v = v->next)  {
        fprintf(o->f, "%s : %s", v->printname, v->v->comment);
        if (v->v->unit)
          fprintf(o->f, " (%s)", v->v->unit);
        fprintf(o->f, "\n");
      }
    }
  }
  return (0);
}

void CloseOutputFiles(void)
{
  OutFile *o;
  for (o = outfile; o; o = o->next)
    if (o->prtype == CDF_PRINT)
      nc_close(o->fid);
    else
      fclose(o->f);
}

double *GetVal(VarDesc *v, int i, int j, int k, BOOL cache)
{
  static double resval;
  if (k == -1)  k = ground[i*row+j].firstabove;
  switch (v->storetype)  {
    case GRID_VAL	:
#ifdef PARALLEL
       if (parallel)  {
         resval = GetValFromWorker(v, i, j, k, cache);
         return (&resval);
       }
       else
#endif
       return (&g[v->v.et][k*layer+i*row+j]);
    case DOUBLE_PTR	:
#ifdef PARALLEL
       if (parallel)  {
         resval = GetValFromWorker(v, i, j, k, cache);
         return (&resval);
       }
       else
#endif
         switch (v->dims)  {
            case ALL_DIM   : return (&v->v.d[k*layer+i*row+j]);
            case X_DIM | Y_DIM : return (&v->v.d[i*row+j]);
            case Y_DIM | Z_DIM : return (&v->v.d[k*row+j]);
            case X_DIM | Z_DIM : return (&v->v.d[k*xrow+i]);
            case X_DIM	: return (&v->v.d[i]);
            case Y_DIM	: return (&v->v.d[j]);
            case Z_DIM	: return (&v->v.d[k]);
          }
       break;
    case GROUND_PARAM	:
#ifdef PARALLEL
       if (parallel && cache)  {
         GetGroundFromWorker();
       }
#endif
       return (&v->v.g[i*row+j].Tg[0]);
    case PROC_VAL       :
#ifdef PARALLEL
       if (parallel)  {
         resval = GetValFromWorker(v, i, j, k, cache);
       }
       else
#endif
         resval = v->v.proc(k, i, j, v);
       return (&resval);
    default : fprintf(stderr, "FATAL: Internal error in mcprint\n");
    	      exit (1);
  }
  return (&resval);
}

#define MAX(a,b) ((a) < (b) ? (b) : (a))

int WriteOutData(long actime, BOOL ignoretime)
{
  OutFile *o;
  OutVar *v;
  size_t coord[4], *cptr[3];
  int ncoord, i, j, k, dim, vardim;
  for (o = outfile; o; o = o->next)  {
    if (actime >= o->nextprint || ignoretime)  {
      if (o->prtype == CDF_PRINT)  {
        coord[0] = o->ntimerec++;
        nc_put_var1_long(o->fid, o->tvarid, coord, &actime);
        ncoord = 1;
        for (v = o->vptr; v; v = v->next)  {
          ncoord = 1;
          vardim = v->v->dims & ~v->has_range;
#ifdef FERRET_FORMAT
          for (dim = Z_DIM, i = 0; dim >= X_DIM; dim >>= 1, i++)
            if (vardim & dim)  cptr[i] = coord + ncoord++;
          for (dim = Z_DIM, i = 0; dim >= X_DIM; dim >>= 1, i++)
            if (~vardim & dim)  cptr[i] = coord + ncoord++;
          for (k = v->zmin; k <= v->zmax; k++)  {
            *cptr[0] = k;
            for (j = v->ymin; j <= v->ymax; j++)  {
              *cptr[1] = j-!printwithborder;
              for (i = v->xmin; i <= v->xmax; i++)  {
                *cptr[2] = i-!printwithborder;
                if (v->v->dims != ALL_DIM || k < 0 || k >= ground[i*row+j].firstabove)
                  nc_put_var1_double(o->fid, v->datid, coord, GetVal(v->v, i, j, k, TRUE));
/*                else
                  nc_put_var1_double(o->fid, v->datid, coord, &FillValue);  */
              }
            }
          }
#else
          for (dim = X_DIM, i = 0; dim <= Z_DIM; dim <<= 1, i++)
            if (vardim & dim)  cptr[i] = coord + ncoord++;
          for (dim = X_DIM, i = 0; dim <= Z_DIM; dim <<= 1, i++)
            if (~vardim & dim)  cptr[i] = coord + ncoord++;
          for (k = v->zmin; k <= v->zmax; k++)  {
            *cptr[2] = k;
            for (i = v->xmin; i <= v->xmax; i++)  {
              *cptr[0] = i-!printwithborder;
              for (j = v->ymin; j <= v->ymax; j++)  {
                *cptr[1] = j-!printwithborder;
                nc_put_var1_double(o->fid, v->datid, coord, GetVal(v->v, i, j, k, TRUE));
              }
            }
          }
#endif
        }
        nc_sync(o->fid);
      }
      else  {
        fprintf(o->f, "\nTIME : %li\n\n           ", actime);
        for (v = o->vptr; v; v = v->next)  {
          fprintf(o->f, "%15s", v->v->name, v->v->comment);
        }
        fprintf(o->f, "\n\n");
        v = o->vptr;
        for (k = v->zmin; k <= v->zmax; k++)
          for (i = v->xmin; i <= v->xmax; i++)
            for (j = v->ymin; j <= v->ymax; j++)  {
              fprintf(o->f, "%2i/%2i/%2i : ", k, i, j);
              for ( ; v; v = v->next)  {
                fprintf(o->f, "%15.6le", *GetVal(v->v, i, j, k, FALSE));
              }
              fprintf(o->f, "\n");
              v = o->vptr;
           }
      }
      if (!ignoretime)
        while (o->nextprint <= actime)
          o->nextprint += o->printinterval;
    }    /* of if (actime ...  */
  }     /* of for (o = ....  */
  return (0);
}

void FlushOutputFiles(void)
{
  OutFile *o;
  for (o = outfile; o; o = o->next)
    if (o->prtype == ASCII_PRINT)  fflush(o->f);
}
