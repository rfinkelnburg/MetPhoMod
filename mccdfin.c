/*
   IMPLEMENTATION MODULE mccdfin.c
   Modul zum Einlesen von Daten aus einem CDF-File
*/

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mccdfin.h"
#include "mcemiss.h"
#ifdef PARALLEL
#include <pvm3.h>
#include "mcparallel.h"
#endif

typedef struct BorderTimeDesc  {
  BorderVarType  vartype;
  InputSection section;
  BOOL bigdim;
  int ntime, itime, cdfid, vid;
  long actime;
  VarDesc actualvar, nextvar;
  long *timetable;
  size_t coords[5], *dimptr[5];
  double *actualdata, *nextdata;
  struct BorderTimeDesc *next;
}  BorderTimeDesc;

BorderTimeDesc *bordertime;

#define F(a,b,c)  ((a)*layer+(b)*row+(c))
#define FH(a,b)   ((a)*row+(b))


#ifdef PARALLEL

void SendBorderTime(BOOL leftest, BOOL rightest)
{
  BorderTimeDesc *bt;
  int nbt = 0, i, dff;
  for (bt = bordertime; bt; bt = bt->next)
    nbt += (bt->section != WEST_BORDER || leftest)  &&
           (bt->section != EAST_BORDER || rightest);
  pvm_pkint(&nbt, 1, 1);
  for (bt = bordertime; bt; bt = bt->next)  {
    if ((bt->section == WEST_BORDER && !leftest)  ||
        (bt->section == EAST_BORDER && !rightest))  continue;
    pvm_pkint((int *)&bt->vartype, 1, 1);
    pvm_pkint((int *)&bt->section, 1, 1);
    pvm_pkint(&bt->bigdim, 1, 1);
    pvm_pkint(&bt->ntime, 1, 1);
    pvm_pkint(&bt->itime, 1, 1);
    pvm_pkint(&bt->cdfid, 1, 1);
    pvm_pkint(&bt->vid, 1, 1);
    pvm_pkstr(bt->actualvar.name);
    pvm_pkint((int *)&bt->actualvar.v.et, 1, 1);
    pvm_pkushort(&bt->actualvar.dims, 1, 1);
    pvm_pkint((int *)&bt->actualvar.storetype, 1, 1);
    pvm_pkint(&bt->actualvar.id, 1, 1);
    pvm_pkint(&bt->actualvar.ncoord, 1, 1);
    pvm_pkint((int *)&bt->nextvar.v.et, 1, 1);
    pvm_pkushort(&bt->nextvar.dims, 1, 1);
    pvm_pkint((int *)&bt->nextvar.storetype, 1, 1);
    pvm_pklong(bt->timetable, bt->ntime, 1);
    pvm_pkuint((unsigned *)bt->coords, 4, 1);
    for (i = 0; i < 4; i++)  {
      dff = bt->dimptr[i] - bt->coords;
      pvm_pkint(&dff, 1, 1);
    }
  }
}

void RecvBorderTime(void)
{
  BorderTimeDesc *bt, *last = NULL;
  char varname[80];
  int nbt, mem, i, dff;
  bordertime = NULL;
  pvm_upkint(&nbt, 1, 1);
  while (nbt--)  {
    bt = (BorderTimeDesc *)malloc(sizeof(BorderTimeDesc));
    pvm_upkint((int *)&bt->vartype, 1, 1);
    pvm_upkint((int *)&bt->section, 1, 1);
    pvm_upkint(&bt->bigdim, 1, 1);
    pvm_upkint(&bt->ntime, 1, 1);
    pvm_upkint(&bt->itime, 1, 1);
    pvm_upkint(&bt->cdfid, 1, 1);
    pvm_upkint(&bt->vid, 1, 1);
    pvm_upkstr(varname);
    pvm_upkint((int *)&bt->actualvar.v.et, 1, 1);
    pvm_upkushort(&bt->actualvar.dims, 1, 1);
    pvm_upkint((int *)&bt->actualvar.storetype, 1, 1);
    pvm_upkint(&bt->actualvar.id, 1, 1);
    pvm_upkint(&bt->actualvar.ncoord, 1, 1);
    pvm_upkint((int *)&bt->nextvar.v.et, 1, 1);
    pvm_upkushort(&bt->nextvar.dims, 1, 1);
    pvm_upkint((int *)&bt->nextvar.storetype, 1, 1);
    bt->timetable = (long *)malloc(bt->ntime * sizeof(long));
    pvm_upklong(bt->timetable, bt->ntime, 1);
    pvm_upkuint(bt->coords, 4, 1);
    for (i = 0; i < 4; i++)  {
      pvm_upkint(&dff, 1, 1);
      bt->dimptr[i] = bt->coords + dff;
    }
    switch (bt->section)  {
      case NORTH_BORDER :
            bt->actualvar.v.d = northborder + GetNamedVar(varname)->v.et*nz*xrow;
	    break;
      case SOUTH_BORDER :
            bt->actualvar.v.d = southborder + GetNamedVar(varname)->v.et*nz*xrow;
	    break;
      case WEST_BORDER  :
            bt->actualvar.v.d = westborder + GetNamedVar(varname)->v.et*nz*row;
	    break;
      case EAST_BORDER  :
            bt->actualvar.v.d = eastborder + GetNamedVar(varname)->v.et*nz*row;
	    break;
      case EMISSIONS    :
            bt->actualvar.v.d = GetEmVarWithID(bt->actualvar.id)->v.d;
	    break;
      default	        :
            bt->actualvar.v = GetNamedVar(varname)->v;
            break;
    }
    switch (bt->vartype)  {
      case XWALL_VAR :
	 mem = xrow*nz;
	 break;
      case WALL_VAR :
	 mem = row*nz;
	 break;
      case LAYER_VAR :
	 mem = layer;
	 break;
      case PROFILE_VAR :
	 mem = nz;
	 break;
      case MESH_VAR :
	 mem = mesh;
	 break;
      case COORD_VAR :
	 mem = bt->actualvar.ncoord;
	 break;
    }
    bt->actualdata = bt->actualvar.v.d;
    bt->nextdata = (double *)calloc(mem, sizeof(double));
    bt->nextvar.v.d = bt->nextdata;
    bt->nextvar.ncoord = bt->actualvar.ncoord;
    *(bordertime ? &last->next : &bordertime) = bt;
    last = bt;
  } 
  if (last)  last->next = NULL;
}

#endif

void TestIfInRange(VarDesc *v, double tmp)
{
  static VarDesc *vstat = NULL;
  static int vcount = 0;
  if (v != vstat)  {
    vstat = v; vcount = 0;
  }
  if (tmp != -999. && (tmp < v->rmin || tmp > v->rmax) && vcount < 12)  {
    if (vcount == 11)  InputError(ERROR, "and more...");
    else
      InputError(ERROR, "Variable %s (%.2lg) in netCDF is out of range (%.2lg - %.2lg).",
         v->name, tmp, v->rmin, v->rmax);
    inputerror = 1;
    vcount++;
  }
}

#define NCGET(iv)  \
  { \
    nc_get_var1_double(cdfid, vid, coords, &tmp);  \
    TestIfInRange(v, tmp);  \
    iv = tmp;  \
  }

int ReadVarFromCDF(int cdfid, int vid, VarDesc *v, size_t **d, size_t *coords, BOOL bigdim)
{
  int lnx, lny, offs, startloop;
  int startoff = ((v->option & STARTING_WITH_ZERO) != 0);
  double tmp;
  nc_type vartype;
  nc_inq_vartype(cdfid, vid, &vartype);
  lnx = nxm + bigdim + startoff; lny = nym + bigdim + startoff;
  offs = !startoff && !bigdim;
  startloop = !startoff && bigdim;
  switch (v->storetype)  {
    case GRID_VAL :
       if (!g[v->v.et])  g[v->v.et] = (double *)calloc(mesh, sizeof(double));
       for (*d[2] = 0; *d[2] < nz; (*d[2])++)
         for (*d[0] = startloop; *d[0] < lnx; (*d[0])++)
           for (*d[1] = startloop; *d[1] < lny; (*d[1])++)
             NCGET(g[v->v.et][F(*d[2],*d[0]+offs,*d[1]+offs)]);
       break;
    case LAYER_PTR :
       if (v->dims != ALL_DIM)  {
         InputError(ERROR, "Internal implementation error (3)");
         return (1);
       }
       for (*d[2] = 0; *d[2] < nz; (*d[2])++)
         for (*d[0] = startloop; *d[0] < lnx; (*d[0])++)
           for (*d[1] = startloop; *d[1] < lny; (*d[1])++)
             NCGET(v->v.d[F(*d[2],*d[0]+offs,*d[1]+offs)]);
       break;
    case FILE_PTR :
       switch (v->dims)  {
         case (X_DIM | Y_DIM) :
            for (*d[0] = startloop; *d[0] < lnx; (*d[0])++)
              for (*d[1] = 0; *d[1] < lny; (*d[1])++)
                NCGET(v->v.d[FH(*d[0]+offs,*d[1]+offs)]);
            break;
         case (Y_DIM | Z_DIM) :
            for (*d[2] = 0; *d[2] < nz; (*d[2])++)
              for (*d[1] = startloop; *d[1] < lny; (*d[1])++)
                NCGET(v->v.d[FH(*d[2],*d[1]+offs)]);
            break;
         default : InputError(ERROR, "Internal implementation error (2)");
                   break;
       }
       break;
    case XFILE_PTR :
       for (*d[2] = 0; *d[2] < nz; (*d[2])++)
	 for (*d[0] = startloop; *d[0] < lnx; (*d[0])++)
	   NCGET(v->v.d[*d[2] * xrow + *d[0]+offs]);
       break;
    case DOUBLE_PTR :
       switch (v->dims)  {
         case X_DIM : for (*d[0] = startloop; *d[0] < lnx; (*d[0])++)
           		NCGET(v->v.d[*d[0]+offs]);
           	      break;
         case Y_DIM : for (*d[1] = startloop; *d[1] < lny; (*d[1])++)
           		NCGET(v->v.d[*d[1]+offs]);
                      break;
         case Z_DIM : for (*d[2] = 0; *d[2] < nz; (*d[2])++)
           		NCGET(v->v.d[*d[2]]);
           	      break;
         case COUNT_DIM :
                      for (*d[4] = 0; *d[4] < v->ncoord; (*d[4])++)
                        NCGET(v->v.d[*d[4]]);
                      break;
         case 0     : NCGET(*v->v.d); break;
         default : InputError(ERROR, "Internal implementation error (1)");
           	   return (1);
       }
       break;
    case GROUND_PARAM :
       for (*d[0] = startloop; *d[0] < lnx; (*d[0])++)
         for (*d[1] = startloop; *d[1] < lny; (*d[1])++)
           NCGET(v->v.g[(*d[0]+offs)*row+*d[1]+offs].Tg[0]);
       break;
    default :
       InputError(ERROR, "Variable %s cannot be set via CDF-File. (Wrong-Type)", v->name);
       return (1);
  }
  return (0);
}

int ReadCDFFile(char *name, InputSection section)
{
  int cdfid, i, j, *dimdim, dimfound, ndim, nvar, natt, recdim, *dim,
    nvarfound, bigdim = 0, coordfound = 0, mem, pointvalues;
  static int ncoord = 0;
  FILE *f;
  BOOL timedvar, noclose = FALSE;
  VarDesc *v, vh, *vo;
  Coordinate *coord = NULL;
  size_t start[1], count[1];
  ptrdiff_t stride[1], imap[1];
  const struct  {
    char name[5];
    Dimension dim;
    int *size;
  }  dimtab[] = {"X", X_DIM, &nxm, "Y", Y_DIM, &nym, "Z", Z_DIM, &nz, "Time", TIME_DIM, NULL,
  	         "CNT",
  	         COUNT_DIM,
  	         &ncoord};
  char txt[MAX_NC_NAME];
  BorderTimeDesc *bt;
  size_t *coords, len, *dimptr[5], lcoord[5], ntimes;
  long *timetable = NULL;
  double *doubtable;
  nc_type datatype;
  if (section == GRID || section == TIME || section == DEPOSITION)  {
    InputError(ERROR, "CDFIMPORT is not allowed in this section.");
    return (1);
  }
  if (!(f = fopen(name, "r")))  {
    InputError(ERROR, "Unable to open CDF-File \"%s\".", name);
    return (1);
  }
  fclose(f);
  if (nc_open(name, NC_NOWRITE, &cdfid))  {
    InputError(ERROR, "Unable to open CDF-File %s.", name);
    return (1);
  }
  nc_inq(cdfid, &ndim, &nvar, &natt, &recdim);
  dim = (int *)calloc(ndim, sizeof(int));
  dimdim = (int *)calloc(ndim, sizeof(int));
  nvarfound = dimfound = 0;
  for (i = ndim; i--; )  {
    nc_inq_dim(cdfid, i, txt, &len);
    for (j = 5; --j >= 0 && strcmp(txt, dimtab[j].name); );
    dimdim[i] = j;
    if (j >= 0 && dimtab[j].size)  {
      if (j < 2)  {
        if (len != *dimtab[j].size && len != *dimtab[j].size+2)  {
          InputError(WARNING, "\"%s\": Dimension %s has wrong size. Is %i should be %i or %i.",
             name, txt, len, *dimtab[j].size, *dimtab[j].size+2);
          dimdim[i] = -1;
        }
        else  {
          dimfound |= 1 << j;
          if (len == *dimtab[j].size+2)  bigdim |= 1 << j;
        }
      }
      else if (j < 4)  {
        if (len != *dimtab[j].size)  {
          InputError(WARNING, "\"%s\": Dimension %s has wrong size. Is %i should be %i.",
             name, txt, len, *dimtab[j].size);
          dimdim[i] = -1;
	}
	else  dimfound |= 1 << j;
      }
      else  {
        *dimtab[j].size = len;
        dimfound |= 1 << j;
      }
    }
    if (j == 3)  ntimes = len;
  }
  if (bigdim && bigdim != ((X_DIM | Y_DIM) & dimfound))  {
    InputError(ERROR, "Mismatch of X/Y-Dimension-sizes.");
    goto error_end;
  }
  bigdim = !!bigdim;
  for (i = nvar; i--; )  {
    nc_inq_var(cdfid, i, txt, &datatype, &ndim, dim, &natt);
    for (v = variable; v && strcmp(txt, v->name); v = v->next);
    if (v)  {
      vo = v;
      if (v->storetype == GRID_VAL &&
	  actualsection >= NORTH_BORDER && actualsection <= EAST_BORDER)
	v = BorderVariable(v, &vh);
      if (v->storetype == GRID_VAL && actualsection == EMISSIONS)  {
        pointvalues = (ndim < 3 && (dimdim[*dim] == 4 || dimdim[dim[1]] == 4) ? *dimtab[4].size : 0);
	if (pointvalues && !coord)  coord = AllocateCoordVar(pointvalues);
        v = EmissionVariable(v, pointvalues, 0, coord);
      }
      if (v->section == section)  {
	if (v->inputtype != NORMAL_NUM)  {
	  InputError(WARNING, "The Variable %s cannot be set via a CDF-File.\n",
	     inplineno, txt);
	  continue;
	}
	if (datatype != NC_DOUBLE && datatype != NC_FLOAT && datatype != NC_LONG)  {
	  InputError(ERROR, "\"%s\": The type of variable %s is not of appropriate type.",
	     name, txt);
	  goto error_end;
	}
	if (v->init == WAS_SET_BY_USER)  {
	  InputError(WARNING, "The variable \"%s\" found in file \"%s\" is already initialized.",
	     txt, name);
	  if (section == EMISSIONS)
	    InputError(WARNING, " ....so summing up emissions.\n");
	  else
	    continue;
	}
	if (v->init == CALCULATED_VALUES)  {
	  InputError(WARNING, "The variable \"%s\" found in file \"%s\" cannot be initialized.\n"
	             "         It's calculated by meteochem.", txt, name);
	  continue;
	}
	if ((v->option & STARTING_WITH_ZERO) && !bigdim)  {
	  InputError(ERROR, "Variable \"%s\" requires Dimensions X and Y to be two fields bigger",
	     v->name);
	  goto error_end;
	}
	memset(dimptr, 0, 5 * sizeof(*dimptr));
	for (j = ndim; --j >= 0 && dimdim[dim[j]] != 3; );
	if (timedvar = j >= 0)  {
	  if (actualsection < ENVIRONMENT || actualsection == INITIAL_DATA)  {
	    InputError(ERROR, "Time-dependent variables are not allowed in this section.");
	    goto error_end;
	  }
	  noclose = TRUE;
	  bt = (BorderTimeDesc *)malloc(sizeof(BorderTimeDesc));
	  if (v->dims == (X_DIM | Z_DIM))
	    bt->vartype = XWALL_VAR;
	  else if (v->dims == (Y_DIM | Z_DIM))
	    bt->vartype = WALL_VAR;
	  else if (v->dims == (X_DIM | Y_DIM))
	    bt->vartype = LAYER_VAR;
	  else if (v->dims == Z_DIM)
	    bt->vartype = PROFILE_VAR;
	  else if (v->dims == ALL_DIM)
	    bt->vartype = MESH_VAR;
	  else if (v->dims == COUNT_DIM)
	    bt->vartype = COORD_VAR;
	  else  {
	    InputError(ERROR, "The variable \"%s\" must not be time-dependent.", v->name);
	    goto error_end;
	  }
	  bt->section = actualsection;
	  bt->bigdim = bigdim;
	  bt->ntime = ntimes;
	  bt->itime = -1;
	  bt->actime = 0;
	  bt->cdfid = cdfid;
	  bt->vid = i;
	  bt->actualvar = *v;
	  bt->nextvar = *v;
	  bt->timetable = timetable;
	  coords = bt->coords;
	  bt->actualdata = v->v.d;
	  switch (bt->vartype)  {
	    case XWALL_VAR :
	       mem = xrow*nz;
	       break;
	    case WALL_VAR :
	       mem = row*nz;
	       break;
	    case LAYER_VAR :
	       mem = layer;
	       break;
	    case PROFILE_VAR :
	       mem = nz;
	       break;
	    case MESH_VAR :
	       mem = mesh;
	       break;
	    case COORD_VAR :
	       mem = v->ncoord;
	       break;
	  }
	  bt->nextdata = (double *)calloc(mem, sizeof(double));
	  if (!(bt->nextvar.v.d = bt->nextdata))  {
	    InputError(ERROR, "Unable to allocate memory in Function \"ReadCDFFile\"");
	    goto error_end;
	  }
	  bt->next = bordertime;
	  bordertime = bt;
	}
	else  coords = lcoord;
	for (j = ndim; j--; )  {
	  if (dimdim[dim[j]] < 0)  {
	    nc_inq_dim(cdfid, dim[j], txt, &len);
	    InputError(ERROR, "\"%s\": Variable %s includes \"%s\", an unkown or unusable dimension.",
	       name, v->name, txt);
	    goto error_end;
	  }
	  dimptr[dimdim[dim[j]]] = coords + j;
	}
	for (j = 5; j--; )
	  if (!dimptr[j])  dimptr[j] = coords + ndim++;
	printf("Reading Variable %s from CDF-File \"%s\"\n", txt, name);
	*dimptr[3] = 0;
	nvarfound++;
	vo->init = WAS_SET_BY_USER;
	if (ReadVarFromCDF(cdfid, i, v, dimptr, coords, !!bigdim))  goto error_end;
	if (timedvar)  memcpy(bt->dimptr, dimptr, 5 * sizeof(long *));
      }
      else if (!strcmp(txt, "Time"))  {
        timetable = (long *)calloc(ntimes, sizeof(long));
        *lcoord = 0;
        if (nc_get_vara_long(cdfid, i, lcoord, &ntimes, timetable))  {
          InputError(ERROR, "A variable called \"Time\" was found, but I couldn't read it.");
        }
        else  {
          if (*timetable > tstart)
            InputError(WARNING, "The first time-slice for value for %s is later (%ld sec) than tstart(%ld sec).",
               v->name, *timetable, tstart);
          for (bt = bordertime; bt; bt = bt->next)
            if (!bt->timetable)  bt->timetable = timetable;
        }
      }
    }
    else if (!strcmp(txt, "XCoord") || !strcmp(txt, "YCoord") || !strcmp(txt, "ZCoord"))  {
      if (ndim == 1 && dimdim[*dim] == 4)  {
	if (!coord)  coord = AllocateCoordVar(*dimtab[4].size);
	*start = 0;
	*count = *dimtab[4].size;
	*stride = 1;
	*imap = &coord[1].x - &coord[0].x;
	nc_get_varm_float(cdfid, i, start, count, stride, imap,
	      &coord[0].x + (*txt - 'X'));
	coordfound |= 1 << (*txt - 'X');
      }
    }
  }
  if (coord)
    if (coordfound != 7)  {
      InputError(ERROR, "Not all necessary coordinates found in file \"%s\".", name);
      goto error_end;
    }
    else
      ConvertEmissionCoords(*dimtab[4].size, coord, 0);
  if (!noclose)  nc_close(cdfid);
  free(dim);
  if (!nvarfound)
    InputError(WARNING, "No usable Variable found in File \"%s\".", name);
  return (0);
error_end :
  nc_close(cdfid);
  free(dim); free(dimdim);
  return (1);
}

void InitCDFInit(void)
{
  bordertime = NULL;
  ncopts = NC_VERBOSE;
}

void ActualiseValuesInTime(long actime)
{
  double ti;
  BorderTimeDesc *bt;
  int i, dtime;
  for (bt = bordertime; bt; bt = bt->next)  {
    if (bt->itime+1 < bt->ntime)  {
      switch (bt->vartype)  {
        case WALL_VAR :
           i = row*nz;
           break;
        case XWALL_VAR :
	   i = xrow*nz;
	   break;
	case LAYER_VAR :
	   i = layer;
	   break;
	case PROFILE_VAR :
	   i = nz;
	   break;
	case MESH_VAR :
	   i = mesh;
	   break;
	case COORD_VAR :
	   i = bt->actualvar.ncoord;
	   break;
      }
      if (actime < bt->timetable[bt->itime+1])  {
        dtime = actime - bt->actime;
        while (i--)
          bt->actualdata[i] += dtime * bt->nextdata[i];
        bt->actime = actime;
      }
      else  {
        while (bt->itime+1 < bt->ntime && actime >= bt->timetable[bt->itime+1])
          *bt->dimptr[3] = ++bt->itime;
#ifdef PARALLEL
        if (parallel)
          if (master)  {
            ReadVarFromCDF(bt->cdfid, bt->vid, &bt->actualvar, bt->dimptr, bt->coords, bt->bigdim);
	    if (bt->vartype == COORD_VAR)
	      PackCoordVariable(&bt->actualvar);
	    else
	      PackVariable(bt->actualvar.v.d, bt->vartype, bt->section);
          }
          else  {
            UnpackVariable(bt->actualvar.v.d);
          }
        else
#endif
          ReadVarFromCDF(bt->cdfid, bt->vid, &bt->actualvar, bt->dimptr, bt->coords, bt->bigdim);
        if (bt->itime + 1 < bt->ntime)  {
          *bt->dimptr[3] = bt->itime+1;
#ifdef PARALLEL
          if (parallel)
            if (master)  {
              ReadVarFromCDF(bt->cdfid, bt->vid, &bt->nextvar, bt->dimptr, bt->coords, bt->bigdim);
              if (bt->vartype == COORD_VAR)
		PackCoordVariable(&bt->nextvar);
	      else
	        PackVariable(bt->nextvar.v.d, bt->vartype, bt->section);
            }
            else  {
              UnpackVariable(bt->nextvar.v.d);
            }
          else
#endif
            ReadVarFromCDF(bt->cdfid, bt->vid, &bt->nextvar, bt->dimptr, bt->coords, bt->bigdim);
          dtime = bt->timetable[bt->itime+1] - bt->timetable[bt->itime];
	  ti = 1. / (double)dtime;
	  dtime = actime - bt->timetable[bt->itime];
	  while (i--)  {
            bt->nextdata[i] = (bt->nextdata[i] - bt->actualdata[i]) * ti;
	    bt->actualdata[i] += dtime * bt->nextdata[i];
	  }
	  bt->actime = actime;
        }
        else  {
#ifdef PARALLEL
          if (master)
#endif
            fprintf(stderr, "WARNING: no further values for variable \"%s\" at time %li\n",
               bt->actualvar.name, actime);
        }
      }
    }
  }
}
