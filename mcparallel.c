/*
   MODULE mcparallel.c
   Die meteochem-Schnittstelle zur PVM3-Bibliothek.
*/


#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <sys/times.h>
#include <limits.h>
#include <pvm3.h>
#include <pvmtev.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mccdfin.h"
#include "mchemparse.h"
#include "mchem.h"
#include "mcemiss.h"
#include "mcdep.h"
#include "mcpress.h"
#include "mcparallel.h"
#include "mc_module.hh"

extern "C" int pvm_catchout(FILE *);

#define MAXWORKER 20
#ifdef IRIX
#define McInPlace PvmDataRaw
#else
#define McInPlace PvmDataInPlace
#endif

/* global */
int master, workers, parallel, myparent, leftbrother, rightbrother,
    mfirstx, mnx, mlastx, xoff;
BOOL leftest, rightest;
char nodename[80];

Worker worker[MAXWORKER];

/* private */
int mytid;
double *meshcache;

extern int worknum; /* Ist in mcglobal.c definiert */
static int isdebug;
static long waitingtime = 0;

void ReceiveWithTimeout(int tid, int tag)
{
  static struct timeval tmout = {120L, 0L};
  struct tms buffer;
  int i, h;
  if (isdebug)  {
    pvm_recv(tid, tag);
    return;
  }
  for (i = 5; i--; )  {
/*    waitingtime -= times(&buffer); */
    h = pvm_trecv(tid, tag, &tmout);
/*    waitingtime += times(&buffer); */
    if (h > 0)  return;
    printf("\nI got no message from worker %x for 2 minutes.\n", tid);
    if (pvm_nrecv(tid, NOTIFYDATA) > 0)  {
      pvm_upkint(&h, 1, 1);
      printf("It seems, that worker %x died without known reason, so\n", h);
      goto terminate;
    }
    if (i)  printf("I didn't find any reason, so wait another 2 minutes...\n");
  }
  printf("Still no message from worker %x. I don't know what happened!\n", tid);
terminate :
  printf("Stopping calculation, by sending TERM to everybody...\n");
  SendSignal(SIGTERM);
  printf("and exiting with error code.\n");
  exit (2);
}

void IsMaster(void)
{
  int xpvm_tid;
  Pvmtmask trace_mask;
  parallel = (mytid = pvm_mytid()) >= 0;
  if (!parallel)  {
    printf("PVM not found or error in initializing machine.\n");
    master = TRUE;
    return;
  }
  master = (myparent = pvm_parent()) < 0;
#ifdef WITHXPVM
  if (master && (xpvm_tid = pvm_gettid( "xpvm", 0 )) > 0)  {
    printf("Connecting to XPVM...\n");
    /* Set Self Trace & Output Destinations & Message Codes */
    pvm_setopt( PvmSelfTraceCode, 666 );
    pvm_setopt( PvmSelfTraceTid, xpvm_tid );
    pvm_setopt( PvmSelfOutputCode, 667 );
    pvm_setopt( PvmSelfOutputTid, xpvm_tid );
    /* Set Future Children's Trace & Output Dests & Codes */
    pvm_setopt( PvmTraceTid, xpvm_tid );
    pvm_setopt( PvmTraceCode, 666 );
    pvm_setopt( PvmOutputTid, xpvm_tid );
    pvm_setopt( PvmOutputCode, 667 );
    /* Generate Default Trace Mask */
    TEV_INIT_MASK( trace_mask );
    TEV_SET_MASK( trace_mask, TEV_MCAST0 );
    TEV_SET_MASK( trace_mask, TEV_SEND0 );
    TEV_SET_MASK( trace_mask, TEV_RECV0 );
    TEV_SET_MASK( trace_mask, TEV_NRECV0 );
    /* Add Other Desired Events Here */
    /* Set Self Trace Mask */
    pvm_settmask( PvmTaskSelf, trace_mask );
    /* Set Future Children's Trace Mask */
    pvm_settmask( PvmTaskChild, trace_mask );
  }
  else
#endif
    pvm_setopt(PvmRoute, PvmRouteDirect);
  gethostname(nodename, 79);
  nodename[79] = 0;
  if (myparent < 0 && myparent != PvmNoParent)  {
    pvm_perror(NULL);
    pvm_exit();
    return;
  }
}

void InitializeWorkers(char *progname, BOOL debug)
{
  int i, j, tid[MAXWORKER], totx, totspeed, nodes, sumspeed;
  struct pvmhostinfo *info, myinfo[MAXWORKER];
  char childname[80];
  isdebug = debug;
  if (master)  {
    if (parallel)  {
      pvm_config(&workers, &i, &info);
      printf("%d %s found...", workers, workers > 1 ? "hosts" : "host");
#ifdef IRIX
      if (!workerswanted)  workerswanted = 7;
      workers = (nsubs ? workerswanted : 4); nodes = 1;
      printf("since this seems to be an SGIMP, I will use %d workers however!\n", workers);
#endif
      if (workers > MAXWORKER)  workers = MAXWORKER;
      i = nx / (OVERLAP + 3);
      if (i < workers)  workers = i;
#ifndef IRIX
      nodes = workers;
#endif
      if (workerswanted)  workers = workerswanted;
    }
    else  workers = nodes = 0;
    if (workers < 2)  {
      parallel = FALSE;
      printf("using sequential mode!\n");
    }
    else  {
      pvm_catchout(stderr);
      for (i = j = 0; i < workers; i++)
        if (workeronthishost || strcmp(info[i].hi_name, nodename))
          if (pvm_spawn(progname, NULL,
              (debug ? PvmTaskDebug : PvmTaskDefault) | PvmTaskHost,
              info[i % nodes].hi_name, 1, tid+j))  {
            myinfo[j] = info[i % nodes];
            j++;
          }
          else  {
            fprintf(stderr, "WARNING: Unable to start \"%s\" on host \"%s\"\n",
            	    progname, info[i % nodes].hi_name);
          }
      workers = j;
      if (workers <= 0)  {
        fprintf(stderr, "ERROR spawning workers!\n");
        pvm_exit();
        exit (0);
      }
      else  {
        printf("Using %d workers!\n", workers);
        /* Divide Grid */
        totx = nxm /*+ (workers-1) * 2 * OVERLAP */;
        totspeed = 0;
        for (i = workers; i--; )  totspeed += myinfo[i].hi_speed;
        j = sumspeed = 0;
        for (i = 0; i < workers; i++)  {
          worker[i].wfirstx = j;
          worker[i].firstx = j + (i > 0) * OVERLAP + 1;
          sumspeed += myinfo[i].hi_speed;
          worker[i].nx = (int)((double)sumspeed / (double)totspeed * totx + 0.5) - worker[i].firstx + 1;
          worker[i].wnx = worker[i].nx + OVERLAP * ((i > 0)+1) + 2;
          j = worker[i].firstx + worker[i].nx - OVERLAP - 1;
        }
/*        j = nxm / workers;
        for (i = workers; i--; )  {
          worker[i].firstx = j*i+1;
          worker[i].nx = j;
          worker[i].wfirstx = j*i - OVERLAP * !!i;
          worker[i].wnx = j + OVERLAP * ((i > 0)+1) + 2;
        }  */
        worker[workers-1].nx = nx - worker[workers-1].firstx;
        worker[workers-1].wnx = worker[workers-1].nx + OVERLAP + 2;
        for (i = 0; i < workers; i++)  {
          worker[i].tid = tid[i];
          pvm_recv(tid[i], NODEINIT);
          pvm_upkstr(childname);
          printf("Worker %d on %s is ready (tid =%5x). It will work on %d lines!\n",
                 i+1, childname, tid[i], worker[i].nx);
        }
        pvm_notify(PvmTaskExit, NOTIFYDATA, workers, tid);
      }
    }
  }
  else  {
    pvm_initsend(PvmDataRaw);
    pvm_pkstr(nodename);
    pvm_send(myparent, NODEINIT);
  }
}

typedef enum {SUBSTANCE, EDUCT, PRODUCT, REACTION, END}  PlaceType;

void SendChemistry(void)
{
  int i, j, type;
  PlaceType placetype;
  pvm_initsend(PvmDataRaw);
  for (i = 0; i < nsubst; i++)
    if (subst[i].inert)  {
      placetype = SUBSTANCE;
      pvm_pkint((int *)&placetype, 1, 1);
      pvm_pkstr(subst[i].name);
    }
  for (i = 0; i < nreact; i++)  {
    for (j = 0; j < reaction[i].neduct; j++)  {
      placetype = EDUCT;
      pvm_pkint((int *)&placetype, 1, 1);
      pvm_pkstr(reaction[i].educt[j].comp->name);
      pvm_pkint(&reaction[i].educt[j].prefact, 1, 1);
    }
    for (j = 0; j < reaction[i].nproduct; j++)  {
      placetype = PRODUCT;
      pvm_pkint((int *)&placetype, 1, 1);
      pvm_pkstr(reaction[i].product[j].comp->name);
      pvm_pkdouble(&reaction[i].product[j].prefact, 1, 1);
    }
    placetype = REACTION;
    pvm_pkint((int *)&placetype, 1, 1);
    pvm_pkdouble(&reaction[i].K, 1, 1);
    pvm_pkdouble(&reaction[i].EoR, 1, 1);
    pvm_pkdouble(&reaction[i].N, 1, 1);
    pvm_pkdouble(&reaction[i].M, 1, 1);
    pvm_pkdouble(&reaction[i].k0, 1, 1);
    pvm_pkdouble(&reaction[i].kinf, 1, 1);
    pvm_pkint((int *)&reaction[i].type, 1, 1);
    pvm_pkint(&reaction[i].netot, 1, 1);
    pvm_pkint(&reaction[i].specidx, 1, 1);
  }
  placetype = END;
  pvm_pkint((int *)&placetype, 1, 1);
  for (i = 0; i < workers; i++)
    pvm_send(worker[i].tid, CHEMISTRYDATA);
}  

void GetChemistry(void)
{
  PlaceType cmd;
  Token name;
  int epr;
  double ppr;
  neduct = netot = oldneduct = nproduct = oldnproduct = nreact = 
  nfast = nslow = 0;
  lineno = nsubst = 1;
  pvm_recv(myparent, CHEMISTRYDATA);
  do  {
    pvm_upkint((int *)&cmd, 1, 1);
    switch (cmd)  {
      case SUBSTANCE :
         pvm_upkstr(name);
         PlaceSubstance(name, TRUE);
         break;
      case EDUCT :
         pvm_upkstr(name);
         pvm_upkint(&epr, 1, 1);
         PlaceEduct(name, epr);
         break;
      case PRODUCT :
         pvm_upkstr(name);
         pvm_upkdouble(&ppr, 1, 1);
         PlaceProduct(name, ppr);
         break;
      case REACTION :
         pvm_upkdouble(&reaction[nreact].K, 1, 1);
         pvm_upkdouble(&reaction[nreact].EoR, 1, 1);
         pvm_upkdouble(&reaction[nreact].N, 1, 1);
         pvm_upkdouble(&reaction[nreact].M, 1, 1);
         pvm_upkdouble(&reaction[nreact].k0, 1, 1);
         pvm_upkdouble(&reaction[nreact].kinf, 1, 1);
         pvm_upkint((int *)&reaction[nreact].type, 1, 1);
         pvm_upkint(&netot, 1, 1);
         pvm_upkint(&reaction[nreact].specidx, 1, 1);
         reaction[nreact].educt = educt + oldneduct;
         reaction[nreact].product = product + oldnproduct;
         reaction[nreact].neduct = neduct - oldneduct;
         reaction[nreact].nproduct = nproduct - oldnproduct;
         reaction[nreact].netot = netot;
         oldneduct = neduct; oldnproduct = nproduct;
         netot = 0; fixedconc = 1.;
         nreact++;
         break;
    }
  }  while (cmd != END);
  AssignReactToSubst();
  CreateChemTable();
}

void SendEmissions(void)
{
  int w, i;
  EmissionDesc *e;
  int newcoord;
  Coordinate *acoord = NULL;
  if (firstemiss)
    printf("  emission data\n");
  for (e = firstemiss; e; e = e->next)  {
    pvm_initsend(PvmDataRaw);
    pvm_pkstr((char *)e->v->name);
    pvm_pkdouble(e->reduction, 1, 1);
    pvm_pkint(&e->v->id, 1, 1);
    pvm_pkint(&e->cem.n, 1, 1);
    newcoord = e->cem.c != acoord;
    pvm_pkint(&newcoord, 1, 1);
    for (w = 0; w < workers; w++)
      pvm_send(worker[w].tid, EMISSIONDATA);
    if (e->cem.n)  {  /* Point-Source entry */
      if (newcoord)  {
        pvm_initsend(PvmDataRaw);
        pvm_pkfloat(&e->cem.c[0].x, e->cem.n, &acoord[1].x - &acoord[0].x);
        pvm_pkfloat(&e->cem.c[0].y, e->cem.n, &acoord[1].y - &acoord[0].y);
        pvm_pkfloat(&e->cem.c[0].z, e->cem.n, &acoord[1].z - &acoord[0].z);
        pvm_pklong(&e->cem.c[0].loc, e->cem.n, &acoord[1].loc - &acoord[0].loc);
	for (w = 0; w < workers; w++)
	  pvm_send(worker[w].tid, EMISSIONDATA);
        acoord = e->cem.c;
      }
      pvm_initsend(PvmDataRaw);
      pvm_pkdouble(e->em, e->cem.n, 1);
      for (w = 0; w < workers; w++)
	pvm_send(worker[w].tid, EMISSIONDATA);
    }
    else
      PackVariable(e->em, LAYER_VAR, EMISSIONS);
  }
  pvm_initsend(PvmDataRaw);
  pvm_pkstr("");
  for (w = 0; w < workers; w++)
    pvm_send(worker[w].tid, EMISSIONDATA);
}

void GetEmissions(void)
{
  int w, newcoord, nc, eid;
  char name[80];
  EmissionDesc *e;
  double *data, reduction;
  Coordinate *coord;
  pvm_recv(myparent, EMISSIONDATA);
  pvm_upkstr(name);
  actualsection = McInterface::mcsect.FindSection(EMISSIONS);
  while (*name)  {
    pvm_upkdouble(&reduction, 1, 1);
    pvm_upkint(&eid, 1, 1);
    pvm_upkint(&nc, 1, 1);
    pvm_upkint(&newcoord, 1, 1);
    if (nc && newcoord)  {
      coord = NEW(nc, Coordinate);
      if (!coord)  {
        fprintf(stderr, "ERROR: Unable to allocate memory in \"GetEmissions\" (mcparallel.c)\n");
        exit (3);
      }
      pvm_recv(myparent, EMISSIONDATA);
      pvm_upkfloat(&coord[0].x, nc, &coord[1].x - &coord[0].x);
      pvm_upkfloat(&coord[0].y, nc, &coord[1].y - &coord[0].y);
      pvm_upkfloat(&coord[0].z, nc, &coord[1].z - &coord[0].z);
      pvm_upklong(&coord[0].loc, nc, &coord[1].loc - &coord[0].loc);
      ConvertEmissionCoords(nc, coord, xoff);
    }
    data = EmissionVariable(GetNamedVar(name), nc, eid, coord)->v.d;
    /*    data = PlaceEmissionVariable(GetNamedVar(name), nc, coord);  affects firstemiss */
    *firstemiss->reduction = reduction;
    if (nc)  {
      pvm_recv(myparent, EMISSIONDATA);
      pvm_upkdouble(data, nc, 1);
    }
    else
      UnpackVariable(data);
    pvm_recv(myparent, EMISSIONDATA);
    pvm_upkstr(name);
  }
}

void SendDataToWorkers(void)
{
  int i, j, k, net, l, r, w;
  BorderType bt;
  GroundParam *gr;
  /* TIME */
  printf("  core modules...\n");
  pvm_initsend(PvmDataRaw);
  pvm_pklong(&tstart, 1, 1);
  pvm_pklong(&tend, 1, 1);
  pvm_pklong(&tincmax, 1, 1);
  pvm_pklong(&irtinc, 1, 1);
  pvm_pkint(&nta, 1, 1);
  pvm_pklong(tincarray, nta, 1);
  pvm_pkdouble(&timeofday, 1, 1);
  pvm_pklong(&dayofyear, 1, 1);
  pvm_pklong(&tchem, 1, 1);
  pvm_pklong(&dumptime, 1, 1);
  pvm_pkstr(dumppath);
  pvm_pkint(&timeddump, 1, 1);
  for (i = 0; i < workers; i++)  {
    pvm_send(worker[i].tid, TIMEDATA);
  }
  /* GRID */
  for (i = 0; i < workers; i++)  {
    pvm_initsend(PvmDataRaw);
    l = (i ? worker[i-1].tid : 
         (westbordertype == CYCLIC && workers > 1 ? worker[workers-1].tid : 0));
    r = (i+1 < workers ? worker[i+1].tid : 
         (eastbordertype == CYCLIC && workers > 1 ? worker[0].tid : 0));
    leftest = !i; rightest = i+1 == workers;
    pvm_pkint(&l, 1, 1);
    pvm_pkint(&r, 1, 1);
    pvm_pkint(&leftest, 1, 1);
    pvm_pkint(&rightest, 1, 1);
    i++; pvm_pkint(&i, 1, 1); i--;
    worknum = 0;
    j = worker[i].wnx - 2;
    pvm_pkint(&j, 1, 1);
    pvm_pkint(&ny, 1, 1);
    pvm_pkint(&nz, 1, 1);
    pvm_pkint(&nsubs, 1, 1);
    pvm_pkint(&worker[i].wfirstx, 1, 1);
    pvm_pkdouble(&dx, 1, 1);
    pvm_pkdouble(&dy, 1, 1);
    pvm_send(worker[i].tid, GRIDDATA);
  }
  /* OPTIONS */
  pvm_initsend(PvmDataRaw);
  pvm_pkint(&groundinterface, 1, 1);
  pvm_pkint(&shortwaveradiation, 1, 1);
  pvm_pkint(&radiationtype, 1, 1);
  pvm_pkint(&calcsoiltemp, 1, 1);
  pvm_pkint((int *)&turbtype, 1, 1);
  pvm_pkint((int *)&advectiontype, 1, 1);
  pvm_pkint(&advection, 1, 1);
  pvm_pkint(&windadvection, 1, 1);
  pvm_pkint(&cloudwater, 1, 1);
  pvm_pkint(&cloudice, 1, 1);
  pvm_pkint(&isdebug, 1, 1);
  pvm_pkint((int *)&coriolistype, 1, 1);
  pvm_pkint(&smiord, 1, 1);
  pvm_pkint(&smnonos, 1, 1);
  pvm_pkint((int *)&pressuretype, 1, 1);
  pvm_pkint((int *)&filtertype, 1, 1);
  pvm_pkint((int *)&radiationtype, 1, 1);
  pvm_pkdouble(&spatialfilter, 1, 1);
  for (i = 0; i < workers; i++)  {
    pvm_send(worker[i].tid, OPTIONDATA);
  }
  /* ENVIRONMENT */
  printf("  environment\n");
  PackVariable(topo, LAYER_VAR, ENVIRONMENT);
  pvm_initsend(PvmDataRaw);
  pvm_pkdouble(zfaces, nz, 1);
/*  pvm_pkdouble(leveli, nz, 1);
  pvm_pkdouble(ldiff, nz, 1); */
  pvm_pkdouble(&reflevel, 1, 1);
  pvm_pkdouble(&Xlong, 1, 1);
  pvm_pkdouble(&Xlat, 1, 1);
  pvm_pkdouble(&timezonediff, 1, 1);
  pvm_pkdouble(&turnmapangle, 1, 1);
  pvm_pkdouble(&stratozone, 1, 1);
  pvm_pkdouble(&topprecipwater, 1, 1);
  pvm_pkdouble(&turbidity, 1, 1);
  pvm_pkdouble(&toppp, 1, 1);
  pvm_pkdouble(&photofact, 1, 1);
  for (i = 0; i < workers; i++)  {
    pvm_send(worker[i].tid, ENVDATA);
  }
  /* INITIAL-DATA */
  printf("  Initial-data\n");
  net = 0;
  for (j = maxentity; j--; )
    net += !!g[j];
  for (i = 0; i < workers; i++)  {
    pvm_initsend(PvmDataRaw);
    pvm_pkint(&net, 1, 1);
    pvm_pkint(gactive, maxentity, 1);
    pvm_send(worker[i].tid, INITDATA);
    for (j = UWIND; j < maxentity; j++)  {
      if (g[j])  {
        pvm_initsend(PvmDataRaw);
        pvm_pkint(&j, 1, 1);
        for (k = 0; k < nz; k++)
          pvm_pkdouble(g[j]+k*layer + worker[i].wfirstx*row,
                       row*worker[i].wnx, 1);
        pvm_send(worker[i].tid, INITDATA);
      }
    }
  }
  for (j = maxentity; j--; )  {
    if (g[j] && j != TEMP)
      free(g[j]);
    if (gactive[j])
      g[j] = g[TEMP];
    else
      g[j] = NULL;
  }
  if (!cloudwater)  g[CLOUDWATER] = g[RAINWATER] = NULL;
  if (!cloudice)    g[CLOUDICE] = NULL;
  meshcache = g[TEMP];
  /* BORDER */
  printf("  border-data\n");
  for (i = 0; i < workers; i++)  {
    pvm_initsend(PvmDataRaw);
    bt = (i || (westbordertype == CYCLIC && workers > 1) ? NEIGHBOUR : westbordertype);
    pvm_pkint((int *)&bt, 1, 1);
    bt = (i+1 < workers || (eastbordertype == CYCLIC && workers > 1) ? NEIGHBOUR : eastbordertype);
    pvm_pkint((int *)&bt, 1, 1);
    pvm_pkint((int *)&southbordertype, 1, 1);
    pvm_pkint((int *)&northbordertype, 1, 1);
    pvm_pkint(northborderval, maxentity, 1);
    pvm_pkint(southborderval, maxentity, 1);
    pvm_pkint(westborderval, maxentity, 1);
    pvm_pkint(eastborderval, maxentity, 1);
    pvm_send(worker[i].tid, BORDERDATA);
  }
  PackVariable(toppress, LAYER_VAR, TOP_BORDER);
  PackVariable(toptemp, LAYER_VAR, TOP_BORDER);
  PackVariable(ugeos, PROFILE_VAR, TOP_BORDER);
  PackVariable(vgeos, PROFILE_VAR, TOP_BORDER);
  for (i = 0; i < maxentity; i++)  {
    PackVariable(westborder+i*row*nz, WALL_VAR, WEST_BORDER);
    PackVariable(eastborder+i*row*nz, WALL_VAR, EAST_BORDER);
    PackVariable(northborder+i*xrow*nz, XWALL_VAR, NORTH_BORDER);
    PackVariable(southborder+i*xrow*nz, XWALL_VAR, SOUTH_BORDER);
  }
  printf("  ground data\n");
  SendGroundToWorkers();
  pvm_initsend(PvmDataRaw);
#ifdef RADFACT
  pvm_pkdouble(Rdirfact, mesh, 1);
  pvm_pkdouble(Rdifffact, mesh, 1);
#else
  pvm_pkdouble(Rdirfact, nz, 1);
  pvm_pkdouble(Rdifffact, nz, 1);
#endif
  for (w = 0; w < workers; w++)  {
    pvm_send(worker[w].tid, GROUNDDATA);
  }
  if (nsubs)  {
    printf("  chemistry\n");
    SendChemistry();
  }
  SendEmissions();
  pvm_initsend(PvmDataRaw);
  SendDepositionData();
  for (i = 0; i < workers; i++)
    pvm_send(worker[i].tid, DEPOSITIONDATA);
  for (i = 0; i < workers; i++)  {
    pvm_initsend(PvmDataRaw);
    SendBorderTime(!i, i+1 == workers);
    pvm_send(worker[i].tid, BORDERDATA);
  }
  McInterface::modmanager.SendToWorkers();
  printf("Done.\n");
}

void ReadDataFromMaster(void)
{
  Entity et;
  int i, j;
  GroundParam *gr;
  /* TIME */
  pvm_recv(myparent, TIMEDATA);
  pvm_upklong(&tstart, 1, 1);
  pvm_upklong(&tend, 1, 1);
  pvm_upklong(&tincmax, 1, 1);
  pvm_upklong(&irtinc, 1, 1);
  pvm_upkint(&nta, 1, 1);
  pvm_upklong(tincarray, nta, 1);
  pvm_upkdouble(&timeofday, 1, 1);
  pvm_upklong(&dayofyear, 1, 1);
  pvm_upklong(&tchem, 1, 1);
  pvm_upklong(&dumptime, 1, 1);
  pvm_upkstr(dumppath);
  pvm_upkint(&timeddump, 1, 1);
  /* GRID */
  pvm_recv(myparent, GRIDDATA);
  pvm_upkint(&leftbrother, 1, 1);
  pvm_upkint(&rightbrother, 1, 1);
  pvm_upkint(&leftest, 1, 1);
  pvm_upkint(&rightest, 1, 1);
  pvm_upkint(&worknum, 1, 1);
  pvm_upkint(&nx, 1, 1);
  pvm_upkint(&ny, 1, 1);
  pvm_upkint(&nz, 1, 1);
  pvm_upkint(&nsubs, 1, 1);
  pvm_upkint(&xoff, 1, 1);
  pvm_upkdouble(&dx, 1, 1);
  pvm_upkdouble(&dy, 1, 1);
  mfirstx = (1 + OVERLAP) * !leftest;
  mnx = nx + 2 - (1 + OVERLAP) * (!leftest + !rightest);
  mlastx = mfirstx + mnx;
  ny--;
  EndOfGrid();
  /* OPTIONS */
  pvm_recv(myparent, OPTIONDATA);
  pvm_upkint(&groundinterface, 1, 1);
  pvm_upkint(&shortwaveradiation, 1, 1);
  pvm_upkint(&radiationtype, 1, 1);
  pvm_upkint(&calcsoiltemp, 1, 1);
  pvm_upkint((int *)&turbtype, 1, 1);
  pvm_upkint((int *)&advectiontype, 1, 1);
  pvm_upkint(&advection, 1, 1);
  pvm_upkint(&windadvection, 1, 1);
  pvm_upkint(&cloudwater, 1, 1);
  pvm_upkint(&cloudice, 1, 1);
  pvm_upkint(&isdebug, 1, 1);
  pvm_upkint((int *)&coriolistype, 1, 1);
  pvm_upkint(&smiord, 1, 1);
  pvm_upkint(&smnonos, 1, 1);
  pvm_upkint((int *)&pressuretype, 1, 1);
  pvm_upkint((int *)&filtertype, 1, 1);
  pvm_upkint((int *)&radiationtype, 1, 1);
  pvm_upkdouble(&spatialfilter, 1, 1);
  EndOfOptions();
  /* GROUND */
  UnpackVariable(topo);
  pvm_recv(myparent, ENVDATA);
  pvm_upkdouble(zfaces, nz, 1);
/*  pvm_upkdouble(leveli, nz, 1);
  pvm_upkdouble(ldiff, nz, 1);  */
  pvm_upkdouble(&reflevel, 1, 1);
  pvm_upkdouble(&Xlong, 1, 1);
  pvm_upkdouble(&Xlat, 1, 1);
  pvm_upkdouble(&timezonediff, 1, 1);
  pvm_upkdouble(&turnmapangle, 1, 1);
  pvm_upkdouble(&stratozone, 1, 1);
  pvm_upkdouble(&topprecipwater, 1, 1);
  pvm_upkdouble(&turbidity, 1, 1);
  pvm_upkdouble(&toppp, 1, 1);
  pvm_upkdouble(&photofact, 1, 1);
  EndOfEnvironment();
  /* INITIAL-DATA */
  pvm_recv(myparent, INITDATA);
  pvm_upkint(&i, 1, 1);
  pvm_upkint(gactive, maxentity, 1);
  while (i--)  {
    pvm_recv(myparent, INITDATA);
    pvm_upkint((int *)&et, 1, 1);
    pvm_upkdouble(g[et], mesh, 1);
  }
  pvm_recv(myparent, BORDERDATA);
  pvm_upkint((int *)&westbordertype, 1, 1);
  pvm_upkint((int *)&eastbordertype, 1, 1);
  pvm_upkint((int *)&southbordertype, 1, 1);
  pvm_upkint((int *)&northbordertype, 1, 1);
  pvm_upkint(northborderval, maxentity, 1);
  pvm_upkint(southborderval, maxentity, 1);
  pvm_upkint(westborderval, maxentity, 1);
  pvm_upkint(eastborderval, maxentity, 1);
  UnpackVariable(toppress);
  UnpackVariable(toptemp);
  UnpackVariable(ugeos);
  UnpackVariable(vgeos);
  for (i = 0; i < maxentity; i++)  {
    if (leftest)   UnpackVariable(westborder+i*row*nz);
    if (rightest)  UnpackVariable(eastborder+i*row*nz);
    UnpackVariable(northborder+i*xrow*nz);
    UnpackVariable(southborder+i*xrow*nz);
  }
  GetGroundFromMaster();
  pvm_recv(myparent, GROUNDDATA);
  pvm_upkdouble(Rdirfact, nz, 1);
  pvm_upkdouble(Rdifffact, nz, 1);
  if (nsubs)  GetChemistry();
  GetEmissions();
  pvm_recv(myparent, DEPOSITIONDATA);
  GetDepositionData();
  pvm_recv(myparent, BORDERDATA);
  RecvBorderTime();
  McInterface::modmanager.GetFromMaster();
}

void SendCommand(void)
{
  int i;
  for (i = 0; i < workers; i++)  {
    pvm_send(worker[i].tid, COMMAND);
  }
}

void SendStatus(int i)
{
  double energy;
  pvm_initsend(PvmDataRaw);
  pvm_pkint(&i, 1, 1);
  pvm_pklong(&actime, 1, 1);
  energy = CalcTopEnergy();
  pvm_pkdouble(&energy, 1, 1);
  pvm_pklong(&waitingtime, 1, 1);
  if (chemtinc)
    waitingtime = 0;
  pvm_pkint(&worknum, 1, 1);
  pvm_pkstr(nodename);
  pvm_send(myparent, STATUS);
}

int IsPlausible(double *energy)
{
  int i, s, wn;
  static long mark = 0;
  long wt;
  long act;
  char hostname[30];
  struct tms buffer;
  double e;
  plausible = TRUE;
  *energy = 0.;
/*  printf("\n"); */
  for (i = 0; i < workers; i++)  {
    ReceiveWithTimeout(-1, STATUS);
    pvm_upkint(&s, 1, 1);
    plausible &= s;
    pvm_upklong(&act, 1, 1);
    pvm_upkdouble(&e, 1, 1);
    pvm_upklong(&wt, 1, 1);
    pvm_upkint(&wn, 1, 1);
    if (chemtinc)
      worker[wn-1].waitingtime = times(&buffer) - mark - wt;
    *energy += e;
    pvm_upkstr(hostname);
    if (act != actime)  {
      printf("ERROR: Worker \"%s\" is at %ld seconds.\n", hostname, act);
      plausible = FALSE;
    }
  }
/*
  if (chemtinc)  {
    for (i = 0; i < workers; i++)
      printf("Worker #%d waited %ld seconds\n", i+1, worker[i].waitingtime);
    mark = times(&buffer);
  }
*/
  return (plausible);
}

void ChangeTai(int *tai)
{
  int i, ata;
  if (master)  {
    *tai = 0;
    for (i = 0; i < workers; i++)  {
      ReceiveWithTimeout(-1, STATUS);
      pvm_upkint(&ata, 1, 1);
      if (ata > *tai)  *tai = ata;
    }
    pvm_initsend(PvmDataRaw);
    pvm_pkint(tai, 1, 1);
    for (i = 0; i < workers; i++)  pvm_send(worker[i].tid, STATUS);
  }
  else  {
    pvm_initsend(PvmDataRaw);
    pvm_pkint(tai, 1, 1);
    pvm_send(myparent, STATUS);
    pvm_recv(myparent, STATUS);
    pvm_upkint(tai, 1, 1);
  }
}

void ChangeVectMult(double *sum)
{
  double wsum;
  int i, nit;
  if (master)  {
    for (nit = 0; 1; nit++)  {
      *sum = 0.;
      for (i = 0; i < workers; i++)  {
        ReceiveWithTimeout(-1, STATUS);
        pvm_upkdouble(&wsum, 1, 1);
        *sum += wsum;
      }
      pvm_initsend(PvmDataRaw);
      pvm_pkdouble(sum, 1, 1);
      for (i = 0; i < workers; i++)
        pvm_send(worker[i].tid, STATUS);
      if (wsum == -9999.)  {
        printf("cj-gr it: %d  ", nit / 2);
        return;
      }
    }
  }
  else  {
    pvm_initsend(PvmDataRaw);
    pvm_pkdouble(sum, 1, 1);
    pvm_send(myparent, STATUS);
    pvm_recv(myparent, STATUS);
    pvm_upkdouble(sum, 1, 1);
  }
}

ParCommand GetCommand(void)
{
  ParCommand c;
  pvm_recv(myparent, COMMAND);
  pvm_upkint((int *)&c, 1, 1);
  return (c);
}

void GetAveragesFromWorkers(double *sum)
{
  int i, k;
  double bar[MAXZGRID];
  memset(sum, 0, nz*sizeof(double));
  for (i = 0; i < workers; i++)  {
    ReceiveWithTimeout(-1, AVERAGES);
    pvm_upkdouble(bar, nz, 1);
    for (k = nz; k--; )  sum[k] += bar[k];
  }
}

void SendAveragesToWorkers(double *bar)
{
  int i;
  pvm_initsend(PvmDataRaw);
  pvm_pkdouble(bar, nz, 1);
  for (i = 0; i < workers; i++)
    pvm_send(worker[i].tid, AVERAGES);
}

void ChangeAveragesWithMaster(double *bar)
{
  pvm_initsend(PvmDataRaw);
  pvm_pkdouble(bar, nz, 1);
  pvm_send(myparent, AVERAGES);
  pvm_recv(myparent, AVERAGES);
  pvm_upkdouble(bar, nz, 1);
}

void ReadMeshFromWorkers(double *m, int z0, int dims)
{
  int k, w, loc, n;
  for (w = 0; w < workers; w++)  {
    ReceiveWithTimeout(worker[w].tid, VALUES);
    loc = (worker[w].firstx - (!w)) * row;
    n = (worker[w].nx + (!w) + (w+1 == workers)) * row;
    for (k = z0; k < (dims & Z_DIM ? nz : 1); k++)
      pvm_upkdouble(m+loc+k*layer, n, 1);
  }
}

void SendMeshToWorkers(double *m, int z0)
{
  int k, w, loc, n;
  for (w = 0; w < workers; w++)  {
    pvm_initsend(McInPlace);
    loc = worker[w].wfirstx * row;
    n = worker[w].wnx * row;
    for (k = z0; k < nz; k++)
      pvm_pkdouble(m+loc+k*layer, n, 1);
    pvm_send(worker[w].tid, VALUES);
  }
}

double *GetGridVarFromWorker(int et)
{
  int cmd, w;
  pvm_initsend(PvmDataRaw);
  cmd = SENDGRID;
  pvm_pkint(&cmd, 1, 1);
  pvm_pkint(&et, 1, 1);
  for (w = 0; w < workers; w++)
    pvm_send(worker[w].tid, COMMAND);
  ReadMeshFromWorkers(meshcache, 0, ALL_DIM);
  return (meshcache);
}

double GetValFromWorker(VarDesc *v, int i, int j, int k, BOOL cache)
{
  static VarDesc *cachedvar = NULL;
  static long cachetime = -1;
  Entity et;
  int cmd, w, pos;
  double result;
  if (cache)  {
    if (cachedvar != v || cachetime != actime)  {
      if (v->storetype == GRID_VAL)  {
        GetGridVarFromWorker(v->v.et);
      }
      else  {
        pvm_initsend(PvmDataRaw);
        cmd = SENDMESH;
        pvm_pkint(&cmd, 1, 1);
        pvm_pkstr((char *)v->name);
	for (w = 0; w < workers; w++)
          pvm_send(worker[w].tid, COMMAND);
        ReadMeshFromWorkers(meshcache, 0, v->dims);
      }
      cachedvar = v; cachetime = actime;
    }
    return (meshcache[k*layer+i*row+j]);
  }
  else  {
    pvm_initsend(PvmDataRaw);
    if (v->storetype == GRID_VAL)  {
      cmd = SENDGRIDVAL;
      pvm_pkint(&cmd, 1, 1);
      et = (int)v->v.et;
      pvm_pkint(&et, 1, 1);
    }
    else  {
      cmd = SENDMESHVAL;
      pvm_pkint(&cmd, 1, 1);
      pvm_pkstr((char *)v->name);
    }
    for (w = 0; w+1 < workers && i >= worker[w+1].firstx; w++);
    pos = i - worker[w].wfirstx;
    pvm_pkint(&pos, 1, 1);
    pvm_pkint(&j, 1, 1);
    k *= !!(v->dims & Z_DIM);
    pvm_pkint(&k, 1, 1);
    pvm_send(worker[w].tid, COMMAND);
    ReceiveWithTimeout(worker[w].tid, VALUES);
    pvm_upkdouble(&result, 1, 1);
    return (result);
  }
}

void ReadGroundFromWorker(void)
{
  int w, i, j, n;
  GroundParam *gr;
  for (w = 0; w < workers; w++) {
    ReceiveWithTimeout(worker[w].tid, VALUES);
    for (i = worker[w].firstx - !w, n = worker[w].nx + !w + (w+1 == workers); n--; i++)
      for (j = 1; j < ny; j++)  {
	gr = ground+i*row+j;
	pvm_upkdouble(&gr->Tg[0], (&gr->slope.z - &gr->Tg[0]) + 1, 1);
	pvm_upkshort(&gr->firstabove, 1, 1);
      }
  }
}

void GetGroundFromWorker(void)
{
  static long cachedtime = -1;
  const ParCommand cmd = SENDGROUND;
  int w;
  GroundParam *gr;
  if (cachedtime != actime)  {
    for (w = 0; w < workers; w++)  {
      pvm_initsend(PvmDataRaw);
      pvm_pkint((int *)&cmd, 1, 1);
      pvm_send(worker[w].tid, COMMAND);
    }
    ReadGroundFromWorker();
    cachedtime = actime;
  }
}

void SendGroundToMaster(void)
{
  int i, j;
  GroundParam *gr;
  pvm_initsend(PvmDataRaw);
  for (i = mfirstx; i < mlastx; i++)
    for (j = 1; j < ny; j++)  {
      gr = ground+i*row+j;
      pvm_pkdouble(&gr->Tg[0], (&gr->slope.z - &gr->Tg[0]) + 1, 1);
      pvm_pkshort(&gr->firstabove, 1, 1); 
    }
  pvm_send(myparent, VALUES);
}

void SendGroundToWorkers(void)
{
  GroundParam *gr;
  int w, i, j, k;
  for (w = 0; w < workers; w++)  {
    pvm_initsend(PvmDataRaw);
    for (i = worker[w].wfirstx, k = worker[w].wnx; k--; i++)
      for (j = 0; j <= ny; j++)  {
        gr = ground+i*row+j;
        pvm_pkdouble(&gr->Tg[0], (&gr->slope.z - &gr->Tg[0]) + 1, 1);
        pvm_pkshort(&gr->firstabove, 1, 1);
      }
    pvm_send(worker[w].tid, GROUNDDATA);
  }
}

void GetGroundFromMaster(void)
{
  GroundParam *gr;
  int i, j;
  pvm_recv(myparent, GROUNDDATA);
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)  {
      gr = ground+i*row+j;
      pvm_upkdouble(&gr->Tg[0], (&gr->slope.z - &gr->Tg[0]) + 1, 1);
      pvm_upkshort(&gr->firstabove, 1, 1);
    }
}

void SendMeshToMaster(double *m, int z0, int dims)
{
  int k, loc, n;
  loc = mfirstx * row;
  n = mnx * row;
  pvm_initsend(McInPlace);
  for (k = z0; k < (dims & Z_DIM ? nz : 1); k++)
    pvm_pkdouble(m+loc+k*layer, n, 1);
  pvm_send(myparent, VALUES);
}

void GetMeshFromMaster(double *m, int z0)
{
  int k;
  pvm_recv(myparent, VALUES);
  pvm_upkdouble(m+z0*layer, mesh-z0*layer, 1);
}

void SendWallToWorker(int net, int x, int tid, CommunicationClass c)
{
  Entity et;
  int k;
  x *= row;
  for (et = 0; et < net; et++)  {
    pvm_initsend(McInPlace);
    for (k = 0; k < nz; k++)
      if (g[et])  pvm_pkdouble(g[et]+k*layer+x, (OVERLAP+1)*row, 1);
    if (groundinterface)
      pvm_pkbyte((char *)(ground+x), (OVERLAP+1)*row*sizeof(GroundParam), 1);
    pvm_send(tid, c);
  }
}

void SendPressureWall(double *x, int xl, int tid, CommunicationClass c)
{
  int k;
  xl *= row;
  pvm_initsend(McInPlace);
  for (k = 0; k < nz; k++)
    pvm_pkdouble(x+k*layer+xl, row, 1);
  pvm_send(tid, c);
}

void GetWallFromWorker(int net, int x, int tid, CommunicationClass c)
{
  Entity et;
  int k;
  x *= row;
  for (et = 0; et < net; et++)  {
    ReceiveWithTimeout(tid, c);
    for (k = 0; k < nz; k++)
      if (g[et])  pvm_upkdouble(g[et]+k*layer+x, (OVERLAP+1)*row, 1);
    if (groundinterface)
      pvm_upkbyte((char *)(ground+x), (OVERLAP+1)*row*sizeof(GroundParam), 1);
  }
}

void GetPressureWall(double *x, int xl, int tid, CommunicationClass c)
{
  int k;
  xl *= row;
  pvm_recv(tid, c);
  for (k = 0; k < nz; k++)
    pvm_upkdouble(x+k*layer+xl, row, 1);
}

/*
void GetLayerFromWorkers(double *l)
{
  int w;
  for (w = 0; w < workers; w++)  {
    pvm_recv(worker[w].tid, FILTERDATA);
    pvm_upkdouble(l + (worker[w].firstx - !w) * row,
    		  (worker[w].nx + !w + (w+1 == workers)) * row, 1);
  }
}

void SendLayerToWorkers(double *l)
{
  int w;
  for (w = 0; w < workers; w++)  {
    pvm_initsend(McInPlace);
    pvm_pkdouble(l + (worker[w].firstx - !w) * row,
    		 (worker[w].nx + !w + (w+1 == workers)) * row, 1);
    pvm_send(worker[w].tid, FILTERDATA);
  }
}

void ChangeLayerWithMaster(double *l)
{
  int f, n;
  f = !leftest * row;
  n = (nxm + leftest + rightest) * row;
  pvm_initsend(McInPlace);
  pvm_pkdouble(l + f, n, 1);
  pvm_send(myparent, FILTERDATA);
  pvm_recv(myparent, FILTERDATA);
  pvm_upkdouble(l + f, n, 1);
}

*/

void PackVariable(double *d, BorderVarType t, InputSection s)
{
  int k, n, loc, mem, w;
  switch (t)  {
  case WALL_VAR :
    pvm_initsend(McInPlace);
    mem = row*nz;
    pvm_pkint(&mem, 1, 1);
    pvm_pkdouble(d, row*nz, 1);
    if (s == WEST_BORDER)  pvm_send(worker[0].tid, BORDERVALUES);
    else if (s == EAST_BORDER)  pvm_send(worker[workers-1].tid, BORDERVALUES);
    else {
      fprintf(stderr, "ERROR: Internal implementation error in \"PackVariable\"\n");
      plausible = FALSE;
    }
    break;
  case XWALL_VAR :
    for (w = 0; w < workers; w++)  {
      pvm_initsend(PvmDataRaw);
      n = worker[w].wnx;
      loc = worker[w].wfirstx;
      mem = nz*n;
      pvm_pkint(&mem, 1, 1);
      for (k = 0; k < nz; k++)
	pvm_pkdouble(d+loc+k*xrow, n, 1);
      pvm_send(worker[w].tid, BORDERVALUES);
    }
    break;
  case LAYER_VAR :
    for (w = 0; w < workers; w++)  {
      pvm_initsend(PvmDataRaw);
      mem = worker[w].wnx * row;
      pvm_pkint(&mem, 1, 1);
      pvm_pkdouble(d+worker[w].wfirstx * row, mem, 1);
      pvm_send(worker[w].tid, BORDERVALUES);
    }
    break;
  case PROFILE_VAR :
    pvm_initsend(McInPlace);
    pvm_pkint(&nz, 1, 1);
    pvm_pkdouble(d, nz, 1);
    for (w = 0; w < workers; w++)  {
      pvm_send(worker[w].tid, BORDERVALUES);
    }
    break;
  case MESH_VAR :
    for (w = 0; w < workers; w++)  {
      pvm_initsend(McInPlace);
      n = worker[w].wnx * row;
      loc = worker[w].wfirstx * row;
      mem = nz*n;
      pvm_pkint(&mem, 1, 1);
      for (k = 0; k < nz; k++)
	pvm_pkdouble(d+loc+k*layer, n, 1);
      pvm_send(worker[w].tid, BORDERVALUES);
    }
    break;
  }
}

void PackCoordVariable(VarDesc *v)
{
  int w;
  pvm_initsend(McInPlace);
  pvm_pkint(&v->ncoord, 1, 1);
  pvm_pkdouble(v->v.d, v->ncoord, 1);
  for (w = 0; w < workers; w++)  {
    pvm_send(worker[w].tid, BORDERVALUES);
  }
}

void SendRestartName(char *name)
{
  int w;
  pvm_initsend(PvmDataRaw);
  w = !!name;
  pvm_pkint(&w, 1, 1);
  if (w)  pvm_pkstr(name);
  for (w = 0; w < workers; w++)  {
    pvm_send(worker[w].tid, RESTARTDATA);
  }
}

BOOL GetRestartName(char *name)
{
  int ok;
  pvm_recv(myparent, RESTARTDATA);
  pvm_upkint(&ok, 1, 1);
  if (ok)  pvm_upkstr(name);
  return (ok);
}

void UnpackVariable(double *d)
{
  int n;
  pvm_recv(myparent, BORDERVALUES);
  pvm_upkint(&n, 1, 1);
  pvm_upkdouble(d, n, 1);
}

void SendToAll(int id)
{
  int w;
  for (w = 0; w < workers; w++)  {
    pvm_send(worker[w].tid, id);
  }
}

void SendSignal(int sig)
{
  int w;
  for (w = 0; w < workers; w++)
    pvm_sendsig(worker[w].tid, sig);
}
