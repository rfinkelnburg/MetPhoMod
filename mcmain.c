
/**************************************************************
   MODULE MCMain 
   Das Hauptmodul des MetPhoMod-Programmes.
***************************************************************/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mcglobal.h"
#include "mcpress.h"
#include "mcdynamics.h"
#include "mcprint.h"
#include "mcclouds.h"
#include "mcttt.h"
#include "mckeps.h"
#include "mcfloat.h"
#include "mcfilter.h"
#include "mchelp.h"
#include "mcemiss.h"
#include "mcadvect.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcparse.h"
#include "mccdfin.h"
#include "mcdep.h"
#include "mchemparse.h"
#include "mchem.h"
#include "mcnest.h"
#include "mcppm.h"
#ifdef PARALLEL
#include <pvm3.h>
#include "mcparallel.h"
#endif
#include "mc_module.hh"

typedef enum {RESTART = 1, NOEXCPT = 2, SYNTAX_ONLY = 4, REOPEN = 8,
              NO_OUTPUT = 16, DEBUG_WORKERS = 32, SEQUENTIAL = 64}  RunningFlags;

void TraceHumid(char *txt)
{
  int i, et;
  for (et = UWIND; et <= HUMIDITY; et++)
    for (i = mesh; i--; )
      if (fabs(g[3][i]) > 1000.)  goto humid_error;
  return;
humid_error :
  printf("et = %d error at point %d (%d/%d/%d) %s\n",
     et, i, (i % layer) / row, (i % layer) % row, i / layer, txt);
  exit (1);
}

long int AnalizeOptions(int *argc, char **argv, char **arg)
{
  int i, j, k, result = 0;
  char *p;
  static char *optiontxt[] = {"-RESTART",  "-noexcpt", "-syntax", "-REOPEN",
  			      "-nooutput", "-debug", "-sequential"};
  for (i = 1; i < *argc; )
    if (*argv[i] == '-')  {
      for (j = 6; j >= 0 && strcmp(argv[i], optiontxt[j]); j--);
      if (j >= 0)  {
        result |= 1 << j;
        (*argc)--;
        for (k = i; k < *argc; k++)  argv[k] = argv[k+1];
        if (j == 0)  {
          if (i < *argc)  {
            *arg = argv[i];
            (*argc)--;
            for (k = i; k < *argc; k++)  argv[k] = argv[k+1];
          }
          else  {
            printf("Error in command line: -RESTART needs an argument\n");
            *argc = 0;
          }
        }
      }
      else  *argc = 0;
    }
    else  i++;
  return (result);
}

/*
void PrintEnergy(char *cmmt)
{
  double sum, tsum;
  int i;
  static double prev = 0.;
  sum = tsum = 0.;
  for (i = nz; i-- > 1; )  {
    sum += level[i] * Cp * g[TEMP][i][2][2]*density[i];      
    tsum += level[i] * Cp * g[TEMP][i][2][2];
  }
  sum += Cp * ground[2][2].a[TEMP] * ground[2][2].z * density[1];
  tsum += Cp * ground[2][2].a[TEMP] * ground[2][2].z;
  printf("%s : %lf %lf  %lf\n", cmmt, sum, sum - prev, tsum);
  prev = sum;
}
*/

void TestChemicals(char *location)
{
  int i, j, k;
  Entity et;
  for (et = SUBS; et < maxentity; et++)
    for (i = mesh; --i; )
      if (!pstat && g[et][i] < 0.)  {
        printf("\"%s\" is %le at %i/%i/%i, %s\n",
               subst[et - HUMIDITY].name, g[et][i],
               (i / row) % xrow, i % row, i / layer, location);
        plausible = FALSE;
      }
}

/*
FILE *PrintHeader(char *fname)
{
  FILE *f;
  int i;
  char *cat[] = {"Begin", "Advection", "TimeStep", "Ground-Turb", "Trans-Turb",
  		 "Filter", "Emit", "Deposit", "Chemistry"};
  f = fopen(fname, "w");
  fprintf(f, "");
  for (i = 0; i < 9; i++)  fprintf(f, "%-15s", cat[i]);
  fprintf(f, "\n");
  while (i--)  fprintf(f, "             O3            NO2             NO");
  fprintf(f, "\n");
  return (f);
}

void PrintConz(FILE *f, int i, int j, int k)
{
  int loc;
  loc = k*layer+j*row+i;
  fprintf(f, "%15.8e%15.8e%15.8e", g[5][loc], g[4][loc], g[6][loc]);
}

void PrintEnergy()
{
  double sum = 0.;
  int i;
  for (i = nzm*layer; i < mesh; i++)
    sum += sqr(flux[WWIND][i]);
  printf("Energy = %lf\n", sum);
}
*/

void SequentialLoop(int startupflags)
{
  int i;
  double theta;
  long daytime, newtinc;
  Entity et;
  Vector sunpos;
  while (plausible && actime < tend)  {
    daytime = ((int)(timeofday*3600.)+actime) % 86400;
    printf("Now at %li seconds = %02li:%02li:%02li  ",
           actime, daytime / 3600,
           (daytime % 3600)/60,
           daytime % 60);
    ActualiseValuesInTime(actime);	/* Actualise Border values */
    for (et = maxentity; et--; )
      CalculateLayerAverage(g[et], pstat, avg+et*nz);
    SetAvgMountains();
    CalcMeanHydrostaticPressure();
    InterpolateToFaces();			/* Interpolate wind speeds from center to faces */
    InterpolateToCenter(-1., pressuretype == NONHYDROSTATIC);
    if (tstart == actime)  Continuity();
    CheckTimeStep(&tinc, &chemtinc, 0.8);
    if (pressuretype == NONHYDROSTATIC)
      ApplyBuoyancy(tinc);			/* Calculate Buoyancy */
    if (pressuretype != NOPRESSURE)  {		/* Calc wind acceleration */
      switch (pressuretype)  {			/* Calculate Pressure field */
         case HYDROSTATIC    : CalcHydrostaticPressure(); break;
         case NONHYDROSTATIC : SolveForPressure(tinc); break;
      }
      ApplyPressure(tinc);
    }
    Continuity();				/* Mass conservation */
    InterpolateToCenter(1., pressuretype == NONHYDROSTATIC);
    if (coriolistype != NOCORIOLIS)  ApplyCoriolis(tinc);
    if (filtertype != NO_FILTER)  {
      SetBoundary(WWIND);
      switch (filtertype)  {
        case PEPPER_FILTER :
           for (i = nz; i--; )
             for (et = HUMIDITY; et--; )
               ApplyFilter(g[et]+i*layer, pstat+i*layer, spatialfilter);
           break;
        case SHAPIRO_FILTER :
           for (i = nz; i--; )  {
             for (et = HUMIDITY; et--; )
               ShapiroFilter(g[et]+i*layer, pstat+i*layer, 3, spatialfilter);
/*             if (turbtype == KEPS_T)  {
               ShapiroFilter(g[TKE]+i*layer, pstat+i*layer, 3, spatialfilter);
               ShapiroFilter(g[EPS]+i*layer, pstat+i*layer, 3, spatialfilter);
             } */
           }
           break;
        case DIFFUSION_FILTER :
           for (i = nz; i--; )
             for (et = HUMIDITY; et--; )
               ShapiroFilter(g[et]+i*layer, pstat+i*layer, 1, spatialfilter);
           break;
      }
    }
    SetBoundary(WWIND+1);
    CheckTimeStep(&newtinc, &chemtinc, 1.);
    if (newtinc < tinc)  {
      fprintf(stderr, "Warning: Time step had to be changed in the middle of iteration\n");
      tinc = newtinc;
    }
    chemtinc = (chemtime <= actime ? chemtinc : 0);
    if (nsubs)  printf("chemtinc = %3ld  ", chemtinc);
    if (advection || windadvection)
      switch (advectiontype)  {
        case MPDATA_A : Advect(tinc, chemtinc); break;
        case PPM_A    : PPM_Transport(tinc, chemtinc, smiord); break;
      }
    if (dampinglayer)  DampTop();
    if (groundinterface || shortwaveradiation)  {
      GroundInterface(timeofday + (double)actime / 3600., dayofyear,
		      tinc, &sunpos);
      sunelevation = (sunpos.z > 0. ? 57.29578 * asin(sunpos.z) : 0.);
    }
    if (cloudwater)  CloudPhysics();
    if (turbtype == TTTT)  {
      if (groundinterface)
        CalcGroundTurbulence(tinc);
      CalcTransilientTurbulence(tinc);
    }
    else if (turbtype == KEPS_T)  {
      CalcKEpsilon(tinc);
      if (groundinterface)
        CalcGroundTurbulence(tinc);
      CalcKTurb(tinc, chemtinc);
    }
    Emit(tinc);
    Deposit(tinc);
    if (!(radiationtype && shortwaveradiation) &&
	// When using the two-stream module, ChemicalTimeStep
	// must be called from within ShortWaveRadiation!
	nsubs && actime % tchem == 0)  ChemicalTimeStep(tchem, &sunpos);
    SetBoundary(chemtinc ? maxentity : SUBS);
    
    // Call ModuleManager
    McInterface::modmanager.DoCalc(actime+tinc);
/*    TestChemicals("Boundary");  */
    actime += tinc; chemtime += chemtinc;
    if (!(startupflags & NO_OUTPUT))  {
      WriteOutData(actime, FALSE);
      WriteToDomainFiles(actime);
    }
    if (dumpnow || (dumptime && (actime % dumptime == 0)))  {
      MakeAFullDump();
      dumpnow = 0;
    }
    if (mailtime)  {
      mailtime = 0;
      MailTheTime();
    }
    printf("Energy = %lf\n", CalcTopEnergy());
  }
}


#ifdef PARALLEL

void MastersLoop(int startupflags)
{
  int i;
  double energy;
  long daytime, newtinc;
  Entity et;
  Vector sunpos;
  const ParCommand timestep = DO_TIME_STEP;
  while (plausible && actime < tend)  {
    daytime = ((int)(timeofday*3600.)+actime) % 86400;
    printf("Now at %li seconds = %02li:%02li:%02li  ",
           actime, daytime / 3600, (daytime % 3600)/60, daytime % 60);
    pvm_initsend(PvmDataRaw);
    pvm_pkint((int *)&timestep, 1, 1);
    SendCommand();
    ActualiseValuesInTime(actime);
    for (et = maxentity; et--; )
      CalculateLayerAverage(g[et], pstat, avg+et*nz);
    CalcMeanHydrostaticPressure();
    CheckTimeStep(&tinc, &chemtinc, 0.8);
    if (pressuretype == NONHYDROSTATIC)  {		/* Calc wind acceleration */
      SolveForPressure(tinc);
    }
    CheckTimeStep(&newtinc, &chemtinc, 1.);
    chemtinc = (chemtime <= actime ? chemtinc : 0);
    if (nsubs)  printf("chemtinc = %3ld  ", chemtinc);
    if (newtinc < tinc)  {
      fprintf(stderr, "Warning: Time step had to be changed in the middle of iteration\n");
      tinc = newtinc;
    }
    
    // Call ModuleManager
    McInterface::modmanager.DoCalc(actime+tinc);

    actime += tinc; chemtime += chemtinc;
    if (dumpnow || (dumptime && (actime % dumptime == 0)))  {
      MakeAFullDump();
      dumpnow = 0;
    }
    plausible = IsPlausible(&energy);
    printf("Energy = %lf\n", energy);
    if (!(startupflags & NO_OUTPUT))  {
      WriteOutData(actime, FALSE);
      WriteToDomainFiles(actime);
    }
    if (mailtime)  {
      mailtime = 0;
      MailTheTime();
    }
/*    if (nt != tinc)  {
      tinc = nt;
      printf("Changing \"tinc\"; now %ld seconds\n", tinc);
    } */
  }
}

void WorkersLoop(int startupflags)
{
  double theta;
  long daytime, newtinc;
  Entity et;
  Vector sunpos;
  VarDesc *v;
  int i, j, k;
  char varname[80];
  while (1)  {
    switch (GetCommand())  {
    case DO_TIME_STEP :
      if (plausible)  {
	daytime = ((int)(timeofday*3600.)+actime) % 86400;
	ActualiseValuesInTime(actime);	/* Actualise Border values */
	for (et = maxentity; et--; )
	  CalculateLayerAverage(g[et], pstat, avg+et*nz);
	SetAvgMountains();
	CalcMeanHydrostaticPressure();
	InterpolateToFaces();			/* Interpolate wind speeds from center to faces */
	InterpolateToCenter(-1., pressuretype == NONHYDROSTATIC);
	if (tstart == actime)  Continuity();
	CheckTimeStep(&tinc, &chemtinc, 0.8);
	if (pressuretype == NONHYDROSTATIC)
	  ApplyBuoyancy(tinc);			/* Calculate Buoyancy */
	if (pressuretype != NOPRESSURE)  {		/* Calc wind acceleration */
	  switch (pressuretype)  {			/* Calculate Pressure field */
	  case HYDROSTATIC    : CalcHydrostaticPressure(); break;
	  case NONHYDROSTATIC : SolveForPressure(tinc); break;
	  }
	  ApplyPressure(tinc);
	}
	Continuity();				/* Mass conservation */
	InterpolateToCenter(1., pressuretype == NONHYDROSTATIC);
	if (coriolistype != NOCORIOLIS)  ApplyCoriolis(tinc);
	if (filtertype != NO_FILTER)  {
	  SetBoundary(WWIND);   /* Set Boundary for Wind only */
	  switch (filtertype)  {
	  case PEPPER_FILTER :
	    for (i = nz; i--; )
	      for (et = HUMIDITY; et--; )
		ApplyFilter(g[et]+i*layer, pstat+i*layer, spatialfilter);
	    break;
	  case SHAPIRO_FILTER :
	    for (i = nz; i--; )
	      for (et = HUMIDITY; et--; )
		ShapiroFilter(g[et]+i*layer, pstat+i*layer, 3, spatialfilter);
	    break;
	  case DIFFUSION_FILTER :
	    for (i = nz; i--; )
	      for (et = HUMIDITY; et--; )
		ShapiroFilter(g[et]+i*layer, pstat+i*layer, 1, spatialfilter);
	    break;
	  }
	}
	SetBoundary(WWIND+1);
	CheckTimeStep(&newtinc, &chemtinc, 1.);
	if (newtinc < tinc)  {
	  tinc = newtinc;
	}
	chemtinc = (chemtime <= actime ? chemtinc : 0);
	if (advection || windadvection)
	  switch (advectiontype)  {
	  case MPDATA_A : Advect(tinc, chemtinc); break;
	  case PPM_A    : PPM_Transport(tinc, chemtinc, smiord); break;
	  }
	if (dampinglayer)  DampTop();
	if (groundinterface || shortwaveradiation)  {
	  GroundInterface(timeofday + (double)actime / 3600., dayofyear, tinc, &sunpos);
	  sunelevation = (sunpos.z > 0. ? 57.29578 * asin(sunpos.z) : 0.);
	}
	if (cloudwater)  CloudPhysics();
	if (turbtype == TTTT)  {
	  if (groundinterface)
	    CalcGroundTurbulence(tinc);
	  CalcTransilientTurbulence(tinc);
	}
	else if (turbtype == KEPS_T)  {
	  CalcKEpsilon(tinc);
	  if (groundinterface)
	    CalcGroundTurbulence(tinc);
	  CalcKTurb(tinc, chemtinc);
	}
	Emit(tinc);
	Deposit(tinc);
	if (!(radiationtype && shortwaveradiation) &&
	    // When using the two-stream module, ChemicalTimeStep
	    // mus be called from within ShortWaveRadiation!
	    nsubs && actime % tchem == 0)  ChemicalTimeStep(tchem, &sunpos);
	SetBoundary(chemtinc ? maxentity : SUBS);
    
	// Call ModuleManager
	McInterface::modmanager.DoCalc(actime+tinc);

	actime += tinc; chemtime += chemtinc;
	if (dumpnow || (dumptime && (actime % dumptime == 0)))  {
	  MakeAFullDump();
	  dumpnow = 0;
	}
      }  /* if (plausible) */
      SendStatus(plausible);
      break;
    case SENDGRIDVAL :
      pvm_upkint((int *)&et, 1, 1);
      pvm_upkint(&i, 1, 1);
      pvm_upkint(&j, 1, 1);
      pvm_upkint(&k, 1, 1);
      pvm_initsend(PvmDataRaw);
      pvm_pkdouble(g[et]+i*row+j+k*layer, 1, 1);
      pvm_send(myparent, VALUES);
      break;
    case SENDMESHVAL :
      pvm_upkstr(varname);
      v = GetNamedVar(varname);
      pvm_upkint(&i, 1, 1);
      pvm_upkint(&j, 1, 1);
      pvm_upkint(&k, 1, 1);
      pvm_initsend(PvmDataRaw);
      if (v->storetype == PROC_VAL)  {
	theta = v->v.proc(k, i, j, v);
	pvm_pkdouble(&theta, 1, 1);
      }
      else
	pvm_pkdouble(GetNamedVar(varname)->v.d+i*row+j+k*layer, 1, 1);
      pvm_send(myparent, VALUES);
      break;
    case SENDGRID :
      pvm_upkint((int *)&et, 1, 1);
      SendMeshToMaster(g[et], 0, ALL_DIM);
      break;
    case SENDMESH :
      pvm_upkstr(varname);
      v = GetNamedVar(varname);
      if (v->storetype == PROC_VAL)  {
	for (k = nz; k--; )
	  for (i = nx+1; i--; )
	    for (j = ny+1; j--; )
	      flux[0][i*row+j+k*layer] = v->v.proc(k, i, j, v);
	SendMeshToMaster(flux[0], 0, v->dims);
	memset(flux[0], 0, mesh * sizeof(double));
      }
      else
	SendMeshToMaster(v->v.d, 0, v->dims);
      break;
    case SENDGROUND :
      SendGroundToMaster();
      break;
    case EXIT :
      pvm_exit();
      exit (0);
    }
  }
}

#endif

void GenerateSpecialChem(char *chemname)
{
  FILE *out;
  const char *outname = "mcspecial.c";
  if (!(out = fopen(outname, "w")))  {
    fprintf(stderr, "Error opening \"%s\".\n", outname);
    exit (1);
  }
  if (GenerateSpecial(out, chemname))
    printf("Errors while building \"%s\".\n", outname);
  else
    printf("\"%s\" was successfully created. Compile and link it now!\n",
       outname);
}

int main(int argc, char *argv[])
{
  Entity et;
  int error, startupflags, i;
  char *restartname, name[128];
#ifdef PARALLEL
  const ParCommand end = EXIT;
#endif
  InitVarTable();
#ifdef PARALLEL
  IsMaster();
  if (master) {
    printf("MetPhoMod   Rel. Pre-2.1-0\n" SYSTEM " parallel Version - " DATE "\n\n");
#else
  printf("MetPhoMod   Rel. Pre-2.1-0\n" SYSTEM " sequential Version - " DATE "\n\n");
#endif
    if (argc > 1 && (!strcmp(argv[1], "-help") || !strcmp(argv[1], "-h")))  {
      McInterface::modmanager.Init1();
      PrintHelp(argc - 2, argv + 2);
      exit (0);
    }
    if (argc > 1 && !strcmp(argv[1], "-genchem"))  {
      if (argc != 3)  goto instant_help;
      GenerateSpecialChem(argv[2]);
      exit (0);
    }
    startupflags = AnalizeOptions(&argc, argv, &restartname);
    if (argc != 2)  {
instant_help :
      printf("Usage : meteochem [-RESTART restart-file [-REOPEN]] [-noexcpt] [-syntax] \\\n"
             "                  [-nooutput] [-debug] [-sequential] inpfile\n"
             "or    : meteochem (-h | -help) [keyword]\n"
             "or    : meteochem -genchem chemfile\n");
      return (1);
    }
    McInterface::modmanager.Init1();
    if (startupflags & SYNTAX_ONLY)  {
      if (ParseInput(argv[1]))  {
        printf("Errors in Input-File\n");
        return (1);
      }
      else  {
        printf("the input-file seems to be OK!\n");
        return (0);
      }
    }
#ifdef PARALLEL
    if (startupflags & SEQUENTIAL)  {
      parallel = FALSE; master = TRUE;
      leftest = rightest = TRUE;
    }
  }
  if (master)  {
#endif
    if (ParseInput(argv[1]))  {
      fprintf(stderr, "Errors parsing Input-File. I don't start the calculation.\n");
      return (1);
    }
    plausible = TRUE;
    if (!(startupflags & NOEXCPT))  InstallFpeHandler();
    tinc = chemtinc = tincmax;
    actime = chemtime = tstart;
#ifdef PARALLEL
    if (!parallel)  printf("\n\n!!! Using sequential mode !!!\n\n");
    else
      McInterface::modmanager.CheckIfReadyForParallel();
  }
  InitializeWorkers(argv[0], startupflags & DEBUG_WORKERS);
  if (parallel && !master)  {
    if (!(startupflags & NOEXCPT))  InstallFpeHandler();
    McInterface::modmanager.Init1();
    ReadDataFromMaster();
    InitializeData();
    plausible = TRUE;
    tinc = chemtinc = tincmax;
    actime = chemtime = tstart;
  }
  if (parallel && master)  {
    printf("Sending Data to workers...\n");
    SendDataToWorkers();
    InitializeData();
  }
  if (!parallel)
#endif
    InitializeData();
  for (et = maxentity; et--; )
    if (g[et])
      CalculateLayerAverage(g[et], pstat, avg+et*nz);
    else
      memset(avg+et*nz, 0, nz*sizeof(double));
  CalcMeanHydrostaticPressure();
#ifdef PARALLEL
  if (!parallel || !master)  {
#endif
    SetAvgMountains();
    SetBoundary(maxentity);
/*    CalcAllKm();  */
#ifdef PARALLEL
  }
  if (parallel && !master)  {
    if (GetRestartName(name))  {
      if (!(ReadFullDump(name, (startupflags & REOPEN)) && actime <= tend))  {
	printf("Error : Restart-file not found\n");
	plausible = FALSE;
      }
    }
  }
  else  {
    if (parallel)
      SendRestartName((startupflags & RESTART) ? restartname : NULL);
#endif
    if (startupflags & RESTART)  {
      if (!ReadFullDump(restartname, (startupflags & REOPEN)))  {
        printf("Error : Restart-file not found\n");
        return (1);
      }
    }
    else if (startupflags & REOPEN)  {
      printf("Error: -REOPEN can only be used together with -RESTART\n");
      return (1);
    }
#ifdef PARALLEL
  }
  if (master)  {
#endif
    if (!(startupflags & NO_OUTPUT))  {
      printf("Opening output-Files...");
      if (OpenOutputFiles(argv[1], (startupflags & RESTART) && (startupflags & REOPEN)))  {
        printf("\nExecution abortet\n");
        exit (1);
      }
      CreateDomainFiles((startupflags & RESTART) && (startupflags & REOPEN));
    }
    InstallSignalHandler();
    printf("\n");
#ifdef PARALLEL
    if (parallel)  printf("Starting calculation...\n");
  }
#else  
  printf("Starting calculation...\n");
#endif
#if PARALLEL
  if (!parallel || master)  
#endif
    if (!(startupflags & NO_OUTPUT))  {
      WriteOutData(tstart, FALSE);
      WriteToDomainFiles(tstart);
    }
/*  TestChemicals("Start");  */
#if PARALLEL
  if (parallel)
    if (master)  MastersLoop(startupflags);
    else	 WorkersLoop(startupflags);
  else
#endif
    SequentialLoop(startupflags);
  printf("Closing output-Files...\n");
  if (plausible)  printf("Calculation terminated successfully\n\n");
  else  {
    if (!(startupflags & NO_OUTPUT))
      WriteOutData(actime, TRUE);
    printf("Calculation stopped because of inplausibilities!\n\n");
  }
  if (!(startupflags & NO_OUTPUT))  {
    CloseOutputFiles();
    CloseDomainFiles();
  }
#ifdef PARALLEL
  if (parallel)  {
    pvm_initsend(PvmDataRaw);
    pvm_pkint((int *)&end, 1, 1);
    SendCommand();
    pvm_exit();
  }
#endif
  return (0);
}
