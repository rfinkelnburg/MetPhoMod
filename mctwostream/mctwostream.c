/*
   MODULE mc2stream
   A general two-stream solver. This solver has been derived from
   "ps2str.f" in the Madronich-Paciage
*/

#include "mctwostream.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef PARALLEL
#include <pvm3.h>
#endif
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mchemparse.h"
#include "mchem.h"
#include "mcparse.h"
#include "mcphotodiss.h"
#include "mcabsorb.h"
#include "mc_module.hh"
#include "mcextrad.h"
#ifdef PARALLEL
#include "mcparallel.h"
#endif

static double *tz, *temp, *airconz, *tlay, *cz, *o3, *so2, *no2;
static int tnz = 0;
static Entity o3et, so2et, no2et;
static double *dto3, *dtso2, *dtno2,
  *dtcld, *omcld, *gcld, *dtaer, *omaer, *gaer, *extrad, *waveleft;
static double *edir, *edn, *eup, *fdir, *fdn, *fup, *ftot;
static BOOL radresol;

Entity SearchVariable(char *name)
{
  VarDesc *v;
  for (v = variable; v && strcmp(v->name, name); v = v->next);
  if (v)  return v->v.et;
  else    return Entity(0);
}

void TwoStreamSetup(int nz, double *zcenter, int nwl)
{
  int i, j;
  if (tnz) {
    fprintf(stderr, "Implementation error: TwoStreamSetup called twice for the same object.\n");
    exit (3);
  }
  tnz = nz + int((50000. - zcenter[nz-1] - reflevel) / 1000.);
  o3et = SearchVariable("O3");
  no2et = SearchVariable("NO2");
  so2et = SearchVariable("SO2");

  // dynamically allocate necessary arrays

  if (!(tz = NEW(tnz, double)) ||
      !(temp = NEW(tnz, double)) ||
      !(tlay = NEW(tnz, double)) ||
      !(airconz = NEW(tnz, double)) ||
      !(cz = NEW(tnz, double)) ||
      !(dto3 = NEW(tnz, double)) ||
      !(dtso2 = NEW(tnz, double)) ||
      !(dtno2 = NEW(tnz, double)) ||
      !(dtcld = NEW(tnz, double)) ||
      !(omcld = NEW(tnz, double)) ||
      !(gcld = NEW(tnz, double)) ||
      !(dtaer = NEW(tnz, double)) ||
      !(omaer = NEW(tnz, double)) ||
      !(gaer = NEW(tnz, double)) ||
      !(o3 = NEW(tnz, double)) ||
      !(edir = NEW(tnz, double)) ||
      !(edn = NEW(tnz, double)) ||
      !(eup = NEW(tnz, double)) ||
      !(fdir = NEW(tnz, double)) ||
      !(fdn = NEW(tnz, double)) ||
      !(fup = NEW(tnz, double)) ||
      !(!nsubs || (ftot = NEW(nz * nwl, double))) ||
      (no2et && !(no2 = NEW(tnz, double))) ||
      (so2et && !(so2 = NEW(tnz, double))) ||
      !(extrad = new double[nwl])) {
    fprintf(stderr, "FATAL ERROR: Error allocating memory in \"TwoStream Setup\"\n");
    exit (4);
  }
}


int SetupLevels(void)
{
  int i;
  // Setup vertical levels, use mc-grid in the lower half and complete
  // with values from the standard atmosphere.
  for (i = 0; i < nz; i++) {
    tz[i] = zcenter[i] + reflevel;
  }
  for (double h = tz[nz-1] + 1000.; h < 50000.; h += 1000., i++) {
    tz[i] = h;
  }
  return i;
}

void TempNAirConz(double *temp, double *airconz,
		  int o3et, double *o3,
		  int so2et, double *so2,
		  int no2et, double *no2)
{
  int i;
  for (i = 0; i < nz; i++) {
    temp[i] = AbsoluteTemp(avg[TEMP*nz+i], pp[i], avg[HUMIDITY*nz+i]);
    airconz[i] = pp[i] / (temp[i] * kBoltz) * 1.e-6;
    if (o3et)  o3[i] = avg[o3et*nz+i] * airconz[i] * 1.e-9;
    else       o3[i] = 3.e-8 * airconz[i]; /* take normal background concentration */
    if (so2et)  so2[i] = avg[so2et*nz+i] * airconz[i] * 1.e-9;
    if (no2et)  no2[i] = avg[no2et*nz+i] * airconz[i] * 1.e-9;
  }

  /* Take the rest from US standard atmosphere */
  const static float stddensity[] = {
    2.55E+19, 2.31E+19, 2.09E+19, 1.89E+19, 1.70E+19, 1.53E+19, 1.37E+19,
    1.23E+19, 1.09E+19, 9.71E+18, 8.60E+18, 7.59E+18, 6.49E+18, 5.54E+18,
    4.74E+18, 4.05E+18, 3.46E+18, 2.96E+18, 2.53E+18, 2.16E+18, 1.85E+18,
    1.57E+18, 1.34E+18, 1.14E+18, 9.76E+17, 8.33E+17, 7.12E+17, 6.09E+17,
    5.21E+17, 4.47E+17, 3.83E+17, 3.28E+17, 2.82E+17, 2.41E+17, 2.06E+17,
    1.76E+17, 1.51E+17, 1.30E+17, 1.12E+17, 9.62E+16, 8.31E+16, 7.19E+16,
    6.23E+16, 5.40E+16, 4.70E+16, 4.09E+16, 3.56E+16, 3.11E+16, 2.74E+16,
    2.42E+16, 2.14E+16 };
  const static float stdtemp[] = {
    288.150, 281.651, 275.154, 268.659, 262.166, 255.676, 249.187, 242.700,
    236.215, 229.733, 223.252, 216.774, 216.650, 216.650, 216.650, 216.650,
    216.650, 216.650, 216.650, 216.650, 216.650, 217.581, 218.574, 219.567,
    220.560, 221.552, 222.544, 223.536, 224.527, 225.518, 226.509, 227.500,
    228.490, 230.973, 233.743, 236.513, 239.282, 242.050, 244.818, 247.584,
    250.350, 253.114, 255.878, 258.641, 261.403, 264.164, 266.925, 269.684,
    270.650, 270.650, 270.650 };
  const static float stdo3[] = {
    1.02e+12, 9.7e+11, 9.2e+11, 8e+11, 6.8e+11, 6.3e+11, 5.8e+11, 5.75e+11,
    5.7e+11, 6.1e+11, 6.5e+11, 8.9e+11, 1.13e+12, 1.575e+12, 2.02e+12,
    2.185e+12, 2.35e+12, 2.65e+12, 2.95e+12, 3.495e+12, 4.04e+12, 4.405e+12,
    4.77e+12, 4.815e+12, 4.86e+12, 4.7e+12, 4.54e+12, 4.285e+12, 4.03e+12,
    3.635e+12, 3.24e+12, 2.88e+12, 2.52e+12, 2.275e+12, 2.03e+12, 1.805e+12,
    1.58e+12, 1.4e+12, 1.22e+12, 1.0465e+12, 8.73e+11, 7.4e+11, 6.07e+11,
    5.025e+11, 3.98e+11, 3.36e+11, 2.74e+11, 2.215e+11, 1.69e+11, 1.36e+11,
    1.03e+11, 8.47e+10, 6.64e+10, 5.24e+10, 3.84e+10, 3.195e+10, 2.55e+10,
    2.08e+10, 1.61e+10, 1.365e+10, 1.12e+10, 9.265e+09, 7.33e+09, 6.07e+09,
    4.81e+09, 3.99e+09, 3.17e+09, 2.445e+09, 1.72e+09, 1.235e+09, 7.5e+08,
    6.45e+08, 5.4e+08, 3.8e+08, 2.2e+08, 1.95e+08, 1.7e+08 };

  int j;
  for ( ; i < tnz; i++) {
    j = int(tz[i]/1000.+0.5);
    temp[i] = stdtemp[j];
    airconz[i] = stddensity[j];
    o3[i] = stdo3[j];
  }
  for (i = tnz-1 ; i--; ) {
    tlay[i] = 0.5 * (temp[i+1] + temp[i]);
  }
}

static int SetConz(int i, int j,
		   const double *airconz,
		   int o3et, double *o3,
		   int so2et, double *so2,
		   int no2et, double *no2)
{
  const double conv = 1.e-9;
  int loc, hloc, fa, k, k2;
  hloc = row*i + j;
  fa = ground[hloc].firstabove;
  for (k = fa, k2 = 0; k < nz; k++, k2++) {
    loc = layer * k + hloc;
    double tmp = airconz[k] * conv;
    if (o3et)  o3[k2] = g[o3et][loc] * tmp;
    else       o3[k2] = 3.e-8 * airconz[k];
    if (so2et) so2[k2] = g[so2et][loc] * tmp;
    if (no2et) no2[k2] = g[no2et][loc] * tmp;
  }
  return (fa);
}

class RaleighBins {
  int nw;
  double *srayl;
public :
  RaleighBins(int nw, const double *wl);
  ~RaleighBins(void);
  void CalcRaleigh(double *tz, double *airconz, double *cz, AllZ *dtrl);
};

static RaleighBins *ralbin;
static OzoneBands *o3bin;
static SO2Bands *so2bin;
static NO2Bands *no2bin;

RaleighBins::RaleighBins(int nw_in, const double *wl)
{
  nw = nw_in;
  if (!(srayl = new double[nw])) {
    InputError(ERROR, "Unable to allocate memory in RaleighBins\n");
    exit (3);
  }
  double wc, xx;
  for (int i = nw; i--; ) {
    wc = 5.e-4 * (wl[i] + wl[i+1]);
    if (wc <= 0.55)
      xx = 3.677 + 0.389*wc + 0.09426/wc;
    else
      xx = 4.04;
    srayl[i] = 4.02e-28 / pow(wc, xx);
  }
}

RaleighBins::~RaleighBins(void)
{
  delete srayl;
}

void RaleighBins::CalcRaleigh(double *tz, double *airconz,
			      double *cz, AllZ *dtrl)
{
  double deltaz;
  int i, iw;
  for (i = tnz-1; i--; ) {
    deltaz = 100. * (tz[i+1] - tz[i]);  /* in cm !? */
    cz[i] = (airconz[i+1] - airconz[i]) / log(airconz[i+1] / airconz[i]) * deltaz;
  }
  cz[tnz-2] += 8.05e5 * airconz[tnz-1];
  for (i = tnz-1; i--; )
    for (iw = nw; iw--; )
      dtrl[iw][i] = srayl[iw] * cz[i];
}

TwoStreamModule::TwoStreamModule(void) :
  Module(&McInterface::modmanager, "TwoStream", DIAGNOSTIC,
	 "TwoStream radiation module multi wave radiation module\n"
         "according to Madronich",
	 "alpha")
{}

char *OnTwoStream(RelyCmmd cmmd, VarDesc *v)
{
  switch (cmmd)  {
    case IS_RELEVANT :
       return ((char *)radiationtype);
    case NAME_OF_RELVAR :
       return ("Radiation-Alg");
    case VAL_OF_RELVAR :
       return (char *)(radiationtype ? "twostream" : "clear-sky");
  }
  return (NULL);
}

void TwoStreamModule::CalcTwoStream(const Vector &sunpos, const double sundistance)
{
  static AllZ dtrl[KW], dto2[KW], dto3[KW],
    dtso2[KW], dtno2[KW], dtcld[KW], omcld[KW], gcld[KW],
    dtaer[KW], omaer[KW], gaer[KW];
  static AllZP dsdh[NZ];
  static double topfdr[KW], topedn[KW];
  static int nid[NZ+1];
  static BOOL first = TRUE, chemstep;
  double zen, totrad, totdir;
  int xs, xe;
  if (actime >= mtime + dtime) {
#ifdef PARALLEL
    if (parallel && master) return; // Nothing to do in this case
    if (parallel && !master)  {
      xs = mfirstx + !mfirstx;
      xe = mlastx - rightest;
    }
    else  {
#endif
      xs = 1; xe = nx;
#ifdef PARALLEL
    }
#endif
    chemstep = nsubs && actime % tchem == 0;
    while (actime > mtime)
      mtime += dtime;
    zen = acos(sunpos.z) / w0;
    /* Calculate and interpolate air conzentration */
    if (first) {
      first = FALSE;
      memset(dtrl, 0, KW * NZ * sizeof(double));
      memset(dto2, 0, KW * NZ * sizeof(double));
      memset(dto3, 0, KW * NZ * sizeof(double));
      memset(dtso2, 0, KW * NZ * sizeof(double));
      memset(dtno2, 0, KW * NZ * sizeof(double));
      memset(dtcld, 0, KW * NZ * sizeof(double));
      memset(omcld, 0, KW * NZ * sizeof(double));
      memset(gcld, 0, KW * NZ * sizeof(double));
      memset(dtaer, 0, KW * NZ * sizeof(double));
      memset(omaer, 0, KW * NZ * sizeof(double));
      memset(gaer, 0, KW * NZ * sizeof(double));
    }
    tnz = SetupLevels();               /* Setup levels */
    sphers_(&tnz, tz, &zen, dsdh, nid);
    if (nid[1] >= 0) {
      TempNAirConz(temp, airconz, o3et, o3, so2et, so2, no2et, no2);
      ralbin->CalcRaleigh(tz, airconz, cz, dtrl);
      nwl++;
      seto2_(&tnz, &nwl, wl, cz, &zen, dto2);  /* a fortran routine! */
      nwl--;
      o3bin->OzoneAbsorption(tnz, tz, o3, tlay, stratozone * 1000., dto3);
      if (so2et)
	so2bin->SO2Absorption(tnz, tz, so2, tlay, dtso2);
      if (no2et)
	no2bin->NO2Absorption(tnz, tz, no2, tlay, dtno2);
      if (nsubs) apd.SetCol(temp, airconz);
      totrad = totdir = 0.;
      {
	int iw;
	const double null = 0.;
	double cosradsoil;
	for (iw = 0; iw < nwl; iw++) {
	  rtlink_(&tnz, tz, &mean_albedo, &zen, dsdh, nid,
		  dtrl[iw], dto3[iw], dto2[iw], dtso2[iw], dtno2[iw],
		  dtcld[iw], omcld[iw], gcld[iw],
		  dtaer[iw], omaer[iw], gaer[iw],
		  edir, edn, eup, fdir, fdn, fup, &null, &sunpos.z);
	  double ff;
	  ff = extrad[iw] * sundistance * (waveleft[iw+1] - waveleft[iw]);
	  if (edn[nzm] > 1.e-20 && fdir[nzm] > 1.e-20) {
	    topfdr[iw] = fdir[nzm] * ff;
	    topedn[iw] = edn[nzm] / fdir[nzm];
	  }
	  else {
	    topfdr[iw] = topedn[iw] = 0.;
	  }
	  totrad += ff * edn[0];
	  totdir += ff * fdir[0];
	  if (nsubs && actime % tchem == 0)
	    for (int k = 0; k < nz; k++) {
	      ftot[k*nwl+iw] = ff * (fdir[k] + fdn[k] + fup[k]) * 
		5.039e11 * wc[iw];
	      //	      if (!k)
	      //		printf("%d %d %lg\n", iw, k, ftot[k*nwl+iw]);
	    }
	}
      }
      int fa, nl, hloc;
      double cosradsoil;
      if (radresol) {
	sphers_(&nz, tz, &zen, dsdh, nid);
	for (int ix = xs; ix < xe; ix++)
	  for (int iy = ny; --iy; ) {
	    GroundParam *gr;
	    hloc = ix*row + iy;
	    gr = &ground[hloc];
	    cosradsoil = AngleBetween(&sunpos, &gr->slope);
	    if (cosradsoil < 0.)  cosradsoil = 0.;
	    fa = SetConz(ix, iy, airconz, o3et, o3, so2et, so2, no2et, no2);
	    nl = nz - fa;
	    o3bin->OzoneAbsorption(nl, tz+fa, o3, tlay+fa, 0., dto3);
	    if (so2et)
	      so2bin->SO2Absorption(nl, tz+fa, so2, tlay+fa, dtso2);
	    if (no2et)
	      no2bin->NO2Absorption(nl, tz+fa, no2, tlay+fa, dtno2);
	    totrad = totdir = 0.;
	    for (int iw = nwl; iw--; ) {
	      double ff = 1. - gr->absorb;
	      rtlink_(&nl, tz+fa, &ff, &zen, dsdh, nid,
		      dtrl[iw]+fa, dto3[iw], dto2[iw]+fa, dtso2[iw], dtno2[iw],
		      dtcld[iw], omcld[iw], gcld[iw],
		      dtaer[iw], omaer[iw], gaer[iw],
		      edir, edn, eup, fdir, fdn, fup, &topedn[iw], &cosradsoil);
	      ff = topfdr[iw];
	      totrad += ff * edn[0];
	      if (cosradsoil > 0.)
		totdir += ff * fdir[0] * cosradsoil;
	      if (nsubs && actime % tchem == 0)
		for (int k = nl; k--; )
		  ftot[k*nwl+iw] = ff * (fdir[k] + fdn[k] + fup[k]) *
		    5.039e11 * wc[iw];
	    }
	    if (chemstep)
	      for (int k = nl; k--; ) {
		apd.CalcRates(ftot+k*nwl, k+fa, photorate);
		BoxChemStep(ix, iy, k+fa, tchem, &sunpos, photorate);
	      }
	    gr->Rdir = totdir;
	    gr->Rg = totrad + gr->Rdir;
	    gr->Rdiff = totrad;	  
	  }
      }
      else {  // if (!radresol) ...
	for (int ix = xs; ix < xe; ix++)
	  for (int iy = ny; --iy; ) {
	    GroundParam *gr;
	    gr = &ground[ix*row+iy];
	    fa = gr->firstabove;
	    cosradsoil = AngleBetween(&gr->slope, &sunpos);
	    if (cosradsoil < 0.) cosradsoil = 0;
	    gr->Rdir = totdir * cosradsoil;
	    gr->Rg = totrad + gr->Rdir;
	    gr->Rdiff = totrad;
	    if (chemstep)
	      for (int k = fa; k < nz; k++) {
		apd.CalcRates(ftot+k*nwl, k, photorate);
		BoxChemStep(ix, iy, k, tchem, &sunpos, photorate);
	      }
	  }
      }
    }
    else {
      if (chemstep)
	memset(photorate, 0, (apd.NUsed()+1) * sizeof(double));
      for (int ix = xs; ix < xe; ix++)
	for (int iy = ny; --iy; ) {
	  GroundParam *gr;
	  gr = &ground[ix*row+iy];
	  int fa = gr->firstabove;
	  gr->Rdir = gr->Rg = gr->Rdiff = 0.;
	  if (chemstep)
	    for (int k = fa; k < nz; k++)
	      BoxChemStep(ix, iy, k, tchem, &sunpos, photorate);
	}
    }
  }   // if (mtime   ...etc
}

void TwoStreamModule::Init1(void)
{
  wl = new double[2000];
  static VectPtr vp;

  // Define Input-sections
  InputSection mysection =
    InputSection((sectdesc = McInterface::mcsect.AddASection(
      "RADIATION", McInterface::mcsect.FindSection(ENVIRONMENT),
      NULL, NULL, CAN_BE_OMITTED))->id);
  vp.nval = &nwl;
  vp.v = wl;

  // Define new Input-Variables
  VarDesc *nv;
  nv = new VarDesc("dt-twostream",
		   "The time step for the two-stream module",
		   "sec", INT_PTR, &dtime, 0, mysection, MUST_BE_SET,
		   120, 1, 86400, INT_NUM, 0, OnTwoStream);
  McInterface::mcvar.AddNewVariable(nv);
  nv = new VarDesc("horizontal-rad-resolution",
		   "Calculate two-stream radiation in every model column",
		   NULL, INT_PTR, &radresol, 0, mysection, SET_DEFAULT,
		   1, 0, 1, ON_OFF, 0, OnTwoStream);
  McInterface::mcvar.AddNewVariable(nv);
  nv = new VarDesc("waveleft", "Left edges of the wavebins of the two stream module",
		   "nm", VECT_PTR, &vp, 2000, mysection, MUST_BE_SET,
		   0, 0, 5000, VECTOR_VAL, 0, OnTwoStream);
  McInterface::mcvar.AddNewVariable(nv);

  // Instantiate the PhotoDiss-Module
  InstantiatePhotoDiss(apd);  // This function is defined in "mcphotodiss.h"
}

void TwoStreamModule::EndOfSection(SectionDesc *section)
{
  if (radiationtype) {
    if (section == sectdesc) {
      if (nsubs && !nsubst) {
	InputError(ERROR, "When using the TwoStream and CHEMISTRY modules together, RADIATION has to be defined after CHEMISTRY.");
	return;
      }
      if (dtime % *tincarray != 0) {
	InputError(ERROR, "\"dt-twostream\" must be a multiple of biggest "
		   "possible time step (%d).", *tincarray);
	return;
      }
      if (chemtinc % dtime != 0) {
	InputError(ERROR, "dt-chem must be a multiple of dt-twostream");
	return;
      }
      if (nwl < 2) {
	InputError(ERROR, "I need at least two values for waveleft.");
	return;
      }
      int i;
      for (i = nwl-1; i--; )
	if (wl[i] >= wl[i+1]) {
	  InputError(ERROR, "waveleft values must be sorted in ascending order.\n"
		     "I have found: %g, %g.\n", wl[i], wl[i+1]);
	  inputerror = 1;
	  return;
	}
      wl = (double *)realloc(wl, nwl * sizeof(double));
      waveleft = wl;
      nwl--;  /* The rightmost wl-value serves as the right boundary.
		 Only nwl-1 wavelength-bins are valid, though */
      wc = new double[nwl];
      for (i = nwl; i--; )
	wc[i] = 0.5 * (wl[i] + wl[i+1]);
      TwoStreamSetup(nz, zcenter, nwl);
      SetupExtRadiation(nwl, wl, extrad);
      ralbin = new RaleighBins(nwl, wl);
      o3bin = new OzoneBands(nwl, wl);
      if (so2et)  so2bin = new SO2Bands(nwl, wl);
      if (no2et)  no2bin = new NO2Bands(nwl, wl);
      mean_albedo = 0.;
      // Setup photodissociation reactions
      apd.SetupReactions(nwl, wl, nz);
    }
    else if (section == McInterface::mcsect.FindSection("OPTIONS")) {
      sectdesc->calltype = MUST_BE_CALLED;
      InputError(WARNING, "The module \"%s\" is still very alpha. Handle results with care.", Name());
    }
  }
}

void TwoStreamModule::Init2(void)
{
  // Calculate mean albedo
  for (int i = nx; --i; ) 
    for (int j = ny; --j; )
      mean_albedo += ground[i*row+j].absorb;
  mean_albedo = 1. - mean_albedo / (nxm * nym);
  mtime = -dtime;
  if (nsubs)
    photorate = (new double[apd.NUsed()]) - 1;
}

#ifdef PARALLEL

#define TSDATA 90  // Would be nice if ID is unique, but not really necessary

char TwoStreamModule::ReadyForParallel(void) const
{
  return (TRUE);
}

void TwoStreamModule::SendToWorkers(void)
{
  if (radiationtype) {
    pvm_initsend(PvmDataRaw);
    pvm_pkint(&radresol, 1, 1);
    pvm_pklong(&dtime, 1, 1);
    pvm_pkint(&nwl, 1, 1);
    pvm_pkdouble(wl, nwl+1, 1);
    for (int w = 0; w < workers; w++)
      pvm_send(worker[w].tid, TSDATA);
  }
}

void TwoStreamModule::GetFromMaster(void)
{
  if (radiationtype) {
    pvm_recv(myparent, TSDATA);
    pvm_upkint(&radresol, 1, 1);
    pvm_upklong(&dtime, 1, 1);
    pvm_upkint(&nwl, 1, 1);
    nwl++;
    pvm_upkdouble(wl, nwl, 1);
    EndOfSection(sectdesc);
  }
}

#endif

long TwoStreamModule::DoCalc(long mintime, long)
{
  return (mtime);
}

long TwoStreamModule::Time(void) const
{
  return mtime;
}

void TwoStreamModule::PrintPhotoDissReact(void) const
{
  apd.PrintReactions();
}

int TwoStreamModule::RegisterReaction(int id)
{
  return apd.UseReaction(id);
}

TwoStreamModule twostreammodule;
