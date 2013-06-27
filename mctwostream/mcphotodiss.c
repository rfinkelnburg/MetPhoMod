/*
   MODULE MCPhotodiss
   Calculates cross sections and quantum yields for several
   photodissoziation reactions
*/

#include "mcphotodiss.h"
#include "mcapd.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mc_group.hh>
#include <mc_group.t>
#include "mcextrad.h"
#include "mcglobal.h"

#ifdef GCC

template class Group<PhotoDissReact>;
template class GroupEnumerator<PhotoDissReact>;

#endif

class PhotoDiss1 : public PhotoDissReact {
  double *f0, *f1;
  const double *wl, *tabs;
  int nwl, from, to;
public :
  PhotoDiss1(AllPhotoDiss *apd) : PhotoDissReact(apd, "NO2 -> NO + O(3P)")
    {
      f0 = f1 = NULL;
    }
  ~PhotoDiss1(void);
  virtual void Setup(int nwl, const double *wl, int nz);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;
};

PhotoDiss1::~PhotoDiss1(void)
{
  if (f0)  delete f0;
  if (f1)  delete f1;
}

#define SIG1 57
#define QY1 65

void PhotoDiss1::Setup(int nw_in, const double *wl_in, int nz)
     /* NO2 -> NO + O(3P) */ 
{ 
  const float swltpl[SIG1+1] = {
    202.02, 204.08, 206.19, 208.33, 210.53, 212.77, 215.09, 217.39, 
    219.78, 222.22, 224.72, 227.27, 229.89, 232.56, 235.29, 238.09, 
    240.96, 243.90, 246.91, 250.00, 253.17, 256.41, 259.74, 263.16, 
    266.67, 270.27, 273.97, 277.78, 281.69, 285.71, 289.85, 294.12, 
    298.51, 303.03, 307.69, 312.50, 317.5, 322.5, 327.5, 332.5, 
    337.5, 342.5, 347.5, 352.5, 357.5, 362.5, 367.5, 372.5, 
    377.5, 382.5, 387.5, 392.5, 397.5, 402.5, 407.5, 412.5, 
    417.5, 422.5 };
  const float sigt[SIG1] = {
    4.145e-19, 4.478e-19, 4.454e-19, 4.641e-19, 4.866e-19, 4.818e-19, 
    5.022e-19, 4.441e-19, 4.713e-19, 3.772e-19, 3.929e-19, 2.74e-19, 
    2.778e-19, 1.689e-19, 1.618e-19, 8.812e-20, 7.472e-20, 3.909e-20, 
    2.753e-20, 2.007e-20, 1.973e-20, 2.111e-20, 2.357e-20, 2.698e-20, 
    3.247e-20, 3.785e-20, 5.03e-20, 5.88e-20, 7e-20, 8.15e-20, 
    9.72e-20, 1.154e-19, 1.344e-19, 1.589e-19, 1.867e-19, 2.153e-19, 
    2.477e-19, 2.807e-19, 3.133e-19, 3.425e-19, 3.798e-19, 4.065e-19, 
    4.313e-19, 4.717e-19, 4.833e-19, 5.166e-19, 5.315e-19, 5.508e-19, 
    5.644e-19, 5.757e-19, 5.927e-19, 5.845e-19, 6.021e-19, 5.781e-19, 
    5.999e-19, 5.651e-19, 5.812e-19 };
  const float sigslope[SIG1] = {
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 
    0, 0, 7.5e-24, 8.2e-24, -5.3e-24, -4.3e-24, 
    -3.1e-24, -1.62e-23, -2.84e-23, -3.57e-23, -5.36e-23, -6.86e-23, 
    -7.86e-23, -1.105e-22, -1.355e-22, -1.277e-22, -1.612e-22, -1.89e-22, 
    -1.219e-22, -1.921e-22, -1.095e-22, -1.322e-22, -1.102e-22, -8.06e-23, 
    -8.67e-23, -9.45e-23, -9.23e-23, -7.38e-23, -5.99e-23, -5.45e-23, 
    -1.129e-22, 1e-25, -1.208e-22 };
  const float qywltpl[QY1+1] = {
    0, 285, 290, 295, 300, 305, 310, 315, 320, 
    325, 330, 335, 340, 345, 350, 355, 360, 
    365, 370, 375, 380, 381, 382, 383, 384, 
    385, 386, 387, 388, 389, 390, 391, 392, 
    393, 394, 395, 396, 397, 398, 399, 400, 
    401, 402, 403, 404, 405, 406, 407, 408, 
    409, 410, 411, 412, 413, 414, 415, 416, 
    417, 418, 419, 420, 421, 422, 423, 424, 
    450};
  const float qytpl[QY1] = {
    1., 1.00, .999, .998, .997, .996, .995, .994, .993, 
    .992, .991, .990, .989, .988, .987, .986, .984, 
    .983, .981, .979, .975, .974, .973, .972, .971, 
    .969, .967, .966, .964, .962, .960, .959, .957, 
    .953, .950, .942, .922, .870, .820, .760, .695, 
    .635, .560, .485, .425, .350, .290, .225, .185, 
    .153, .130, .110, .094, .083, .070, .059, .048, 
    .039, .030, .023, .018, .012, .008, .004, .000 };
  nwl = nw_in;
  wl = wl_in;
  f0 = new double[nwl];
  f1 = new double[nwl];
  double *sig = new double[nwl];
  if (!f0 || !f1 || !sig) {
    fprintf(stderr, "FATAL: Error allocating memory in \"PhotoDiss1::Setup\"\n");
    exit (3);
  }
  Rebin(SIG1, swltpl, sigt, nwl, wl, f0);
  Rebin(SIG1, swltpl, sigslope, nwl, wl, f1);
  Rebin(QY1, qywltpl, qytpl, nwl, wl, sig);
  for (int i = nwl; i--; ) {
    f0[i] *= sig[i];
    f1[i] *= sig[i];
  }
  for (from = 0; from < nwl && f0[from] == 0. && f1[from] == 0.; from++);
  for (to = nwl-1; to >= 0 && f0[to] == 0. && f1[to] == 0.; to--);
  delete sig;
}

void PhotoDiss1::SetCol(int serialid, const double *tabs_in, const double *airc)
{
  tabs = tabs_in;
}

double PhotoDiss1::CalcRate(const double *flux, int k) const
{
  int i;
  double rate = 0.;
  for (i = from; i <= to; i++)
    rate += flux[i] * (f0[i] + tabs[k] * f1[i]);
  return (rate);
}

class PhotoDiss2a : public PhotoDissReact {
  double *s226, *s263, *s298, *qy, *qexp;
  double *xso3;
  double *so3, *q;
  const double *wl;
  int nwl, sid, setup;
public :
  PhotoDiss2a(AllPhotoDiss *apd) : PhotoDissReact(apd, "O3 -> O2 + O(1D)")
    {
      wl = xso3 = s226 = s263 = s298 = qy = qexp = so3 = q = NULL;
      setup = 0;
    }
  ~PhotoDiss2a(void);
  virtual void Setup(int nwl_in, const double *wl_in, int nz);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;
  double CalcRate2(const double *flux, int k) const;
};

void PhotoDiss2a::Setup(int nwl_in, const double *wl_in, int nz)
{
#include "ph_o3data.c"
#include "ph2data.c"
  if (setup) return;
  nwl = nwl_in; wl = wl_in;
  sid = -1;
  xso3 = new double[nwl];
  s226 = new double[nwl];
  s263 = new double[nwl];
  s298 = new double[nwl];
  qy = new double[nwl];
  qexp = new double[nwl];
  so3 = new double[nwl*nz];
  q = new double[nwl*nz];
  double *p1, *p2;
  p1 = new double[nwl];
  p2 = new double[nwl];
  if (!s226 || !s263 || !s298 || !qy || !qexp || !p1 || !p2 ||
      !xso3 || !so3 || !q) {
    fprintf(stderr, "FATAL: Unable to allocate memory in PhotoDiss2a::PhotoDiss2a\n");
    exit (3);
  }
  Rebin(NSO3, swl, s226tpl, nwl, wl, s226);
  Rebin(NSO3, swl, s263tpl, nwl, wl, s263);
  Rebin(NSO3, swl, s298tpl, nwl, wl, s298);
  int i;
  for (i = nwl; i--; ) {
    const double sfct = 1.e-20;
    s226[i] *= sfct;
    s263[i] *= sfct;
    s298[i] *= sfct;
  }
  Rebin(NXO3, xwl, xtpl, nwl, wl, xso3);
  Interpolate(QY2, qwl, phi298, nwl, wl, p1);
  Interpolate(QY2, qwl, phi230, nwl, wl, p2);
  for (i = nwl; i--; ) {
    if (p1[i] != 0.) {
      qexp[i] = log(p1[i] / p2[i]) / (1./230. - 1./298.);
      qy[i] = exp(log(p1[i]) + qexp[i] / 298.);
    }
    else
      qy[i] = qexp[i] = 0.;
  }
  setup = 1;
  delete p1, p2;
}

PhotoDiss2a::~PhotoDiss2a(void)
{
  if (xso3)  delete xso3;
  if (s226)  delete s226;
  if (s263)  delete s263;
  if (s298)  delete s298;
  if (qy)  delete qy;
  if (qexp)  delete qexp;
  if (so3)  delete so3;
  if (q)  delete q;
}

void PhotoDiss2a::SetCol(int serialid, const double *tabs, const double *airc)
{
  int i, k, loc;
  if (serialid != sid) {
    for (k = nz; k--; )
      for (i = 0; i < nwl; i++) {
	loc = k*nwl+i;
	if (wl[i] > 240.5 && wl[i+1] < 350.)
	  if (tabs[k] < 263)
	    so3[loc] = s226[i] + (s263[i] - s226[i]) / (263.-226.) *
	      (tabs[k]-226.);
	  else
	    so3[loc] = s263[i] + (s298[i] - s263[i]) / (298.-263.) *
	      (tabs[k]-263.);
	else
	  so3[loc] = xso3[i];
	q[loc] = qy[i] * exp(-qexp[i] / tabs[k]);
      }
    sid = serialid;
  }
}

double PhotoDiss2a::CalcRate(const double *flux, int k) const
{
  double res = 0.;
  for (int iw = nwl; iw--; )
    res += so3[iw] * q[iw] * flux[iw];
  return (res);
}

double PhotoDiss2a::CalcRate2(const double *flux, int k) const
{
  double res = 0.;
  for (int iw = nwl; iw--; )
    res += so3[iw] * (1. - q[iw]) * flux[iw];
  return (res);
}

class PhotoDiss2b : public PhotoDissReact {
  const double *wl;
  int nwl;
  PhotoDiss2a *o3b;
public :
  PhotoDiss2b(AllPhotoDiss *apd, PhotoDiss2a *o3b_in) : 
    PhotoDissReact(apd, "O3 -> O2 + O(3P)")
    { o3b = o3b_in; }
  virtual void Setup(int nwl, const double *wl, int nz_in);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;
};

void PhotoDiss2b::Setup(int nwl, const double *wl, int nz)
{
  o3b->Setup(nwl, wl, nz);
}

void PhotoDiss2b::SetCol(int serialid, const double *tabs, const double *airc)
{
  o3b->SetCol(serialid, tabs, airc);
}

double PhotoDiss2b::CalcRate(const double *flux, int k) const
{
  return (o3b->CalcRate2(flux, k));
}

class SimplePhotoDiss : public PhotoDissReact {
  int first, last;
  const double *f0;
public :
  SimplePhotoDiss(AllPhotoDiss *apd, const char *name) : 
    PhotoDissReact(apd, name)
    { }
  virtual void Setup(int nwl, const double *f0, int nz);
  virtual double CalcRate(const double *flux, int k) const;
};

void SimplePhotoDiss::Setup(int nwl, const double *f0_in, int nz)
{
  f0 = f0_in;
  for (first = 0; first < nwl && f0[first] == 0.; first++);
  for (last = first; last < nwl && f0[last] != 0.; last++);
}

double SimplePhotoDiss::CalcRate(const double *flux, int k) const
{
  double rate = 0.;
  for (int i = first; i <= last; i++)
    rate += flux[i] * f0[i];
  return rate;
}

class PhotoDiss3 : public SimplePhotoDiss {
  double *f0;
public :
  PhotoDiss3(AllPhotoDiss *apd) : SimplePhotoDiss(apd, "NO3 -> NO + O2")
    { }
  ~PhotoDiss3(void);
  virtual void Setup(int nwl, const double *wl, int nz);
};

PhotoDiss3::~PhotoDiss3(void)
{
  delete f0;
}

void PhotoDiss3::Setup(int nwl, const double *wl, int nz)
{
#include "ph3data.c"  // contains the spectra of the reaction
  double q;
  f0 = new double[nwl];
  if (!f0) {
    fprintf(stderr, "Error allocating memory in PhotoDiss3::Setup.\n");
    exit (3);
  }
  Rebin(NWL3, swl3, s0_3, nwl, wl, f0);
  for (int iw = nwl; iw--; ) {
    double wc;
    wc = 0.5 * (wl[iw] + wl[iw+1]);
    if (wc < 584. || wc >= 640.)  q = 0.;
    else if (wc >= 595.)  q = 0.35 * (1. - (wc - 595.) / 45.);
    else  q = 0.35 * (wc - 584.) / 11.;
    f0[iw] *= q;
  }
  SimplePhotoDiss::Setup(nwl, f0, nz);
}

class PhotoDiss4 : public SimplePhotoDiss {
  double *f0;
public :
  PhotoDiss4(AllPhotoDiss *apd) : SimplePhotoDiss(apd, "NO3 -> NO2 + O(3P)")
    { }
  ~PhotoDiss4(void);
  virtual void Setup(int nwl, const double *wl, int nz);
};

PhotoDiss4::~PhotoDiss4(void)
{
  delete f0;
}

void PhotoDiss4::Setup(int nwl, const double *wl, int nz)
{
#include "ph4data.c"  // contains the spectra of the reaction
  double q;
  f0 = new double[nwl];
  if (!f0) {
    fprintf(stderr, "Error allocating memory in PhotoDiss4::Setup.\n");
    exit (3);
  }
  Rebin(NWL4, swl4, s0_4, nwl, wl, f0);
  for (int iw = nwl; iw--; ) {
    double wc;
    wc = 0.5 * (wl[iw] + wl[iw+1]);
    if (wc < 584.) q = 1.;
    else if (wc >= 640.)  q = 0.;
    else if (wc >= 595.)  q = 0.65 * (1. - (wc - 595.) / 45.);
    else  q = 1. - 0.35 * (wc - 584.) / 11.;
    f0[iw] *= q;
  }
  SimplePhotoDiss::Setup(nwl, f0, nz);
}

class PhotoDiss5a : public PhotoDissReact {
  const double *wl;
  double *sstd, *sn2o5, *qy;
  int nwl, sid, setup;
public :
  PhotoDiss5a(AllPhotoDiss *apd) : PhotoDissReact(apd, "N2O5 -> NO3 + NO + O(3P)")
    {
      wl = sstd = sn2o5 = qy = NULL;
      setup = 0;
    }
  ~PhotoDiss5a(void);
  virtual void Setup(int nwl_in, const double *wl_in, int nz);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;
  double CalcRate2(const double *flux, int k) const;
};

PhotoDiss5a::~PhotoDiss5a(void)
{
  if (sstd)  delete sstd;
  if (sn2o5)  delete sn2o5;
  if (qy)  delete qy;
}

#define NWL5 17

void PhotoDiss5a::Setup(int nwl_in, const double *wl_in, int nz)
{
  const float wltpl[NWL5] = {
    200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255,
    260, 265, 270, 275, 280};
  const float stpl[NWL5] = {
    920e-20, 820e-20, 560e-20, 370e-20, 220e-20, 144e-20, 99e-20,
    77e-20, 62e-20, 52e-20, 40e-20, 32e-20, 26e-20, 20e-20, 16.1e-20,
    13.0e-20, 11.7e-20};
  if (setup) return;
  nwl = nwl_in; wl = wl_in;
  sstd = new double[nwl];
  sn2o5 = new double[nwl*nz];
  qy = new double[nwl];
  if (!sstd || ! sn2o5 || !qy) {
    fprintf(stderr, "Unable to allocate memory in PhotoDiss5a::Setup\n");
    exit (2);
  }
  Interpolate(NWL5, wltpl, stpl, nwl, wl, sstd);
  for (int iw = nwl; iw--; ) {
    double wc = 0.5 * (wl[iw]+wl[iw+1]);
    if (wl[iw] > 400.)  sstd[iw] = 0.;
    else if (wl[iw] >= 280.)
      sstd[iw] = 1.e-20 * exp(2.735 + (4728.5-17.127*wc)/262.5);
    double tmp = 3.832441 - 0.012809638 * wc;
    if (tmp > 1.)  tmp = 1.;
    if (tmp < 0.)  tmp = 0.;
    qy[iw] = tmp;
  }
  setup = 1;
  sid = -1;
}

void PhotoDiss5a::SetCol(int serialid, const double *tabs, const double *airc)
{
  if (serialid != sid) {
    for (int k = nz; k--; )
      for (int iw = nwl; iw--; ) {
	int loc = k*nwl+iw;
	double t, wc = (wl[iw] + wl[iw+1]) * 0.5;
	if (wc >= 285. && wc >= 380.) {
	  t = tabs[k];
	  if (t < 225.)  t = 225.;
	  else if (t > 300.)  t = 300.;
	  sn2o5[loc] = 1.e-20 * exp(2.735 + (4728.5-17.127*wc) / t);
	}
	else
	  sn2o5[loc] = sstd[iw];
      }
    sid = serialid;
  }
}

double PhotoDiss5a::CalcRate(const double *flux, int k) const
{
  double rate = 0.;
  for (int iw = nwl; iw--; )
    rate += qy[iw] * sn2o5[k*nwl+iw] * flux[iw];
  return (rate);
}

double PhotoDiss5a::CalcRate2(const double *flux, int k) const
{
  double rate = 0.;
  for (int iw = nwl; iw--; )
    rate += (1. - qy[iw]) * sn2o5[k*nwl+iw] * flux[iw];
  return (rate);
}

class PhotoDiss5b : public PhotoDissReact {
  const double *wl;
  int nwl;
  PhotoDiss5a *n2o5a;
public :
  PhotoDiss5b(AllPhotoDiss *apd, PhotoDiss5a *n2o5a_in) : 
    PhotoDissReact(apd, "N2O5 -> NO3 + NO2")
    { n2o5a = n2o5a_in; }
  virtual void Setup(int nwl, const double *wl, int nz_in);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;
};

void PhotoDiss5b::Setup(int nwl, const double *wl, int nz)
{
  n2o5a->Setup(nwl, wl, nz);
}

void PhotoDiss5b::SetCol(int serialid, const double *tabs, const double *airc)
{
  n2o5a->SetCol(serialid, tabs, airc);
}

double PhotoDiss5b::CalcRate(const double *flux, int k) const
{
  return (n2o5a->CalcRate2(flux, k));
}

class PhotoDiss6 : public SimplePhotoDiss {
  double *f0;
public :
  PhotoDiss6(AllPhotoDiss *apd) : SimplePhotoDiss(apd, "HNO2 -> OH + NO")
    { f0 = NULL; }
  ~PhotoDiss6(void);
  virtual void Setup(int nwl, const double *wl, int nz);
};

PhotoDiss6::~PhotoDiss6(void)
{
  if (f0)  delete f0;
}

void PhotoDiss6::Setup(int nwl, const double *wl, int nz)
{
#include "ph6data.c"  // contains the spectra of the reaction
  f0 = new double[nwl];
  if (!f0) {
    fprintf(stderr, "Error allocating memory in PhotoDiss4::Setup.\n");
    exit (3);
  }
  Interpolate(NWL6, swl6, s0_6, nwl, wl, f0);
  SimplePhotoDiss::Setup(nwl, f0, nz);
}

class PhotoDiss7 : public PhotoDissReact {
  double *qa, *qb, *sq;
  int first, last, nwl;
public :
  PhotoDiss7(AllPhotoDiss *apd) : PhotoDissReact(apd, "HNO3 -> OH + NO2")
    { qa = qb = sq = NULL; }
  ~PhotoDiss7(void);
  virtual void Setup(int nwl, const double *wl, int nz);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;
};

PhotoDiss7::~PhotoDiss7(void)
{
  if (qa)  delete qa;
  if (qb)  delete qb;
  if (sq)  delete sq;
}

void PhotoDiss7::Setup(int nwl_in, const double *wl, int nz)
{
#include "ph7data.c"
  nwl = nwl_in;
  qa = new double[nwl];
  qb = new double[nwl];
  sq = new double[nwl*nz];
  if (!qa || ! qb || !sq) {
    fprintf(stderr, "Error allocating memory in \"PhotoDiss7::PhotoDiss7\".\n");
    exit (3);
  }
  Interpolate(NWL7, swl7, s0_7, nwl, wl, qa);
  Interpolate(NWL7, swl7, sexp7, nwl, wl, qb);
  for (first = 0; first < nwl && qa[first] == 0.; first++);
  for (last = nwl-1; last > 0 && qa[last] == 0.; last--);
}

void PhotoDiss7::SetCol(int serialid, const double *tabs, const double *airc)
{
  for (int k = nz; k--; )
    for (int iw = nwl; iw--; )
      sq[k*nwl+iw] = qa[iw] * exp(qb[iw] * (tabs[k] - 298.));
}

double PhotoDiss7::CalcRate(const double *flux, int k) const
{
  double res = 0.;
  double *s = sq + k*nwl;
  for (int iw = nwl; iw--; )
    res += flux[iw] * s[iw];
  return res;
}

class PhotoDiss8 : public SimplePhotoDiss {
  double *f0;
public :
  PhotoDiss8(AllPhotoDiss *apd) : SimplePhotoDiss(apd, "HNO4 -> HO2 + NO2")
    { f0 = NULL; }
  ~PhotoDiss8(void);
  virtual void Setup(int nwl, const double *wl, int nz);
};

PhotoDiss8::~PhotoDiss8(void)
{
  if (f0)  delete f0;
}

void PhotoDiss8::Setup(int nwl, const double *wl, int nz)
{
#define NWL8 28
  const float wltpl[NWL8] = {
    190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 
    250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 
    310, 315, 320, 325 };
  const float stpl[NWL8] = {
    1.01e-17, 8.16e-18, 5.63e-18, 3.67e-18, 2.39e-18, 1.61e-18, 1.18e-18, 9.35e-19, 
    7.92e-19, 6.82e-19, 5.81e-19, 4.89e-19, 4.12e-19, 3.5e-19, 2.85e-19, 2.3e-19, 
    1.81e-19, 1.34e-19, 9.3e-20, 6.2e-20, 3.9e-20, 2.4e-20, 1.4e-20, 9e-21, 
    5e-21, 3e-21, 2e-21, 1e-21 };
  f0 = new double[nwl];
  if (!f0) {
    fprintf(stderr, "Error allocating memory in PhotoDiss4::Setup.\n");
    exit (3);
  }
  Interpolate(NWL8, wltpl, stpl, nwl, wl, f0);
  SimplePhotoDiss::Setup(nwl, f0, nz);
}

class PhotoDiss9 : public PhotoDissReact {
  int nwl, first, last;
  const double *wl;
  double *qs, *s;
public :
  PhotoDiss9(AllPhotoDiss *apd) : PhotoDissReact(apd, "H2O2 -> 2 OH")
    {  }
  ~PhotoDiss9(void);
  virtual void Setup(int nwl, const double *wl, int nz);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;
};

void PhotoDiss9::Setup(int nwl_in, const double *wl_in, int nz)
{
  nwl = nwl_in;
  wl = wl_in;
#define NWL9 33
  const float wltpl[NWL9] = {
    190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 
    250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 
    310, 315, 320, 325, 330, 335, 340, 345, 350 };
  const float qstpl[NWL9] = {
    6.72e-19, 5.64e-19, 4.75e-19, 4.08e-19, 3.57e-19, 3.07e-19, 2.58e-19, 2.17e-19, 
    1.82e-19, 1.5e-19, 1.24e-19, 1.02e-19, 8.3e-20, 6.7e-20, 5.3e-20, 4.2e-20, 
    3.3e-20, 2.6e-20, 2e-20, 1.5e-20, 1.2e-20, 9e-21, 6.8e-21, 5.1e-21, 
    3.9e-21, 2.9e-21, 2.2e-21, 1.6e-21, 1.3e-21, 1e-21, 7e-22, 5e-22, 
    4e-22 };
  qs = new double[nwl];
  s = new double[nwl*nz];
  if (!qs || !s) {
    fprintf(stderr, "Unable to allocate memory in PhotoDiss9::Setup.\n");
    exit (3);
  }
  Interpolate(NWL9, wltpl, qstpl, nwl, wl, qs);
  for (first = 0; first < nwl && qs[first] == 0.; first++);
  for (last = nwl-1; last > 0 && qs[last] == 0.; last--);
}

void PhotoDiss9::SetCol(int serialid, const double *tabs, const double *airc)
{
  const double
    A0 = 6.4761E+04,
    A1 = -9.2170972E+02,
    A2 = 4.535649,
    A3 = -4.4589016E-03,
    A4 = -4.035101E-05,
    A5 = 1.6878206E-07,
    A6 = -2.652014E-10,
    A7 = 1.5534675E-13,
    B0 = 6.8123E+03,
    B1 = -5.1351E+01,
    B2 = 1.1522E-01,
    B3 = -3.0493E-05,
    B4 = -1.0924E-07;
  for (int iw = nwl; iw--; ) {
    if (wl[iw] >= 260. && wl[iw] < 350.) {
      double wc = 0.5 * (wl[iw] + wl[iw+1]);
      double sumA, sumB;
      sumA = ((((((A7*wc + A6)*wc + A5)*wc + 
		 A4)*wc +A3)*wc + A2)*wc + 
	      A1)*wc + A0;
      sumB = (((B4*wc + B3)*wc + B2)*wc + 
	      B1)*wc + B0;
      for (int k = nz; k--; ) {
	double t = tabs[k];
	if (t < 200.) t = 200.;
	if (t > 400.) t = 400.;
	double chi = 1. / (1 + exp(-1265. / t));
	s[k*nwl+iw] = (chi * sumA + (1 - chi)* sumB) * 1.e-21;
      }
    }
    else
      for (int k = nz; k--; )
	s[k*nwl+iw] = qs[iw];
  }
}

double PhotoDiss9::CalcRate(const double *flux, int k) const
{
  double res = 0.;
  for (int iw = first; iw <= last; iw++)
    res += flux[iw] * s[k*nwl+iw];
  return (res);
}

class PhotoDiss10a : public PhotoDissReact {
  const double *wl;
  double *sa, *sb, *s1, *s2, *y1, *y2;
  int nwl, nz, first, last, sid, setup;
public :
  PhotoDiss10a(AllPhotoDiss *apd) : PhotoDissReact(apd, "CH2O -> H + HCO")
    {
      sa = sb = s1 = s2 = y1 = y2 = NULL;
      setup = 0;
    }
  ~PhotoDiss10a(void);
  virtual void Setup(int nwl_in, const double *wl_in, int nz);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;
  double CalcRate2(const double *flux, int k) const;
};

PhotoDiss10a::~PhotoDiss10a(void)
{
  if (sa)  delete sa;
  if (sb)  delete sb;
  if (s1)  delete s1;
  if (s2)  delete s2;
  if (y1)  delete y1;
  if (y2)  delete y2;
}

void PhotoDiss10a::Setup(int nwl_in, const double *wl_in, int nz_in)
{
#define NWL10 23
  const float wltpl[NWL10+1] = {
    301.25, 303.75, 306.25, 308.75, 311.25, 313.75, 316.25, 318.75, 
    321.25, 323.75, 326.25, 328.75, 331.25, 333.75, 336.25, 338.75, 
    341.25, 343.75, 346.25, 348.75, 351.25, 353.75, 356.25, 358.75 };
  const float satpl[NWL10] = {
    1.37e-20, 4.43e-20, 3.27e-20, 2.24e-20, 8.82e-21, 3.47e-20, 
    3.94e-20, 1.69e-20, 1.16e-20, 4.71e-21, 4.61e-20, 2.34e-20, 
    1.31e-20, 1.14e-21, 1.3e-21, 3.54e-20, 8.98e-21, 1.31e-20, 
    5.2e-22, 3.1e-22, 2.16e-21, 1.63e-20, 9.9e-22 };
  const float sbtpl[NWL10] = {
    -2.1e-24, -4.73e-23, -1.06e-23, -7.24e-24, 2.48e-23, -3.64e-23, 
    -2.3e-23, 6.59e-24, -1.52e-23, 1.18e-24, -8.86e-23, -2.15e-23, 
    -1.53e-23, 4.32e-24, 5e-25, -8.96e-23, 1.86e-23, -2.64e-23, 
    9.57e-24, 4.38e-24, 9.48e-24, -4.05e-23, 1.27e-23 };
#define NY1WL 19
  const float y1wltpl[NY1WL] = {
    255., 260., 265., 270., 275., 280., 285., 288., 
    290., 295., 300., 305., 310., 312., 315., 320., 
    325., 330., 335.};
  const float y1tpl[NY1WL] = {
    .32, .315, .33, .38, .44, .525, .635, .685, 
    .714, .74, .75, .755, .750, .74, .715, .595, 
    .46, .310, .12 };
#define NY2WL 35
  const float y2wltpl[NY2WL] = {
    0., 255., 260., 262, 265., 270., 272., 275., 278., 
    280., 285., 290., 295., 300., 305., 310., 314., 
    315., 318., 320., 325., 330., 332., 333., 333.5, 
    334., 334.5, 335., 336., 340., 345., 350., 355., 
    360., 10000. };
  const float y2tpl[NY2WL] = {
    .49, .49, .497, .499, .495, .48, .465, .4, .35, 
    .34, .32, .295, .275, .255, .245, .26, .278, 
    .29, .35, .395, .5, .67, .76, .775, .778, 
    .78, .778, .775, .75, .635, .5, .377, .225, 
    .105, .105 };
  if (setup) return;
  wl = wl_in; nwl = nwl_in; nz = nz_in;
  sa = new double[nwl];
  sb = new double[nwl];
  s1 = new double[nwl*nz];
  s2 = new double[nwl*nz];
  y1 = new double[nwl];
  y2 = new double[nwl];
  if (!sa || !sb || !s1 || !s2 || !y1 || !y2) {
    fprintf(stderr, "Unable to allocate memory in \"PhotoDiss10a::Setup\".\n");
    exit (3);
  }
  Rebin(NWL10, wltpl, satpl, nwl, wl, sa);
  Rebin(NWL10, wltpl, sbtpl, nwl, wl, sb);
  for (first = 0; first < nwl && sa[first] == 0.; first++);
  for (last = nwl-1; last > 0 && sb[last] == 0.; last--);
  Interpolate(NY1WL, y1wltpl, y1tpl, nwl, wl, y1);
  Interpolate(NY2WL, y2wltpl, y2tpl, nwl, wl, y2);
  setup = 1;
  sid = -1;
}

void PhotoDiss10a::SetCol(int serialid, const double *tabs, const double *airc)
{
  if (serialid != sid) {
    for (int iw = first; iw <= last; iw++) {
      double wc = 0.5 * (wl[iw] + wl[iw+1]);
      for (int k = nz; k--; ) {
	double xs, t = tabs[k] - 273.15;
	if (t < -50.)  t = -50.;
	if (t > 20.)  t = 20.;
	xs = sa[iw] + t*sb[iw];
	if (xs < 0.)  xs = 0.;
	s1[k*nwl+iw] = xs * y1[iw];
	double qy;
	if (wc >= 300. && y2[iw] > 0.) {
	  double phi1, phi2, phi20, ak300, akt;
	  phi1 = y1[iw];
	  phi2 = y2[iw];
	  phi20 = 1. - phi1;
	  ak300 = ((1./phi2)-(1./phi20))/2.54E+19;
	  akt =ak300*(1.+61.69*(1.-tabs[k]/300.)*(wc/329.-1.));
	  qy = 1. / ( (1./phi20) + airc[k]*akt);
	}
	else {
	  qy = y2[iw];
	}
	s2[k*nwl+iw] = xs * qy;
      }
    }
    sid = serialid;
  }
}

double PhotoDiss10a::CalcRate(const double *flux, int k) const
{
  double res = 0.;
  for (int iw = first; iw <= last; iw++)
    res += flux[iw] * s1[k*nwl+iw];
  return res;
}

double PhotoDiss10a::CalcRate2(const double *flux, int k) const
{
  double res = 0.;
  for (int iw = first; iw <= last; iw++)
    res += flux[iw] * s2[k*nwl+iw];
  return res;
}

class PhotoDiss10b : public PhotoDissReact {
  const double *wl;
  int nwl;
  PhotoDiss10a *ph10a;
public :
  PhotoDiss10b(AllPhotoDiss *apd, PhotoDiss10a *ph10a_in) : 
    PhotoDissReact(apd, "CH2O -> H2 + CO")
    { ph10a = ph10a_in; }
  virtual void Setup(int nwl, const double *wl, int nz_in);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;
};

void PhotoDiss10b::Setup(int nwl, const double *wl, int nz)
{
  ph10a->Setup(nwl, wl, nz);
}

void PhotoDiss10b::SetCol(int serialid, const double *tabs, const double *airc)
{
  ph10a->SetCol(serialid, tabs, airc);
}

double PhotoDiss10b::CalcRate(const double *flux, int k) const
{
  return (ph10a->CalcRate2(flux, k));
}

class PhotoDiss11 : public SimplePhotoDiss {
  double *f0;
public :
  PhotoDiss11(AllPhotoDiss *apd) : SimplePhotoDiss(apd, "CH3OOH -> CH3O + NO2")
    { f0 = NULL; }
  ~PhotoDiss11(void);
  virtual void Setup(int nwl, const double *wl, int nz);
};

PhotoDiss11::~PhotoDiss11(void)
{
  if (f0)  delete f0;
}

void PhotoDiss11::Setup(int nwl, const double *wl, int nz)
{
#include "ph11data.c"  // contains the spectra of the reaction
  f0 = new double[nwl];
  if (!f0) {
    fprintf(stderr, "Error allocating memory in PhotoDiss4::Setup.\n");
    exit (3);
  }
  Interpolate(NWL11, swl11, s0_11, nwl, wl, f0);
  SimplePhotoDiss::Setup(nwl, f0, nz);
}

class PhotoDiss12 : public SimplePhotoDiss {
  double *f0;
public :
  PhotoDiss12(AllPhotoDiss *apd) : SimplePhotoDiss(apd, "CH3ONO2 -> CH3O+NO2")
    { f0 = NULL; }
  ~PhotoDiss12(void);
  virtual void Setup(int nwl, const double *wl, int nz);
};

PhotoDiss12::~PhotoDiss12(void)
{
  if (f0)  delete f0;
}

void PhotoDiss12::Setup(int nwl, const double *wl, int nz)
{
#include "ph12data.c"  // contains the spectra of the reaction
  f0 = new double[nwl];
  if (!f0) {
    fprintf(stderr, "Error allocating memory in PhotoDiss4::Setup.\n");
    exit (3);
  }
  Interpolate(NWL12, swl12, s0_12, nwl, wl, f0);
  SimplePhotoDiss::Setup(nwl, f0, nz);
}

class PhotoDiss13 : public PhotoDissReact {
  const double *wl;
  double *s0, *stpl, *qy;
  int nwl, nz;
public : 
  PhotoDiss13(AllPhotoDiss *apd) : PhotoDissReact(apd, "CH3CHO -> CH3+CHO")
    {
      s0 = stpl = qy = NULL;
    }
  ~PhotoDiss13(void);
  virtual void Setup(int nwl_in, const double *wl_in, int nz);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;  
};

PhotoDiss13::~PhotoDiss13(void)
{
  if (s0) delete s0;
  if (stpl) delete stpl;
  if (qy) delete qy;
}

void PhotoDiss13::Setup(int nwl_in, const double *wl_in, int nz_in)
{
#include "ph13data.c"
  nwl = nwl_in; wl = wl_in; nz = nz_in;
  stpl = new double[nwl];
  qy = new double[nwl];
  s0 = new double[nwl*nz];
  if (!stpl || !qy || !s0) {
    fprintf(stderr, "Unable to allocate memory in PhotoDiss13::Setup\n");
    exit (2);
  }
  Interpolate(NSWL13, swl13, s0_13, nwl, wl, stpl);
  Interpolate(NQWL13, qwl13, q13, nwl, wl, qy);
}

void PhotoDiss13::SetCol(int serialid, const double *tabs, const double *airc)
{
  for (int k = nz; k--; )
    for (int iw = nwl; iw--; ) {
      int loc = k*nwl+iw;
      s0[loc] = stpl[iw] *
	(qy[iw] ? 1. / (1. + (1./qy[iw] - 1.)*airc[k]/2.465e19) : 0);
    }
}

double PhotoDiss13::CalcRate(const double *flux, int k) const
{
  double rate = 0.;
  for (int iw = nwl; iw--; )
    rate += s0[k*nwl+iw] * flux[iw];
  return (rate);
}

class PhotoDiss14 : public PhotoDissReact {
  const double *wl;
  double *qy, *stpl;
  int nwl, nz;
public :
  PhotoDiss14(AllPhotoDiss *apd) : PhotoDissReact(apd, "CH3COCH3 ->")
    {
      qy = stpl = NULL;
    }
  ~PhotoDiss14(void);
  virtual void Setup(int nwl_in, const double *wl_in, int nz);
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
  virtual double CalcRate(const double *flux, int k) const;
};

PhotoDiss14::~PhotoDiss14(void)
{
  if (qy) delete qy;
  if (stpl) delete stpl;
}

void PhotoDiss14::Setup(int nwl_in, const double *wl_in, int nz_in)
{
#include "ph14data.c"
  nwl = nwl_in; wl = wl_in; nz = nz_in;
  stpl = new double[nwl];
  qy = new double[nz];
  if (!stpl || !qy) {
    fprintf(stderr, "Unable to allocate memory in PhotoDiss14::Setup\n");
    exit (2);
  }
  Interpolate(NSWL14, swl14, s0_14, nwl, wl, stpl);
}

void PhotoDiss14::SetCol(int serialid, const double *tabs, const double *airc)
{
  for (int k = nz; k--; )
    qy[k] = 0.0766 + 0.09415*exp(-airc[k]/3.222e18);
}

double PhotoDiss14::CalcRate(const double *flux, int k) const
{
  double rate = 0.;
  for (int iw = nwl; iw--; )
    rate += stpl[iw] * qy[k] * flux[iw];
  return rate;
}

// Instantiation of PhotoDiss-Objects

void InstantiatePhotoDiss(AllPhotoDiss &apd)
{
  static BOOL init = FALSE;
  if (!init) {
    new PhotoDiss1(&apd);
    new PhotoDiss2b(&apd, new PhotoDiss2a(&apd));
    new PhotoDiss3(&apd);
    new PhotoDiss4(&apd);
    new PhotoDiss5b(&apd, new PhotoDiss5a(&apd));
    new PhotoDiss6(&apd);
    new PhotoDiss7(&apd);
    new PhotoDiss8(&apd);
    new PhotoDiss9(&apd);
    new PhotoDiss10b(&apd, new PhotoDiss10a(&apd));
    new PhotoDiss11(&apd);
    new PhotoDiss12(&apd);
    new PhotoDiss13(&apd);
    new PhotoDiss14(&apd);
    init = TRUE;
  }
}
