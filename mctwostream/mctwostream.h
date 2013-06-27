/*
   MODULE mc2stream
   A general two-stream solver. This solver was directly derived from
   "ps2str.f" in the Madronich-Paciage
*/

#ifndef INCLUDE_MCTWOSTREAM
#define INCLUDE_MCTWOSTREAM

#ifndef INCLUDE_MCGROUND
#include "mcground.h"
#endif

#ifndef INCLUDE_MC_MODULE
#include "mc_module.hh"
#endif

#ifndef INCLUDE_MC_SECTION
#include "mc_section.hh"
#endif

#ifndef INCLUDE_MCPHOTODISS
#include "mcphotodiss.h"
#endif

class TwoStreamModule : public Module {
  long mtime, dtime;
  double *wl, mean_albedo, *photorate, *wc;
  int nwl;
  AllPhotoDiss apd;
  SectionDesc *sectdesc;
public :
  TwoStreamModule(void);
  virtual void Init1(void); 
  virtual void EndOfSection(SectionDesc *section);
  virtual void Init2(void);
  virtual long DoCalc(long mintime, long);
#ifdef PARALLEL
  virtual char ReadyForParallel(void) const;
  virtual void SendToWorkers(void);
  virtual void GetFromMaster(void);  
#endif
  virtual long Time(void) const;
  void CalcTwoStream(const Vector &sunpos, const double sundistance);
  void PrintPhotoDissReact(void) const;
     // Prints all available photodissoziation reactions
  int RegisterReaction(int id);
    // Registers reaction #id. If successful, this function returns
    // an pos number > 0. Otherwise it returns 0.
};

double *UpdateJValues(int nreact, int *rid);

void CalcTwoStream(const Vector &sunpos, const double sundistance);

extern TwoStreamModule twostreammodule;

#endif
