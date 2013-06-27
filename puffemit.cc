/*
   MODULE puffemit
   Allows Puff-Emissions to meteochem.
   This Module is at the same time a test case and an example to the new
   standard way to add modules to MetPhoMod
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcglobal.h"
#include "mcground.h"
#include "mcparse.h"
#include "mc_section.hh"
#include "mc_group.hh"
#include "mc_group.t"
#include "mc_commands.hh"
#include "mc_variable.hh"
#include "mc_module.hh"


class Puff {
  Group<VarDesc> var;
  VarDesc *grvar;
  int cx, cy, size, time;
  short triggered;
  double conz, *data;
public :
  const char *name;
  Puff(const Ident &ident, InputSection section);
  void Emit(long actime);
  const char *Name() { return name; }
};

class PuffEmit : public Module, public Command {
  Group<Puff> puff;
  int mytime;
  InputSection mysection;
public :
  // Constructors
  PuffEmit(void);
  
  // These functions are required by the Module-Base-Type
  virtual void Init1(void);
  virtual long DoCalc(long mintime, long maxtime = 0);
  virtual long Time(void)  const;

  virtual short Exec(Ident &);
  virtual short Exec(Group<Ident> &);
  long DoCalc(void);
#ifdef PARALLEL
  virtual char ReadyForParallel(void) const;
#endif
};

double Concentration(int k, int i, int j, VarDesc *v)
{
  return (g[v->userdata][i*row+j+k*layer] * density[k]);
}

Puff::Puff(const Ident &ident, InputSection section)
{
  VarDesc *v;
  v = McInterface::mcvar.FindVariable(ident.Name(), CHEMISTRY);
  if (!v) {
    InputError(ERROR, "species %s was not found.", ident.Name());
    inputerror = 1;
    return;
  }
  data = g[v->v.et];
  var.Register(new VarDesc("cx", "x-coord of center of puffdomain",
			   "cells", INT_PTR, &cx, 0,
			   section, MUST_BE_SET, 0, 1, 1000,
			   INT_NUM, 0, NULL));
  var.Register(new VarDesc("cy", "y-coord of center of puffdomain",
			   "cells", INT_PTR, &cy, 0,
			   section, MUST_BE_SET, 0, 1, 1000,
			   INT_NUM, 0, NULL));
  var.Register(new VarDesc("size", "size of puffdomain",
			   "cells", INT_PTR, &size, 0,
			   section, MUST_BE_SET, 0, 1, 50,
			   INT_NUM, 0, NULL));
  var.Register(new VarDesc("time", "Purr release time", "sec",
			   INT_PTR, &time, 0, InputSection(section),
			   MUST_BE_SET, 0, 0, 8000000, INT_NUM, 0, NULL));
  var.Register(new VarDesc("conz", "conzentration to establish in puffdomain",
			   "microgram/m^3", DOUBLE_PTR, &conz, 0,
			   section, MUST_BE_SET, 0, 0, 10000,
			   INT_NUM, 0, NULL));
  grvar = McInterface::mcvar.AddNewVariable(
             new VarDesc((char *)ident.Name(), "Puff-emission group variable",
			 NULL, GROUP_VAL, (void *)&var, 0, section,
			 MUST_BE_SET,
			 0, 0, 0, VAR_GROUP, 0, NULL, v->v.et));
  name = grvar->name;
  {
    char newname[32];
    sprintf(newname, "%s_C", name);
    McInterface::mcvar.AddNewVariable(
      new VarDesc(strdup(newname), "Absolute concentration of species",
		  "microgram/m^3", PROC_VAL, (void *)Concentration,
		  ALL_DIM, section, CALCULATED_VALUES, 0, 0, 0,
		  NORMAL_NUM, 0, NULL, v->v.et));
    triggered = 0;
  }
}

void Puff::Emit(long actime)
{
  int fa, loc, k;
  if (!triggered && actime >= time) {
    loc = cx*row + cy;
    k = fa = ground[loc].firstabove;
    //    for (k = fa; k < fa+3; k++) {
    data[loc + k*layer] = conz / density[k];
    //}
    triggered = 1;
  }
}

PuffEmit::PuffEmit(void) :
  Module(&McInterface::modmanager, "PUFFEMIT", Module::DIAGNOSTIC,
	 "Puff Emission Module", "well tested, not parallelized"),
  Command("PUFF")
{
}

void PuffEmit::Init1(void)  // Not needed
{
  // Establish a section
  mysection = InputSection(
    McInterface::mcsect.AddASection(
      "PUFFEMIT", McInterface::mcsect.FindSection(CHEMISTRY),
      NULL, NULL, CAN_BE_OMITTED)->id);
  // Add a command
  McInterface::mccommand.Register(this);

  mytime = 0;
}

long PuffEmit::DoCalc(long mintime, long maxtime)
{
  mytime = mintime;
  GroupEnumerator<Puff> p(puff);
  for ( ; !p; ++p)
    p()->Emit(mytime);
  return mytime;
}

long PuffEmit::Time(void) const
{
  return mytime;
}

short PuffEmit::Exec(Ident &ident)
{
  puff.Register(new Puff(ident, mysection));
  return (0);
}

short PuffEmit::Exec(Group<Ident> &gi)
{
  for (GroupEnumerator<Ident> e(gi); !e; ++e)
    puff.Register(new Puff(*e(), mysection));
  return (0);
}

#ifdef PARALLEL
char PuffEmit::ReadyForParallel(void) const
{
  return (!puff.N());
}
#endif

PuffEmit puffemit;
// Instantiate one module instance

/*
template class Group<Puff>;
template class GroupEnumerator<Puff>;
*/
