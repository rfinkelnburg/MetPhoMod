/*
  MODULE MCAPD (AllPhotoDiss)
  Contains a group class for all the available photodissociations
*/

#include "mcapd.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "../mc_group.t"

PhotoDissReact::PhotoDissReact(AllPhotoDiss *apd, const char *name_in)
{
  name = name_in;
  id = apd->RegisterReaction(this);
}

void PhotoDissReact::SetCol(int serialid, 
			    const double *tabs, const double *airc)
{}

AllPhotoDiss::AllPhotoDiss(void)
{
  serialid = 0;
}

int AllPhotoDiss::RegisterReaction(PhotoDissReact *react_in)
{
  return (react.Register(react_in));
}

int AllPhotoDiss::UseReaction(const int id)
{
  PhotoDissReact *r;
  if (r = react.Find(id)) {
    return used.Register(r);
  }
  else
    return 0;
}

void AllPhotoDiss::SetupReactions(int nwl, const double *wl, int nz)
{
  GroupEnumerator<PhotoDissReact> e(used);
  for ( ; !e; ++e)
    e()->Setup(nwl, wl, nz);
}

void AllPhotoDiss::SetCol(const double *tabs, const double *airc)
{
  GroupEnumerator<PhotoDissReact> e(used);
  for ( ; !e; ++e)
    e()->SetCol(serialid++, tabs, airc);  
}

void AllPhotoDiss::CalcRates(const double *flux, int k, double *rates)
{
  GroupEnumerator<PhotoDissReact> e(used);
  for ( ; !e; ++e) {
    *(++rates) = e()->CalcRate(flux, k);
//    if (!k) printf("%18s = %.3le\n", e()->name, *rates);
  }
}

void AllPhotoDiss::PrintReactions(void) const
{
  GroupEnumerator<PhotoDissReact> e(react);
  printf("Available photodissociation reactions:\n\n"
	 "--ID-- -REACTION------\n");
  for ( ; !e; ++e)
    printf("%6d %s\n", e.ID(), e()->name);
}

int AllPhotoDiss::NUsed(void) const
{
  return used.N();
}
