/*
  MODULE MCAPD (AllPhotoDiss)
  Contains a group class for all the available photodissociations
*/

#ifndef INCLUDE_MCAPD
#define INCLUDE_MCAPD

#include <stddef.h>

#ifndef INCLUDE_MC_GROUP
#include "mc_group.hh"
#endif

class AllPhotoDiss;

class PhotoDissReact {
  int id;
public :
  const char *name;  // Must be global for group-template
  // creator
  PhotoDissReact(AllPhotoDiss *apd, const char *name_in);
  virtual void Setup(int nwl, const double *wl, int nz) = 0;
  virtual void SetCol(int serialid, const double *tabs, const double *airc);
    // Default is: Do nothing!
  virtual double CalcRate(const double *flux, int k) const = 0;
};

class AllPhotoDiss {
  Group<PhotoDissReact> react, used;
  int serialid;
public :
  // Constructors:
  AllPhotoDiss(void);

  // Modifiers
  int RegisterReaction(PhotoDissReact *react);
    // Registers a reaction in the PhotoDiss Database. Returns
    // the reaction ID
  int UseReaction(const int id);
    // Enters a mark, which says that this reaction is used in the
    // actual simulation. The function returns the register ID
    // of the Reaction, or 0 if reaction wth id has not been found.
  void SetupReactions(int nwl, const double *wl, int nz);
    // initializes internal data for all used reactions.
  void SetCol(const double *tabs, const double *airc);
    // Set and prepare column values.
  void CalcRates(const double *flux, int k, double *rates);
    // Calculates reaction rates for all involved reactions using
    // radiation flux "flux". The results will be written into the
    // "rates" array.
  // Accessors
  void PrintReactions(void) const;
    // Prints a description of all reactions to the screen.
  int NUsed(void) const;
    // return the number of used reactions.
};

#endif
