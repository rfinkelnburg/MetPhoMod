/*
   MODULE mc_commands
   The Command class enables the User to define additional input commands
*/

#ifndef INCLUDE_MC_COMMANDS
#define INCLUDE_MC_COMMANDS

#ifndef INCLUDE_MC_GROUP
#include "mc_group.hh"
#endif

class Ident {
  friend class Group<Ident>;
protected :
  char name[32];
public :
  Ident(const char *name_in);
  
  const char *Name(void) const { return name; }
};

class Command : public Ident {
  friend class Group<Command>;
protected :
  Command(const char *name_in) : Ident(name_in) {}

public :
  const char *Name(void) const
    { return Ident::Name(); }

  virtual short Exec(void);
  virtual short Exec(const Ident &);
  virtual short Exec(Group<Ident> &);
  // These routines are called, when the Usercommand is encountered
  // in the input-file. If one of these routines are not overloaded
  // the interpretor assumes, that the command can not be called
  // with the according number of arguments.
};

#endif
