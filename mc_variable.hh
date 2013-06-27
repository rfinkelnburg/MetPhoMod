/*
   MODULE mc_variable
   A for the management of the MC-Variable data base.
   mc_variable is a wrapper module which enables to manage the
   C-structures in a C++ style
*/


#ifndef INCLUDE_MC_VARIABLE
#define INCLUDE_MC_VARIABLE

#ifndef INCLUDE_MCGLOBAL
#include "mcglobal.h"
#endif

#ifndef INCLUDE_MCPARSE
#include "mcparse.h"      // Contains the definition of VarDesc, and variable
#endif

typedef VarDesc Variable;  // This could change to something else some day...

class VariableManager {
private :
  VarDesc *root;
public :
  void Init(VarDesc *variable);
  // Initializes the variable manager
  Variable *AddNewVariable(VarDesc *newvar);
  // Establishes and allocates a new variable
  Variable *FindVariable(const char *name, InputSection section = InputSection(-1));
  Variable *FindVariable(const int id);
  int NewID(void);
  // returns a new, unused variable-ID
};

#endif
