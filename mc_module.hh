/*
   MODULE mc_module
   Implements a standardized way of adding modules to the MetPhoMod
   Program. It includes two major parts:

     1) The module managing class: This class manages modules. It keeps
        a list of active modules, and it includes a timer, which
	automaticly calls the modules then necessary.
     2) An abstract template module. All user implemented modules should
        inherit from this one.
*/

#ifndef INCLUDE_MODULE
#define INCLUDE_MODULE

#ifndef INCLUDE_MC_VARIABLE
#include "mc_variable.hh"
#endif

#ifndef INCLUDE_MC_SECTION
#include "mc_section.hh"
#endif

#ifndef INCLUDE_MC_COMMANDS
#include "mc_commands.hh"
#endif

class ModuleManager;

class McInterface {
  static short initialize;
public :
  McInterface();
  static ModuleManager modmanager;
  static VariableManager mcvar;
  static SectionManager mcsect;
  static Group<Command> mccommand;
};

class Module {
public :
  typedef enum {DIAGNOSTIC, PROGNOSTIC}  ModType;
  McInterface interface;
private :
  ModuleManager *mymanager;
  Module *next;
  ModType modtype;
  const char *name, *comment, *implementation_status;
  friend class ModuleManager;
protected :
  // Standard constructor.
  Module(ModuleManager *mm, char *modulename, ModType modtype_in,
	 const char *comment_in, const char *istatus);
  // Registers the actual module automatically at mm.
public :
  ~Module(void);

  // Modifiers
  virtual void Init1(void);
  virtual void Init2(void);
  virtual void EndOfSection(SectionDesc *section);
  // These three commands can be overloaded, to provide module
  // initialisation routines. The three routines are called in a
  // specific initialisation situation.
  //
  // Init1: Before reading the input file
  //    This phase is normally used, to define new sections, variables
  //    and commands.
  // Init2: After reading the input file, before starting the simulation
  //    Normally used for allocation of multidimensional arrays, etc.
  // EndOfSection: Called after reading of section "section" is completed.
#ifdef PARALLEL
  virtual char ReadyForParallel(void) const = 0;
  // True if the Module is ready for a parallel run. This function
  // must be overloaded by derived classes. When a module is not
  // active in a calculation, the answer may be yes even if no
  // implementation exists.
  virtual void SendToWorkers(void);
  virtual void GetFromMaster(void);
  // This function pair allows to initialize the data in parallel mode.
  // The master process calls "SendToWorkers", the workers call
  // "GetFromMaster". Both procedures are called *before* Init2.
  // These Routines do nothing by default.
#endif

  virtual long DoCalc(long mintime, long maxtime = 0) = 0;
  // Asks the Module to do a time iteration step. If maxtimest is
  // unequal to zero, the module should not do a timestep, which carries
  // the module time further than maxtime. mintime should be reached
  // in any case. This routine must be overloaded!

  // Accessors
  virtual long Time(void)  const = 0;
  // Returns the time of the module
  ModType Type(void)  { return modtype; }
  // Is this a diagnostic or a prognostic module?
  //
  // Prognostic modules perform an integration in time.
  // Diagnostic modules are not really time dependent. They calculate
  // their result statically, based on the state of other variables.
  
  const char *Name()
  {
    return name;
  }
};

class ModuleManager {
private :
  Module *firstprog, *lastprog;
  Module *firstdiag, *lastdiag;
  long actime;

  void RegisterModule(Module *m, char *modulename);
  void DeleteModule(Module *m);
  friend class Module;
public :
  // Constructors
  void Init(void);
  // Creates a module manager

  // Modifiers
  void Init1(void);
  void Init2(void);
  void EndOfSection(SectionDesc *section);
  // Initializes all registered Modules
#ifdef PARALLEL
  void CheckIfReadyForParallel(void) const;
  void SendToWorkers(void);
  void GetFromMaster(void);
#endif

  long DoCalc(long mintime, long maxtime = 0);
  // Asks the ModuleManager to do a certain time integration step.
  // The ModuleManager will do a time step between mintimestep and
  // maxtimestep. If maxtimestep is not given, the time step can
  // have any length.

  // Accessors
  int Time(void)  const { return actime; }
  // Returns the time of the modul manager
};

#endif
