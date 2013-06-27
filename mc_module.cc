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

#include "mc_module.hh"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include "mcparse.h"

#include "mc_group.t"

// Implementation of class Module

Module::Module(ModuleManager *mm, char *modulename, ModType modtype_in,
	       const char *comment_in, const char *istatus)
// Registers the actual module automaticly at mm.
{
  mm->RegisterModule(this, modulename);
  mymanager = mm;
  modtype = modtype_in;
  name = modulename;
  comment = comment_in;
  implementation_status = istatus;
}

Module::~Module(void)
{
  mymanager->DeleteModule(this);
}

// Dummy routines for virtual initialisation routines are provided.
void Module::Init1(void)
{}

void Module::Init2(void)
{}

void Module::EndOfSection(SectionDesc *section)
{}

#ifdef PARALLEL

void Module::SendToWorkers(void)
{}

void Module::GetFromMaster(void)
{}

#endif

// Implementation of class ModuleManager

void ModuleManager::Init(void)
{
  firstprog = lastprog = firstdiag = lastdiag = NULL;
}

void ModuleManager::RegisterModule(Module *m, char *modulename)
{
  m->next = NULL;
  if (m->Type() == Module::PROGNOSTIC) {
    if (lastprog) {
      lastprog->next = m;
      lastprog = m;
    }
    else
      lastprog = firstprog = m;
  }
  else {
    if (lastdiag) {
      lastdiag->next = m;
      lastdiag = m;
    }
    else
      lastdiag = firstdiag = m;
  }
}

void ModuleManager::DeleteModule(Module *m)
{
  Module *i;
  if (firstprog) {
    if (firstprog == m) {
      firstprog = m->next;
      if (!firstprog)  lastprog = NULL;
      return;
    }
    for (i = firstprog; i->next && i->next != m; i = i->next);
    if (i) {
      i->next = m->next;
      if (!i->next)  lastprog = i;
      return;
    }
  }
  if (firstdiag) {
    if (firstdiag == m) {
      firstdiag = m->next;
      if (!firstdiag)  firstdiag = NULL;
      return;
    }
    for (i = firstdiag; i->next && i->next != m; i = i->next);
    if (i) {
      i->next = m->next;
      if (!i->next)  lastdiag = i;
      return;
    }
  }
}

void ModuleManager::Init1(void)
{
  Module *m;
  for (m = firstprog; m; m = m->next)
    m->Init1();
  for (m = firstdiag; m; m = m->next)
    m->Init1();
}

void ModuleManager::Init2(void)
{
  Module *m;
  for (m = firstprog; m; m = m->next)
    m->Init2();
  for (m = firstdiag; m; m = m->next)
    m->Init2();
}

void ModuleManager::EndOfSection(SectionDesc *section)
{
  Module *m;
  for (m = firstprog; m; m = m->next)
    m->EndOfSection(section);
  for (m = firstdiag; m; m = m->next)
    m->EndOfSection(section);
}

#ifdef PARALLEL

void ModuleManager::CheckIfReadyForParallel(void) const
{
  Module *m;
  char ready = 1;
  for (m = firstprog; m; m = m->next)
    if (!m->ReadyForParallel()) {
      InputError(ERROR, "Module \"%s\" can not be run in parallel mode. (not implemented)",
		 m->Name());
      ready = 0;
    }
  for (m = firstdiag; m; m = m->next)
    if (!m->ReadyForParallel()) {
      InputError(ERROR, "Module \"%s\" can not be run in parallel mode. (not implemented)",
		 m->Name());
      ready = 0;
    }
  if (!ready)  exit (2);
}  

void ModuleManager::SendToWorkers(void)
{
  Module *m;
  for (m = firstprog; m; m = m->next) {
    printf("  %s\n", m->Name());
    m->SendToWorkers();
  }
  for (m = firstdiag; m; m = m->next) {
    printf("  %s\n", m->Name());
    m->SendToWorkers();
  }
}

void ModuleManager::GetFromMaster(void)
{
  Module *m;
  for (m = firstprog; m; m = m->next)
    m->GetFromMaster();
  for (m = firstdiag; m; m = m->next)
    m->GetFromMaster();
}

#endif

long ModuleManager::DoCalc(long mintime, long maxtime)
{
  long firstime;
  Module *m, *firstm;

  // Call diagnostic modules
  for (Module *dm = firstdiag; dm; dm = dm->next)
    dm->DoCalc(mintime);

  // Call prognostic modules, if any
  if (firstprog) {
    do {
      // Search the module with the backmost time
      firstime = firstprog->Time(); firstm = firstprog;
      for (m = firstprog; m; m = m->next) {
	if (m->Time() < firstime)  {
	  firstime = m->Time();
	  firstm = m;
	}
      }

      // If another integration step is necessary, do it!
      if (firstime < mintime) {
	firstm->DoCalc(mintime, maxtime);
      }
    }  while (firstime < mintime);
    return (actime = firstime);
  }
  else
    return (actime = mintime);
}

short McInterface::initialize = 1;
SectionManager McInterface::mcsect;
VariableManager McInterface::mcvar;
ModuleManager McInterface::modmanager;
Group<Command> McInterface::mccommand;

McInterface::McInterface(void)
{
  if (initialize) {
    mcsect.Init();
    mcvar.Init(variable);
    modmanager.Init();
    initialize = 0;
  }
}
