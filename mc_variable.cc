/*
   MODULE mc_variable
   A for the management of the MC-Variable data base.
   mc_variable is a wrapper module which enables to manage the
   C-structures in a C++ style
*/

#include "mc_variable.hh"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "mcground.h"
#include "mc_group.hh"
#include "mcparse.h"

void VariableManager::Init(VarDesc *variable)
{
  root = variable;
}

VarDesc *VariableManager::AddNewVariable(VarDesc *newvar)
{
  VarDesc *last;
  for (last = variable; last->next; last = last->next);
  last->next = newvar;
  newvar->next = NULL;
  newvar->id = NewID();
  return (newvar);
}

Variable *VariableManager::FindVariable(const char *name, InputSection section)
{
  VarDesc *v;
  for (v = variable; v && (strcmp(v->name, name) || (section >= 0 && section != v->origsection)); v = v->next);
  return (v);
}

Variable *VariableManager::FindVariable(const int id)
{
  VarDesc *v;
  for (v = variable; v && v->id != id; v = v->next);
  return (v);
}

int VariableManager::NewID(void)
{
  int id = 0;
  VarDesc *v;
  for (v = variable; v && v->id != id; v = v->next)
    if (v->id > id)  id = v->id;
  return (id+1);
}

#include "mc_group.t"

// template class Group<VarDesc>;
