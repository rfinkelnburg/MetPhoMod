/*
   MODULE mc_commands
   The Command class enables the User to define additional input commands
*/


#include "mc_commands.hh"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include "mc_group.t"
#include "mcparse.h"

Ident::Ident(const char *name_in)
{
  strncpy(name, name_in, 31);
}

static char wrongarg[] = "Command \"%s\" can not be called with this number of arguments.";

short Command::Exec(void)
{
  InputError(ERROR, wrongarg, Name());
  return 1;
}

short Command::Exec(const Ident &i)
{
  InputError(ERROR, wrongarg, Name());
  return 1;
}

short Command::Exec(Group<Ident> &i)
{
  InputError(ERROR, wrongarg, Name());
  return 1;
}

/*
template class Group<Ident>;
template class Group<Command>;
template class GroupEnumerator<Ident>;
template class GroupEnumerator<Command>;
template class Holder<Ident>;
template class Holder<Command>;
*/
