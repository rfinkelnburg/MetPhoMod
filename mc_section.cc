/*
   MODULE mc_section
   Manages input sections
*/

#include "mc_section.hh"
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "mcglobal.h"
#include "mcrely.h"
#include "mcnest.h"

#define MAXSECTION 18

static SectionDesc section[MAXSECTION] =
{"TIME", TIME, NULL, TRUE, FALSE, NULL, EndOfTime, MUST_BE_CALLED, NULL,
 "GRID",  GRID, &section[TIME], TRUE, FALSE, NULL, EndOfGrid, MUST_BE_CALLED, NULL,
 "OPTIONS", OPTIONS, &section[GRID], TRUE, FALSE, NULL, EndOfOptions, MUST_BE_CALLED, NULL,
 "ENVIRONMENT", ENVIRONMENT, &section[OPTIONS], TRUE, FALSE, NULL, EndOfEnvironment, MUST_BE_CALLED, NULL,
 "INITIAL-DATA", INITIAL_DATA, &section[GRID], TRUE, FALSE, NULL, NULL, MUST_BE_CALLED, NULL,
 "BORDER-DATA", BORDER_DATA, &section[INITIAL_DATA], TRUE, FALSE, NULL, NULL, MUST_BE_CALLED, NULL,
 "GROUND", GROUND_BORDER, &section[ENVIRONMENT], TRUE, FALSE, TestOnSection, NULL, MUST_BE_CALLED, NULL,
 "TOP", TOP_BORDER, &section[ENVIRONMENT], TRUE, FALSE, TestOnSection, NULL, MUST_BE_CALLED, NULL,
 "NORTH", NORTH_BORDER, &section[BORDER_DATA], TRUE, FALSE, TestOnSection, NULL, MUST_BE_CALLED, NULL,
 "SOUTH", SOUTH_BORDER, &section[BORDER_DATA], TRUE, FALSE, TestOnSection, NULL, MUST_BE_CALLED, NULL,
 "WEST", WEST_BORDER, &section[BORDER_DATA], TRUE, FALSE, TestOnSection, NULL, MUST_BE_CALLED, NULL,
 "EAST", EAST_BORDER, &section[BORDER_DATA], TRUE, FALSE, TestOnSection, NULL, MUST_BE_CALLED, NULL,
 "CHEMISTRY", CHEMISTRY, &section[INITIAL_DATA], TRUE, FALSE, NULL, NULL, CALL_WHEN_SUBSTANCES, NULL,
 "EMISSIONS", EMISSIONS, &section[CHEMISTRY], FALSE, FALSE, SetGridToEmission, NULL, CAN_BE_OMITTED, NULL,
 "REDUCTION", REDUCTION, &section[EMISSIONS], FALSE, FALSE, NULL, NULL, CAN_BE_OMITTED, NULL,
 "DEPOSITION", DEPOSITION, &section[CHEMISTRY], FALSE, FALSE, SetGridToEmission, NULL, CAN_BE_OMITTED, NULL,
 "NESTING", NESTING, &section[GRID], FALSE, FALSE, NULL, EndOfNesting, CAN_BE_OMITTED, NULL,
 "OUTPUT", OUTPUT, &section[GRID], FALSE, FALSE, NULL, NULL, CAN_BE_OMITTED, NULL};
 

void SectionManager::Init(void)
{
  int i;
  for (i = 0; i < MAXSECTION-1; i++) {
    section[i].next = &section[i+1];
  }
}

SectionDesc *SectionManager::AddASection(const char *name, const SectionDesc *dependson,
					 StartprocTemplate startproc,
					 EndprocTemplate endproc,
					 CallType calltype,
					 BOOL initialize)
{
  SectionDesc *ns, *ls;
  int maxid = 0;
  ns = new SectionDesc;
  ns->name = name;
  ns->dependson = dependson;
  ns->startproc = startproc;
  ns->endproc = endproc;
  ns->calltype = calltype;
  for (ls = section; ls->next; ls = ls->next) {
    if (ls->id > maxid)  maxid = ls->id;
  }
  if (ls->id > maxid)  maxid = ls->id;
  ns->initialize = initialize;
  ns->found = FALSE;
  ls->next = ns;
  ns->id = maxid+1;
  ns->next = NULL;
  return ns;
}

void SectionManager::TestOnSections(void)
{
  SectionDesc *s;
  for (s = section; s; s = s->next)
    if (s->calltype == MUST_BE_CALLED)  {
      fprintf(stderr, "ERROR: Unable to find section %s.\n", s->name);
      inputerror = 1;
    }
    else if (s->calltype == CALL_WHEN_SUBSTANCES && nsubs > 0)  {
      fprintf(stderr, "ERROR: Unable to find section %s, although there are %i Substances defined\n",
         s->name, nsubs);
      inputerror = 1;
    }
  if (section[CHEMISTRY].calltype == WAS_CALLED && tchem % tinc)  {
    fprintf(stderr, "ERROR: tchem must be a multiple of tinc.");
    inputerror = 1;
  }
}

SectionDesc *SectionManager::FindSection(const char *name) const
{
  SectionDesc *s;
  for (s = section; s && strcmp(s->name, name); s = s->next);
  return (s);
}

SectionDesc *SectionManager::FindSection(const int id) const
{
  SectionDesc *s;
  for (s = section; s && s->id != id; s = s->next);
  return (s);
}
