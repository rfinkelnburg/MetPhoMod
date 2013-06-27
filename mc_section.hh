/*
   MODULE mc_section
   Manages input sections
*/

#ifndef INCLUDE_MC_SECTION
#define INCLUDE_MC_SECTION

#ifndef BOOL
#define BOOL int
#define TRUE 1
#define FALSE 0
#endif

#ifndef INCLUDE_MCPARSE
#include "mcparse.h"
#endif

class SectionManager {
  SectionDesc *root;
public :
  void Init(void);
  SectionDesc *AddASection(const char *name, const SectionDesc *dependson,
			   StartprocTemplate startproc,
			   EndprocTemplate endproc,
			   CallType calltype,
			   BOOL initialize = TRUE);
  void TestOnSections(void);
  SectionDesc *FindSection(const char *name) const;
  SectionDesc *FindSection(const int id) const;
};

#endif
