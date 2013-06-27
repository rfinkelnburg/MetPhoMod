

#ifndef INCLUDE_MCRELY
#define INCLUDE_MCRELY

#ifndef INCLUDE_MCPARSE
#include "mcparse.h"
#endif

int SetGridToSect(InputSection oldsec, InputSection newsec);
int SetGridToEmission(InputSection oldsec, InputSection newsec);
int TestOnSection(InputSection oldsec, InputSection newsec);
char *OnGroundInterface(RelyCmmd cmmd, VarDesc *v);
char *OnCloudWater(RelyCmmd cmmd, VarDesc *v);
char *OnCloudIce(RelyCmmd cmmd, VarDesc *v);
char *OnBorderType(RelyCmmd cmmd, VarDesc *v);
char *OnCoriolis(RelyCmmd cmmd, VarDesc *v);
char *OnAdvection(RelyCmmd cmmd, VarDesc *v);
char *OnChemistry(RelyCmmd cmmd, VarDesc *v);
char *OnKEpsTurb(RelyCmmd cmmd, VarDesc *v);

#endif
