/*
   DEFINITION MODULE mchem.h
   Strukturdefinitionen fuer die Chemieberechnungen.
*/

extern int allfast;
extern ProductionDesc *prodd;

double AirConz(int k, int i, int j);

void SetSubstKonz(char *name, double conz);

void NoteTransport(BOOL transport);

#ifdef MCGROUND

void ChemicalTimeStep(long dt, Vector *sunpos);

#endif
