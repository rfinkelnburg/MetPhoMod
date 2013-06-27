/*
   DEFINITION MODULE mchem.h
   Strukturdefinitionen fuer die Chemieberechnungen.
*/

extern int allfast;
extern ProductionDesc *prodd;

double AirConz(int k, int i, int j, struct VarDesc *v);

void SetSubstKonz(char *name, double conz);

void NoteTransport(BOOL transport);

#ifdef MCGROUND

void BoxChemStep(int i, int j, int k, long dt,
		 const Vector *sunpos, double *photorate);

void ChemicalTimeStep(long dt, const Vector *sunpos);

#endif
