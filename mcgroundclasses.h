/*
   MODULE McGroundClasses
   Dieses Modul ermöglicht eine Behandlung der Bodenparameter mit Hilfe
   von Klassen.
*/

int InitializeGroundClasses(void);

int AddClassVar(char *name);

int AddClassIdx(int idx);

int AddClassVal(double val);

int EndOfClassLine(void);

void FillGroundValues(void);
