/*
   DEFINITION MODULE mccdfin.h
   Modul zum Einlesen von Daten aus einem CDF-File
*/

#define MCCDFIN

typedef enum
   {WALL_VAR, XWALL_VAR, LAYER_VAR, PROFILE_VAR, MESH_VAR, COORD_VAR, GROUND_VAR}  BorderVarType;

int ReadCDFFile(char *name, InputSection section);

void ActualiseValuesInTime(long actime);

#ifdef PARALLEL

void SendBorderTime(BOOL leftest, BOOL rightest);

void RecvBorderTime(void);

#endif
