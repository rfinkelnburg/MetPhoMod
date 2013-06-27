/*
   IMPLEMENTATION MODULE mpdata.c
   Dieses Modul wurde direkt aus einer originalen Version von 
   Smolarkievicz und Schaer abgeleitet.
*/

void MinMax(double *min, double *max, ...);

void MpData1(double *u, double *x, int n, BOOL iord, BOOL nonos);

void NpData1(double *u, double *x, int n, BOOL iord, BOOL nonos);

void MpData2(double *ux, double *uy, double *x, int nx, int ny, int iord, BOOL nonos);

void NpData2(double *ux, double *uy, double *x, int nx, int ny, int iord, BOOL nonos);

void MpAdvect(int tinc, double *ux, double *uy, double *x, double *xp, double *wwind, double *dcm,
              double zfact, double dens, double denskp,
              int nx, int ny, int iord, BOOL nonos);
