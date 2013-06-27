/*
   DEFINITION MODULE mcfilter.c
   Dieses Modul implementiert ein raeumliches Filter. Es handelt sich um das
   auf Filter nach Pepper et al. (1979) - gemaess Pielke p. 329
*/

void ApplyFilter(double *layer, PointStatus *pstat, double delta);

void ShapiroFilter(double *x, PointStatus *pstat, int ord, double filtval);
