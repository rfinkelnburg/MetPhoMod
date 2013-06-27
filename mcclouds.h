/*
   DEFINITION MODULE mcclouds.c
   Implementiert Wolkenphysik
*/


double TotalClouds(int k, int i, int j, struct VarDesc *v);

double TotalWaterClouds(int k, int i, int j, struct VarDesc *v);

double TotalIceClouds(int k, int i, int j, struct VarDesc *v);

double TerminalVelocityOfRain(int k, int i, int j, struct VarDesc *v);

void CloudPhysics(void);
