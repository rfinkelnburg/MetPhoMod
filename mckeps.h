/*
   DEFINITION MODULE McKEps
   Implementiert Boundary-Layer Turbulenz nach dem k-epsilon Schema.
*/

void CalcAllKm(void);
/* Berechnet die K-Faktoren aufgrund von TKE und EPS */

void CalcKEpsilon(long tinc);
/* Berechnet die neuen K-Parameter, sowie neue Werte fuer TKE, EPS
   und THV aufgrund der K-Eps Parametrisierung. */

void CalcKTurb(long tinc, long chemtinc);
/* Rechnet die vertikale Turbulenz aufgrund der vorhandenen K-Parameter. */
