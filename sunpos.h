/*
   DEFINITION MODULE SunPos
   Ein Modul zur Berechnung des Sonnenstandes aufgrund der Uhrzeit und dem
   Tag des Jahres.
   Originaler Autor: W.Eugster
   Adaptiert fuer meteochem
*/

#ifndef INCLUDE_SUNPOS
#define INCLUDE_SUNPOS

extern double S0;    /* Solarkonstante */

struct Vector  {
  double x, y, z;
};

void SunVector(Vector *s, int day,
               double localtime, double longitude, double latitude,
               double timezonediff);

double AngleBetween(const Vector *s, const Vector *t);
/* Berechnet den Winkel zwischen den zwei Vektoren s und t */

double CosZOfTheSun(int day, double localtime, double longitude,
                    double latitude, double timezonediff);
/* Berechnet die Hoehe der Sonne über Boden. */

int JulianDay(int day, int month, int year);

double RatioDistMeanDist(int day);

#endif
