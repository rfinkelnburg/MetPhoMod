/*

   This file is a part of

   MetPhoMod v2.0, the comprehensive three dimensional, prognostic
		   mesoscale atmospheric summer smog model.

   Copyright (C) 1996-1999, Silvan Perego


   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, refer to the WEB page

   http://www.gnu.org/copyleft/gpl.html

*/

/*
   DEFINITION MODULE SunPos
   Ein Modul zur Berechnung des Sonnenstandes aufgrund der Uhrzeit und dem
   Tag des Jahres.
   Originaler Autor: W.Eugster
   Adaptiert fuer meteochem
*/

extern double S0;    /* Solarkonstante */

typedef struct  {
  double x, y, z;
}  Vector;

double Declination(int day);
/* Berechnet die Deklination der Sonne aufgrund des Julianischen Datums
   innerhalb des Jahres. */

double SolarTime(int day, double localtime, double longitude,
                 double timezonediff);
/* Berechnet die die Sonnenzeit. */

double HourAngle(int day, double localtime, double longitude,
                 double timezonediff);

void SunVector(Vector *s, int day,
               double localtime, double longitude, double latitude,
               double timezonediff);

double AngleBetween(Vector *s, Vector *t);
/* Berechnet den Winkel zwischen den zwei Vektoren s und t */

double CosZOfTheSun(int day, double localtime, double longitude,
                    double latitude, double timezonediff);
/* Berechnet die Hoehe der Sonne über Boden. */

int JulianDay(int day, int month, int year);

double RatioDistMeanDist(int day);
