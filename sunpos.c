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
   IMPLEMENTATION MODULE Sonne
*/

/* Original:	sonne.pas,	*/
/* 		Version 2.04 	*/

/* Autor:	Werner Eugster	*/
/* Version:	1.03		*/
/* Date:	27. August 1992	*/

/* compile with -lm switch	*/

/* Input needed:					*/
/*	geogr. longitude in degrees			*/
/*	geogr. latitude in degrees			*/
/*	Julian Day of the Year				*/
/*	hour of day (decimal)				*/
/*	difference to the actual timezone in hours	*/
/*							*/
/* The difference to the actual timezone in hours is	*/
/* defined as follows:					*/
/* if the position given by the geogr. longitude and	*/
/* latitude lies within the MET timezone, the		*/
/* difference is 0.0 hours as long as your time in	*/
/* hour of day is in MET. But if you calculate in	*/
/* MEST (MET Summertime), the difference will be	*/
/* +1.0 hour.						*/


#include <stdio.h>
#include <math.h>
#include "sunpos.h"

#define		TRUE 1
#define		FALSE !TRUE

extern double Pi;
double S0 = 1376.;


extern double sqr(double x);

double Declination(int day)
/* day: Julian Day of the year		*/
{
  double de;
  de = (-23.4683 * cos(((0.9856 * (day)) + 9.3)*Pi/180.));
  de *= (Pi/180.);	/* Converting degrees to radians	*/
  return (de);
}

double EquationOfTime(int day)
/* day: Julian Day of the year		*/
{
  double thetanull, zg, argu;
  zg = (0.018-0.413*cos(0.017167173*day)+3.4925*cos(0.034334345*day) +
    0.1019*cos(0.051501519*day)+7.3057*sin(0.017167173*day) +
    9.3014*sin(0.034334345*day)+0.3256*sin(0.051501519*day)) * Pi / 180.;
  return(zg);
}

double SolarTime(int day, double localtime, double longitude, double timezonediff)
/* day: 	Julian Day of the year		*/
/* localtime:	Local time in hours		*/
/* longitude:	Longitude in degrees		*/
/* timezonediff:Difference of time zone in hours*/
{
  double zg, st;
  zg = EquationOfTime(day);
  zg *= 12. / Pi;
  st = localtime - timezonediff + zg + longitude/15.;
  return(st);	/* Solar Time in UTC	*/
}

double HourAngle(int day, double localtime, double longitude, double timezonediff)
/* day: 	Julian Day of the year		*/
/* localtime:	Local time in hours		*/
/* longitude:	Longitude in degrees		*/
/* timezonediff:Difference of time zone in hours*/
{
  double sw;
  sw = (SolarTime(day, localtime, longitude, timezonediff) - 12.) * 15.;
  sw *= (Pi/180.);	/* Converting degrees to radians	*/
			/* The Hour Angle of the Sun is		*/
			/* relative to the South Direction!	*/
  return(sw);
}

void SunVector(Vector *s, int day, double localtime,
               double longitude, double latitude, double timezonediff)
/* s:	resulting Sun Unity Vector	*/
/* day:	Julian Day of the year		*/
{
  double ds, t;
  double lat;
  ds = Declination(day);
  t = HourAngle(day,localtime,longitude,timezonediff);
  lat = latitude * Pi/180.;
  /* local vector to the sun	*/
  s->x = -cos(ds) * sin(t);
  s->y = -cos(ds) * cos(t) * sin(lat) + sin(ds) * cos(lat);
  s->z = cos(ds) * cos(t) * cos(lat) + sin(ds) * sin(lat);
}

void TopoVector(Vector *s, double exposition, double slopeangle)
/* s:		resulting Perpendicular Unity Vector	*/
/* exposition:	Angle of the exposition rel. to North*/
/* slopeangle:	Angle of the slope in radians	*/
{
  s->x = -sin(exposition) * sin(slopeangle);
  s->y = -cos(exposition) * sin(slopeangle);
  s->z = cos(slopeangle);
}

double AngleBetween(Vector *s, Vector *t)
{
  double zenithangle;
  zenithangle = (s->x*t->x + s->y*t->y + s->z*t->z) /
                sqrt((s->x*s->x + s->y*s->y + s->z*s->z) *
                     (t->x*t->x + t->y*t->y + t->z*t->z));
  /* this is the cosine of the angle...	*/
  if (zenithangle < 0.)  return(-1.0);
/*  zenithangle = acos(zenithangle); */
  return (zenithangle);
}


double CosZOfTheSun(int day, double localtime, double longitude,
                    double latitude, double timezonediff)
{
  Vector s;
  double elev;
  SunVector(&s, day, localtime, longitude, latitude, timezonediff);
  elev = s.z / sqrt(sqr(s.x) + sqr(s.y) + sqr(s.z));
		/* this is now the cosine of the zenit angle	*/
  if (elev < 0.)  elev = 0.;
  return (elev);
}

int JulianDay(int day, int month, int year)
{
  static int numdays[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
  int jd;
  jd = numdays[month-1] + day + (month > 2 && !(year & 3));
  return (jd);
}


double RatioDistMeanDist(int day)
{
  double d0;
  d0 = (double)(day - 1) * 0.017214;
  return (1.00011 + 0.034221 * cos(d0) + 0.00128 * sin(d0) +
          0.000719 * cos(2.*d0) + 0.000077 * sin(2.*d0));
  /* Pielke: 8-58a */
}
