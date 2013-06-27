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
#undef FALSE
#define		FALSE !TRUE

extern double Pi;
double S0 = 1376.;


extern double sqr(double x);

void SunVector(Vector *s, int day, double localtime,
               double longitude, double latitude, double timezonediff)
/* s:	resulting Sun Unity Vector	*/
/* day:	Julian Day of the year		*/
{
  double ds, t, sintz, costz, sin2tz, cos2tz, sin3tz, cos3tz;
  double lat, eqr, jd;
  const double dr = Pi / 180.;
  jd = 2.*Pi*(day - 1 + (localtime - timezonediff) / 24.) / 365.;
  sintz = sin(jd);
  costz = cos(jd);
  sin2tz = 2.*sintz*costz;
  cos2tz = costz*costz-sintz*sintz;
  sin3tz = sintz*cos2tz + costz*sin2tz;
  cos3tz = costz*cos2tz - sintz*sin2tz;
  ds = 0.006918 - 0.399912*costz  + 0.070257*sintz -
    0.006758*cos2tz + 0.000907*sin2tz -
    0.002697*cos3tz + 0.001480*sin3tz;

  eqr   = (0.000075 + 0.001868*costz  - 0.032077*sintz
    - 0.014615*cos2tz - 0.040849*sin2tz) * 12 / Pi;

  t = dr * (15 * (localtime - 12. + eqr - timezonediff) + longitude);
  lat = latitude * dr;
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

double AngleBetween(const Vector *s, const Vector *t)
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
