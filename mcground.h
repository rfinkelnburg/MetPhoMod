/*
   DEFINITION MODULE mcground
   In diesem Modul werden Bodenbilanzen gerechnet. 
*/

#ifndef INCLUDE_MCGROUND
#define INCLUDE_MCGROUND

#ifndef INCLUDE_MCGLOBAL
#include "mcglobal.h"
#endif

#ifndef INCLUDE_SUNPOS
#include "sunpos.h"
#endif

#define MCGROUND
#define NGROUNDLAYER 6

typedef struct  {
  double Tg[NGROUNDLAYER], /* Aktuelle Temperatur auf 0cm 1cm 10cm 1m*/
  	 z0,            /* Bodenrauhigkeit */
	 Cg, ks,	/* Waermekapazitaet, thermal diffusivity */
	 X, r,		/* Relative Feuchte des Bodens und Feuchte-Widerstand */
	 albedo,
         absorb,        /* Anteil der einfallenden Strahlung, die absorbiert wird (Boden und Pflanzen kombiniert) */
	 wtheta, wq,	/* Sensibler und latenter Waermefluss in die Atmosphaere */
	 zL,            /* Stabilitaetsparameter */
	 uw, vw,
	 Qg,		/* Bodenwaermefluss */
	 emissity,
	 sigf,		/* Pflanzlicher Bedeckungsgrad */
	 La,		/* Leaf Area Index */
	 Tf,		/* Temperatur der pflanzlichen Oberflaeche */
	 albedof,	/* Albedo des pflanz. Bewuchses */
	 emissityf,	/* Emmissitaet der Pflanzen */
	 rcf,		/* Pflanzlicher Widerstand gegen Verdunstung */
	 fwatercontent,	/* Parameter fuer den Wassergehalt der Pflanze */
	 Rg,		/* Global-Strahlung */
	 Rgcum,		/* Cumulated global radiation over one day */
	 R,		/* Net-Radiation */
	 IRdown,	/* IR-Einstrahlung */
	 Rdir,		/* Direkte Strahlung */
	 Rdiff,		/* Diffuse Strahlung */
    //	 z,		/* Hoehe des ersten Punktes */ removed with v2.2
	 ustar,		/* friction velocity */
	 phim,		/* Similarity function for momentum */
	 phit,		/* Similarity function for scalars */
	 blh,		/* Height of boundary layer */
	 ra,		/* Turbulenter Widerstand gegen Deposition */
	 Drmax,		/* Maximaler Depositionswiderstand */
	 Drmin,		/* Minimaler Depositionswiderstand */
	 Drnight,	/* Depositionswiderstand in der Nacht */
	 Drwet,		/* Depositionswiderstand bei nasser Oberflaeche */
	 rain_intensity,/* Intensity of rain */
	 cum_rain,      /* Cumulated rain */
	 snow_intensity,/* Intensity of snow */
	 cum_snow,      /* Cumulated snow */
	 groundclass;   /* Class of GroundCell */
  //  EntitySet a;  Entity-set has been removed from v2.2
  Vector slope;		/* Vektor senkrecht zur Oberfläche */
  short int firstabove;
}  GroundParam;

extern GroundParam *ground;
extern double grlevel[];

double VirtPot(double T, double q, double p);

double AbsoluteTemp(double theta, double press, double humid);
/* Berechnet die Absolute Temperatur ausgehend von der virtuel-potentiellen-
   der Feuchte und des Druckes. */

void CalcSlope(Vector *v, int i, int j);

void GroundInterface(double timeofday, int dayofyear, double tinc, Vector *sunpos);

#endif
