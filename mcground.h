/*
   DEFINITION MODULE mcground
   In diesem Modul werden Bodenbilanzen gerechnet. 
*/

#define MCGROUND
#define NGROUNDLAYER 6

typedef struct  {
  double Tg[NGROUNDLAYER], /* Aktuelle Temperatur auf 0cm 1cm 10cm 1m*/
  	 z0,            /* Bodenrauhigkeit */
	 Cg, ks,	/* Waermekapazitaet, thermal diffusivity */
	 X, r,		/* Relative Feuchte des Bodens und Feuchte-Widerstand */
	 albedo,
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
	 z,		/* Hoehe des ersten Punktes */
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
  EntitySet a;          /* Zustand der Luft im ersten Punkt oberhalb des Bodens */
  Vector slope;		/* Vektor senkrecht zur Oberfläche */
  short int firstabove;
}  GroundParam;

extern GroundParam *ground;
extern double grlevel[];

double VirtPot(double T, double q, double p);

double AbsoluteTemp(double theta, double press, double humid);
/* Berechnet die Absolute Temperatur ausgehend von der virtuel-potentiellen-
   der Feuchte und des Druckes. */

double AbsolutePointTemp(int k, int i, int j);

double EquivalentPointTemp(int k, int i, int j);

double RelativeHumidity(int k, int i, int j);

void CalcSlope(Vector *v, int i, int j);

void GroundInterface(double timeofday, int dayofyear, double tinc, Vector *sunpos);

