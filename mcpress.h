/* DEFINITION MODULE MCServ
   Diverse Rechenroutienen. */


double CalcPileEnergy(int i, int j);

double SatHumidity(double T, double p);
/* Berechnet die Saettigungsfeuchte in Abhaengigkeit von T und p.
   Das Resultat wird in g/g (resp. kg/kg) Mischungsverhaeltnis an-
   gegeben. */

void CalcMeanHydrostaticPressure(void);
/* Berechnet die mittlere hydrostatische Druckverteilung */

double CalcTopEnergy(void);

void CalculateLayerAverage(double *l, PointStatus *stat, double *bar);
/* Berechnet den Durchschnitt ueber alle nz Layers beginned mit Layer *l.
   Der Status der Punkte wird in *stat abgefragt. Das Resultat wird in
   bar[] gespeichert. */

void CalcHydrostaticPressure(void);
/* Calculates the small scale pressure distibution, with the hydrostatic law. */

void SolveForPressure(int tinc);
/* This routine applies the conjugated gradient method, to solve for the non-hydrostatic
   pressure. The conjugated gradient routine was implemented according to N.R. p. 84. */
