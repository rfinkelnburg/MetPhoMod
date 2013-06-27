/* DEFINITION MODULE MCLoop
   In diesem Modul sind die eigentlichen Itegrationsschritte zusammengefasst. */

/* void CalcTurbulentDerivatives(double tinc);  */

void SetAvgMountains(void);

void SetBoundaryFaces(void);

void InterpolateToFaces(void);
/* Interpolates the Wind-Speeds at the cell-boundaries.
   If a neighbour of a cell is in the mountain, the flux is zero. */

void InterpolateToCenter(double fact, int owfact);

void ApplyCoriolis(int tinc);
/* This procedure will apply the Coriolis-Forces to the Wind-Speeds in the center of the
   Grid-Cells. */

void ApplyBuoyancy(int tinc);
/* This Procedure will add the Bouyancy-Term to the vertical-Wind-Speed. This is normally
   is only necessary in non-hydrostatic mode.
   The equation used is:
      dw     |theta'   Cv p'|
      -- = g |------ - -----|
      dt     |theta0   Cp p0|
*/

void ApplyPressure(int tinc);
/* Will apply the wind acceleration, caused by pressure gradient, to the winds at
   the faces. */

void Continuity(void);

void CheckTimeStep(long *tinc, long *chemtinc, double cfact);
/*  This routine controls, if the actual time-step fullfills the Courant-Condition:
   
    V*dt <= cfact * dx (resp. dy, dz)
*/

void DampTop(void);

void SetBoundary(int net);
