/*
   MODULE mcextrad
   Contains the radiation at the top of the atmosphere and interpolates
   it to the desired grid
*/

#ifndef INCLUDE_MCEXTRAD
#define INCLUDE_MCEXTRAD

void Rebin(int no, const float *wco, const float *vo,
	   int nw, const double *wl, double *v);

void Interpolate(int no, const float *wco, const float *vo,
		 int nw, const double *wl, double *v);

void SetupExtRadiation(int nw, const double *wl, double *v);

#endif
