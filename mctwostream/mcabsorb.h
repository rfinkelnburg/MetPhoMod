/*
   MODULE mcadsorb
   Calculates optical absorption thicknesses for different species.
*/


#ifndef INCLUDE_MCABSORB
#define INCLUDE_MCABSORB

#define KW 300  /* These values must in any case correspond to the */
#define NZ 200  /* according consts in the "param" file */

typedef double AllZ[NZ];
typedef double AllZP[NZ+1];

class OzoneBands {
  int nw;
  double *xso3, *s226, *s263, *s298;
  const double *wl;
public :
  OzoneBands(int nw, const double *wl);
  ~OzoneBands(void);
  void OzoneAbsorption(int nz, const double *z, const double *o3,
		       const double *temp, double dobnew, AllZ *dto3);
     /* Ozone Absorption:
	
	INPUT-VARIABLES:
	  nz      : Number of z-levels
	  z       : height of z-levels
	  o3      : ozone concentrations in molec/cm^3
	  temp    : temperature at each level
	  dobnew  : total ozone (in Dobson-units)
	
	OUTPUT-Variables:
	  dto3    : Optical depth due to ozone absorption for
	            each wavelength on each level.
     */
};

class SO2Bands {
  int nw;
  double *xsso2;
public :
  SO2Bands(int nw, const double *wl);
  ~SO2Bands(void);
  void SO2Absorption(int nz, const double *z, const double *so2,
		     const double *temp, AllZ *dtso2);
};

class NO2Bands {
  int nw;
  double *xsno2;
public :
  NO2Bands(int nw, const double *wl);
  ~NO2Bands(void);
  void NO2Absorption(int nz, const double *z, const double *no2,
		     const double *temp, AllZ *dtno2);
};

/* Interface definitions for included fortran routines. */

extern "C" {

void seto2_(int *nz, int *nw, const double *wl, double *cz,
	    double *zen, AllZ *dto2);

void sphers_(int *nz, double *z, double *zen, AllZP *dsdh, int *nid);

void rtlink_(int *nz, double *z, double *ag, double *zen,
	     AllZP *dsdh, int *nid, double *dtrl, double *dto3,
	     double *dto2, double *dtso2, double *dtno2,
	     double *dtcld, double *omcld, double *gcld,
	     double *dtaer, double *omaer, double *gaer,
	     double *edir, double *edn, double *eup,
	     double *fdir, double *fdn, double *fup,
	     const double *fdn0, const double *cosradsoil);

}

#endif
