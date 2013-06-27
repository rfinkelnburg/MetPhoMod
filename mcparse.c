/*
   IMPLEMENTATION MODULE mcparse.c
   Dieses Modul liest und interpretiert mit Hilfe des Modules mcgrammar.y
   das Input-File des Meteochem-Programmes.
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "mcglobal.h"
#include "sunpos.h"
#include "mcground.h"
#include "mcgroundclasses.h"
#include "mcpress.h"
#include "mcparse.h"
#include "mcemiss.h"
#include "mcfilter.h"
#include "mcdep.h"
#include "mchemparse.h"
#include "mchem.h"
#include "mcclouds.h"
#include "mcnest.h"
#include "mpdata.h"
#include "mcppm.h"
#include "mcrely.h"
#include "mc_module.hh"
#include "mc_section.hh"
#include "mc_variable.hh"
#ifdef PARALLEL
#include <pvm3.h>
#include "mcparallel.h"
#endif
#include "mc_group.t"

#define LAYER_DIM (X_DIM | Y_DIM)
#define NO_DIM 0

typedef struct  {
  BOOL isvalid;
  VarDesc *vd;        /* Index auf das variable-Feld */
  int has_range, xmin, xmax, ymin, ymax, zmin, zmax;
}  SetVarType;

SetVarType setvar;
int inplineno, inputerror;
IntVectPtr taptr = {&nta, tincarray};

OptionDesc option[MAXOPTION] =
{"off", "switch off", ON_OFF, FALSE,
 "on" , "switch on", ON_OFF, TRUE,
 "none", "no coriolis forces are considered", CORIOLIS_OPTION, NOCORIOLIS,
 "full", "all coriolis effects are considered", CORIOLIS_OPTION, FULLCORIOLIS,
 "differential", "coriolis acceleration is assumed to be included in p0",
    CORIOLIS_OPTION, DIFFERENTIALCORIOLIS,
 "geostrophic", "Coriolis forces are calculated relative to the synoptic Wind (Vars: Ugeos, Vgeos)",
    CORIOLIS_OPTION, GEOSTROPHICCORIOLIS,
 "no-pressure", "Pressure field and Wind acceleration will not be calculated", PRESSUREALG_OPTION, NOPRESSURE,
 "hydrostatic", "Calculate hydrostatic pressure-fields", PRESSUREALG_OPTION, HYDROSTATIC,
 "non-hydrostatic", "Calculate nonhydrostatic pressure-fields", PRESSUREALG_OPTION, NONHYDROSTATIC,
 "free"       , "Bordervalues are identical to their nearest neighbour", BORDERTYPE_OPTION, FREEBORDER,
 "cyclic"     , "cyclic border conditions", BORDERTYPE_OPTION, CYCLIC,
 "mirrored"   , "mirror border conditions", BORDERTYPE_OPTION, MIRRORED,
 "constant-inflow", "Border values are given by the user", BORDERTYPE_OPTION, CONSTANTBORDER,
 "auto-constant", "Border values depend on wind direction", BORDERTYPE_OPTION, AUTOCONSTANTBORDER,
 "sponge"     , "Border values are fixed, but effects are damped", BORDERTYPE_OPTION, SPONGEBORDER,
 "no-filter", "No spatial filtering is applied", FILTERTYPE_OPTION, NO_FILTER,
 "pepper",    "A Pepper-Filter is applied", FILTERTYPE_OPTION, PEPPER_FILTER,
 "shapiro",   "A Shapiro-Filter is applied", FILTERTYPE_OPTION, SHAPIRO_FILTER,
 "diffusion", "A diffusive filter is applied (= Shapiro first order)", FILTERTYPE_OPTION, DIFFUSION_FILTER,
 "no-turbulence","No boundary layer turbulence is considered", TURBTYPE_OPTION, NO_TURB,
 "TTTT", "Transilient Turbulence Theory Turbulence Parametrisation", TURBTYPE_OPTION, TTTT,
 "k-eps", "K-epsilon (E-epsilon) Turbulence parametrisation", TURBTYPE_OPTION, KEPS_T,
 "mpdata", "mpdata advection algorithm (Smolarkievicz, 1990)", ADVECTIONTYPE_OPTION, MPDATA_A,
 "ppm", "PPM advection algorithm (Colella and Woodward, 1983), Clappier (1998)", ADVECTIONTYPE_OPTION, PPM_A,
 "calculate", "Calculate the soil temperatures using MethPhomod's soil module", SOILTEMP_OPTION, TRUE,
 "preset", "Use preset soil temperatures", SOILTEMP_OPTION, FALSE,
 "clear-sky", "Use simple clear-sky radiation module", RADIATION_OPTION, FALSE,
 "twostream", "Use twostream multiscattering radiation module, based on Madronich's TUV", RADIATION_OPTION, TRUE
};


typedef enum  {
  LENGTH_U, TIME_U, SPEED_U, TEMP_U, PRESS_U,
  POWER_P_AREA_U, CHEMISTRY_U, MIXING_U, ARC_U
}  UnitCat;
 
char unit[][6] = {
"m", "sec", "m/s", "K", "Pa", "W/m^2", "ppb", "kg/kg", "deg"};

#define MAXVARIABLE 144

GroundParam grmark;

VarDesc variable[MAXVARIABLE] =
{VarDesc("U", "Wind in WE-Direction", unit[SPEED_U], GRID_VAL, (void *)UWIND, ALL_DIM, INITIAL_DATA,
	 SET_DEFAULT, 0, -120., 120., NORMAL_NUM, 0, NULL),
 VarDesc("Ground-Interface", "Calculate Ground-Interactions? (on/off)", NULL, INT_PTR, &groundinterface,
	 NO_DIM, OPTIONS,SET_DEFAULT, 1, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, NULL),
 VarDesc("Turbulence-Type", "Type of turbulence parametrisation", NULL, INT_PTR, &turbtype, NO_DIM,
	 OPTIONS, SET_DEFAULT, 2, 0, 2, TURBTYPE_OPTION, 0, NULL),
 VarDesc("Shielding_Factor", "Shielding factor of plants", NULL, GROUND_PARAM, &grmark.sigf, LAYER_DIM,
	 GROUND_BORDER, MUST_BE_SET, 0, 0, 1, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Scalar-Advection", "Calculate Advection of Scalars? (on/off)", NULL, INT_PTR, &advection, NO_DIM,
	 OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, NULL),
 VarDesc("Print-Borders", "Print values at the border in the output-File? (on/off)", NULL, INT_PTR,
	 &printwithborder, NO_DIM, OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, 0, NULL),
 VarDesc("Wind-Advection", "Calculate Advection of Wind Speed? (on/off)", NULL, INT_PTR, &windadvection,
	 NO_DIM, OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, NULL),
 VarDesc("Solar-Radiation", "Calculate the solar radiation (on/off)", NULL, INT_PTR,
	 &shortwaveradiation, NO_DIM, OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, NO_ERRORS_ALLOWED,
	 NULL),
 VarDesc("Clouds", "Consider Clouds? (on/off)", NULL, INT_PTR, &cloudwater, NO_DIM,
	 OPTIONS, SET_DEFAULT, 0, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, NULL),
 VarDesc("Cloud-Ice", "Consider Cloud Ice and Snow? (on/off)", NULL, INT_PTR, &cloudice, NO_DIM,
	 OPTIONS, SET_DEFAULT, 0, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, NULL),
 VarDesc("Timed-Dump", "Make multiple dump-files with time stamp.", NULL, INT_PTR, &timeddump, NO_DIM,
	 OPTIONS, SET_DEFAULT, 0, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, NULL),
 VarDesc("V", "Wind in SN-Direction", unit[SPEED_U], GRID_VAL, (long *)VWIND, ALL_DIM, INITIAL_DATA,
	 SET_DEFAULT, 0, -120., 120., NORMAL_NUM, 0, NULL),
 VarDesc("W", "Wind in vertical direction", unit[SPEED_U], GRID_VAL, (long *)WWIND, ALL_DIM, INITIAL_DATA,
	 CALCULATED_VALUES, 0, -30., 30., NORMAL_NUM, 0, NULL),
 VarDesc("Uflux", "Uflux", unit[SPEED_U], DOUBLE_PTR, NULL, ALL_DIM, INITIAL_DATA,
	 CALCULATED_VALUES, 0, -120., 120., NORMAL_NUM, 0, NULL),
 VarDesc("Vflux", "Wind in SN-Direction", unit[SPEED_U], DOUBLE_PTR, NULL, ALL_DIM, INITIAL_DATA,
	 CALCULATED_VALUES, 0, -120., 120., NORMAL_NUM, 0, NULL),
 VarDesc("Wflux", "Wind in vertical direction", unit[SPEED_U], DOUBLE_PTR, NULL, ALL_DIM, INITIAL_DATA,
	 CALCULATED_VALUES, 0, -30., 30., NORMAL_NUM, 0, NULL),
 VarDesc("TEMP", "virt. pot. Air Temperature", unit[TEMP_U], GRID_VAL, (long *)TEMP, ALL_DIM, INITIAL_DATA,
	 MUST_BE_SET, 290., 240., 400., NORMAL_NUM, 0, NULL),
 VarDesc("Q", "Humidity of air", unit[MIXING_U], GRID_VAL, (long *)HUMIDITY, ALL_DIM, INITIAL_DATA,
	 MUST_BE_SET, 0.005, 0., 0.03, NORMAL_NUM, 0, NULL),
 VarDesc("Qwat", "Cloud water", unit[MIXING_U], GRID_VAL, (long *)CLOUDWATER, ALL_DIM, INITIAL_DATA,
	 SET_DEFAULT, 0., 0., 0.03, NORMAL_NUM, 0, OnCloudWater),
 VarDesc("Qrain", "Concentration of rain drops", unit[MIXING_U], GRID_VAL, (long *)RAINWATER, ALL_DIM,
	 INITIAL_DATA, SET_DEFAULT, 0., 0., 0.03, NORMAL_NUM, 0, OnCloudWater),
 VarDesc("Qice", "Cloud Ice/Cloud Snow", unit[MIXING_U], GRID_VAL, (long *)CLOUDICE, ALL_DIM, INITIAL_DATA,
	 SET_DEFAULT, 0., 0., 0.03, NORMAL_NUM, 0, OnCloudIce),
 VarDesc("TKE", "Turbulent kintetic energy", "Pa/m^3", GRID_VAL, (long *)TKE, ALL_DIM, INITIAL_DATA,
	 SET_DEFAULT, 0., 0., 1.e5, NORMAL_NUM, 0, OnKEpsTurb),
 VarDesc("eps", "Dissipation", "???", GRID_VAL, (long *)EPS, ALL_DIM, INITIAL_DATA, 
	 SET_DEFAULT, 0., 0., 1.e5, NORMAL_NUM, 0, OnKEpsTurb),
 VarDesc("Tabs", "Absolute Air-Temperature", unit[TEMP_U], PROC_VAL, (void *)AbsolutePointTemp, ALL_DIM,
	 INITIAL_DATA, CALCULATED_VALUES, 290., 240., 340., NORMAL_NUM, 0, NULL),
 VarDesc("Te", "equivalent potential Air-Temperature", unit[TEMP_U], PROC_VAL, (void *)EquivalentPointTemp,
	 ALL_DIM, INITIAL_DATA, CALCULATED_VALUES, 290., 240., 340., NORMAL_NUM, 0, NULL),
 VarDesc("H", "Relative Humidity", "%", PROC_VAL, (void *)RelativeHumidity, ALL_DIM, INITIAL_DATA,
	 CALCULATED_VALUES, 290., 240., 340., NORMAL_NUM, 0, NULL),
 VarDesc("AirConz", "Concentration of air", "molec/cm^3", PROC_VAL, (void *)AirConz, ALL_DIM, INITIAL_DATA,
	 CALCULATED_VALUES, 290., 240., 340., NORMAL_NUM, 0, NULL),
 VarDesc("Vtr", "Terminal velocity of rain", unit[SPEED_U], PROC_VAL, (void *)TerminalVelocityOfRain,
	 ALL_DIM, INITIAL_DATA, CALCULATED_VALUES, 0., 0., 50., NORMAL_NUM, 0, NULL),
 VarDesc("Start", "Time to start simulation", unit[TIME_U], LONG_PTR, &tstart, NO_DIM, TIME,
	 SET_DEFAULT, 0, 0, 1000000, INT_NUM, 0, NULL),
 VarDesc("End", "Time to end simulation", unit[TIME_U], LONG_PTR, &tend, NO_DIM, TIME, MUST_BE_SET,
	 0, 0, 100000000, INT_NUM, 0, NULL),
 VarDesc("dt", "Time-Step of simulation", unit[TIME_U], LONG_PTR, &tincmax, NO_DIM, TIME, MUST_BE_SET,
	 20, 1, 3600, INT_NUM, 0, NULL),
 VarDesc("dt-gallery", "Possible-Values for dt", unit[TIME_U], INT_VECT_PTR, &taptr, MAXTINC,
	 TIME, MUST_BE_SET, 20, 1, 200, INT_VECTOR_VAL, 0, NULL),
 VarDesc("dt-IR", "Time-Step for atmospheric IR-Radiation", unit[TIME_U], LONG_PTR, &irtinc, NO_DIM,
	 GROUND_BORDER, MUST_BE_SET, 60, 0, 600, INT_NUM, 0, OnGroundInterface),
 VarDesc("dt-chem", "Time-Step of chemical species", unit[TIME_U], LONG_PTR, &tchem, NO_DIM, CHEMISTRY,
	 MUST_BE_SET, 20, 1, 600, INT_NUM, 0, OnChemistry),
 VarDesc("photofact", "Factor to multiply to the photochemical rate constants", NULL, DOUBLE_PTR,
	 &photofact, NO_DIM, CHEMISTRY, SET_DEFAULT, 1., 0., 100., NORMAL_NUM, 0, OnChemistry),
 VarDesc("fully-implicit", "Calculate Chemistry fully implicit? (on/off)", NULL, INT_PTR, &allfast, NO_DIM,
	 CHEMISTRY, SET_DEFAULT, 0, 0, 1, ON_OFF, 0, OnChemistry),
 VarDesc("dumptime", "Interval to make a memory dump", NULL, LONG_PTR, &dumptime, NO_DIM, TIME, SET_DEFAULT,
	 0, 0, 1000000, INT_NUM, 0, NULL),
 VarDesc("dumppath", "Path of directory for dump-files", NULL, CHAR_PTR, dumppath, NO_DIM, TIME,
	 SET_DEFAULT, 0, 0, 0, TEXT_VAL, 0, NULL),
 VarDesc("Time", "Day-Time of start of simulation", unit[TIME_U], DOUBLE_PTR, &timeofday, NO_DIM,
	 TIME, MUST_BE_SET, 0, 0, 24, TIME_VAL, 0, NULL),
 VarDesc("Date", "Date of start of simulation", "Day of year (M-D-Y)", LONG_PTR, &dayofyear, NO_DIM, TIME,
	 MUST_BE_SET, 0, 0, 366, DATE_VAL, 0, NULL),
 VarDesc("Km", "Turbulence Parameter", NULL, DOUBLE_PTR, NULL, ALL_DIM, INITIAL_DATA,
	 CALCULATED_VALUES, 0, 0., 100., NORMAL_NUM, 0, NULL),
 VarDesc("IRheat", "atmospheric radiative heating(+)/cooling(-)", "K/sec", DOUBLE_PTR, NULL, ALL_DIM,
	 INITIAL_DATA, CALCULATED_VALUES, 0, 0., 100., NORMAL_NUM, 0, NULL),
 VarDesc("nx", "Number of Nodes in x (WE)-Direction", NULL, INT_PTR, &nx, NO_DIM, GRID, MUST_BE_SET,
	 0, 1, 500, INT_NUM, NO_ERRORS_ALLOWED, NULL),
 VarDesc("ny", "Number of Nodes in y (SN)-Direction", NULL, INT_PTR, &ny, NO_DIM, GRID, MUST_BE_SET,
	 0, 1, 500, INT_NUM, NO_ERRORS_ALLOWED, NULL),
 VarDesc("nz", "Number of Nodes in z (vertical)-Direction", NULL, INT_PTR, &nz, NO_DIM, GRID,
	 MUST_BE_SET, 0, 4, 100, INT_NUM, NO_ERRORS_ALLOWED, NULL),
 VarDesc("nsubs", "Number of chemical species", NULL, INT_PTR, &nsubs, NO_DIM, GRID, SET_DEFAULT,
	 0, 0, 2000, INT_NUM, 0, NULL),
 VarDesc("dx", "Grid Distance in the x-Direction", unit[LENGTH_U], DOUBLE_PTR, &dx, NO_DIM, GRID,
	 MUST_BE_SET, 0, 300, 40000, SINGLE_VAL, 0, NULL),
 VarDesc("dy", "Grid Distance in the y-Direction", unit[LENGTH_U], DOUBLE_PTR, &dy, NO_DIM, GRID,
	 MUST_BE_SET, 0, 300, 40000, SINGLE_VAL, 0, NULL),
 VarDesc("workers", "Number of workers in a parallel run", NULL, INT_PTR, &workerswanted, NO_DIM,
	 GRID, SET_DEFAULT, 0, 2, 64, SINGLE_VAL, 0, NULL),
 VarDesc("topo", "topography. Normally height above sea level", unit[LENGTH_U], DOUBLE_PTR, NULL,
	 LAYER_DIM, ENVIRONMENT, SET_DEFAULT, 0, 0, 8000, NORMAL_NUM, 0, NULL),
 VarDesc("reflevel", "Deepest ground-level in model domain.", unit[LENGTH_U], DOUBLE_PTR, &reflevel,
	 NO_DIM, ENVIRONMENT, MUST_BE_SET, 0, 0, 8000, NORMAL_NUM, 0, NULL),
 VarDesc("level", "Height of model layers above \"reflevel\"", unit[LENGTH_U], DOUBLE_PTR, NULL,
	 Z_DIM, ENVIRONMENT, MUST_BE_SET, 0, 0, 15000, NORMAL_NUM, 0, NULL),
 VarDesc("Longitude", "geographical Longitude (East of Greenwich)", unit[ARC_U], DOUBLE_PTR, &Xlong, NO_DIM, 
	 ENVIRONMENT, MUST_BE_SET, 0, -180., 180., SINGLE_VAL, 0, OnGroundInterface),
 VarDesc("Latitude", "Geographical Latitude", unit[ARC_U], DOUBLE_PTR, &Xlat,
	 NO_DIM, ENVIRONMENT, MUST_BE_SET, 0, -90, 90, SINGLE_VAL, 0, OnCoriolis),
 VarDesc("Zonediff", "Local-Time - UTC", "h", DOUBLE_PTR, &timezonediff, NO_DIM, ENVIRONMENT,
	 MUST_BE_SET, 0, -24, 25, SINGLE_VAL, 0, OnGroundInterface),
 VarDesc("MapAngle", "Angle between model-y and north (positive when model domain was turned to the left)",
	 "deg", DOUBLE_PTR, &turnmapangle, NO_DIM, ENVIRONMENT,
	 SET_DEFAULT, 0, -180., 180., SINGLE_VAL, 0, OnGroundInterface),
 VarDesc("Stratospheric-Ozone", "Stratospheric ozone", "m", DOUBLE_PTR, &stratozone, NO_DIM,
	 ENVIRONMENT, SET_DEFAULT, 0.4, 0.05, 0.8, SINGLE_VAL, 0, OnGroundInterface),
 VarDesc("Precipitable-Water", "Precipitable Water above modelling Domain", "g/cm^2", DOUBLE_PTR,
	 &topprecipwater, NO_DIM, ENVIRONMENT, SET_DEFAULT, 1., 0, 5, SINGLE_VAL, 0,
	 OnGroundInterface),
 VarDesc("Turbidity", "Turbidity-Parameter (optimize manually!)", NULL, DOUBLE_PTR, &turbidity, NO_DIM,
	 ENVIRONMENT, SET_DEFAULT, 0.08, 0, 0.8, SINGLE_VAL, 0, OnGroundInterface),
 VarDesc("TopPressure", "Pressure on top of modelling domain", unit[PRESS_U], DOUBLE_PTR, &toppp, NO_DIM,
	 ENVIRONMENT, SET_DEFAULT, 79000, 1000, 100000, SINGLE_VAL, 0, NULL),
 VarDesc("Print-Residual", "Print convergence information of SOR-Algorithm (only in non-hydrostatic mode)",
	 NULL, INT_PTR, &printresid, NO_DIM, OPTIONS, SET_DEFAULT, 0, 0, 1, ON_OFF, 0,  NULL),
 VarDesc("Main-works-too", "In parallel mode: one worker is placed on the main host", NULL, INT_PTR,
	 &workeronthishost, NO_DIM, OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, 0, NULL),
 VarDesc("Damping-Layer", "Use a four damping layers at the top", NULL, INT_PTR,
	 &dampinglayer, NO_DIM, OPTIONS, SET_DEFAULT, 0, 0, 1, ON_OFF, 0, NULL),
 VarDesc("Coriolis", "How to consider Coriolis Forces (none/full/differential)", NULL, INT_PTR,
	 &coriolistype, NO_DIM,
	 OPTIONS, SET_DEFAULT, (double)GEOSTROPHICCORIOLIS, 0, 8, CORIOLIS_OPTION, 0, NULL),
 VarDesc("Advection-Type", "How to calculate Advection (mpdata/ppm)", NULL, INT_PTR, &advectiontype, NO_DIM,
	 OPTIONS, SET_DEFAULT, (double)MPDATA_A, 0, 8, ADVECTIONTYPE_OPTION, 0, NULL),
 VarDesc("Soil-Temperature", "How to treat soil temperature (calculate/preset)", NULL, INT_PTR, &calcsoiltemp,
	 NO_DIM, OPTIONS, SET_DEFAULT, (double)TRUE, 0, 1, SOILTEMP_OPTION, 0, NULL),
 VarDesc("Radiation-Alg", "How to calculate radiation", NULL, INT_PTR, &radiationtype,
	 NO_DIM, OPTIONS, SET_DEFAULT, (double)FALSE, 0, 1, RADIATION_OPTION, 0, NULL),
 VarDesc("iord", "Iteration order of Smolarkievicz Advection algorithm", NULL, INT_PTR, &smiord, NO_DIM,
	 OPTIONS, SET_DEFAULT, 1., 0, 5, INT_NUM, 0, OnAdvection),
 VarDesc("No-Oscillating", "Use Non-oscillating option in Smolarkievicz Advection algorithm (on/off)",
	 NULL, INT_PTR, &smnonos,
	 NO_DIM, OPTIONS, SET_DEFAULT, 0., 0, 1, ON_OFF, 0, OnAdvection),
 VarDesc("PressureAlg", "Indicates the method of pressure evaluation", NULL, INT_PTR,
	 &pressuretype, NO_DIM, OPTIONS, SET_DEFAULT, (double)HYDROSTATIC, 0, 2, PRESSUREALG_OPTION,
	 0, NULL),
 VarDesc("Spatial-Filter-Type", "The type of the spatial low-pass filter to apply", NULL, INT_PTR, &filtertype,
	 NO_DIM, OPTIONS, SET_DEFAULT, (double)NO_FILTER, 0, 3, FILTERTYPE_OPTION, 0, NULL),
 VarDesc("Spatial-Filter-Val", "Applies a spatial filter after each time-step with damping value", NULL,
	 DOUBLE_PTR, &spatialfilter, NO_DIM, OPTIONS, SET_DEFAULT, 0.001, 0, 1, SINGLE_VAL, 0, NULL),
 VarDesc("dPress", "p' = (p - p0)", unit[PRESS_U], DOUBLE_PTR, NULL, ALL_DIM, INITIAL_DATA, 
	 CALCULATED_VALUES, 0, -100., 100., NORMAL_NUM, 0, NULL),
 VarDesc("density", "density(z)", "kg/m^3", DOUBLE_PTR, NULL, Z_DIM, INITIAL_DATA,
	 CALCULATED_VALUES, 0, 0, 100, NORMAL_NUM, 0, NULL),
 VarDesc("P0", "middle air pressure", "Pa", DOUBLE_PTR, NULL, Z_DIM, INITIAL_DATA,
	 CALCULATED_VALUES, 0, 60000, 130000, NORMAL_NUM, 0, NULL),
 VarDesc("delta_TopPressure", "p' in the top Layer (Pa)", unit[PRESS_U], DOUBLE_PTR, NULL, LAYER_DIM,
	 TOP_BORDER, SET_DEFAULT, 0, -2000, 2000, NORMAL_NUM, STARTING_WITH_ZERO, NULL),
 VarDesc("TopTemp", "virtual-potential Temparature above domain", unit[TEMP_U], DOUBLE_PTR, NULL, LAYER_DIM,
	 TOP_BORDER, MUST_BE_SET, 0, 0., 50., NORMAL_NUM, STARTING_WITH_ZERO, NULL),
 VarDesc("Ugeos", "Geostrophic Wind (U)", unit[SPEED_U], DOUBLE_PTR, NULL, Z_DIM,
	 TOP_BORDER, MUST_BE_SET, 0, -100, 100, NORMAL_NUM, 0, OnCoriolis),
 VarDesc("Vgeos", "Geostrophic Wind (V)", unit[SPEED_U], DOUBLE_PTR, NULL, Z_DIM,
	 TOP_BORDER, MUST_BE_SET, 0, -100, 100, NORMAL_NUM, 0, OnCoriolis),
 VarDesc("T_Surface", "Temperature at the ground surface (K)", unit[TEMP_U], GROUND_PARAM, &grmark.Tg[0],
	 LAYER_DIM, GROUND_BORDER, MUST_BE_SET, 0, 230, 370, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("T_Gr1", "Temperature at 1st Ground-Level (4cm below Surface)", unit[TEMP_U], GROUND_PARAM,
	 &grmark.Tg[1], LAYER_DIM, GROUND_BORDER,
	 CALCULATED_VALUES, 0, 230, 370, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("T_Gr2", "Temperature at 2nd Ground-Level (16cm below Surface)", unit[TEMP_U], GROUND_PARAM,
	 &grmark.Tg[2], LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 230, 370, NORMAL_NUM, 0,
	 OnGroundInterface),
 VarDesc("T_Gr3", "Temperature at 3rd Ground-Level (36cm below Surface)", unit[TEMP_U], GROUND_PARAM,
	 &grmark.Tg[3], LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 230, 370, NORMAL_NUM, 0, 
	 OnGroundInterface),
 VarDesc("T_Gr4", "Temperature at 4th Ground-Level (64cm below Surface)", unit[TEMP_U], GROUND_PARAM,
	 &grmark.Tg[4], LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 230, 370, NORMAL_NUM, 0,
	 OnGroundInterface),
 VarDesc("T_Ground", "Temperature at 1m below ground surface (K)", unit[TEMP_U], GROUND_PARAM, &grmark.
	 Tg[NGROUNDLAYER-1], LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 285, 240, 350, NORMAL_NUM, 0,
	 OnGroundInterface),
 VarDesc("Ground_U", "U-Wind in Ground-Layer", unit[SPEED_U], GROUND_PARAM, &grmark.a[UWIND],
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, -50., 50., NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Ground_V", "V-Wind in Ground-Layer", unit[SPEED_U], GROUND_PARAM, &grmark.a[VWIND],
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, -50., 50., NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Ground_TEMP", "Temperature in Ground-Layer", unit[TEMP_U], GROUND_PARAM, &grmark.a[TEMP],
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 230., 330., NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Ground_Q", "Humidity in Ground-Layer", unit[MIXING_U], GROUND_PARAM, &grmark.a[HUMIDITY],
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Rain_Intensity", "Rain Intensity at Ground", "mm/sec", GROUND_PARAM, &grmark.rain_intensity,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, OnCloudWater),
 VarDesc("Cumulated_Rain", "Cumulated rain since start of model run", "mm", GROUND_PARAM, &grmark.cum_rain,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, OnCloudWater),
 VarDesc("Snow_Intensity", "Rain Intensity at Ground", "mm/sec", GROUND_PARAM, &grmark.snow_intensity,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, OnCloudWater),
 VarDesc("Cumulated_Snow", "Cumulated snow since start of model run", "mm (Water equivalent)", GROUND_PARAM,
	 &grmark.cum_snow, LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, 
	 NORMAL_NUM, 0, OnCloudWater),
 VarDesc("TotalClouds", "Total cloud water and ice content in vertical pile", "mm", PROC_VAL, (void *)TotalClouds,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, OnCloudIce),
 VarDesc("TotalCloudWater", "Total cloud water content in vertical pile", "mm", PROC_VAL, (void *)TotalWaterClouds,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, OnCloudWater),
 VarDesc("TotalCloudIce", "Total cloud ice content in vertical pile", "mm", PROC_VAL, (void *)TotalIceClouds,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, OnCloudIce),
 VarDesc("T_Foliage", "Temperature of plant surface", unit[TEMP_U], GROUND_PARAM, &grmark.Tf, LAYER_DIM,
	 GROUND_BORDER, CALCULATED_VALUES, 0., 220., 350., NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Qg", "Ground-Flux", unit[POWER_P_AREA_U], GROUND_PARAM, &grmark.Qg, LAYER_DIM, GROUND_BORDER,
	 CALCULATED_VALUES, 0., -500., 500., NORMAL_NUM, 0, NULL),
 VarDesc("Rel_Humidity", "Relative Humidity of the air in the ground-pores", unit[MIXING_U], GROUND_PARAM, 
	 &grmark.X, LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 1, 0, 1, NORMAL_NUM, 0,
	 OnGroundInterface),
 VarDesc("Evaporation_Resistance", "Evaporation-Resistance of soil", NULL, GROUND_PARAM, &grmark.r,
	 LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 1.e-3, 0., 1.e-2, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("GroundClass", "Landusage class at cell, according to classlist in input file", NULL, GROUND_PARAM,
	 &grmark.groundclass,
	 LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 0, 0, 1e4, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Ground_Heat_Capacity", "Ground-Heat-Capacity", "J/(m^3*)", GROUND_PARAM, &grmark.Cg,
	 LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 2.5e6, 1e4, 1e8, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Ground_Diffusivity", "ks (=Ground-Diffusivity)", NULL, GROUND_PARAM, &grmark.ks,
	 LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 0.53e-6, 1e-8, 5e-3, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Albedo", "Albedo of soil", NULL, GROUND_PARAM, &grmark.albedo, LAYER_DIM, GROUND_BORDER,
	 SET_DEFAULT, 0.2, 0, 1, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Ground_Emissivity", "Ground-Emissivity", NULL, GROUND_PARAM, &grmark.emissity, LAYER_DIM,
	 GROUND_BORDER, SET_DEFAULT, 0.98, 0.4, 1, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Z0", "Ground roughness", unit[LENGTH_U], GROUND_PARAM, &grmark.z0, LAYER_DIM, GROUND_BORDER,
	 SET_DEFAULT, 1, 1e-6, 30, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Zt", "Height of first Point above Ground", unit[LENGTH_U], GROUND_PARAM, &grmark.z, LAYER_DIM,
	 GROUND_BORDER, SET_DEFAULT, 10, 3, 25, NORMAL_NUM, 0, NULL),
 VarDesc("Leaf_Area_Index", "Leaf-Area-Index", NULL, GROUND_PARAM, &grmark.La, LAYER_DIM, GROUND_BORDER,
	 SET_DEFAULT, 4, 0., 20, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Plant_Albedo", "Albedo of plants", NULL, GROUND_PARAM, &grmark.albedof, LAYER_DIM,
	 GROUND_BORDER, SET_DEFAULT, 0.15, 0.01, 0.8, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Plant_Emissivity", "Emissivity of plants", NULL, GROUND_PARAM, &grmark.emissityf, LAYER_DIM,
	 GROUND_BORDER, SET_DEFAULT, 0.95, 0.4, 1., NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Plant_Resistance", "Resistance of plants against evapotranspiration", NULL, GROUND_PARAM, 
	 &grmark.rcf, LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 150., 5., 500., NORMAL_NUM,
	 0, OnGroundInterface),
 VarDesc("Plant_Water_Content", "Water content of plants", NULL, GROUND_PARAM, &grmark.fwatercontent,
	 LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 0.6, 0.01, 1, NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("Rglobal", "Global Radiation", unit[POWER_P_AREA_U], GROUND_PARAM, &grmark.Rg, LAYER_DIM,
	 GROUND_BORDER, CALCULATED_VALUES, 0., 0., 1200., NORMAL_NUM, 0, NULL),
 VarDesc("Rglobal_cumulated", "Global Radiation integrated for the actual day", "J/m^2", GROUND_PARAM,
	 &grmark.Rgcum, LAYER_DIM,
	 GROUND_BORDER, CALCULATED_VALUES, 0., 0., 1200000., NORMAL_NUM, 0, NULL),
 VarDesc("Rnet", "Net radiation", unit[POWER_P_AREA_U], GROUND_PARAM, &grmark.R, LAYER_DIM,
	 GROUND_BORDER, CALCULATED_VALUES, 0., 0., 1000., NORMAL_NUM, 0, NULL),
 VarDesc("IRdown", "IR irradiation", unit[POWER_P_AREA_U], GROUND_PARAM, &grmark.IRdown, LAYER_DIM,
	 GROUND_BORDER, CALCULATED_VALUES, 0., 0., 1000., NORMAL_NUM, 0, NULL),
 VarDesc("Sun_Elevation_Angle", "Elevation angle of the sun (0 at night-time)", "deg", GROUND_PARAM,
	 &sunelevation, TIME_DIM, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 90.,
	 NORMAL_NUM, 0, OnGroundInterface),
 VarDesc("wtheta", "sensible heat flux", "m*K/sec", GROUND_PARAM, &grmark.wtheta, LAYER_DIM,
	 GROUND_BORDER, CALCULATED_VALUES, 0., -100., 100., NORMAL_NUM, 0, NULL),
 VarDesc("wq", "latent heat flux", "m*kg/(sec*kg)", GROUND_PARAM, &grmark.wq, LAYER_DIM, GROUND_BORDER,
	 CALCULATED_VALUES, 0., -1., 1., NORMAL_NUM, 0, NULL),
 VarDesc("zL", "z/L, where L = monin Obukhov Length.", unit[LENGTH_U], GROUND_PARAM, &grmark.zL,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0., -100., 100., NORMAL_NUM, 0, NULL),
 VarDesc("ustar", "u* the friction velocity", unit[SPEED_U], GROUND_PARAM, &grmark.ustar,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0., -100., 100., NORMAL_NUM, 0, NULL),
 VarDesc("phim", "a Businger Turbulence Parameter", NULL, GROUND_PARAM, &grmark.phim,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0., -100., 100., NORMAL_NUM, 0, NULL),
 VarDesc("phit", "a Businger Turbulence Parameter", NULL, GROUND_PARAM, &grmark.phit,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0., -100., 100., NORMAL_NUM, 0, NULL),
 VarDesc("BLH", "Height of boundary layer", NULL, GROUND_PARAM, &grmark.blh,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 10000., NORMAL_NUM, 0, NULL),
 VarDesc("R0_max", "Deposition resistance maximum", "s/m", GROUND_PARAM, &grmark.Drmax,
	 LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 400., 0., 5000., NORMAL_NUM, 0, NULL),
 VarDesc("R0_min", "Deposition resistance minimum (must be 0 on water surfaces!)", "s/m", GROUND_PARAM,
	 &grmark.Drmin,
	 LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 300., 0., 5000., NORMAL_NUM, 0, NULL),
 VarDesc("R0_night", "Deposition resistance in the night", "s/m", GROUND_PARAM, &grmark.Drnight,
	 LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 300., 0., 5000., NORMAL_NUM, 0, NULL),
 VarDesc("R0_wet", "Deposition resistance on a wet surface", "s/m", GROUND_PARAM, &grmark.Drnight,
	 LAYER_DIM, GROUND_BORDER, SET_DEFAULT, 300., 0., 5000., NORMAL_NUM, 0, NULL),
#ifdef RADFACT
 VarDesc("Rdir_Fact", "Factor for direkt Short-Wave Radiation", NULL, DOUBLE_PTR, NULL, ALL_DIM,
	 GROUND_BORDER, SET_DEFAULT, 1., 0., 3., NORMAL_NUM, 0, NULL),
 VarDesc("Rdiff_Fact", "Factor for diffuse Short-Wave Radiation", NULL, DOUBLE_PTR, NULL, ALL_DIM,
	 GROUND_BORDER, SET_DEFAULT, 1., 0., 6., NORMAL_NUM, 0, NULL),
#else
 VarDesc("Rdir_Fact", "Factor for direkt Short-Wave Radiation", NULL, DOUBLE_PTR, NULL, Z_DIM,
	 GROUND_BORDER, SET_DEFAULT, 1., 0., 3., NORMAL_NUM, 0, NULL),
 VarDesc("Rdiff_Fact", "Factor for diffuse Short-Wave Radiation", NULL, DOUBLE_PTR, NULL, Z_DIM,
	 GROUND_BORDER, SET_DEFAULT, 1., 0., 6., NORMAL_NUM, 0,  NULL),
#endif
 VarDesc("Rdir", "Direct Short-wave Radiation", unit[POWER_P_AREA_U], GROUND_PARAM, &grmark.Rdir,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 0., NORMAL_NUM, 0, NULL),
 VarDesc("Rdiff", "Diffuse Short-wave Radiation", unit[POWER_P_AREA_U], GROUND_PARAM, &grmark.Rdiff,
	 LAYER_DIM, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 0., NORMAL_NUM, 0, NULL),
 VarDesc("BorderType-West", "Type of Borderconditions", NULL, INT_PTR, &westbordertype, NO_DIM, BORDER_DATA,
	 SET_DEFAULT, (double)SPONGEBORDER, 0, 64, BORDERTYPE_OPTION, NO_ERRORS_ALLOWED, NULL),
 VarDesc("BorderType-East", "Type of Borderconditions", NULL, INT_PTR, &eastbordertype, NO_DIM, BORDER_DATA,
	 SET_DEFAULT, (double)SPONGEBORDER, 0, 64, BORDERTYPE_OPTION, NO_ERRORS_ALLOWED, NULL),
 VarDesc("BorderType-North", "Type of Borderconditions", NULL, INT_PTR, &northbordertype, NO_DIM, BORDER_DATA,
	 SET_DEFAULT, (double)SPONGEBORDER, 0, 64, BORDERTYPE_OPTION, NO_ERRORS_ALLOWED, NULL),
 VarDesc("BorderType-South", "Type of Borderconditions", NULL, INT_PTR, &southbordertype, NO_DIM, BORDER_DATA,
	 SET_DEFAULT, (double)SPONGEBORDER, 0, 64, BORDERTYPE_OPTION, NO_ERRORS_ALLOWED, NULL),
 VarDesc(nestvarname[0], "x-coordinates of origin of nested subdomain", unit[LENGTH_U], LONG_PTR,
	 &dummynestvar, NO_DIM, NESTING, MUST_BE_SET, 0, 0, 500000, INT_NUM, 0,
	 OnNestDomain),
 VarDesc(nestvarname[1], "y-coordinates of origin of nested subdomain", unit[LENGTH_U], LONG_PTR,
	 &dummynestvar, NO_DIM, NESTING, MUST_BE_SET, 0, 0, 500000, INT_NUM, 0,
	 OnNestDomain),
 VarDesc(nestvarname[2], "number of points along x of nested subdomain", NULL, LONG_PTR,
	 &dummynestvar, NO_DIM, NESTING, MUST_BE_SET, 0, 0, 1000, INT_NUM, 0,
	 OnNestDomain),
 VarDesc(nestvarname[3], "number of points along y of nested subdomain", NULL, LONG_PTR,
	 &dummynestvar, NO_DIM, NESTING, MUST_BE_SET, 0, 0, 1000, INT_NUM, 0,
	 OnNestDomain),
 VarDesc(nestvarname[4], "grid cell size along x of nested subdomain", unit[LENGTH_U], LONG_PTR,
	 &dummynestvar, NO_DIM, NESTING, MUST_BE_SET, 0, 0, 10000, INT_NUM, 0,
	 OnNestDomain),
 VarDesc(nestvarname[5], "grid cell size along y of nested subdomain", unit[LENGTH_U], LONG_PTR,
	 &dummynestvar, NO_DIM, NESTING, MUST_BE_SET, 0, 0, 10000, INT_NUM, 0,
	 OnNestDomain),
 VarDesc(nestvarname[6], "turn-angle of nested subdomain", "deg", DOUBLE_PTR,
	 &dummynestvar, NO_DIM, NESTING, MUST_BE_SET, 0, -360, 360, NORMAL_NUM, 0,
	 OnNestDomain),
 VarDesc(nestvarname[7], "printing interval", unit[TIME_U], LONG_PTR,
	 &dummynestvar, NO_DIM, NESTING, MUST_BE_SET, 0, 0, 28800, INT_NUM, 0,
	 OnNestDomain)
};

ValueDesc val;
SectionDesc *actualsection;
Date startdate;
Time starttime;

void mcierror(char *s)
{
  extern char *mcitext;
  fprintf(stderr, "Line %i. %s\n", inplineno, s);
  fflush(stderr);
}

void InputError(ErrorLevel level, const char *s, ...)
{
  static BOOL sectionprinted;
  va_list ap;
  char txt[256], txt2[256];
  if (s)  {
    va_start(ap, s);
    vsprintf(txt, s, ap);
    va_end(ap);
    if (!sectionprinted)  {
      printf("----------- %s -----------\n", actualsection->name);
      fflush(stdout);
      sectionprinted = TRUE;
    }
    switch (level)  {
      case WARNING :
         printf("Line %i. WARNING: %s\n", inplineno, txt);
         fflush(stdout);
         break;
      case ERROR   :
         sprintf(txt2, "ERROR: %s", txt);
         mcierror(txt2);
         break; 
    }
  }
  else  sectionprinted = FALSE;
}

char *OnGroundInterface(RelyCmmd cmmd, VarDesc *v)
{
  switch (cmmd)  {
    case IS_RELEVANT :
       return ((char *)groundinterface);
    case NAME_OF_RELVAR :
       return ("Ground-Interface");
    case VAL_OF_RELVAR :
       return (char *)(groundinterface ? "on" : "off");
  }
}

char *OnCloudWater(RelyCmmd cmmd, VarDesc *v)
{
  switch (cmmd)  {
    case IS_RELEVANT :
       return ((char *)cloudwater);
    case NAME_OF_RELVAR :
       return ("Clouds");
    case VAL_OF_RELVAR :
       return (char *)(cloudwater ? "on" : "off");
  }
}

char *OnCloudIce(RelyCmmd cmmd, VarDesc *v)
{
  switch (cmmd)  {
    case IS_RELEVANT :
       return ((char *)cloudice);
    case NAME_OF_RELVAR :
       return ("Cloud-Ice");
    case VAL_OF_RELVAR :
       return (char *)(cloudice ? "on" : "off");
  }
}

char *OnAdvection(RelyCmmd cmmd, VarDesc *v)
{
  switch (cmmd)  {
    case IS_RELEVANT :
       return ((char *)(advection || windadvection));
    case NAME_OF_RELVAR :
       return ("?Advection");
    case VAL_OF_RELVAR :
       return (char *)((advection && windadvection) ? "on" : "off");
  }
}

char *OnCoriolis(RelyCmmd cmmd, VarDesc *v)
{
  static char *coption[] = {
    "none", "full", "differential", "geostrophic"
  };
  if (*v->name == 'L')
    switch (cmmd)  {
      case IS_RELEVANT :
         return ((char *)(coriolistype != NOCORIOLIS || groundinterface));
      case NAME_OF_RELVAR :
         return ("Coriolis/Ground-Interface");
      case VAL_OF_RELVAR :
         return (coption[coriolistype]);
    }
  else  {
    switch (cmmd)  {
      case IS_RELEVANT :
         return ((char *)(coriolistype == GEOSTROPHICCORIOLIS));
      case NAME_OF_RELVAR :
         return ("Coriolis");
      case VAL_OF_RELVAR :
         return (coption[coriolistype]);
    }
  }
}

char *OnChemistry(RelyCmmd cmmd, VarDesc *v)
{
  switch (cmmd)  {
    case IS_RELEVANT :
       return ((char *)(nsubs > 0));
    case NAME_OF_RELVAR :
       return ("Chemistry");
    case VAL_OF_RELVAR :
       return (char *)((nsubs > 0) ? "on" : "off");
  }
}

char *OnKEpsTurb(RelyCmmd cmmd, VarDesc *v)
{
  switch (cmmd)  {
    case IS_RELEVANT :
       return ((char *)(turbtype == KEPS_T));
    case NAME_OF_RELVAR :
       return ("Turbulence-Type");
    case VAL_OF_RELVAR :
       return (char *)((nsubs > 0) ? "on" : "off");
  }
}

char *OnBorderType(RelyCmmd cmmd, VarDesc *v)
{
  int i, val;
  switch (cmmd)  {
    case IS_RELEVANT :
       switch (v->section)  {
         case NORTH_BORDER :
            return ((char *)((northbordertype & DATABORDER) > 0));
         case SOUTH_BORDER :
            return ((char *)((southbordertype & DATABORDER) > 0));
         case WEST_BORDER :
            return ((char *)((westbordertype & DATABORDER) > 0));
         case EAST_BORDER :
            return ((char *)((eastbordertype & DATABORDER) > 0));
       }
    case NAME_OF_RELVAR :
       switch (v->section)  {
         case NORTH_BORDER :
            return ("BorderType-North");
         case SOUTH_BORDER :
            return ("BorderType-South");
         case WEST_BORDER :
            return ("BorderType-West");
         case EAST_BORDER :
            return ("BorderType-East");
       }
    case VAL_OF_RELVAR :
       switch (v->section)  {
         case NORTH_BORDER :
            val = northbordertype;
         case SOUTH_BORDER :
            val = southbordertype;
         case WEST_BORDER :
            val = westbordertype;
         case EAST_BORDER :
            val = eastbordertype;
       }
       for (i = MAXOPTION;
            --i >= 0 && (option[i].type != BORDERTYPE_OPTION ||
            		 option[i].value != val); );
       if (i >= 0)  return (option[i].name);
       else         return ("noval (Error in implementation)");
  }
}

int FillArrayWithNumber(VarDesc *v, double dval, long ival, SetVarType *setv)
{
  int i, j, k;
  switch (v->storetype)  {
    case GRID_VAL :
       if (g[v->v.et])
         for (k = setv->zmin; k <= setv->zmax; k++)
           for (i = setv->xmin; i <= setv->xmax; i++)
             for (j = setv->ymin; j <= setv->ymax; j++)
               g[v->v.et][k*layer+i*row+j] = dval;
       break;
    case DOUBLE_PTR :
       switch (v->dims)  {
         case (ALL_DIM) :
            for (k = setv->zmin; k <= setv->zmax; k++)
              for (i = setv->xmin; i <= setv->xmax; i++)
                for (j = setv->ymin; j <= setv->ymax; j++)
                  v->v.d[k*layer+i*row+j] = dval;
            break;
         case (X_DIM | Y_DIM) :
            for (i = setv->xmin; i <= setv->xmax; i++)
              for (j = setv->ymin; j <= setv->ymax; j++)
                v->v.d[i*row+j] = dval;
            break;
         case (Y_DIM | Z_DIM) :
            for (k = setv->zmin; k <= setv->zmax; k++)
              for (j = setv->ymin; j <= setv->ymax; j++)
                v->v.d[k*row+j] = dval;
            break;
         case (X_DIM | Z_DIM) :
            for (k = setv->zmin; k <= setv->zmax; k++)
              for (i = setv->xmin; i <= setv->xmax; i++)
                v->v.d[k*xrow+i] = dval;
            break;
         case X_DIM : for (i = setv->xmin; i <= setv->xmax; i++)
           		v->v.d[i] = dval;
           	      break;
         case Y_DIM : for (j = setv->ymin; j <= setv->ymax; j++)
           		v->v.d[j] = dval;
                      break;
         case Z_DIM : for (k = setv->zmin; k <= setv->zmax; k++)
           		v->v.d[k] = dval;
           	      break;
         case 0     : *v->v.d = dval; break;
         default : InputError(ERROR, "Internal implementation error (2)");
                   break;
       }
       break;
    case LONG_PTR :
       *v->v.li = (long)ival;
       break;
    case INT_PTR :
       *v->v.i = ival;
       break;
    case GROUND_PARAM :
       for (i = setv->xmin; i <= setv->xmax; i++)
         for (j = setv->ymin; j <= setv->ymax; j++)
           v->v.g[i*row+j].Tg[0] = dval;
       break;
  }
  return (0);
}

void SetDefaultVal(VarDesc *v)
{
  SetVarType sv;
  int startoff;
  startoff = ((v->option & STARTING_WITH_ZERO) == 0);
  sv.vd = v;
  sv.has_range = 0;
  sv.xmin = sv.xmax = sv.ymin = sv.ymax = startoff;
  sv.zmin = sv.zmax = 0;
  if (v->dims & X_DIM)  sv.xmax = nx - startoff;
  if (v->dims & Y_DIM)  sv.ymax = ny - startoff;
  if (v->dims & Z_DIM)  sv.zmax = nzm;
  FillArrayWithNumber(v, v->defval, (long)v->defval, &sv);
}

VarDesc *BorderVariable(VarDesc *v, VarDesc *vh)
{
  *vh = *v;
  vh->storetype = DOUBLE_PTR;
  switch (actualsection->id)  {
    case NORTH_BORDER : vh->dims = X_DIM | Z_DIM;
			vh->v.d = northborder+vh->v.et*nz*xrow;
			vh->init = InitType(vh->option >> 8);
			break;
    case SOUTH_BORDER : vh->dims = X_DIM | Z_DIM;
			vh->v.d = southborder+vh->v.et*nz*xrow;
			vh->init = InitType(vh->option >> 8);
			break;
    case WEST_BORDER  : vh->dims = Y_DIM | Z_DIM;
			vh->v.d = westborder+vh->v.et*nz*row;
			vh->init = InitType(vh->option >> 8);
			break;
    case EAST_BORDER  : vh->dims = Y_DIM | Z_DIM;
			vh->v.d = eastborder+vh->v.et*nz*row;
			vh->init = InitType(vh->option >> 8);
			break;
    default	      : return (v);
  }
  return (vh);
}

int ChangeSectionTo(char *sname)
{
  VarDesc *v, vh;
  SectionDesc *sect;
  if (!sname)  {
    actualsection = NULL;
    return (0);
  }
  sect = McInterface::mcsect.FindSection(sname);
  if (sect->dependson && !sect->dependson->found)  {
    InputError(ERROR, "Section %s depends on section %s, which must be defined first.",
       sect->name, sect->dependson->name);
    return (1);
  }
  if (sect->startproc && sect->startproc(InputSection(actualsection->id), InputSection(sect->id)))  return (1);
  actualsection = sect;
  sect->found = TRUE;
  sect->calltype = WAS_CALLED;
  InputError(ERROR, NULL);
  if (sect->initialize)
    for (v = variable; v; v = v->next)  {
      if (v->section == actualsection->id && v->init == SET_DEFAULT)
	SetDefaultVal(BorderVariable(v, &vh));
  }
  return (0);
}

int SetGridToSect(InputSection oldsec, InputSection newsec)
{
  int i;
  VarDesc *v;
  for (v = variable; v; v = v->next)
    if (v->storetype == GRID_VAL)  {
      v->section = newsec;
      v->init = InitType(v->option >> 8);
      v->relieson = (newsec >= NORTH_BORDER && newsec <= EAST_BORDER ? OnBorderType : NULL);
    }
  return (0);
}

int SetGridToEmission(InputSection oldsec, InputSection newsec)
{
  int i;
  VarDesc *v;
  for (v = variable; v; v = v->next)
    if (v->storetype == GRID_VAL && v->origsection == CHEMISTRY)  {
      v->section = newsec;
      v->init = InitType(v->option >> 8);
      v->relieson = NULL;
    }
  return (0);
}

int TestOnSection(InputSection oldsec, InputSection newsec)
{
  if (oldsec < BORDER_DATA || oldsec > EAST_BORDER)  {
    InputError(ERROR, "Section \"%s\" must be called within \"BORDER-DATA\"", 
	       McInterface::mcsect.FindSection(newsec)->name);
    return (1);
  }
  if (newsec >= NORTH_BORDER && newsec <= EAST_BORDER)
    SetGridToSect(InputSection(oldsec), InputSection(newsec));
  return (0);
}

int FindOption(char *oname)
{
  int i;
  for (i = MAXOPTION; --i >= 0 && strcmp(oname, option[i].name); );
  if (i < 0)  {
    InputError(ERROR, "%s is not a valid option", oname);
    return (1);
  }
  val.type = option[i].type;
  val.v.ival = option[i].value;
  return (0);
}

void InitArrayVar(void)
{
  val.type = ARRAY_VAL;
  val.v.a.nx = val.v.a.ny = val.v.a.nz = val.v.a.has_dimension = val.v.a.nval = val.v.a.has_range = 0;
  if (setvar.isvalid)  {
    val.v.a.xmin = val.v.a.xmax = setvar.xmin;
    val.v.a.ymin = val.v.a.ymax = setvar.ymin;
    val.v.a.zmin = val.v.a.zmax = setvar.zmin;
  }
  else  {
    val.v.a.xmin = val.v.a.xmax = val.v.a.ymin = val.v.a.ymax = 1;
    val.v.a.zmin = val.v.a.zmax = 0;
  }
}

int SetDimension(Dimension d)
{
  if (setvar.isvalid && !(setvar.vd->dims & d))  {
    InputError(ERROR, "Illegal dimension for this Variable.");
    return (1);
  }
  val.v.a.has_dimension |= d;
  if (setvar.isvalid)
    switch (d)  {
      case X_DIM : val.v.a.xmin = setvar.xmin;
      		   val.v.a.xmax = setvar.xmax;
      		   break;
      case Y_DIM : val.v.a.ymin = setvar.ymin;
      		   val.v.a.ymax = setvar.ymax;
      		   break;
      case Z_DIM : val.v.a.zmin = setvar.zmin;
      		   val.v.a.zmax = setvar.zmax;
      		   break;
    }
  else
    switch (d)  {
      case X_DIM : val.v.a.xmax = (val.v.a.nx = nx) - 1; break;
      case Y_DIM : val.v.a.ymax = (val.v.a.ny = ny) - 1; break;
      case Z_DIM : val.v.a.zmax = (val.v.a.nz = nz) - 1; break;
    }
  return (0);
}

int SetRange(Dimension d, int min, int max)
{
  if (setvar.isvalid && (setvar.has_range & d))  {
    InputError(ERROR, "Range redeclared in value definition.");
    return (1);
  }
  if (val.v.a.has_range & d)  {
    InputError(ERROR, "You cannot define a Dimension twice.");
    return (1);
  }
  val.v.a.has_range |= d;
  switch (d)  {
    case X_DIM : val.v.a.xmin = min; val.v.a.xmax = max; break;
    case Y_DIM : val.v.a.ymin = min; val.v.a.ymax = max; break;
    case Z_DIM : val.v.a.zmin = min; val.v.a.zmax = max; break;
  }
  return (0);
}

void AllocateArraySpace()
{
  val.v.a.maxval = (val.v.a.xmax - val.v.a.xmin + 1)*(val.v.a.ymax - val.v.a.ymin + 1)*(val.v.a.zmax - val.v.a.zmin + 1);
  val.v.a.a = (double *)malloc(sizeof(double) * val.v.a.maxval);
  val.v.a.deletable = TRUE;
  val.vectype = ValueDesc::IS_VECT;
}

int AddArrayVal(double num)
{
  if (val.vectype == ValueDesc::IS_NORMAL) {
    InputError(ERROR, "This is not an array or vector variable!");
    return (1);
  }
  if (val.v.a.nval >= val.v.a.maxval)  {
    InputError(ERROR, "Too many initialisers. You only need %i", val.v.a.maxval);
    return (1);
  }
  if (val.vectype == ValueDesc::IS_INT_VECT) {
    InputError(ERROR, "This is an integer vector; floats are not allowed!");
    return (1);
  }
  else
    val.v.a.a[val.v.a.nval++] = num;
  return (0);
}

int AddIntArrayVal(long num)
{
  if (val.vectype == ValueDesc::IS_NORMAL) {
    InputError(ERROR, "This is not an array or vector variable!\n");
    return (1);
  }
  if (val.v.a.nval >= val.v.a.maxval)  {
    InputError(ERROR, "Too many initialisers. You only need %i", val.v.a.maxval);
    return (1);
  }
  if (val.vectype == ValueDesc::IS_INT_VECT)
    ((long *)val.v.a.a)[val.v.a.nval++] = num;
  else
    val.v.a.a[val.v.a.nval++] = (double)num;
  return (0);
}

#ifdef PARALLEL

VarDesc *GetNamedVar(char *name)
{
  VarDesc *v;
  for (v = variable; v && strcmp(name, v->name); v = v->next);
  return (v);
}

#endif

int SetVariable(const char *name)
{
  VarDesc *v;
  v = McInterface::mcvar.FindVariable(name);
  if (!v)  {
    InputError(ERROR, "%s is not a known variable or keyword.", name);
    return (1);
  }
  return SetVariable(v);
}

int SetVariable(const char *groupname, const char *name)
{
  VarDesc *v;
  v = McInterface::mcvar.FindVariable(groupname, InputSection(actualsection->id));
  if (!v)  {
    InputError(ERROR, "%s is not a known variable or keyword.", groupname);
    return (1);
  }
  if (v->storetype != GROUP_VAL) {
    InputError(ERROR, "%s is not a VariableGroup.", groupname);
    return (1);
  }
  v->init = WAS_SET_BY_USER;
  v = v->v.group->Find(name);
  if (!v) {
    InputError(ERROR, "%s is not a member of %s.", name, groupname);
    return (1);
  }
  return SetVariable(v);
}

int SetVariable(VarDesc *v)
{
  VarDesc *vc;
  static VarDesc vh;
  BOOL rval;
  vc = v;
  if (v->storetype == GRID_VAL)  {
    if (actualsection->id == DEPOSITION)
      vc = DepositionVariable(v, &vh);
    else if (actualsection->id == EMISSIONS)
      vc = EmissionVariable(v, 0, 0, NULL);
    else if (actualsection->id == REDUCTION)
      vc = EmissionVariable(v, 0, 0, NULL);
    else if (actualsection->id >= NORTH_BORDER && actualsection->id <= EAST_BORDER)
      vc = BorderVariable(v, &vh);
  }
  if (vc->section != actualsection->id)  {
    InputError(ERROR, "The variable %s was found in the wrong section. It belongs to %s",
       vc->name, McInterface::mcsect.FindSection(v->section)->name);
    return (1);
  }
  rval = !(v->option & STARTING_WITH_ZERO);
  setvar.vd = v; setvar.has_range = setvar.zmin = setvar.zmax = 0;
  setvar.xmin = setvar.xmax = setvar.ymin = setvar.ymax = rval;
  if (vc->inputtype & (VECTOR_VAL | INT_VECTOR_VAL))  {
    val.v.a.nx = val.v.a.ny = val.v.a.nz = val.v.a.has_dimension = val.v.a.nval = val.v.a.has_range = 0;
    val.v.a.maxval = vc->dims;   /* dims bezeichnet bei VECTOR_VAL die Anzahl der Elemente. */
    val.v.a.a = (vc->storetype == INT_VECT_PTR || vc->storetype == VECT_PTR
		 ? vc->v.vp->v : vc->v.d);
    val.v.a.deletable = FALSE;
    val.vectype = (vc->inputtype == VECTOR_VAL ?
		   ValueDesc::IS_VECT : ValueDesc::IS_INT_VECT);
    val.type = vc->inputtype;
  }
  else  {
    if (v->dims & X_DIM)  setvar.xmax = nx - rval;
    if (v->dims & Y_DIM)  setvar.ymax = ny - rval;
    if (v->dims & Z_DIM)  setvar.zmax = nzm;
    val.vectype = ValueDesc::IS_NORMAL;
  }
  setvar.isvalid = TRUE;
  return (0);
}

int RangeError(void)
{
  const char err1[] = "\"GROUNDLEVEL\" is not allowed in this context.";
  if (setvar.has_range & X_DIM)  {
    if (setvar.xmin == -1)  {
      InputError(ERROR, err1);
      return (1);
    }
    if ((val.v.a.xmin < setvar.xmin) || (val.v.a.xmax > setvar.xmax))  {
      InputError(ERROR, "X-Range (%i..%i) does not suit variable (%i..%i).",
         val.v.a.xmin, val.v.a.xmax, setvar.xmin, setvar.xmax);
      return (1);
    }
  }
  if (setvar.has_range & Y_DIM)  {
    if (setvar.xmin == -1)  {
      InputError(ERROR, err1);
      return (1);
    }
    if ((val.v.a.ymin < setvar.ymin) || (val.v.a.ymax > setvar.ymax))  {
      InputError(ERROR, "Y-Range (%i..%i) does not suit variable (%i..%i).",
         val.v.a.ymin, val.v.a.ymax, setvar.ymin, setvar.ymax);
      return (1);
    }
  }
  if (setvar.has_range & Z_DIM)  {
    if (setvar.xmin == -1)  {
      InputError(ERROR, err1);
      return (1);
    }
    if ((val.v.a.zmin < setvar.zmin) || (val.v.a.zmax > setvar.zmax))  {
      InputError(ERROR, "Z-Range (%i..%i) does not suit variable (%i..%i).",
         val.v.a.zmin, val.v.a.zmax, setvar.zmin, setvar.zmax);
      return (1);
    }
  }
  return (0);
}

int CopyRangeToVar()
{
  setvar.has_range = val.v.a.has_range;
  if (RangeError())  return (1);
  if (val.v.a.has_range & X_DIM)  {
    setvar.xmin = val.v.a.xmin;
    setvar.xmax = val.v.a.xmax;
  }
  if (val.v.a.has_range & Y_DIM)  {
    setvar.ymin = val.v.a.ymin;
    setvar.ymax = val.v.a.ymax;
  }
  if (val.v.a.has_range & Z_DIM)  {
    setvar.zmin = val.v.a.zmin;
    setvar.zmax = val.v.a.zmax;
  }
  if (setvar.has_range & ~setvar.vd->dims)  {
    InputError(ERROR, "Variable \"%s\" has invalid Dimensions.", setvar.vd->name);
    return (1);
  }
  return (0);
}

int AssignVariable()
{
  VarDesc *v;
  int i, j, k, ioff, joff, koff, ifa, jfa, kfa;
  long ival;
  double dval;
  VarDesc vh;
  v = setvar.vd;
  setvar.isvalid = FALSE;
  if (inputerror)  goto error_end;
  if (v->storetype == GRID_VAL && actualsection->id >= NORTH_BORDER &&
      actualsection->id <= EAST_BORDER)  {
    if (val.type == AVG_VAL)  {
      if (val.v.ival > nz || val.v.ival < 1)  {
        InputError(ERROR, "Value %d is out of range (1-%d)", val.v.ival, nz);
        inputerror = TRUE;
      }
      switch (actualsection->id)  {
	case NORTH_BORDER : northborderval[v->v.et] = val.v.ival; break;
	case SOUTH_BORDER : southborderval[v->v.et] = val.v.ival; break;
	case WEST_BORDER  : westborderval[v->v.et] = val.v.ival; break;
	case EAST_BORDER  : eastborderval[v->v.et] = val.v.ival; break;
      }
      v->init = WAS_SET_BY_USER;
      return (0);
    }
    v = BorderVariable(v, &vh);
  }
  if (v->storetype == GRID_VAL && actualsection->id == EMISSIONS)
    v = EmissionVariable(v, 0, 0, NULL);
  if (v->storetype == GRID_VAL && actualsection->id == REDUCTION)
    v = EmissionVariable(v, 0, 0, NULL);
  if (v->storetype == GRID_VAL && actualsection->id == DEPOSITION)
    v = DepositionVariable(v, &vh);
#ifdef PARALLEL
  if (parallel && v->storetype == GRID_VAL && gactive[v->v.et] && !g[v->v.et] &&
      !(g[v->v.et] = NEW(mesh, double)))  {
    InputError(ERROR, "Error allocating memory!");
    pvm_exit();
    exit (1);
  }
#endif
  if (v->storetype == GRID_VAL && !gactive[v->v.et])  {
    InputError(ERROR, "Variable \"%s\" is not active in the current calculation.",
       v->name);
    goto error_end;
  }
  if (!(val.type & v->inputtype))  {
    InputError(ERROR, "Illegal Type for this variable");
    goto error_end;
  }
/*  if (v->init == CALCULATED_VALUES)  {
    InputError(ERROR,
       "The variable \"%s\" will be calculated by the program and cannot be set by the user.", v->name);
    return (1);
  }  */
  if (v->relieson && !v->relieson(IS_RELEVANT, v))  {
    InputError(WARNING, "The variable \"%s\" has no effect when option \"%s\" is set to \"%s\".",
       v->name, v->relieson(NAME_OF_RELVAR, v), v->relieson(VAL_OF_RELVAR, v));
  }
  if (val.type == ARRAY_VAL)  {
    if (val.v.a.has_dimension & ~v->dims)  {
      InputError(ERROR, "Illegal dimension for this variable.");
      goto error_end;
    }
    if (val.v.a.has_range & setvar.has_range)  {
      InputError(ERROR, "Range defined twice. (In Variable and in Data)");
      goto error_end;
    }
    if (RangeError())  goto error_end;
    if (val.v.a.has_range & X_DIM)  {
      setvar.xmin = val.v.a.xmin;
      setvar.xmax = val.v.a.xmax;
    }
    if (val.v.a.has_range & Y_DIM)  {
      setvar.ymin = val.v.a.ymin;
      setvar.ymax = val.v.a.ymax;
    }
    if (val.v.a.has_range & Z_DIM)  {
      setvar.zmin = val.v.a.zmin;
      setvar.zmax = val.v.a.zmax;
    }
    for (i = val.v.a.nval; i--; )
      if (val.v.a.a[i] < v->rmin || val.v.a.a[i] > v->rmax)  {
        InputError(ERROR, "This Value (%le) is not within the permitted range (%le..%le)",
           val.v.a.a[i], v->rmin, v->rmax);
        goto error_end;
      }
    setvar.has_range |= val.v.a.has_range;
    ifa = (val.v.a.has_dimension & X_DIM ? setvar.xmax - setvar.xmin + 1 : 1);
    jfa = (val.v.a.has_dimension & Y_DIM ? setvar.ymax - setvar.ymin + 1 : 1);
    kfa = (val.v.a.has_dimension & Z_DIM ? setvar.zmax - setvar.zmin + 1 : 1);
    if (inputerror)  goto error_end;
    switch (v->storetype)  {
      case GRID_VAL :
         for (k = setvar.zmin; k <= setvar.zmax; k++)  {
           koff = ifa * jfa * (k - setvar.zmin) * ((val.v.a.has_dimension & Z_DIM) > 0);
           for (i = setvar.xmin; i <= setvar.xmax; i++)  {
             ioff = koff + jfa * (i - setvar.xmin) * ((val.v.a.has_dimension & X_DIM) > 0);
             for (j = setvar.ymin; j <= setvar.ymax; j++)  {
               g[v->v.et][k*layer+i*row+j] =
                 val.v.a.a[ioff + (j - setvar.ymin) * ((val.v.a.has_dimension & Y_DIM) > 0)];
             }
           }
         }
         break;
      case DOUBLE_PTR :
         switch (v->dims)  {
	   case (ALL_DIM) :
              for (k = setvar.zmin; k <= setvar.zmax; k++)  {
                koff = ifa * jfa * (k - setvar.zmin) * ((val.v.a.has_dimension & Z_DIM) > 0);
                for (i = setvar.xmin; i <= setvar.xmax; i++)  {
                  ioff = koff + jfa * (i - setvar.xmin) * ((val.v.a.has_dimension & X_DIM) > 0);
                  for (j = setvar.ymin; j <= setvar.ymax; j++)  {
                    v->v.d[k*layer+i*row+j] =
                      val.v.a.a[ioff + (j - setvar.ymin) * ((val.v.a.has_dimension & Y_DIM) > 0)];
                  }
                }
              }
              break;
           case (X_DIM | Y_DIM) :
              for (i = setvar.xmin; i <= setvar.xmax; i++)  {
                ioff = jfa * (i - setvar.xmin) * ((val.v.a.has_dimension & X_DIM) > 0);
                for (j = setvar.ymin; j <= setvar.ymax; j++)  {
                  v->v.d[i*row+j] =
                    val.v.a.a[ioff + (j - setvar.ymin) * ((val.v.a.has_dimension & Y_DIM) > 0)];
                }
              }
              break;
           case (Y_DIM | Z_DIM) :
              for (k = setvar.zmin; k <= setvar.zmax; k++)  {
                koff = jfa * (k - setvar.zmin) * ((val.v.a.has_dimension & Z_DIM) > 0);
                for (j = setvar.ymin; j <= setvar.ymax; j++)  {
                  v->v.d[k*row+j] =
                    val.v.a.a[koff + (j - setvar.ymin) * ((val.v.a.has_dimension & Y_DIM) > 0)];
                }
              }
              break;
	   case (X_DIM | Z_DIM) :
              for (k = setvar.zmin; k <= setvar.zmax; k++)  {
                koff = ifa * (k - setvar.zmin) * ((val.v.a.has_dimension & Z_DIM) > 0);
                for (i = setvar.xmin; i <= setvar.xmax; i++)  {
                  v->v.d[k*xrow+i] =
                    val.v.a.a[koff + (i - setvar.xmin) * ((val.v.a.has_dimension & X_DIM) > 0)];
                }
              }
              break;
           case X_DIM : for (i = setvar.xmin; i <= setvar.xmax; i++)
           		  v->v.d[i] = val.v.a.a[i - setvar.xmin];
           		break;
           case Y_DIM : for (j = setvar.ymin; j <= setvar.ymax; j++)
           		  v->v.d[j] = val.v.a.a[j - setvar.ymin];
           		break;
           case Z_DIM : for (k = setvar.zmin; k <= setvar.zmax; k++)
           		  v->v.d[k] = val.v.a.a[k - setvar.zmin];
           		break;
           default : InputError(ERROR, "Internal implementation error (2)");
                     break;
         }
         break;
      case GROUND_PARAM :
         for (i = setvar.xmin; i <= setvar.xmax; i++)  {
           ioff = jfa * (i - setvar.xmin) * ((val.v.a.has_dimension & X_DIM) > 0);
           for (j = setvar.ymin; j <= setvar.ymax; j++)  {
             v->v.g[i*row+j].Tg[0] =
               val.v.a.a[ioff + (j - setvar.ymin) * ((val.v.a.has_dimension & Y_DIM) > 0)];
           }
         }
         break;
      default :
         InputError(ERROR, "Internal implementation error.");
         goto error_end;
    }
    if (val.v.a.deletable)  free(val.v.a.a);
  }
  else  {     /* kein Array */
    switch (val.type)  {
      case INT_NUM		: ival = val.v.ival; dval = (double)val.v.ival; break;
      case DOUBLE_NUM		: dval = val.v.dval; break;
      case ON_OFF		:
      case CORIOLIS_OPTION	:
      case BORDERTYPE_OPTION    :
      case FILTERTYPE_OPTION    :
      case TURBTYPE_OPTION      :
      case ADVECTIONTYPE_OPTION :
      case PRESSUREALG_OPTION	:
      case RADIATION_OPTION     :
      case SOILTEMP_OPTION      : ival = val.v.ival; dval = (double)val.v.ival; break;
      case TIME_VAL		: dval = val.v.time.hour + val.v.time.minute / 60.; 
      				  starttime = val.v.time;
      				  break;
      case DATE_VAL		: ival = JulianDay(val.v.date.day, val.v.date.month, val.v.date.year);
      				  startdate = val.v.date;
     				  dval = (double)ival; break;
      case AVG_VAL		: break;
      case VECTOR_VAL		:
	if (val.v.a.nval < val.v.a.maxval)  {
	  if (v->storetype == VECT_PTR)
	    *setvar.vd->v.vp->nval = val.v.a.nval;
	  else {
	    InputError(ERROR, "Not enough vector elements (%i). Expected %i.",
		       val.v.a.nval, val.v.a.maxval);
	    return (1);
	  }
	}
	goto normal_end;
      case INT_VECTOR_VAL       : *setvar.vd->v.ivp->nval = val.v.a.nval;
      				  goto normal_end;
      case TEXT_VAL		: strcpy(setvar.vd->v.txt, val.v.txt);
      				  goto normal_end;
      default : InputError(ERROR, "Internal implementation error (4)");
      		return (1);
    }
    if (val.type != AVG_VAL && (dval < v->rmin || dval > v->rmax))  {
      InputError(ERROR, "This Value (%le) is not within the permitted range (%le..%le)", dval, v->rmin, v->rmax);
      goto error_end;
    }
    if (FillArrayWithNumber(v, dval, ival, &setvar))  goto error_end;
  }
normal_end :
  if (!setvar.has_range)  setvar.vd->init = WAS_SET_BY_USER;
  else if (v->init == MUST_BE_SET)  setvar.vd->init = WAS_PARTITALLY_SET;
  return (0);
error_end :
  if (val.type == ARRAY_VAL && val.v.a.deletable)  free(val.v.a.a);
  setvar.vd->init = MADE_AN_ERROR;
  return (1);
}

typedef struct UVarDesc  {
  char name[128];
  ValueDesc vd;
  UVarDesc *next;
} UVarDesc;

UVarDesc *uvar;

int DefineAVariable(char *name)
{
  UVarDesc *newv;
  VarDesc *v;
  for (v = variable; v && strcmp(v->name, name); v = v->next);
  if (v)  {
    InputError(ERROR, "\"%s\" is fixed. It can not be redefined by the user.", name);
    return (1);
  }
  for (newv = uvar; newv && strcmp(name, newv->name); newv = newv->next);
  if (newv)  {
    InputError(ERROR, "There is a User-Variable \"%s\" already defined", name);
    return (1);
  }
  newv = (UVarDesc *)malloc(sizeof(UVarDesc));
  memcpy(&newv->vd, &val, sizeof(ValueDesc));
  strcpy(newv->name, name);
  if (val.type == ARRAY_VAL)
    newv->vd.v.a.deletable = FALSE;
  newv->next = uvar;
  uvar = newv;
  return (0);
}

int CopyVarToVal(char *name)
{
  UVarDesc *v;
  for (v = uvar; v && strcmp(name, v->name); v = v->next);
  if (!v)  {
    InputError(ERROR, "User-Variable \"%s\" not found.", name);
    return (1);
  }
  memcpy(&val, &v->vd, sizeof(ValueDesc));
  return (0);
}

void ClearAllUserVars()
{
  UVarDesc *v, *vn;
  for (v = uvar; v; v = vn)  {
    vn = v->next;
    free(v);
  }
}

int EndOfTime(void)
{
  int i;
  BOOL wrong = FALSE;
  tinc = tincmax;
  for (tai = nta-1; tai >= 0 && tinc != tincarray[tai]; tai--);
  if (wrong = (tai < 0))
    InputError(ERROR, "dt (%d) was not found in dt-gallery.", tinc);
  else  {
    for (i = nta; --i; )
      wrong |= tincarray[i-1] <= tincarray[i];
    if (wrong)
      InputError(ERROR, "dt-gallery values must be indicated in descendent order.");
  }
  return (wrong);
}

void PutToVarTable(void *p, char *name)
{
  VarDesc *v;
  for (v = variable; v && strcmp(v->name, name); v = v->next);
  if (!v)  {
    InputError(ERROR, "Internal Error in mcparse::PutToVarTable: Variable %s unkown!\n", name);
    exit (1);
  }
  v->v.li = (long *)p;
}


void AdjustGridsInVarTable(void)
{
  VarDesc *v;
  for (v = variable; v; v = v->next)
    if (v->storetype == GROUND_PARAM)
      v->v.d = (double *)ground + (v->v.d - (double *)&grmark);
}

int EndOfGrid(void)
{
  int i;
  short allocok;
  if (inputerror)  return (1);
  if (nz > MAXZGRID)  {
    InputError(ERROR, "Vertical levels are restricted to %d points in the actual Version.\n"
       "(You tried to allocate %d)", MAXZGRID, nz);
    return (1);
  }
  nxm = nx++; nym = ny++; nzm = nz - 1; maxentity = SUBS + nsubs;
  mesh = nz * (layer = (nx+1)*(row = ny+1));
  xrow = nx+1;
  allocok = !!(g = (double **)calloc(maxentity, sizeof(double *))) &&
            !!(gactive = (BOOL *)calloc(maxentity, sizeof(BOOL)));
  for (i = 0; i < maxentity && allocok; i++)
    if (i < CLOUDWATER || i > EPS)
      gactive[i] = TRUE;
#ifdef PARALLEL
  if (parallel && master)
    parallel = nx / (OVERLAP+3) > 1;
  if (!parallel || !master)  {
#endif
    for (i = 0; i < maxentity && allocok; i++)
      if (i < CLOUDWATER || i > EPS)
        allocok = !!(g[i] = (double *)calloc(mesh, sizeof(double)));
      allocok &= 
         !!(flux[UWIND] = (double *)calloc(mesh, sizeof(double))) &&
         !!(flux[VWIND] = (double *)calloc(mesh, sizeof(double))) &&
         !!(flux[WWIND] = (double *)calloc(mesh, sizeof(double)));
#ifdef PARALLEL
   }
#endif
  if (!allocok ||
      !(avg = (double *)calloc(maxentity * nz, sizeof(double))) ||
      !(topo = (double *)calloc(layer, sizeof(double))) ||
      !(level = (double *)calloc(nz, sizeof(double))) ||
      !(leveli = (double *)calloc(nz, sizeof(double))) ||
      !(ldiff = (double *)calloc(nz, sizeof(double))) ||
      !(zcenter = (double *)calloc(nz, sizeof(double))) ||
      !(zfaces = (double *)calloc(nz, sizeof(double))) ||
      !(toppress = (double *)calloc(layer, sizeof(double))) ||
      !(toptemp = (double *)calloc(layer, sizeof(double))) ||
      !(ugeos = (double *)calloc(nz, sizeof(double))) ||
      !(vgeos = (double *)calloc(nz, sizeof(double))) ||
      !(Km = (double *)calloc(mesh, sizeof(double))) ||
      !(press = (double *)calloc(mesh, sizeof(double))) ||
#ifdef RADFACT
      !(Rdirfact = (double *)calloc(mesh, sizeof(double))) ||
      !(Rdifffact = (double *)calloc(mesh, sizeof(double))) ||
#else
      !(Rdirfact = (double *)calloc(nz, sizeof(double))) ||
      !(Rdifffact = (double *)calloc(nz, sizeof(double))) ||
#endif
      !(tmplayer = (double *)calloc(layer, sizeof(double))) ||
      !(tmplayer2 = (double *)calloc(layer, sizeof(double))) ||
      !(nlevelpt = (int *)calloc(nz, sizeof(int))) ||
      !(pstat = (PointStatus *)calloc(mesh, sizeof(PointStatus))) ||
      !(ground = (GroundParam *)calloc(layer, sizeof(GroundParam))) ||
      !(pp = (double *)calloc(nz+1, sizeof(double))+1) ||
      !(density = (double *)calloc(nz+2, sizeof(double))+1) ||
      !(northborder = (double *)calloc(maxentity * nz * xrow, sizeof(double))) ||
      !(southborder = (double *)calloc(maxentity * nz * xrow, sizeof(double))) ||
      !(westborder = (double *)calloc(maxentity * nz * row, sizeof(double))) ||
      !(eastborder = (double *)calloc(maxentity * nz * row, sizeof(double))) ||
      !(northborderval = (int *)calloc(maxentity, sizeof(int))) ||
      !(southborderval = (int *)calloc(maxentity, sizeof(int))) ||
      !(westborderval = (int *)calloc(maxentity, sizeof(int))) ||
      !(eastborderval = (int *)calloc(maxentity, sizeof(int))))  {
    InputError(ERROR, "ERROR in allocating Memory.\nProbably not enough memory on this Computer for the posed problem.");
    return (1);
  }
  dxi = 1. / dx; dyi = 1. / dy;
  dx2 = 2. * dx; dy2 = 2. * dy;
  dx2i = 1. / dx2; dy2i = 1. / dy2;
  ixx = 1. / (dx * dx); iyy = 1. / (dy * dy);
  PutToVarTable(Km, "Km");
  PutToVarTable(topo, "topo");
  PutToVarTable(press, "dPress");
  PutToVarTable(toppress, "delta_TopPressure");
  PutToVarTable(toptemp, "TopTemp");
  PutToVarTable(ugeos, "Ugeos");
  PutToVarTable(vgeos, "Vgeos");
  PutToVarTable(flux[UWIND], "Uflux");
  PutToVarTable(flux[VWIND], "Vflux");
  PutToVarTable(flux[WWIND], "Wflux");
  PutToVarTable(Rdirfact, "Rdir_Fact");
  PutToVarTable(Rdifffact, "Rdiff_Fact");
  PutToVarTable(pp, "P0");
  PutToVarTable(density, "density");
  PutToVarTable(zfaces, "level");
  AdjustGridsInVarTable();
  return (0);
}

int EndOfOptions(void)
{
  if (cloudwater)
    gactive[CLOUDWATER] = gactive[RAINWATER] = TRUE;
  if (cloudice)
    gactive[CLOUDICE] = TRUE;
  if (turbtype == KEPS_T)
    gactive[TKE] = gactive[EPS] = TRUE;
#ifdef PARALLEL
  if (!parallel || !master)  {
#endif
    if (cloudwater)  {
      if (!(g[CLOUDWATER] = (double *)calloc(mesh, sizeof(double))) ||
          !(g[RAINWATER] = (double *)calloc(mesh, sizeof(double))))  goto no_memory;
    }
    if (cloudice)  {
      if (!(g[CLOUDICE] = (double *)calloc(mesh, sizeof(double))))  goto no_memory;
    }
    if (turbtype == KEPS_T)  {
      if (!(g[TKE] = (double *)calloc(mesh+layer, sizeof(double))) ||
	  !(g[EPS] = (double *)calloc(mesh+layer, sizeof(double))))  goto no_memory;
    }
    if (groundinterface)  {
      if (!(irdiv = (double *)calloc(mesh, sizeof(double))))  goto no_memory;
      PutToVarTable(irdiv, "IRheat");
    }
#ifdef PARALLEL
  }
#endif
  return (0);
no_memory :
  InputError(ERROR, "ERROR in allocating Memory.\nProbably not enough memory on this Computer for the posed problem.");
  return (1);
}

void AdjustTopoBorders(void)
{
  int i, j, loc;
  double minh;
  if (ny > 6)  {
    if (northbordertype & (AUTOCONSTANTBORDER | FREEBORDER))
      for (i = nx; --i; )  {
	minh = topo[i*row+ny-2];
	for (j = ny-1; j <= ny; j++)
	  if (topo[i*row+j] < minh)
	    topo[i*row+j] = minh;
	  else
	    minh = topo[i*row+j];
      }
    if (southbordertype & (AUTOCONSTANTBORDER | FREEBORDER))
      for (i = nx; --i; )  {
	minh = topo[i*row+2];
	for (j = 1; j; j--)
	  if (topo[i*row+j] < minh)
	    topo[i*row+j] = minh;
	  else
	    minh = topo[i*row+j];
      }
  }
  if (nx > 6)  {
    if (eastbordertype & (AUTOCONSTANTBORDER | FREEBORDER))
      for (j = ny+1; j--; )  {
	minh = topo[(nx-2)*row+j];
	for (i = nx-1; i <= nx; i++)
	  if (topo[i*row+j] < minh)
	    topo[i*row+j] = minh;
	  else
	    minh = topo[i*row+j];
      }
    if (westbordertype & (AUTOCONSTANTBORDER | FREEBORDER))
      for (j = ny+1; j--; )  {
	minh = topo[2*row+j];
	for (i = 1; i; i--)
	  if (topo[i*row+j] < minh)
	    topo[i*row+j] = minh;
	  else
	    minh = topo[i*row+j];
      }
  }
  if (ny > 6)  {
    if (northbordertype & (AUTOCONSTANTBORDER | FREEBORDER))
      for (i = nx+1; i--; )  {
	minh = topo[i*row+ny-2];
	for (j = ny-1; j <= ny; j++)
	  if (topo[i*row+j] < minh)
	    topo[i*row+j] = minh;
	  else
	    minh = topo[i*row+j];
      }
    if (southbordertype & (AUTOCONSTANTBORDER | FREEBORDER))
      for (i = nx+1; i--; )  {
	minh = topo[i*row+2];
	for (j = 1; j; j--)
	  if (topo[i*row+j] < minh)
	    topo[i*row+j] = minh;
	  else
	    minh = topo[i*row+j];
      }
  }
}

int EndOfEnvironment()
{
  int i, j, k, loc;
  double height, min, max;
  Entity et;
  /* Check plausibility of levels */
  //  printf("in EndOfEnvironment\n");
  // printf(""); // Some problems with AIX
  for (i = 0; i < nzm; i++)  {
    if (zfaces[i] > zfaces[i+1])  {
      InputError(ERROR, "levels must be given in increasing order!\n"
      	   "(I found %lf - %lf)", zfaces[i], zfaces[i+1]);
      inputerror = 1;
      return (1);
    }
    if (zfaces[i] + 5. > zfaces[i+1])  {
      InputError(ERROR, "Difference between levels should at least be 5m!\n"
      	   "(I found %lf - %lf = %lfm)", zfaces[i], zfaces[i+1], zfaces[i+1] - zfaces[i]);
      inputerror = 1;
      return (1);
    }
  }
  toplevel = zfaces[nzm];
  for (i = nz; i--; )
    level[i] = zfaces[i] - (i ? zfaces[i-1] : 0);
  for (i = nz; i--; )  {
    zcenter[i] = zfaces[i] - level[i] * 0.5;
    leveli[i] = 1. / level[i];
    ldiff[i] = 0.5 * (level[i] + (i ? level[i-1] : 0.));
  }
/* Round topography to quads */
/*  ShapiroFilter(topo, pstat+nzm*layer, 3, 1.);  */
  for (i = layer; --i; )
    if (topo[i] < reflevel)  {
      topo[i] = reflevel;
      ground[i].firstabove = 0;
    }
    else  {
      if (topo[i] > zfaces[nzm-1] + reflevel)  {
	topo[i] = zfaces[nzm-1] + reflevel;
	InputError(WARNING, "Point %i/%i is higher than second highest level, I need one free line above topography!. I set point to %lf\n", i / row, i % row, topo[i]);
      }
      height = topo[i] - reflevel;
      j = 0;
      if (2 * height < zfaces[0])  topo[i] = reflevel;
      else  {
        for ( ; j < nzm && fabs(zfaces[j] - height) > fabs(zfaces[j+1] - height); j++);
        topo[i] = zfaces[j] + reflevel;
        j++;
      }
      ground[i].firstabove = j;
    }
  for (i = nx; --i; )
    for (j = ny; --j; )  {
      loc = i*row+j;
      MinMax(&min, &max, &topo[loc-1], &topo[loc+1], &topo[loc-row], &topo[loc+row], NULL);
      if (topo[loc] < min)  topo[loc] = min;
    }
  for (k = 0; k < nz; k++)  {
    height = zfaces[k] + reflevel - 1.;
    nlevelpt[k] = 0;
    for (i = layer; i--; )
      if (height > topo[i])
	pstat[k*layer+i] = NORMALPOINT;
      else
	pstat[k*layer+i] = UNDERGROUND;
    for (i = nx; --i; )
      for (j = ny; --j; )
        nlevelpt[k] += (height > topo[i*row+j]);
#ifdef PARALLEL
    if (master)
#endif
      if (!nlevelpt[k])  {
        InputError(ERROR, "No points outside topography on level %d.", k);
        inputerror = TRUE;
      }
  }
  for (i = nx; --i; )
    for (j = ny; --j; )  {
      loc = i*row + j;
      k = ground[loc].firstabove;
      for (et = HUMIDITY+1; et--; )
        ground[loc].a[et] = (g[et] ? g[et][k*layer+loc] : 0.);
      CalcSlope(&ground[loc].slope, i, j);
    }
  Coriol3 = 1.45444e-4 * sin(Pi * Xlat / 180.);
  Coriol2 = 1.45444e-4 * cos(Pi * Xlat / 180.);
  return (0);
}

int EndOfSection()
{
  int error;
  VarDesc *v;
  error = 0;
  if (actualsection->initialize)
    for (v = variable; v; v = v->next)
      if (v->section == actualsection->id)  {
	if (!v->relieson || v->relieson(IS_RELEVANT, v))  {
	  if (v->init == MUST_BE_SET)  {
	    InputError(ERROR, "The variable %s *must* be set in this section.", v->name);
	    error = 1;
	  }
	  else if (v->init == WAS_PARTITALLY_SET)  {
	    InputError(ERROR, "I don't know if variable %s is totally initialized.\n         You are using subranges only!",
	       v->name);
	    error = 1;
	  }
/*
	  else if (v->init == SET_DEFAULT)  {
	    InputError(WARNING, "No value for variable %s.\n         I will use a default of %le",
	       v->name, v->defval);
	    fflush(stdout);
	  }
*/
	}
	if (v->init == MADE_AN_ERROR && (v->option & NO_ERRORS_ALLOWED))  {
          InputError(ERROR, "FATAL: Cannot recover from previous errors!");
          exit (1);
	}
      }
  //  printf("in section %s\n", actualsection->name);
  if (!error) McInterface::modmanager.EndOfSection(actualsection);
  return (error || actualsection->endproc && actualsection->endproc());
}

int EndOfGround()
{
  int i, j, k;
  GroundParam *gr;
  AdjustTopoBorders();
  EndOfEnvironment();
  for (j = 1; j < ny; j++)
    for (i = 1; i < nx; i++)  {
      gr = &ground[i*row+j];
      gr->Tf = gr->Tg[0];
      for (k = NGROUNDLAYER-1; --k; )
        gr->Tg[k] = gr->Tg[NGROUNDLAYER-1] + (gr->Tg[0] - gr->Tg[NGROUNDLAYER-1]) *
	  exp(-sqrt(3.63610e-5 / gr->ks) * grlevel[k]);
      gr->absorb = 1. - (1. - gr->sigf) * gr->albedo - gr->sigf * gr->albedof;
    }
  return (0);
}

void InitializeData(void)
{
  int i, j, k, im, ip, jm, jp, loc;
  Entity et;
  EndOfGround();
  FillGroundValues();
  /* Setze alle Werte innerhalb der Berge 0 */
  for (i = nx+1; i--; )
    for (j = ny+1; j--; )
      for (k = ground[i*row+j].firstabove; k--; )
        for (et = maxentity; et--; )
          if (g[et])  g[et][k*layer+i*row+j] = 0.;
  /* Setze die Ground-Werte und pstat am Rand gemaess den Randbedingungen */
  switch (northbordertype)  {
    case FREEBORDER : case CONSTANTBORDER : case SPONGEBORDER :
    case AUTOCONSTANTBORDER :
       jm = nym; break;
    case CYCLIC :
       jm = 1; break;
    case MIRRORED :
       jm = nym - 1; break;
    case NEIGHBOUR :
       jm = ny; break;
  }
  switch (southbordertype)  {
    case FREEBORDER : case CONSTANTBORDER : case SPONGEBORDER :
    case AUTOCONSTANTBORDER :
       jp = 1; break;
    case CYCLIC :
       jp = nym; break;
    case MIRRORED :
       jp = 2; break;
    case NEIGHBOUR :
       jp = 0; break;
  }
  switch (westbordertype)  {
    case FREEBORDER : case CONSTANTBORDER : case SPONGEBORDER :
    case AUTOCONSTANTBORDER :
       ip = 1; break;
    case CYCLIC :
       ip = nxm; break;
    case MIRRORED :
       ip = 2; break;
    case NEIGHBOUR :
       ip = 0; break;
  }
  switch (eastbordertype)  {
    case FREEBORDER : case CONSTANTBORDER : case SPONGEBORDER :
    case AUTOCONSTANTBORDER :
       im = nxm; break;
    case CYCLIC :
       im = 1; break;
    case MIRRORED :
       im = nxm-1; break;
    case NEIGHBOUR :
       im = nx; break;
  }
  /* Setze Ground-Luft-Werte gleich dem ersten Punkt darueber. */
  for (i = nx; --i; )
    for (j = ny; --j; )
      for (et = HUMIDITY+1; et--; )  {
        loc = i*row+j;
        ground[loc].a[et] = g[et][ground[loc].firstabove*layer+loc];
      }
  for (i = nx+1; i--; ) {
    loc = i*row;
    topo[loc+ny] = topo[loc+jm];
    topo[loc] = topo[loc+jp];
    ground[loc+ny] = ground[loc+jm];
    ground[loc] = ground[loc+jp];
    for (k = nz; k--; )
      pstat[k*layer+loc] = pstat[k*layer+loc+jp];
    for (k = nz; k--; )
      pstat[k*layer+loc+ny] = pstat[k*layer+loc+jm];
  }
  for (j = ny+1; j--; )  {
    topo[nx*row+j] = topo[im*row+j];
    topo[j] = topo[ip*row+j];
    ground[nx*row+j] = ground[im*row+j];
    ground[j] = ground[ip*row+j];
    for (k = nz; k--; )
      pstat[k*layer+j] = pstat[k*layer+ip*row+j];
    for (k = nz; k--; )
      pstat[k*layer+nx*row+j] = pstat[k*layer+im*row+j];
  }
  turnmapangle *= Pi / 180.;
  InitDeposition();
  NoteTransport(TRUE);
  if (advectiontype == PPM_A)
    InitializePPM();
  McInterface::modmanager.Init2();
}

extern "C" int mciwrap()
{
  return (1);
}

void InitVarTable(void)
{
  int i;
  for (i = MAXVARIABLE-1; i--; )
    variable[i].next = variable + i + 1;
  variable[MAXVARIABLE-1].next = NULL;
}

void ShowTitle(char *title)
{
  int i;
  for (i = strlen(title)+4; i--; ) printf("*"); printf("\n");
  printf("* %s *\n", title);
  for (i = strlen(title)+4; i--; ) printf("*"); printf("\n\n");
  worktitle = strdup(title);
}

void TestDataConsistency(void)
{
  if (nx < 7 && (westbordertype == SPONGEBORDER || eastbordertype == SPONGEBORDER))  {
    fprintf(stderr, "ERROR: nx must at least be 6 (is actually %i) for Sponge-Border\n", nx-1);
    inputerror = 1;
  }
  if (ny < 7 && (northbordertype == SPONGEBORDER || southbordertype == SPONGEBORDER))  {
    fprintf(stderr, "ERROR: ny must at least be 6 (is actually %i) for Sponge-Border\n", ny-1);
    inputerror = 1;
  }
}

int mciparse(void);  // mciparse is function generated by yacc (or bison...)

// extern int mcidebug;

int ParseInput(char *fname)
{
  extern FILE *mciin;
  if (!(mciin = fopen(fname, "r")))  {
    fprintf(stderr, "FATAL ERROR: Unable to open File \"%s\"\n", fname);
    return (2);
  }
  inputerror = 0;
  inplineno = 1;
//  mcidebug = 1;
  uvar = NULL;
  setvar.isvalid = FALSE;
//  McInterface::modmanager.Init1();
  ChangeSectionTo(NULL);
  mciparse();
  fclose(mciin);
  ClearAllUserVars();
  McInterface::mcsect.TestOnSections();
  TestDataConsistency();
  return (inputerror);
}
