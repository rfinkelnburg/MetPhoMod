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
#ifdef PARALLEL
#include <pvm3.h>
#include "mcparallel.h"
#endif

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
 "ppm", "PPM advection algorithm (Colella and Woodward, 1983), Clappier (1998)", ADVECTIONTYPE_OPTION, PPM_A
};

int SetGridToSect(InputSection oldsec, InputSection newsec);
int SetGridToEmission(InputSection oldsec, InputSection newsec);
int TestOnSection(InputSection oldsec, InputSection newsec);
char *OnGroundInterface(RelyCmmd cmmd, VarDesc *v);
char *OnCloudWater(RelyCmmd cmmd, VarDesc *v);
char *OnCloudIce(RelyCmmd cmmd, VarDesc *v);
char *OnBorderType(RelyCmmd cmmd, VarDesc *v);
char *OnCoriolis(RelyCmmd cmmd, VarDesc *v);
char *OnAdvection(RelyCmmd cmmd, VarDesc *v);
char *OnChemistry(RelyCmmd cmmd, VarDesc *v);
char *OnKEpsTurb(RelyCmmd cmmd, VarDesc *v);

SectionDesc section[] =
{"TIME",  NULL_SECTION, TRUE, NULL, EndOfTime, MUST_BE_CALLED,
 "GRID",  TIME, TRUE, NULL, EndOfGrid, MUST_BE_CALLED,
 "OPTIONS", GRID, TRUE, NULL, EndOfOptions, MUST_BE_CALLED,
 "ENVIRONMENT", OPTIONS, TRUE, NULL, EndOfEnvironment, MUST_BE_CALLED,
 "INITIAL-DATA", GRID, TRUE, NULL, NULL, MUST_BE_CALLED,
 "BORDER-DATA", INITIAL_DATA, TRUE, NULL, NULL, MUST_BE_CALLED,
 "GROUND", ENVIRONMENT, TRUE, TestOnSection, NULL, MUST_BE_CALLED,
 "TOP", ENVIRONMENT, TRUE, TestOnSection, NULL, MUST_BE_CALLED,
 "NORTH", BORDER_DATA, TRUE, TestOnSection, NULL, MUST_BE_CALLED,
 "SOUTH", BORDER_DATA, TRUE, TestOnSection, NULL, MUST_BE_CALLED,
 "WEST", BORDER_DATA, TRUE, TestOnSection, NULL, MUST_BE_CALLED,
 "EAST", BORDER_DATA, TRUE, TestOnSection, NULL, MUST_BE_CALLED,
 "CHEMISTRY", INITIAL_DATA, TRUE, NULL, NULL, CALL_WHEN_SUBSTANCES,
 "EMISSIONS", CHEMISTRY, FALSE, SetGridToEmission, NULL, CAN_BE_OMITTED,
 "REDUCTION", EMISSIONS, FALSE, NULL, NULL, CAN_BE_OMITTED,
 "DEPOSITION", CHEMISTRY, FALSE, SetGridToEmission, NULL, CAN_BE_OMITTED,
 "NESTING", GRID, FALSE, NULL, EndOfNesting, CAN_BE_OMITTED,
 "OUTPUT", GRID, FALSE, NULL, NULL, CAN_BE_OMITTED};
 
typedef enum  {
  LENGTH_U, TIME_U, SPEED_U, TEMP_U, PRESS_U,
  POWER_P_AREA_U, CHEMISTRY_U, MIXING_U, ARC_U
}  UnitCat;
 
char unit[][6] = {
"m", "sec", "m/s", "K", "Pa", "W/m^2", "ppb", "kg/kg", "deg"};

#define MAXVARIABLE 142

GroundParam grmark;

VarDesc variable[MAXVARIABLE] =
{"U", "Wind in WE-Direction", unit[SPEED_U], GRID_VAL, (long *)UWIND, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
     SET_DEFAULT, 0, -120., 120., NORMAL_NUM, SET_DEFAULT << 8, 0, 0, NULL, NULL,
 "Ground-Interface", "Calculate Ground-Interactions? (on/off)", NULL, INT_PTR, (long *)&groundinterface,
     NO_DIM, OPTIONS, OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "Turbulence-Type", "Type of turbulence parametrisation", NULL, INT_PTR, (long *)&turbtype, NO_DIM,
     OPTIONS, OPTIONS, SET_DEFAULT, 2, 0, 2, TURBTYPE_OPTION, 0, 0, 0, NULL, NULL,
 "Shielding_Factor", "Shielding factor of plants", NULL, GROUND_PARAM, (long *)&grmark.sigf, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, MUST_BE_SET, 0, 0, 1, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Scalar-Advection", "Calculate Advection of Scalars? (on/off)", NULL, INT_PTR, (long *)&advection, NO_DIM,
     OPTIONS, OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "Print-Borders", "Print values at the border in the output-File? (on/off)", NULL, INT_PTR, (long *)&printwithborder,
     NO_DIM, OPTIONS, OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, 0, 0, 0, NULL, NULL,
 "Wind-Advection", "Calculate Advection of Wind Speed? (on/off)", NULL, INT_PTR, (long *)&windadvection, NO_DIM,
     OPTIONS, OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "Solar-Radiation", "Calculate the solar radiation (on/off)", NULL, INT_PTR, (long *)&shortwaveradiation, NO_DIM,
     OPTIONS, OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "Clouds", "Consider Clouds? (on/off)", NULL, INT_PTR, (long *)&cloudwater, NO_DIM,
     OPTIONS, OPTIONS, SET_DEFAULT, 0, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "Cloud-Ice", "Consider Cloud Ice and Snow? (on/off)", NULL, INT_PTR, (long *)&cloudice, NO_DIM,
     OPTIONS, OPTIONS, SET_DEFAULT, 0, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "Timed-Dump", "Make multiple dump-files with time stamp.", NULL, INT_PTR, (long *)&timeddump, NO_DIM,
     OPTIONS, OPTIONS, SET_DEFAULT, 0, 0, 1, ON_OFF, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "V", "Wind in SN-Direction", unit[SPEED_U], GRID_VAL, (long *)VWIND, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
     SET_DEFAULT, 0, -120., 120., NORMAL_NUM, SET_DEFAULT << 8, 0, 12, NULL, NULL,
 "W", "Wind in vertical direction", unit[SPEED_U], GRID_VAL, (long *)WWIND, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
     CALCULATED_VALUES, 0, -30., 30., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Uflux", "Uflux", unit[SPEED_U], LAYER_PTR, NULL, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
     CALCULATED_VALUES, 0, -120., 120., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Vflux", "Wind in SN-Direction", unit[SPEED_U], LAYER_PTR, NULL, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
     CALCULATED_VALUES, 0, -120., 120., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Wflux", "Wind in vertical direction", unit[SPEED_U], LAYER_PTR, NULL, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
     CALCULATED_VALUES, 0, -30., 30., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "TEMP", "virt. pot. Air Temperature", unit[TEMP_U], GRID_VAL, (long *)TEMP, ALL_DIM, INITIAL_DATA, INITIAL_DATA, MUST_BE_SET,
     290., 240., 400., NORMAL_NUM, MUST_BE_SET << 8, 0, 0, NULL, NULL,
 "Q", "Humidity of air", unit[MIXING_U], GRID_VAL, (long *)HUMIDITY, ALL_DIM, INITIAL_DATA, INITIAL_DATA, MUST_BE_SET,
    0.005, 0., 0.03, NORMAL_NUM, MUST_BE_SET << 8, 0, 0, NULL, NULL,
 "Qwat", "Cloud water", unit[MIXING_U], GRID_VAL, (long *)CLOUDWATER, ALL_DIM, INITIAL_DATA, INITIAL_DATA, SET_DEFAULT,
    0., 0., 0.03, NORMAL_NUM, SET_DEFAULT << 8, 0, 0, OnCloudWater, NULL,
 "Qrain", "Concentration of rain drops", unit[MIXING_U], GRID_VAL, (long *)RAINWATER, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
    SET_DEFAULT, 0., 0., 0.03, NORMAL_NUM, SET_DEFAULT << 8, 0, 0, OnCloudWater, NULL,
 "Qice", "Cloud Ice/Cloud Snow", unit[MIXING_U], GRID_VAL, (long *)CLOUDICE, ALL_DIM, INITIAL_DATA, INITIAL_DATA, SET_DEFAULT,
    0., 0., 0.03, NORMAL_NUM, SET_DEFAULT << 8, 0, 0, OnCloudIce, NULL,
 "TKE", "Turbulent kintetic energy", "Pa/m^3", GRID_VAL, (long *)TKE, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
    SET_DEFAULT, 0., 0., 1.e5, NORMAL_NUM, SET_DEFAULT << 8, 0, 0, OnKEpsTurb, NULL,
 "eps", "Dissipation", "???", GRID_VAL, (long *)EPS, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
    SET_DEFAULT, 0., 0., 1.e5, NORMAL_NUM, SET_DEFAULT << 8, 0, 0, OnKEpsTurb, NULL,
 "Tabs", "Absolute Air-Temperature", unit[TEMP_U], PROC_VAL, (long *)AbsolutePointTemp, ALL_DIM, INITIAL_DATA,
     INITIAL_DATA, CALCULATED_VALUES, 290., 240., 340., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Te", "equivalent potential Air-Temperature", unit[TEMP_U], PROC_VAL, (long *)EquivalentPointTemp, ALL_DIM, INITIAL_DATA,
     INITIAL_DATA, CALCULATED_VALUES, 290., 240., 340., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "H", "Relative Humidity", "%", PROC_VAL, (long *)RelativeHumidity, ALL_DIM, INITIAL_DATA,
     INITIAL_DATA, CALCULATED_VALUES, 290., 240., 340., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "AirConz", "Concentration of air", "molec/cm^3", PROC_VAL, (long *)AirConz, ALL_DIM, INITIAL_DATA,
     INITIAL_DATA, CALCULATED_VALUES, 290., 240., 340., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Vtr", "Terminal velocity of rain", unit[SPEED_U], PROC_VAL, (long *)TerminalVelocityOfRain, ALL_DIM, INITIAL_DATA,
     INITIAL_DATA, CALCULATED_VALUES, 0., 0., 50., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Start", "Time to start simulation", unit[TIME_U], LONG_PTR, &tstart, NO_DIM, TIME, TIME, SET_DEFAULT,
     0, 0, 1000000, INT_NUM, 0, 0, 0, NULL, NULL,
 "End", "Time to end simulation", unit[TIME_U], LONG_PTR, &tend, NO_DIM, TIME, TIME, MUST_BE_SET,
     0, 0, 100000000, INT_NUM, 0, 0, 0, NULL, NULL,
 "dt", "Time-Step of simulation", unit[TIME_U], LONG_PTR, &tincmax, NO_DIM, TIME, TIME, MUST_BE_SET,
     20, 1, 3600, INT_NUM, 0, 0, 0, NULL, NULL,
 "dt-gallery", "Possible-Values for dt", unit[TIME_U], INT_VECT_PTR, (long *)&taptr, MAXTINC, TIME, TIME, MUST_BE_SET,
     20, 1, 200, INT_VECTOR_VAL, 0, 0, 0, NULL, NULL,
 "dt-IR", "Time-Step for atmospheric IR-Radiation", unit[TIME_U], LONG_PTR, &irtinc, NO_DIM,
     GROUND_BORDER, GROUND_BORDER, MUST_BE_SET, 60, 0, 600, INT_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "dt-chem", "Time-Step of chemical species", unit[TIME_U], LONG_PTR, &tchem, NO_DIM, CHEMISTRY, CHEMISTRY,
     MUST_BE_SET, 20, 1, 600, INT_NUM, 0, 0, 0, OnChemistry, NULL,
 "photofact", "Factor to multiply to the photochemical rate constants", NULL, DOUBLE_PTR,
     (long *)&photofact, NO_DIM, CHEMISTRY, CHEMISTRY, SET_DEFAULT, 1., 0., 100., NORMAL_NUM, 0,
     0, 0, OnChemistry, NULL,
 "fully-implicit", "Calculate Chemistry fully implicit? (on/off)", NULL, INT_PTR, (long *)&allfast, NO_DIM,
     CHEMISTRY, CHEMISTRY, SET_DEFAULT, 0, 0, 1, ON_OFF, 0, 0, 0, OnChemistry, NULL,
 "dumptime", "Interval to make a memory dump", NULL, LONG_PTR, &dumptime, NO_DIM, TIME, TIME, SET_DEFAULT, 0,
 0, 1000000, INT_NUM, 0, 0, 0, NULL, NULL,
 "dumppath", "Path of directory for dump-files", NULL, CHAR_PTR, (long *)dumppath, NO_DIM, TIME, TIME, SET_DEFAULT,
     0, 0, 0, TEXT_VAL, 0, 0, 0, NULL, NULL,
 "Time", "Day-Time of start of simulation", unit[TIME_U], DOUBLE_PTR, (long *)&timeofday, NO_DIM, TIME, TIME,
     MUST_BE_SET, 0, 0, 24, TIME_VAL, 0, 0, 0, NULL, NULL,
 "Date", "Date of start of simulation", "Day of year (M-D-Y)", LONG_PTR, &dayofyear, NO_DIM, TIME, TIME,
     MUST_BE_SET, 0, 0, 366, DATE_VAL, 0, 0, 0, NULL, NULL,
 "Km", "Turbulence Parameter", NULL, LAYER_PTR, NULL, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
     CALCULATED_VALUES, 0, 0., 100., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "IRheat", "atmospheric radiative heating(+)/cooling(-)", "K/sec", LAYER_PTR, NULL, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
     CALCULATED_VALUES, 0, 0., 100., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "nx", "Number of Nodes in x (WE)-Direction", NULL, INT_PTR, (long *)&nx, NO_DIM, GRID, GRID, MUST_BE_SET,
     0, 1, 500, INT_NUM, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "ny", "Number of Nodes in y (SN)-Direction", NULL, INT_PTR, (long *)&ny, NO_DIM, GRID, GRID, MUST_BE_SET,
     0, 1, 500, INT_NUM, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "nz", "Number of Nodes in z (vertical)-Direction", NULL, INT_PTR, (long *)&nz, NO_DIM, GRID, GRID, MUST_BE_SET,
     0, 4, 100, INT_NUM, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "nsubs", "Number of chemical species", NULL, INT_PTR, (long *)&nsubs, NO_DIM, GRID, GRID, SET_DEFAULT,
     0, 0, 2000, INT_NUM, 0, 0, 0, NULL, NULL,
 "dx", "Grid Distance in the x-Direction", unit[LENGTH_U], DOUBLE_PTR, (long *)&dx, NO_DIM, GRID, GRID, MUST_BE_SET,
     0, 100, 40000, SINGLE_VAL, 0, 0, 0, NULL, NULL,
 "dy", "Grid Distance in the y-Direction", unit[LENGTH_U], DOUBLE_PTR, (long *)&dy, NO_DIM, GRID, GRID, MUST_BE_SET,
     0, 100, 40000, SINGLE_VAL, 0, 0, 0, NULL, NULL,
 "workers", "Number of workers in a parallel run", NULL, INT_PTR, (long *)&workerswanted, NO_DIM, GRID, GRID, SET_DEFAULT,
     0, 2, 64, SINGLE_VAL, 0, 0, 0, NULL, NULL,
 "topo", "topography. Normally height above sea level", unit[LENGTH_U], FILE_PTR, NULL,
     LAYER_DIM, ENVIRONMENT, ENVIRONMENT, SET_DEFAULT, 0, 0, 8000, NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "reflevel", "Deepest ground-level in model domain.", unit[LENGTH_U], DOUBLE_PTR, (long *)&reflevel,
     NO_DIM, ENVIRONMENT, ENVIRONMENT, MUST_BE_SET, 0, 0, 8000, NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "level", "Height of model layers above \"reflevel\"", unit[LENGTH_U], DOUBLE_PTR, NULL,
    Z_DIM, ENVIRONMENT, ENVIRONMENT, MUST_BE_SET, 0, 0, 15000, NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Longitude", "geographical Longitude (East of Greenwich)", unit[ARC_U], DOUBLE_PTR, (long *)&Xlong, NO_DIM, 
     ENVIRONMENT, ENVIRONMENT, MUST_BE_SET, 0, -180., 180., SINGLE_VAL, 0, 0, 0, OnGroundInterface, NULL,
 "Latitude", "Geographical Latitude", unit[ARC_U], DOUBLE_PTR, (long *)&Xlat,
     NO_DIM, ENVIRONMENT, ENVIRONMENT, MUST_BE_SET, 0, -90, 90, SINGLE_VAL, 0, 0, 0, OnCoriolis, NULL,
 "Zonediff", "Local-Time - UTC", "h", DOUBLE_PTR, (long *)&timezonediff, NO_DIM, ENVIRONMENT, ENVIRONMENT,
     MUST_BE_SET, 0, -24, 25, SINGLE_VAL, 0, 0, 0, OnGroundInterface, NULL,
 "MapAngle", "Angle between model-y and north (positive when model domain was turned to the left)",
     "deg", DOUBLE_PTR, (long *)&turnmapangle, NO_DIM, ENVIRONMENT, ENVIRONMENT,
     SET_DEFAULT, 0, -180., 180., SINGLE_VAL, 0, 0, 0, OnGroundInterface, NULL,
 "Stratospheric-Ozone", "Stratosperic ozone", "cm?", DOUBLE_PTR, (long *)&stratozone, NO_DIM,
     ENVIRONMENT, ENVIRONMENT, SET_DEFAULT, 0.4, 0.05, 0.8, SINGLE_VAL, 0, 0, 0, OnGroundInterface, NULL,
 "Precipitable-Water", "Precipitable Water above modelling Domain", "g/cm^2", DOUBLE_PTR,
     (long *)&topprecipwater, NO_DIM, ENVIRONMENT, ENVIRONMENT, SET_DEFAULT, 1., 0, 5, SINGLE_VAL, 0,
     0, 0, OnGroundInterface, NULL,
 "Turbidity", "Turbidity-Parameter (optimize manually!)", NULL, DOUBLE_PTR, (long *)&turbidity, NO_DIM,
     ENVIRONMENT, ENVIRONMENT, SET_DEFAULT, 0.08, 0, 0.8, SINGLE_VAL, 0, 0, 0, OnGroundInterface, NULL,
 "TopPressure", "Pressure on top of modelling domain", unit[PRESS_U], DOUBLE_PTR, (long *)&toppp, NO_DIM,
     ENVIRONMENT, ENVIRONMENT, SET_DEFAULT, 79000, 1000, 100000, SINGLE_VAL, 0, 0, 0, NULL, NULL,
 "Print-Residual", "Print convergence information of SOR-Algorithm (only in non-hydrostatic mode)", NULL, INT_PTR,
     (long *)&printresid, NO_DIM, OPTIONS, OPTIONS, SET_DEFAULT, 0, 0, 1, ON_OFF, 0, 0, 0, NULL, NULL,
 "Main-works-too", "In parallel mode: one worker is placed on the main host", NULL, INT_PTR,
     (long *)&workeronthishost, NO_DIM, OPTIONS, OPTIONS, SET_DEFAULT, 1, 0, 1, ON_OFF, 0, 0, 0, NULL, NULL,
 "Damping-Layer", "Use a four damping layers at the top", NULL, INT_PTR,
     (long *)&dampinglayer, NO_DIM, OPTIONS, OPTIONS, SET_DEFAULT, 0, 0, 1, ON_OFF, 0, 0, 0, NULL, NULL,
 "Coriolis", "How to consider Coriolis Forces (none/full/differential)", NULL, INT_PTR, (long *)&coriolistype, NO_DIM,
     OPTIONS, OPTIONS, SET_DEFAULT, (double)GEOSTROPHICCORIOLIS, 0, 8, CORIOLIS_OPTION, 0, 0, 0, NULL, NULL,
 "Advection-Type", "How to calculate Advection (mpdata/ppm)", NULL, INT_PTR, (long *)&advectiontype, NO_DIM,
     OPTIONS, OPTIONS, SET_DEFAULT, (double)MPDATA_A, 0, 8, ADVECTIONTYPE_OPTION, 0, 0, 0, NULL, NULL,
 "iord", "Iteration order of Smolarkievicz Advection algorithm", NULL, INT_PTR, (long *)&smiord, NO_DIM,
     OPTIONS, OPTIONS, SET_DEFAULT, 1., 0, 5, INT_NUM, 0, 0, 0, OnAdvection, NULL,
 "No-Oscillating", "Use Non-oscillating option in Smolarkievicz Advection algorithm (on/off)", NULL, INT_PTR, (long *)&smnonos,
     NO_DIM, OPTIONS, OPTIONS, SET_DEFAULT, 0., 0, 1, ON_OFF, 0, 0, 0, OnAdvection, NULL,
 "PressureAlg", "Indicates the method of pressure evaluation", NULL, INT_PTR,
     (long *)&pressuretype, NO_DIM, OPTIONS, OPTIONS, SET_DEFAULT, (double)HYDROSTATIC, 0, 2, PRESSUREALG_OPTION,
     0, 0, 0, NULL, NULL,
 "Spatial-Filter-Type", "The type of the spatial low-pass filter to apply", NULL, INT_PTR, (long *)&filtertype,
     NO_DIM, OPTIONS, OPTIONS, SET_DEFAULT, (double)NO_FILTER, 0, 3, FILTERTYPE_OPTION, 0, 0, 0, NULL, NULL,
 "Spatial-Filter-Val", "Applies a spatial filter after each time-step with damping value", NULL, DOUBLE_PTR,
     (long *)&spatialfilter, NO_DIM, OPTIONS, OPTIONS, SET_DEFAULT, 0.001, 0, 1, SINGLE_VAL, 0, 0, 0, NULL, NULL,
 "dPress", "p' = (p - p0)", unit[PRESS_U], LAYER_PTR, NULL, ALL_DIM, INITIAL_DATA, INITIAL_DATA,
     CALCULATED_VALUES, 0, -100., 100., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "density", "density(z)", "kg/m^3", DOUBLE_PTR, NULL, Z_DIM, INITIAL_DATA, INITIAL_DATA,
     CALCULATED_VALUES, 0, 0, 100, NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "P0", "middle air pressure", "Pa", DOUBLE_PTR, NULL, Z_DIM, INITIAL_DATA, INITIAL_DATA,
     CALCULATED_VALUES, 0, 60000, 130000, NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "delta_TopPressure", "p' in the top Layer (Pa)", unit[PRESS_U], FILE_PTR, NULL, LAYER_DIM,
     TOP_BORDER, TOP_BORDER, SET_DEFAULT, 0, -2000, 2000, NORMAL_NUM, STARTING_WITH_ZERO, 0, 0, NULL, NULL,
 "TopTemp", "virtual-potential Temparature above domain", unit[TEMP_U], FILE_PTR, NULL, LAYER_DIM,
     TOP_BORDER, TOP_BORDER, MUST_BE_SET, 0, 0., 50., NORMAL_NUM, STARTING_WITH_ZERO, 0, 0, NULL, NULL,
 "Ugeos", "Geostrophic Wind (U)", unit[SPEED_U], DOUBLE_PTR, NULL, Z_DIM,
     TOP_BORDER, TOP_BORDER, MUST_BE_SET, 0, -100, 100, NORMAL_NUM, 0, 0, 0, OnCoriolis, NULL,
 "Vgeos", "Geostrophic Wind (V)", unit[SPEED_U], DOUBLE_PTR, NULL, Z_DIM,
     TOP_BORDER, TOP_BORDER, MUST_BE_SET, 0, -100, 100, NORMAL_NUM, 0, 0, 0, OnCoriolis, NULL,
 "T_Surface", "Temperature at the ground surface (K)", unit[TEMP_U], GROUND_PARAM, (long *)&grmark.Tg[0],
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, MUST_BE_SET, 0, 230, 370, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "T_Gr1", "Temperature at 1st Ground-Level (4cm below Surface)", unit[TEMP_U], GROUND_PARAM, (long *)&grmark.Tg[1],
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 230, 370, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "T_Gr2", "Temperature at 2nd Ground-Level (16cm below Surface)", unit[TEMP_U], GROUND_PARAM, (long *)&grmark.Tg[2],
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 230, 370, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "T_Gr3", "Temperature at 3rd Ground-Level (36cm below Surface)", unit[TEMP_U], GROUND_PARAM, (long *)&grmark.Tg[3],
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 230, 370, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "T_Gr4", "Temperature at 4th Ground-Level (64cm below Surface)", unit[TEMP_U], GROUND_PARAM, (long *)&grmark.Tg[4],
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 230, 370, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "T_Ground", "Temperature at 1m below ground surface (K)", unit[TEMP_U], GROUND_PARAM, (long *)&grmark.
     Tg[NGROUNDLAYER-1], LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 285, 240, 350, NORMAL_NUM, 0,
     0, 0, OnGroundInterface, NULL,
 "Ground_U", "U-Wind in Ground-Layer", unit[SPEED_U], GROUND_PARAM, (long *)&grmark.a[UWIND],
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, -50., 50., NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Ground_V", "V-Wind in Ground-Layer", unit[SPEED_U], GROUND_PARAM, (long *)&grmark.a[VWIND],
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, -50., 50., NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Ground_TEMP", "Temperature in Ground-Layer", unit[TEMP_U], GROUND_PARAM, (long *)&grmark.a[TEMP],
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 230., 330., NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Ground_Q", "Humidity in Ground-Layer", unit[MIXING_U], GROUND_PARAM, (long *)&grmark.a[HUMIDITY],
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Rain_Intensity", "Rain Intensity at Ground", "mm/sec", GROUND_PARAM, (long *)&grmark.rain_intensity,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, 0, 0, OnCloudWater, NULL,
 "Cumulated_Rain", "Cumulated rain since start of model run", "mm", GROUND_PARAM, (long *)&grmark.cum_rain,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, 0, 0, OnCloudWater, NULL,
 "Snow_Intensity", "Rain Intensity at Ground", "mm/sec", GROUND_PARAM, (long *)&grmark.snow_intensity,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, 0, 0, OnCloudWater, NULL,
 "Cumulated_Snow", "Cumulated snow since start of model run", "mm (Water equivalent)", GROUND_PARAM,
     (long *)&grmark.cum_snow, LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, 
     NORMAL_NUM, 0, 0, 0, OnCloudWater, NULL,
 "TotalClouds", "Total cloud water and ice content in vertical pile", "mm", PROC_VAL, (long *)TotalClouds,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, 0, 0, OnCloudIce, NULL,
 "TotalCloudWater", "Total cloud water content in vertical pile", "mm", PROC_VAL, (long *)TotalWaterClouds,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, 0, 0, OnCloudWater, NULL,
 "TotalCloudIce", "Total cloud ice content in vertical pile", "mm", PROC_VAL, (long *)TotalIceClouds,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0, 0., 0.02, NORMAL_NUM, 0, 0, 0, OnCloudIce, NULL,
 "T_Foliage", "Temperature of plant surface", unit[TEMP_U], GROUND_PARAM, (long *)&grmark.Tf, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., 220., 350., NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Qg", "Ground-Flux", unit[POWER_P_AREA_U], GROUND_PARAM, (long *)&grmark.Qg, LAYER_DIM, GROUND_BORDER,
     GROUND_BORDER, CALCULATED_VALUES, 0., -500., 500., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Rel_Humidity", "Relative Humidity of the air in the ground-pores", unit[MIXING_U], GROUND_PARAM, (long *)
     &grmark.X, LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 1, 0, 1, NORMAL_NUM, 0,
     0, 0, OnGroundInterface, NULL,
 "Evaporation_Resistance", "Evaporation-Resistance of soil", NULL, GROUND_PARAM, (long *)&grmark.r,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 1.e-3, 0., 1.e-2, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "GroundClass", "Landusage class at cell, according to classlist in input file", NULL, GROUND_PARAM,
     (long *)&grmark.groundclass,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 0, 0, 1e4, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Ground_Heat_Capacity", "Ground-Heat-Capacity", "J/(m^3*)", GROUND_PARAM, (long *)&grmark.Cg,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 2.5e6, 1e4, 1e8, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Ground_Diffusivity", "ks (=Ground-Diffusivity)", NULL, GROUND_PARAM, (long *)&grmark.ks,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 0.53e-6, 1e-8, 5e-3, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Albedo", "Albedo of soil", NULL, GROUND_PARAM, (long *)&grmark.albedo, LAYER_DIM, GROUND_BORDER, GROUND_BORDER,
     SET_DEFAULT, 0.2, 0, 1, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Ground_Emissivity", "Ground-Emissivity", NULL, GROUND_PARAM, (long *)&grmark.emissity, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 0.98, 0.4, 1, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Z0", "Ground roughness", unit[LENGTH_U], GROUND_PARAM, (long *)&grmark.z0, LAYER_DIM, GROUND_BORDER,
     GROUND_BORDER, SET_DEFAULT, 1, 1e-6, 30, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Zt", "Height of first Point above Ground", unit[LENGTH_U], GROUND_PARAM, (long *)&grmark.z, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 10, 3, 25, NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Leaf_Area_Index", "Leaf-Area-Index", NULL, GROUND_PARAM, (long *)&grmark.La, LAYER_DIM, GROUND_BORDER,
     GROUND_BORDER, SET_DEFAULT, 4, 0., 20, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Plant_Albedo", "Albedo of plants", NULL, GROUND_PARAM, (long *)&grmark.albedof, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 0.15, 0.01, 0.8, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Plant_Emissivity", "Emissivity of plants", NULL, GROUND_PARAM, (long *)&grmark.emissityf, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 0.95, 0.4, 1., NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Plant_Resistance", "Resistance of plants against evapotranspiration", NULL, GROUND_PARAM, (long *)
     &grmark.rcf, LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 150., 5., 500., NORMAL_NUM,
     0, 0, 0, OnGroundInterface, NULL,
 "Plant_Water_Content", "Water content of plants", NULL, GROUND_PARAM, (long *)&grmark.fwatercontent,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 0.6, 0.01, 1, NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "Rglobal", "Global Radiation", unit[POWER_P_AREA_U], GROUND_PARAM, (long *)&grmark.Rg, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 1200., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Rglobal_cumulated", "Global Radiation integrated for the actual day", "J/m^2", GROUND_PARAM, (long *)&grmark.Rgcum, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 1200000., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Rnet", "Net radiation", unit[POWER_P_AREA_U], GROUND_PARAM, (long *)&grmark.R, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 1000., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "IRdown", "IR irradiation", unit[POWER_P_AREA_U], GROUND_PARAM, (long *)&grmark.IRdown, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 1000., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Sun_Elevation_Angle", "Elevation angle of the sun (0 at night-time)", "deg", GROUND_PARAM,
     (long *)&sunelevation, TIME_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 90.,
     NORMAL_NUM, 0, 0, 0, OnGroundInterface, NULL,
 "wtheta", "sensible heat flux", "m*K/sec", GROUND_PARAM, (long *)&grmark.wtheta, LAYER_DIM,
     GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., -100., 100., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "wq", "latent heat flux", "m*kg/(sec*kg)", GROUND_PARAM, (long *)&grmark.wq, LAYER_DIM, GROUND_BORDER,
     GROUND_BORDER, CALCULATED_VALUES, 0., -1., 1., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "zL", "z/L, where L = monin Obukhov Length.", unit[LENGTH_U], GROUND_PARAM, (long *)&grmark.zL,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., -100., 100., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "ustar", "u* the friction velocity", unit[SPEED_U], GROUND_PARAM, (long *)&grmark.ustar,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., -100., 100., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "phim", "a Businger Turbulence Parameter", NULL, GROUND_PARAM, (long *)&grmark.phim,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., -100., 100., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "phit", "a Businger Turbulence Parameter", NULL, GROUND_PARAM, (long *)&grmark.phit,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., -100., 100., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "BLH", "Height of boundary layer", NULL, GROUND_PARAM, (long *)&grmark.blh,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 10000., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "R0_max", "Deposition resistance maximum", "s/m", GROUND_PARAM, (long *)&grmark.Drmax,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 400., 0., 5000., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "R0_min", "Deposition resistance minimum (must be 0 on water surfaces!)", "s/m", GROUND_PARAM,
     (long *)&grmark.Drmin,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 300., 0., 5000., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "R0_night", "Deposition resistance in the night", "s/m", GROUND_PARAM, (long *)&grmark.Drnight,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 300., 0., 5000., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "R0_wet", "Deposition resistance on a wet surface", "s/m", GROUND_PARAM, (long *)&grmark.Drnight,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 300., 0., 5000., NORMAL_NUM, 0, 0, 0, NULL, NULL,
#ifdef RADFACT
 "Rdir_Fact", "Factor for direkt Short-Wave Radiation", NULL, LAYER_PTR, NULL, ALL_DIM,
     GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 1., 0., 3., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Rdiff_Fact", "Factor for diffuse Short-Wave Radiation", NULL, LAYER_PTR, NULL, ALL_DIM,
     GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 1., 0., 6., NORMAL_NUM, 0, 0, 0, NULL, NULL,
#else
 "Rdir_Fact", "Factor for direkt Short-Wave Radiation", NULL, DOUBLE_PTR, NULL, Z_DIM,
     GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 1., 0., 3., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Rdiff_Fact", "Factor for diffuse Short-Wave Radiation", NULL, DOUBLE_PTR, NULL, Z_DIM,
     GROUND_BORDER, GROUND_BORDER, SET_DEFAULT, 1., 0., 6., NORMAL_NUM, 0, 0, 0, NULL, NULL,
#endif
 "Rdir", "Direct Short-wave Radiation", unit[POWER_P_AREA_U], GROUND_PARAM, (long *)&grmark.Rdir,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 0., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "Rdiff", "Diffuse Short-wave Radiation", unit[POWER_P_AREA_U], GROUND_PARAM, (long *)&grmark.Rdiff,
     LAYER_DIM, GROUND_BORDER, GROUND_BORDER, CALCULATED_VALUES, 0., 0., 0., NORMAL_NUM, 0, 0, 0, NULL, NULL,
 "BorderType-West", "Type of Borderconditions", NULL, INT_PTR, (long *)&westbordertype, NO_DIM, BORDER_DATA,
     BORDER_DATA, SET_DEFAULT, (double)SPONGEBORDER, 0, 64, BORDERTYPE_OPTION, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "BorderType-East", "Type of Borderconditions", NULL, INT_PTR, (long *)&eastbordertype, NO_DIM, BORDER_DATA,
     BORDER_DATA, SET_DEFAULT, (double)SPONGEBORDER, 0, 64, BORDERTYPE_OPTION, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "BorderType-North", "Type of Borderconditions", NULL, INT_PTR, (long *)&northbordertype, NO_DIM, BORDER_DATA,
     BORDER_DATA, SET_DEFAULT, (double)SPONGEBORDER, 0, 64, BORDERTYPE_OPTION, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 "BorderType-South", "Type of Borderconditions", NULL, INT_PTR, (long *)&southbordertype, NO_DIM, BORDER_DATA,
     BORDER_DATA, SET_DEFAULT, (double)SPONGEBORDER, 0, 64, BORDERTYPE_OPTION, NO_ERRORS_ALLOWED, 0, 0, NULL, NULL,
 nestvarname[0], "x-coordinates of origin of nested subdomain", unit[LENGTH_U], LONG_PTR,
     (long *)&dummynestvar, NO_DIM, NESTING, NESTING, MUST_BE_SET, 0, 0, 500000, INT_NUM, 0,
     0, 0, OnNestDomain, NULL,
 nestvarname[1], "y-coordinates of origin of nested subdomain", unit[LENGTH_U], LONG_PTR,
     (long *)&dummynestvar, NO_DIM, NESTING, NESTING, MUST_BE_SET, 0, 0, 500000, INT_NUM, 0,
     0, 0, OnNestDomain, NULL,
 nestvarname[2], "number of points along x of nested subdomain", NULL, LONG_PTR,
     (long *)&dummynestvar, NO_DIM, NESTING, NESTING, MUST_BE_SET, 0, 0, 1000, INT_NUM, 0,
     0, 0, OnNestDomain, NULL,
 nestvarname[3], "number of points along y of nested subdomain", NULL, LONG_PTR,
     (long *)&dummynestvar, NO_DIM, NESTING, NESTING, MUST_BE_SET, 0, 0, 1000, INT_NUM, 0,
     0, 0, OnNestDomain, NULL,
 nestvarname[4], "grid cell size along x of nested subdomain", unit[LENGTH_U], LONG_PTR,
     (long *)&dummynestvar, NO_DIM, NESTING, NESTING, MUST_BE_SET, 0, 0, 10000, INT_NUM, 0,
     0, 0, OnNestDomain, NULL,
 nestvarname[5], "grid cell size along y of nested subdomain", unit[LENGTH_U], LONG_PTR,
     (long *)&dummynestvar, NO_DIM, NESTING, NESTING, MUST_BE_SET, 0, 0, 10000, INT_NUM, 0,
     0, 0, OnNestDomain, NULL,
 nestvarname[6], "turn-angle of nested subdomain", "deg", DOUBLE_PTR,
     (long *)&dummynestvar, NO_DIM, NESTING, NESTING, MUST_BE_SET, 0, -360, 360, NORMAL_NUM, 0,
     0, 0, OnNestDomain, NULL,
 nestvarname[7], "printing interval", unit[TIME_U], LONG_PTR,
     (long *)&dummynestvar, NO_DIM, NESTING, NESTING, MUST_BE_SET, 0, 0, 28800, INT_NUM, 0,
     0, 0, OnNestDomain, NULL,
};

ValueDesc val;
InputSection actualsection;
Date startdate;
Time starttime;

mcierror(char *s)
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
      printf("----------- %s -----------\n", section[actualsection].name);
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
       return (variable[1].name);
    case VAL_OF_RELVAR :
       return (groundinterface ? "on" : "off");
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
       return (cloudwater ? "on" : "off");
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
       return (cloudice ? "on" : "off");
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
       return ((advection && windadvection) ? "on" : "off");
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
       return ((nsubs > 0) ? "on" : "off");
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
       return ((nsubs > 0) ? "on" : "off");
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
    case LAYER_PTR :
       for (k = setv->zmin; k <= setv->zmax; k++)
         for (i = setv->xmin; i <= setv->xmax; i++)
           for (j = setv->ymin; j <= setv->ymax; j++)
             v->v.d[k*layer+i*row+j] = dval;
       break;
    case FILE_PTR :
       switch (v->dims)  {
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
         default : InputError(ERROR, "Internal implementation error (2)");
                   break;
       }
       break;
    case XFILE_PTR :
       for (k = setv->zmin; k <= setv->zmax; k++)
         for (i = setv->xmin; i <= setv->xmax; i++)
           v->v.d[k*xrow+i] = dval;
       break;
    case DOUBLE_PTR :
       switch (v->dims)  {
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
         default : InputError(ERROR, "Internal Error");
           	   return (1);
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
  vh->storetype = FILE_PTR;
  switch (actualsection)  {
    case NORTH_BORDER : vh->dims = X_DIM | Z_DIM;
			vh->v.d = northborder+vh->v.et*nz*xrow;
			vh->storetype = XFILE_PTR;
			vh->init = vh->option >> 8;
			break;
    case SOUTH_BORDER : vh->dims = X_DIM | Z_DIM;
			vh->v.d = southborder+vh->v.et*nz*xrow;
			vh->storetype = XFILE_PTR;
			vh->init = vh->option >> 8;
			break;
    case WEST_BORDER  : vh->dims = Y_DIM | Z_DIM;
			vh->v.d = westborder+vh->v.et*nz*row;
			vh->storetype = FILE_PTR;
			vh->init = vh->option >> 8;
			break;
    case EAST_BORDER  : vh->dims = Y_DIM | Z_DIM;
			vh->v.d = eastborder+vh->v.et*nz*row;
			vh->storetype = FILE_PTR;
			vh->init = vh->option >> 8;
			break;
    default	      : return (v);
  }
  return (vh);
}

int ChangeSectionTo(char *sname)
{
  InputSection i;
  VarDesc *v, vh;
  static int sections_found;
  if (!sname)  {
    actualsection = NULL_SECTION;
    sections_found = 1<<NULL_SECTION;
    return (0);
  }
  for (i = 0; strcmp(sname, section[i].name); i++);
  if (!(sections_found & (1 << section[i].dependson)))  {
    InputError(ERROR, "Section %s depends on section %s, which must be defined first.",
       section[i].name, section[section[i].dependson].name);
    return (1);
  }
  if (section[i].startproc && section[i].startproc(actualsection, i))  return (1);
  actualsection = i;
  sections_found |= (1 << i);
  section[i].calltype = WAS_CALLED;
  InputError(ERROR, NULL);
  if (section[i].initialize)
    for (v = variable; v; v = v->next)  {
      if (v->section == actualsection && v->init == SET_DEFAULT)
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
      v->init = v->option >> 8;
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
      v->init = v->option >> 8;
      v->relieson = NULL;
    }
  return (0);
}

int TestOnSection(InputSection oldsec, InputSection newsec)
{
  if (oldsec < BORDER_DATA || oldsec > EAST_BORDER)  {
    InputError(ERROR, "Section \"%s\" must be called within \"BORDER-DATA\"", section[newsec]);
    return (1);
  }
  if (newsec >= NORTH_BORDER && newsec <= EAST_BORDER)
    SetGridToSect(oldsec, newsec);
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
}

int AddArrayVal(double num)
{
  if (val.v.a.nval >= val.v.a.maxval)  {
    InputError(ERROR, "Too many initialisers. You only need %i", val.v.a.maxval);
    return (1);
  }
  val.v.a.a[val.v.a.nval++] = num;
  return (0);
}

int AddIntArrayVal(long num)
{
  if (val.v.a.nval >= val.v.a.maxval)  {
    InputError(ERROR, "Too many initialisers. You only need %i", val.v.a.maxval);
    return (1);
  }
  ((long *)val.v.a.a)[val.v.a.nval++] = num;
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

int SetVariable(char *name)
{
  VarDesc *v, *vc;
  static VarDesc vh;
  BOOL rval;
  for (v = variable; v && strcmp(name, v->name); v = v->next);
  if (!v)  {
    InputError(ERROR, "%s is not a known variable or keyword.", name);
    return (1);
  }
  vc = v;
  if (v->storetype == GRID_VAL)  {
    if (actualsection == DEPOSITION)
      vc = DepositionVariable(v, &vh);
    else if (actualsection == EMISSIONS)
      vc = EmissionVariable(v, 0, 0, NULL);
    else if (actualsection == REDUCTION)
      vc = EmissionVariable(v, 0, 0, NULL);
    else if (actualsection >= NORTH_BORDER && actualsection <= EAST_BORDER)
      vc = BorderVariable(v, &vh);
  }
  if (vc->section != actualsection)  {
    InputError(ERROR, "The variable %s was found in the wrong section. It belongs to %s",
       name, section[v->section].name);
    return (1);
  }
  rval = !(v->option & STARTING_WITH_ZERO);
  setvar.vd = v; setvar.has_range = setvar.zmin = setvar.zmax = 0;
  setvar.xmin = setvar.xmax = setvar.ymin = setvar.ymax = rval;
  if (vc->inputtype & (VECTOR_VAL | INT_VECTOR_VAL))  {
    val.v.a.nx = val.v.a.ny = val.v.a.nz = val.v.a.has_dimension = val.v.a.nval = val.v.a.has_range = 0;
    val.v.a.maxval = vc->dims;   /* dims bezeichnet bei VECTOR_VAL die Anzahl der Elemente. */
    val.v.a.a = (vc->storetype == INT_VECT_PTR ? (double *)vc->v.ivp->v : vc->v.d);
    val.v.a.deletable = FALSE;
  }
  else  {
    if (v->dims & X_DIM)  setvar.xmax = nx - rval;
    if (v->dims & Y_DIM)  setvar.ymax = ny - rval;
    if (v->dims & Z_DIM)  setvar.zmax = nzm;
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
  if (v->storetype == GRID_VAL && actualsection >= NORTH_BORDER &&
      actualsection <= EAST_BORDER)  {
    if (val.type == AVG_VAL)  {
      if (val.v.ival > nz || val.v.ival < 1)  {
        InputError(ERROR, "Value %d is out of range (1-%d)", val.v.ival, nz);
        inputerror = TRUE;
      }
      switch (actualsection)  {
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
  if (v->storetype == GRID_VAL && actualsection == EMISSIONS)
    v = EmissionVariable(v, 0, 0, NULL);
  if (v->storetype == GRID_VAL && actualsection == REDUCTION)
    v = EmissionVariable(v, 0, 0, NULL);
  if (v->storetype == GRID_VAL && actualsection == DEPOSITION)
    v = DepositionVariable(v, &vh);
#ifdef PARALLEL
  if (parallel && v->storetype == GRID_VAL && gactive[v->v.et] && !g[v->v.et] &&
      !(g[v->v.et] = calloc(mesh, sizeof(double))))  {
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
      case LAYER_PTR :
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
      case FILE_PTR :
         switch (v->dims)  {
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
           default : InputError(ERROR, "Internal implementation error (2)");
                     break;
         }
         break;
      case XFILE_PTR :
         for (k = setvar.zmin; k <= setvar.zmax; k++)  {
           koff = ifa * (k - setvar.zmin) * ((val.v.a.has_dimension & Z_DIM) > 0);
           for (i = setvar.xmin; i <= setvar.xmax; i++)  {
             v->v.d[k*xrow+i] =
               val.v.a.a[koff + (i - setvar.xmin) * ((val.v.a.has_dimension & X_DIM) > 0)];
           }
         }
         break;
      case DOUBLE_PTR :
         switch (v->dims)  {
           case X_DIM : for (i = setvar.xmin; i <= setvar.xmax; i++)
           		  v->v.d[i] = val.v.a.a[i - setvar.xmin];
           		break;
           case Y_DIM : for (j = setvar.ymin; j <= setvar.ymax; j++)
           		  v->v.d[j] = val.v.a.a[j - setvar.ymin];
           		break;
           case Z_DIM : for (k = setvar.zmin; k <= setvar.zmax; k++)
           		  v->v.d[k] = val.v.a.a[k - setvar.zmin];
           		break;
           default : InputError(ERROR, "Internal Error");
           	      goto error_end;
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
      case PRESSUREALG_OPTION	: ival = val.v.ival; dval = (double)val.v.ival; break;
      case TIME_VAL		: dval = val.v.time.hour + val.v.time.minute / 60.; 
      				  starttime = val.v.time;
      				  break;
      case DATE_VAL		: ival = JulianDay(val.v.date.day, val.v.date.month, val.v.date.year);
      				  startdate = val.v.date;
     				  dval = (double)ival; break;
      case AVG_VAL		: break;
      case VECTOR_VAL		: if (val.v.a.nval < val.v.a.maxval)  {
      				    InputError(ERROR, "Not enough vector elements (%i). Expected %i.",
      				    	       val.v.a.nval, val.v.a.maxval);
      				    return (1);
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
  struct UVarDesc *next;
}  UVarDesc;

UVarDesc *uvar;

int DefineAVariable(char *name)
{
  UVarDesc *new;
  VarDesc *v;
  for (v = variable; v && strcmp(v->name, name); v = v->next);
  if (v)  {
    InputError(ERROR, "\"%s\" is fixed. It can not be redefined by the user.", name);
    return (1);
  }
  for (new = uvar; new && strcmp(name, new->name); new = new->next);
  if (new)  {
    InputError(ERROR, "There is a User-Variable \"%s\" already defined", name);
    return (1);
  }
  new = (UVarDesc *)malloc(sizeof(UVarDesc));
  memcpy(&new->vd, &val, sizeof(ValueDesc));
  strcpy(new->name, name);
  if (val.type == ARRAY_VAL)
    new->vd.v.a.deletable = FALSE;
  new->next = uvar;
  uvar = new;
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
    parallel = nx / (OVERLAP+4) > 1;
  if (!parallel || !master)  {
#endif
    for (i = 0; i < maxentity && allocok; i++)
      if (i < CLOUDWATER || i > EPS)
        allocok = !!(g[i] = (double *)calloc(mesh, sizeof(double)));
      allocok &= 
         !!(flux[UWIND] = (double *)calloc(mesh, sizeof(double))) &&
         !!(flux[VWIND] = (double *)calloc(mesh, sizeof(double))) &&
         !!(flux[WWIND] = (double *)calloc(mesh, sizeof(double)));
      press += layer;
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
      !(press = (double *)calloc(mesh+layer, sizeof(double))) ||
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
  PutToVarTable(level, "level");
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
  double height, toplevel, min, max;
  Entity et;
  /* Check plausibility of levels */
  toplevel = level[nzm];
  for (i = 0; i < nzm; i++)  {
    if (level[i] > level[i+1])  {
      InputError(ERROR, "levels must be given in increasing order!\n"
      	   "(I found %lf - %lf)", level[i], level[i+1]);
      inputerror = 1;
      return (1);
    }
    if (level[i] + 5. > level[i+1])  {
      InputError(ERROR, "Difference between levels should at least be 5m!\n"
      	   "(I found %lf - %lf = %lfm)", level[i], level[i+1], level[i+1] - level[i]);
      inputerror = 1;
      return (1);
    }
  }
  memcpy(zfaces, level, nz * sizeof(double));
/* Round topography to quads */
/*  ShapiroFilter(topo, pstat+nzm*layer, 3, 1.);  */
  for (i = layer; --i; )
    if (topo[i] < reflevel)  {
      topo[i] = reflevel;
      ground[i].firstabove = 0;
    }
    else  {
      if (topo[i] > level[nzm-1] + reflevel)  {
	topo[i] = level[nzm-1] + reflevel;
	InputError(WARNING, "Point %i/%i is higher than second highest level, I need one free line above topography!. I set point to %lf\n", i / row, i % row, topo[i]);
      }
      height = topo[i] - reflevel;
      j = 0;
      if (2 * height < level[0])  topo[i] = reflevel;
      else  {
        for ( ; j < nzm && fabs(level[j] - height) > fabs(level[j+1] - height); j++);
        topo[i] = level[j] + reflevel;
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
    height = level[k] + reflevel - 1.;
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

void TestOnSections()
{
  int i;
  for (i = 0; i < NULL_SECTION; i++)
    if (section[i].calltype == MUST_BE_CALLED)  {
      fprintf(stderr, "ERROR: Unable to find section %s.\n", section[i].name);
      inputerror = 1;
    }
    else if (section[i].calltype == CALL_WHEN_SUBSTANCES && nsubs > 0)  {
      fprintf(stderr, "ERROR: Unable to find section %s, although there are %i Substances defined\n",
         section[i].name, nsubs);
      inputerror = 1;
    }
  if (section[CHEMISTRY].calltype == WAS_CALLED && tchem % tinc)  {
    fprintf(stderr, "ERROR: tchem must be a multiple of tinc.");
    inputerror = 1;
  }
}

int EndOfSection()
{
  int error;
  VarDesc *v;
  error = 0;
  /*  if (actualsection == GROUND_BORDER)
      FillGroundValues();   This seems to be at a dangerous wrong position here! */
  if (section[actualsection].initialize)
    for (v = variable; v; v = v->next)
      if (v->section == actualsection)  {
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
  return (error || section[actualsection].endproc && section[actualsection].endproc());
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
    }
  toplevel = level[nzm];
  memcpy(zfaces, level, nz * sizeof(double));
  for (i = nzm; i--; )
    level[i+1] -= level[i];
  for (i = nz; i--; )  {
    zcenter[i] = zfaces[i] - level[i] * 0.5;
    leveli[i] = 1. / level[i];
    ldiff[i] = 0.5 * (level[i] + (i ? level[i-1] : 0.));
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
}

int mciwrap()
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

int ParseInput(char *fname)
{
  extern FILE *mciin;
  if (!(mciin = fopen(fname, "r")))  {
    fprintf(stderr, "FATAL ERROR: Unable to open File \"%s\"\n", fname);
    return (2);
  }
  inputerror = 0;
  inplineno = 1;
  uvar = NULL;
  setvar.isvalid = FALSE;
  ChangeSectionTo(NULL);
  mciparse();
  fclose(mciin);
  ClearAllUserVars();
  TestOnSections();
  TestDataConsistency();
/*  if (!inputerror)
    InitializeData();  */
  return (inputerror);
}

