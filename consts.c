/* November 28, 4057 B.C.*/
/* #define STARTDATE  239935.5 */

#include "kep.h"
#include "consts.h"

/* Conversion factors between degrees and radians */
double DTR = 1.7453292519943295769e-2;
double RTD = 5.7295779513082320877e1;
double RTS = 2.0626480624709635516e5; /* arc seconds per radian */
double STR = 4.8481368110953599359e-6; /* radians per arc second */
double PI = 3.14159265358979323846;
double TPI = 6.28318530717958647693;

extern double PI;

/* Standard epochs.  Note Julian epochs (J) are measured in
 * years of 365.25 days.
 */
double J2000 = 2451545.0;	/* 2000 January 1.5 */
double B1950 = 2433282.423;	/* 1950 January 0.923 Besselian epoch */
double J1900 = 2415020.0;	/* 1900 January 0, 12h UT */

double sqrt(), asin(), log();

/* coordinates of object
 */
int objnum = 0;	/* I.D. number of object */

/* ecliptic polar coordinates:
 * longitude, latitude in radians
 * radius in au
 */
double FAR obpolar[3];


struct orbit objorb = {
"Object         ",
2446800.5,
0.0,
0.0,
102.884,
0.999999,
0.985611,
0.016713,
1.1791,
2446800.5,
-3.86,
0.0,
0, /* &ear404, */
0.0,
0.0,
0.0
};

double appobj[3];


/* coordinates of Earth
 */
/* Heliocentric rectangular equatorial position
 * of the earth at time TDT re equinox J2000
 */
double FAR rearth[3];
double vearth[3];

/* Corresponding polar coordinates of earth:
 * longitude and latitude in radians, radius in au
 */
double FAR eapolar[3];

extern struct plantbl ear404;

/* Julian date of ephemeris */
double JD;
double TDT;
double UT;

/* flag = 0 if TDT assumed = UT,
 *      = 1 if input time is TDT,
 *      = 2 if input time is UT.
 */
int jdflag = 0;

/* correction vector, saved for display  */
double dp[3];

/* display enable flag */
int prtflg = 0;

double dradt, ddecdt;

double zgetdate(), gethms();
double func(), search();
int apparent();

struct orbit *elobject;
double robject[3] = {0.0, 0.0, 0.0};

double emrat = 81.300585;
