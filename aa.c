/* This program calculates orbits of planetary bodies and reduces
 * the coordinates of planets or stars to geocentric and topocentric
 * place.  An effort has been made to use rigorous methods throughout.
 *
 * References to AA page numbers are to The Astronomical Almanac, 1986
 * published by the U.S. Government Printing Office.
 *
 * The program uses as a default the PLAN404 approximations to DE404
 * for planetary positions.
 *
 * Warning! Your atan2() function may not work the same as the one
 * assumed by this program.
 * atan2(x,y) computes atan(y/x), result between 0 and 2pi.
 *
 * S. L. Moshier, November, 1987
 */



/* Conversion factors between degrees and radians */
double DTR = 1.7453292519943295769e-2;
double RTD = 5.7295779513082320877e1;
double RTS = 2.0626480624709635516e5; /* arc seconds per radian */
double STR = 4.8481368110953599359e-6; /* radians per arc second */
double PI = 3.14159265358979323846;
extern double PI;

/* Standard epochs.  Note Julian epochs (J) are measured in
 * years of 365.25 days.
 */
double J2000 = 2451545.0;	/* 2000 January 1.5 */
double B1950 = 2433282.423;	/* 1950 January 0.923 Besselian epoch */
double J1900 = 2415020.0;	/* 1900 January 0, 12h UT */

/* Data structures containing orbital elements of
 * objects that orbit the sun.  See kep.h for the definition.
 */
#include "kep.h"

#ifdef __BORLANDC__
#include <stdlib.h>
#endif

#ifdef _MSC_VER
#if _MSC_VER >= 1000
#include <stdlib.h>
#include <string.h>
#endif
#endif

/* approximate motion of right ascension and declination
 * of object, in radians per day
 */
double FAR dradt;
double FAR ddecdt;

/* Space for star description read from a disc file.
 */
struct star fstar;

/* Space for orbit read from a disc file. Entering 99 for the
 * planet number yields a prompt for a file name containg ASCII
 * strings specifying the elements.
 */
struct orbit forbit;

/* Orbits for each planet.  The indicated orbital elements are
 * not actually used, since the positions are are now calculated
 * from a formula.  Magnitude and semidiameter are still used.
 */
 /* Programs to compute perturbations. */
extern struct plantbl mer404, ven404, ear404, mar404;
extern struct plantbl jup404, sat404, ura404, nep404, plu404;

struct orbit mercury = {
"Mercury        ",
2446800.5, /* January 5.0, 1987 */
7.0048,
48.177,
29.074,
0.387098,
4.09236,
0.205628,
198.7199,
2446800.5,
-0.42,
3.36,
&mer404,
0.0,
0.0,
0.0
};

struct orbit venus = {
"Venus          ",
2446800.5,
3.3946,
76.561,
54.889,
0.723329,
1.60214,
0.006757,
9.0369,
2446800.5,
/* Note the calculated apparent visual magnitude for Venus
 * is not very accurate.
 */
-4.40,
8.34,
&ven404,
0.0,
0.0,
0.0
};

/* Fixed numerical values will be used for earth if read in from a file
 * named earth.orb.  See kfiles.c, kep.h.
 */
struct orbit earth = {
"Earth          ",
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
&ear404,
0.0,
0.0,
0.0
};
extern struct orbit earth;

struct orbit mars = {
"Mars           ",
2446800.5,
1.8498,
49.457,
286.343,
1.523710,
0.524023,
0.093472,
53.1893,
2446800.5,
-1.52,
4.68,
&mar404,
0.0,
0.0,
0.0
};

struct orbit jupiter = {
"Jupiter        ",
2446800.5,
1.3051,
100.358,
275.129,
5.20265,
0.0830948,
0.048100,
344.5086,
2446800.5,
-9.40,
98.44,
&jup404,
0.0,
0.0,
0.0
};

struct orbit saturn = {
"Saturn         ",
2446800.5,
2.4858,
113.555,
337.969,
9.54050,
0.0334510,
0.052786,
159.6327,
2446800.5,
-8.88,
82.73,
&sat404,
0.0,
0.0,
0.0
};

struct orbit uranus = {
"Uranus         ",
2446800.5,
0.7738,
73.994,
98.746,
19.2233,
0.0116943,
0.045682,
84.8516,
2446800.5,
-7.19,
35.02,
&ura404,
0.0,
0.0,
0.0
};

struct orbit neptune = {
"Neptune        ",
2446800.5,
1.7697,
131.677,
250.623,
30.1631,
0.00594978,
0.009019,
254.2568,
2446800.5,
-6.87,
33.50,
&nep404,
0.0,
0.0,
0.0
};

struct orbit pluto = {
"Pluto          ",
2446640.5,
17.1346,
110.204,
114.21,
39.4633,
0.00397570,
0.248662,
355.0554,
2446640.5,
-1.0,
2.07,
&plu404,
0.0,
0.0,
0.0
};

/*
int otest(), ctest();
*/
struct orbit test = {
"Test orbit     ",
2446800.5,
1.8498,
49.457,
286.343,
1.523710,
0.524023,
0.093472,
53.1893,
2446800.5,
-1.52,
4.68,
0,
0.0,
0.0,
0.0
};


/* coordinates of object
 */
int objnum = 0;	/* I.D. number of object */
double robject[3] = {0.0}; /* position */
/* ecliptic polar coordinates:
 * longitude, latitude in radians
 * radius in au
 */
double FAR obpolar[3];

/* coordinates of Earth
 */
/* Heliocentric rectangular equatorial position
 * of the earth at time TDT re equinox J2000
 */
double FAR rearth[3];
/* Corresponding polar coordinates of earth:
 * longitude and latitude in radians, radius in au
 */
double FAR eapolar[3];

/* Julian date of ephemeris
 */
double JD;
double TDT;
double UT;
extern double deltat_value;

/* flag = 0 if TDT assumed = UT,
 *      = 1 if input time is TDT,
 *      = 2 if input time is UT.
 */
int jdflag = 0;

/* correction vector, saved for display  */
double dp[3];

/* display formats for printf()
 */
extern char *intfmt, *dblfmt;

/* display enable flag
 */
int prtflg = 1;

/* Tabulation parameters
 */
static double djd = 1.0;
static int ntab = 1;

struct orbit *elobject;	/* pointer to orbital elements of object */

/* Main program starts here.
 */
int main()
{
int i;

double zgetdate(), gethms();

kinit();

loop:

prtflg = 1;
printf( "Enter starting date of tabulation\n" );
JD = zgetdate(); /* date */
JD += gethms();	/* time of day */
update(); /* find UT and ET */
printf( "Julian day %.7f, delta T = %.2f s\n", JD, deltat_value );

getnum( "Enter interval between tabulations in days", &djd, dblfmt );
getnum( "Number of tabulations to display", &ntab, intfmt );
if( ntab <= 0 )
	ntab = 1;

loop1:
getnum( "Planet number 0-9 or 88 to read star, 99 to read orbit",
	&objnum, intfmt );

switch(objnum)
	{
	case -1: exit(0);
	case 0: elobject = 0;
		printf( "\n                   The Sun\n" );
		break;
	case 1: elobject = &mercury; break;
	case 2: elobject = &venus; break;
	case 3: elobject = 0;
		printf( "\n                   The Moon\n" );
		break;
	case 4: elobject = &mars; break;
	case 5: elobject = &jupiter; break;
	case 6: elobject = &saturn; break;
	case 7: elobject = &uranus; break;
	case 8: elobject = &neptune; break;
	case 9: elobject = &pluto; break;
	case 10: elobject = &test; break;
	case 88:
	        elobject = (struct orbit *)&fstar;
		i = getstar( (struct star *) elobject );
		if( i == 1 )
			goto loop1;
		if( i == 0 )
			break;
		goto operr;
	case 99:
		elobject = &forbit;
		i = getorbit( elobject );
		if( i == 1 )
			goto loop1;
		if( i == 0 )
			break;
	default:
operr:		printf( "Operator error.\n" );
		goto loop;
	}

if( elobject == (struct orbit *)&fstar )
	showcname( &elobject->obname[0] );
else if( elobject )
	printf( "\n                  %s\n", &elobject->obname[0] );

for( i=0; i<ntab; i++ )
	{
/* print Julian date */
	printf( "\nJD %.2f,  ", JD );
	update();

/* Always calculate heliocentric position of the earth */
	kepler( TDT, &earth, rearth, eapolar );

	switch( objnum )
		{
		case 0: dosun(); iter_trnsit( dosun ); break;
		case 3: domoon(); iter_trnsit( domoon ); break;
		case 88: dostar(); iter_trnsit( dostar ); break;
		default: doplanet(); iter_trnsit( doplanet ); break;
		}
	printf( "\n" );
	JD += djd;
	}
goto loop;

#ifdef _MSC_VER
return 0;
#endif
}

