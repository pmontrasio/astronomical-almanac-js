/* This program is based on aa.c and calculates the distance of a planet from Earth
 */

#include "consts.h"

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

#ifdef __EMSCRIPTEN__
#include <stdlib.h>
#endif

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

extern double deltat_value;

double distance;

/* correction vector, saved for display  */
double dp[3];

/* display formats for printf()
 */
extern char *intfmt, *dblfmt;

/* display enable flag
 */

/* Tabulation parameters
 */
static double djd = 1.0;
static int ntab = 1;

struct orbit *elobject;	/* pointer to orbital elements of object */

double km_to_planet(julian_day, planet)
  double julian_day;
  int planet;
 {
  int i;
  JD = julian_day;
  objnum = planet;

  kinit();

  update();
  djd = 1.0;
  ntab = 1;

  switch(objnum)
    {
    case 1: elobject = &mercury; break;
    case 2: elobject = &venus; break;
    case 4: elobject = &mars; break;
    case 5: elobject = &jupiter; break;
    case 6: elobject = &saturn; break;
    case 7: elobject = &uranus; break;
    case 8: elobject = &neptune; break;
    case 9: elobject = &pluto; break;
    default:
operr: printf( "Operator error.\n" );
      exit(0);
    }

  /*  printf( "\n                  %s\n", &elobject->obname[0] );*/

  /* Always calculate heliocentric position of the earth */
  kepler( TDT, &earth, rearth, eapolar );
  distance = planet_distance(elobject) * 149597871.0; /* km */
  return distance;
}
