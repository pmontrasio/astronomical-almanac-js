/* Search program to find dates of equinox or new moon.  */

/* November 28, 4057 B.C.*/
/* #define STARTDATE  239935.5 */

#define STARTDATE  625379.5

/* January 31, 1982 */
#define ENDDATE    2445000.5

/* Define one of these nonzero, the others zero.  */
#define NEWMOON 1
#define FULLMOON 0
#define FIRST_QUARTER_MOON 0
#define THIRD_QUARTER_MOON 0
#define EQUINOX 0

#include "kep.h"
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

/* display formats for printf() */
extern char *intfmt, *dblfmt;

/* display enable flag */
int prtflg = 0;

double dradt, ddecdt;
extern double FAR moonpol[];
extern double FAR moonpp[];

double zgetdate(), gethms();
double func(), search();
int apparent();

struct orbit *elobject;
double robject[3] = {0.0, 0.0, 0.0};
static int first_search;


int main()
{
  double t, t0;

  kinit();
  objnum = 0;
  t0 = STARTDATE;
  first_search = 0;
  while( t0 <= ENDDATE )
    {
      prtflg = 0;
#if EQUINOX
      t = search( t0, 0.0, 364.0 );
#endif
#if NEWMOON
      t = search( t0, 0.0, 27.0 );
#endif
#if FULLMOON
      t = search( t0, PI, 27.0 );
#endif
#if FIRST_QUARTER_MOON
      t = search( t0, -0.5*PI, 27.0 );
#endif
#if THIRD_QUARTER_MOON
      t = search( t0, 0.5*PI, 27.0 );
#endif
      TDT = t;
      printf("%.4f ", t);
      prtflg = 1;
      jtocal(t);
      t0 = t;
    }
exit(0);
}


/* Search for desired conjunction.
   On the first call to this function, search forward from date T
   in steps of delta / 8 days until the error changes sign.
   On subsequent calls, step forward delta days, then search forward
   in steps of 1 day until the error changes sign.
   Then, in both cases, reduce the error by interval halving
   until the function equals Y with the desired precision.  */

double search(t, y, delta)
double t, y, delta;
{
  double tl, tm, th, yl, ym, yh, el, eh, em, dt;

  if (first_search == 0)
    {
      dt = 0.125 * delta;
      first_search = 1;
      th = t;
    }
  else
    {
      dt = 1.0;
      th = t + delta;
    }
  yh = func(th);
  eh = yh - y;
  if (eh > PI)
    eh -= TPI;
  if (eh <= -PI)
    eh += TPI;
  /* Bracket the solution date.  */
  for (;;)
    {
      tl = th;
      yl = yh;
      el = eh;
      th = th + dt;
      yh = func(th);
      eh = yh - y;
    if(eh > PI)
      eh -= TPI;
    if(eh <= -PI)
      eh += TPI;
    /* Keep searching while error is large.  */
    if (fabs (eh) > 0.5*PI)
	continue;
    /* Look for sign change.  */
    if((el * eh) <= 0.0)
	break;
    }

/* Search by simple interval halving.  */
while( (th - tl) > .00001 )
  {
    tm = 0.5 * (tl + th);
    ym = func(tm);
    em = ym - y;
    if(em > PI)
      em -= TPI;
    if(em <= -PI)
      em += TPI;
    /* Replace the interval boundary that has the same sign as em.  */
     if( em * eh > 0 )
      {
	yh = ym;
	th = tm;
	eh = em;
      }
    else
      {
	yl = ym;
	tl = tm;
	el = em;
      }
  }
 tm = tl + (th - tl) * (-el)/(yh -yl);
return (tm);
}

/* Compute desired relation of apperent ecliptic longitude
   as a function of the ephemeris date.  */
double func(t)
double t;
{
  double p[3], q[3], polar[3];
  double s;
#if NEWMOON || FULLMOON || FIRST_QUARTER_MOON || THIRD_QUARTER_MOON
  double m;
#endif

  TDT = t;
  /* Longitude of the sun.  */
  objnum = 0;
  apparent(p,q);
  lonlat(p,TDT,polar,0);
  s = polar[0];
#if NEWMOON || FULLMOON || FIRST_QUARTER_MOON || THIRD_QUARTER_MOON
  /* Longitude of the moon.  */
  objnum = 3;
  apparent(p,q);
  lonlat(p,TDT,polar,0);
  m = polar[0];
  return (s - m);
#endif
#if EQUINOX
  return s;
#endif
}


/* Find apparent coordinates of object at Julian date TDT.
   Q is heliocentric position of object.
   P is geocentric, object minus earth.
   Both outputs are in equatorial rectangular coordinates
   and are referred to the equinox and ecliptic of date.  */

int apparent( p, q )
double p[], q[];
{
double polar[3];
int i;
static double TDTearth = -1.0e38;


/* Calculate heliocentric position of the earth */
if(TDTearth != TDT)
	{
	kepler( TDT, &earth, rearth, eapolar );
	TDTearth = TDT;
	}

if (objnum == 3)
	{
	moonll(TDT, p, polar);
	for(i=0; i<3; i++)
		{
		q[i] = p[i] + rearth[i];
		}
	return 0;
	}
else if( objnum == 0 )
	{
	for(i=0; i<3; i++)
		{
		q[i] = 0.0;
		p[i] = -rearth[i];
		}
	}
else if ((objnum > 0) && (objnum < 10))
	{
	kepler( TDT, &objorb, q, polar );
	}
else
	exit(1);

/* Adjust for light time (planetary aberration)
*/
if( objnum != 0 )
	lightt( &objorb, q, rearth );

/* Find Euclidean vectors between earth, object, and the sun
 */
for( i=0; i<3; i++ )
	p[i] = q[i] - rearth[i];

angles( p, q, rearth );

/* Find unit vector from earth in direction of object
 */
for( i=0; i<3; i++ )
	p[i] /= EO;

/* Correct position for light deflection
 */
if(objnum != 0)
	relativity( p, q, rearth );

/* Correct for annual aberration
 */
annuab( p );

/* Precession of the equinox and ecliptic
 * from J2000.0 to ephemeris date
 */
precess( p, TDT, -1 );

/* Ajust for nutation
 * at current ecliptic.
 */
epsiln( TDT );
nutate( TDT, p );

/* Return dimensions in au.  */
for(i=0; i<3; i++)
	p[i] *= EO;
return 0;
}
