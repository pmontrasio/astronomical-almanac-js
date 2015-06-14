/* Search program to find dates of equinox or new/full moon.  */

double STARTDATE;
double ENDDATE;

#define NEWMOON 0
#define FULLMOON 1
#define SPRING 2
#define SUMMER 3
#define AUTUMN 4
#define WINTER 5

#define _XOPEN_SOURCE /* glibc2 needs this to declare strptime */
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include "kep.h"
/* Conversion factors between degrees and radians */
double DTR = 1.7453292519943295769e-2;
double RTD = 5.7295779513082320877e1;
double RTS = 2.0626480624709635516e5; /* arc seconds per radian */
double STR = 4.8481368110953599359e-6; /* radians per arc second */
double PI = 3.14159265358979323846;
double TPI = 6.28318530717958647693;

/* unused, but must be defined for rplanet.c and rstar.c */
struct orbit *elobject;	/* pointer to orbital elements of object */

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
static int event=NEWMOON;

double zgetdate(), gethms();
double func(), search();
int apparent();

int main(int argc, char **argv)
{
  double t=0, t0;
  int c;
  time_t now=time(NULL);
  struct tm date=*localtime(&now);

  kinit();

  STARTDATE=caltoj( date.tm_year+1900, date.tm_mon+1, date.tm_mday );
  ENDDATE=caltoj( date.tm_year+1900+1, date.tm_mon+1, date.tm_mday );

  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"start", 1, 0, 's'},
      {"end", 1, 0, 'e'},
      {"spring", 0, 0, 'V'},
      {"vernal", 0, 0, 'V'},
      {"summer", 0, 0, 'S'},
      {"autumn", 0, 0, 'A'},
      {"fall", 0, 0, 'A'},
      {"winter", 0, 0, 'W'},
      {"newmoon", 0, 0, 'N'},
      {"fullmoon", 0, 0, 'F'},
      {"help", 0, 0, 'h'},
      {0, 0, 0, 0}
    };
    c = getopt_long (argc, argv, "s:e:VSAWNFh",
		     long_options, &option_index);
    if (c == -1)
      break;
    switch (c) {
    case 's': 
      if (strptime(optarg, "%Y-%m-%d", &date))
	STARTDATE=caltoj( date.tm_year+1900, date.tm_mon+1, date.tm_mday );
      else
	STARTDATE=atof(optarg); 
      break;
    case 'e':
      if (strptime(optarg, "%Y-%m-%d", &date))
	ENDDATE=caltoj( date.tm_year+1900, date.tm_mon+1, date.tm_mday );
      else
	ENDDATE=atof(optarg); 
      break;
    case 'V': event=SPRING; break;
    case 'S': event=SUMMER; break;
    case 'A': event=AUTUMN; break;
    case 'W': event=WINTER; break;
    case 'N': event=NEWMOON; break;
    case 'F': event=FULLMOON; break;
    case 'h': 
    default:
      printf("conjunct -  find dates of equinox or new/full moon\n"
	     "options:\n"
	     " -s, --start DATE  Julian or ISO 8601 starting date (default today)\n"
	     " -e, --end DATE  ending date (default 1 year from today)\n"
	     " -V, --vernal    find times of spring equinox\n"
	     " -S, --summer    find times of summer solstice\n"
	     " -A, --autumn    find times of autumn equinox\n"
	     " -W, --winter    find times of winter solstice\n"
	     " -F, --fullmoon  find times of full moon\n"
	     " -N, --newmoon   find times of new moon (default)\n");
      exit(0);		     
      break;
    }
  }
  if (STARTDATE >= ENDDATE)
    {
      printf("start date=%f >= end date=%f\n", STARTDATE, ENDDATE);
      exit(1);
    }
  
  objnum = 0;
  //  printf("184:STARTDATE=%6.1f\n", STARTDATE);
  t0 = STARTDATE - (((event==FULLMOON)||(event==NEWMOON))?40:400);
  while( 1 )
    {
      prtflg = 0;
      switch(event)
	{
	case SPRING:
	  t = search( t0, 0.0, 364.0 );
	  break;
	case AUTUMN:
	  t = search( t0, PI, 364.0 );
	  break;
	case WINTER:
	  t = search( t0, -PI/2, 364.0 );
	  break;
	case SUMMER:
	  t = search( t0, PI/2, 364.0 );
	  break;
	case NEWMOON:
	  t = search( t0, 0.0, 27.0 );
	  break;
	case FULLMOON:
	  t = search( t0, PI, 27.0 );
	  break;
	}
      TDT = t;
      if (t > ENDDATE)
	break;
      if (t >= STARTDATE)
	{
	  printf("%.4f ", t);
	  prtflg = 1;
	  jtocal(t);
	}
      t0 = t;
    }
exit(0);
}

/* Search for desired conjunction.
   On the first call to this function, search forward from date T
   in steps of DELTA / 8 days until the error changes sign.
   On subsequent calls, step forward DELTA days from T, then search forward
   in steps of 1 day until the error changes sign.
   Then, in both cases, reduce the error by interval halving
   until the function equals Y with the desired precision.  */

static int first_search = 0;

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
	  /* Keep searching while the error is large.  */
	  if (fabs (eh) > 0.5*PI)
	    continue;
	  /* Look for sign change */
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
	  /* Replace the inteval boundary that has the same sign as em.  */
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
}

/* Compute desired relation of apparent ecliptic longitude
   as a function of the ephemeris date.  */
double func(t)
double t;
{
  double p[3], q[3], polar[3];
  double s;
  double m;
  double val;

  TDT = t;
  /* Longitude of the sun.  */
  objnum = 0;
  apparent(p,q);
  lonlat(p,TDT,polar,0);
  val = s = polar[0];
  switch(event)
    {
    case NEWMOON:
    case FULLMOON:
      /* Longitude of the moon.  */
      objnum = 3;
      apparent(p,q);
      lonlat(p,TDT,polar,0);
      m = polar[0];
      val = s - m;
      break;
    case SPRING:
    case WINTER:
    case AUTUMN:
    case SUMMER:
      val = s;
      break;
    }
  return val;
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

/* Adjust for nutation
 * at current ecliptic.
 */
epsiln( TDT );
nutate( TDT, p );

/* Return dimensions in au.  */
for(i=0; i<3; i++)
	p[i] *= EO;
return 0;
}
