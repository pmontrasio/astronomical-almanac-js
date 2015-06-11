/* Calculate time of transit
 * assuming RA and Dec change uniformly with time
 *
 * -- S. L. Moshier
 */

#include "kep.h"
extern double glat;
extern double azimuth, elevation, refracted_elevation; /* degrees */

/* Earth radii per au */
#define DISFAC 2.3454780e4

/* cosine of 90 degrees 50 minutes: */
#define COSSUN -0.014543897651582657
/* cosine of 90 degrees 34 minutes: */
#define COSZEN -9.8900378587411476e-3

/* Returned transit, rise, and set times in radians (2 pi = 1 day) */
double r_trnsit;
double r_rise;
double r_set;
double elevation_threshold, semidiameter;
int f_trnsit, southern_hemisphere;
int no_rise_set (double, int (*)());

double search_halve (double tl, double yl, double th, double yh, int (* func)());
int iter_trnsit( int (* func)() );
static void iter_func( double t, int (* func)() );



int trnsit(J, lha, dec)
double J; /* Julian date */
double lha; /* local hour angle, radians */
double dec; /* declination, radians */
{
double x, y, z, N, D, TPI;
double lhay, cosdec, sindec, coslat, sinlat;

f_trnsit = 0;
TPI = 2.0*PI;
/* Initialize to no-event flag value. */
r_rise = -10.0;
r_set = -10.0;
/* observer's geodetic latitude, in radians */
x = glat * DTR;
coslat = cos(x);
sinlat = sin(x);

cosdec = cos(dec);
sindec = sin(dec);

if (sinlat < 0)
  southern_hemisphere = 1;
 else
   southern_hemisphere = 0;

/* Refer to same start of date as iter_trnsit,
   so r_trnsit means the same thing in both programs.  */
 x = floor(J - 0.5) + 0.5;
 x = (J - x) * TPI;
/* adjust local hour angle */
y = lha;
/* printf ("%.7f,", lha); */
while( y < -PI )
	y += TPI;
while( y > PI )
	y -= TPI;
lhay = y;
y =  y/( -dradt/TPI + 1.00273790934);
r_trnsit = x - y;
/* printf ("rt %.7f ", r_trnsit); */
/* Ordinarily never print here.  */
if( prtflg > 1 )
	{
	printf( "approx local meridian transit" );
	hms( r_trnsit );
	printf( "UT  (date + %.5f)\n", r_trnsit/TPI );
	}
if( (coslat == 0.0) || (cosdec == 0.0) )
	goto norise; 

/* The time at which the upper limb of the body meets the
 * horizon depends on the body's angular diameter.
 */
switch( objnum )
	{
/* Sun */
	case 0:
	  N = COSSUN;
	  semidiameter = 0.2666666666666667;
	  elevation_threshold = -0.8333333333333333;
	  break;

/* Moon, elevation = -34' - semidiameter + parallax
 * semidiameter = 0.272453 * parallax + 0.0799"
 */
	case 3:
		N = 1.0/(DISFAC*obpolar[2]);
		D = asin( N ); /* the parallax */
		semidiameter = 0.2725076*D + 3.874e-7;
		N =  -9.890199094634534e-3 - semidiameter + D;
		semidiameter *= RTD;
		elevation_threshold = -34.0/60.0 - semidiameter;
		N = sin(N);
		break;

/* Other object */
	default:
	  N = COSZEN;
	  semidiameter = 0.0;
	  elevation_threshold = -0.5666666666666666;
	  break;
	}
y = (N - sinlat*sindec)/(coslat*cosdec);

if( (y < 1.0) && (y > -1.0) )
	{
	f_trnsit = 1;
/* Derivative of y with respect to declination
 * times rate of change of declination:
 */
	z = -ddecdt*(sinlat + COSZEN*sindec);
	z /= TPI*coslat*cosdec*cosdec;
/* Derivative of acos(y): */
	z /= sqrt( 1.0 - y*y);
	y = acos(y);
	D = -dradt/TPI + 1.00273790934;
	r_rise = x - (lhay + y)*(1.0 + z)/D;
	r_set = x - (lhay - y)*(1.0 - z)/D;
	/* Ordinarily never print here.  */
		if( prtflg > 1 )
		{
		printf( "rises approx " );
		hms(r_rise);
		printf( "UT,  sets approx " );
		hms(r_set);
		printf( "UT\n" );
		}
	}
norise:		;
return(0);
}





extern struct orbit earth;

/* Julian dates of rise, transet and set times.  */
double t_rise;
double t_trnsit;
double elevation_trnsit;
double t_set;
/* Compute estimate of lunar rise and set times for iterative solution.  */

static void
iter_func(t, func)
double t;
int (* func)();
{
  int prtsave;

  prtsave = prtflg;
  prtflg = 0;
  JD = t;
  update(); /* Set UT and TDT */
  kepler( TDT, &earth, rearth, eapolar );
  (*func)();
  prtflg = prtsave;
}


/* Iterative computation of rise, transit, and set times.  */

int
iter_trnsit( func )
int (* func)();
{
  double JDsave, TDTsave, UTsave;
  double date, date_trnsit, t0, t1;
  double rise1, set1, trnsit1, loopctr, retry;
  double TPI;
  int prtsave;

  loopctr = 0;
  prtsave = prtflg;
  TPI = PI + PI;
  JDsave = JD;
  TDTsave = TDT;
  UTsave = UT;
  retry = 0;
  /* Start iteration at time given by the user.  */
  t1 = UT;

  /* Find transit time. */
  do
    {
      t0 = t1;
      date = floor(t0 - 0.5) + 0.5;
      iter_func(t0, func);
      t1 = date + r_trnsit / TPI;
      if (++loopctr > 10)
	{
	  printf ("? Transit time did not converge.\n");
	  goto no_trnsit;
	}
    }
  while (fabs(t1 - t0) > .0001);

  t_trnsit = t1;
  elevation_trnsit = elevation;
  trnsit1 = r_trnsit;
  set1 = r_set;


  if (f_trnsit == 0)
    {
      /* Rise or set time not found.  Apply a search technique to
	 check near inferior transit if object is above horizon now.  */
      t_rise = -1.0;
      t_set = -1.0;
      if (elevation > elevation_threshold)
	no_rise_set (t_trnsit, func);
      goto prtrnsit;
    }

  /* Set current date to be that of the transit just found.  */
  date_trnsit = date;
  t1 = date + r_rise / TPI;
  /* Choose rising no later than transit.  */
  if (t1 >= t_trnsit)
    {
      date -= 1.0;
      t1 = date + r_rise / TPI;
    }
  loopctr = 0;
  do
    {
      t0 = t1;
      iter_func(t0, func);
      /* Skip out if no event found.  */

      if (f_trnsit == 0)
	{
	  /* Rise or set time not found.  Apply search technique.  */
	  t_rise = -1.0;
	  t_set = -1.0;
	  no_rise_set (t_trnsit, func);
	  goto prtrnsit;
	}
      if (++loopctr > 10)
	{
	  printf ("? Rise time did not converge.\n");
	  f_trnsit = 0;
	  goto prtrnsit;
	}
      t1 = date + r_rise / TPI;
      if (t1 > t_trnsit)
	{
	  date -= 1;
	  t1 = date + r_rise / TPI;
	}
    }
  while (fabs(t1 - t0) > .0001);
  rise1 = r_rise;
  t_rise = t1;

  /* Set current date to be that of the transit.  */
  date = date_trnsit;
  r_set = set1;
  /* Choose setting no earlier than transit.  */
  t1 = date + r_set / TPI;
  if (t1 <= t_trnsit)
    {
      date += 1.0;
      t1 = date + r_set / TPI;
    }
  loopctr = 0;
  do
    {
      t0 = t1;
      iter_func(t0, func);
      if (f_trnsit == 0)
	{
	  /* Rise or set time not found.  Apply search technique.  */
	  t_rise = -1.0;
	  t_set = -1.0;
	  no_rise_set (t_trnsit, func);
	  goto prtrnsit;
	}
      if (++loopctr > 10)
	{
	  printf ("? Set time did not converge.\n");
	  f_trnsit = 0;
	  goto prtrnsit;
	}
      t1 = date + r_set / TPI;
      if (t1 < t_trnsit)
	{
	  date += 1.0;
	  t1 = date + r_set / TPI;
	}
    }
  while (fabs(t1 - t0) > .0001);

  t_set = t1;
  r_trnsit = trnsit1;
  r_rise = rise1;

prtrnsit:
  prtflg = 1;
  printf( "local meridian transit " );
  UT = t_trnsit;
  jtocal (t_trnsit);
  if (t_rise != -1.0)
    {
      printf( "rises " );
      UT = t_rise;
      jtocal (t_rise);
    }
  if (t_set != -1.0)
    {
      printf( "sets " );
      UT = t_set;
      jtocal (t_set);
      if (t_rise != -1.0)
	{
	  t0 = t_set - t_rise;
	  if ((t0 > 0.0) && (t0 < 1.0))
	    printf("Visible hours %.4f\n", 24.0 * t0);
	}
    }

  if ((fabs(JDsave - t_rise) > 0.5) && (fabs(JDsave - t_trnsit) > 0.5)
      && (fabs(JDsave - t_set) > 0.5))
    printf("? wrong event date.\n");

no_trnsit:
  JD = JDsave;
  TDT = TDTsave;
  UT = UTsave;
  /* Reset to original input date entry.  */
  prtflg = 0;
  update();
  prtflg = prtsave;
  f_trnsit = 1;
  return 0;
}


/* If the initial approximation fails to locate a rise or set time,
   this function steps between the transit time and the previous
   or next inferior transits to find an event more reliably.  */

#define STEP_SCALE 0.5

int no_rise_set (t0, func)
     double t0;
     int (* func)();
{
  double t_trnsit0 = t_trnsit;
  double el_trnsit0 = elevation_trnsit;
  double t, e;
  double t_above, el_above, t_below, el_below;

  /* Step time toward previous inferior transit to find
     whether a rise event was missed.  The step size is a function
     of the azimuth and decreases near the transit time.  */
  t_above = t_trnsit0;
  el_above = el_trnsit0;
  t_below = -1.0;
  el_below = el_above;
  t = t_trnsit0 - 0.25;
  e = 1.0;
  while (e > 0.005)
    {
      iter_func(t, func);
      if (elevation > elevation_threshold)
	{
	  /* Object still above horizon.  */
	  t_above = t;
	  el_above = elevation;
	}
      else
	{
	  /* Object is below horizon.  Rise event is bracketed.
	     Proceed to interval halving search.  */
	  t_below = t;
	  el_below = elevation;
	  goto search_rise;
	}
      /* Step time by an amount proportional to the azimuth deviation.  */
      e = azimuth/360.0;
      if (azimuth < 180.0)
	{
	  if (southern_hemisphere == 0)
	    t -= STEP_SCALE * e;
	  else
	    t += STEP_SCALE * e;
	}
      else
	{
	  e = 1.0 - e;
	  if (southern_hemisphere == 0)
	    t += STEP_SCALE * e;
	  else
	    t -= STEP_SCALE * e;
	}
    };

  /* No rise event detected.  */

  if (elevation > elevation_threshold)
    {
      /* printf ("Previous inferior transit is above horizon.\n"); */
      t_rise = -1.0;
      goto next_midnight;
    }
  /* Find missed rise time. */
search_rise:
  t_rise = search_halve (t_below, el_below,
			 t_above, el_above, func);
  f_trnsit = 1;

next_midnight:

  /* Step forward in time toward the next inferior transit.  */
  t_above = t_trnsit0;
  el_above = el_trnsit0;
  t_below = -1.0;
  el_below = el_above;
  t = t_trnsit0 + 0.25;
  e = 1.0;
  while (e > 0.005)
    {
      iter_func(t, func);
      if (elevation > elevation_threshold)
	{
	  /* Object still above horizon.  */
	  t_above = t;
	  el_above = elevation;
	}
      else
	{
	  /* Object is below horizon.  Event is bracketed.
	     Proceed to interval halving search.  */
	  t_below = t;
	  el_below = elevation;
	  goto search_set;
	}
      /* Step time by an amount proportional to the azimuth deviation.  */
      e = azimuth/360.0;
      if (azimuth < 180.0)
	{
	  if (southern_hemisphere == 0)
	    t -= STEP_SCALE * e;
	  else
	    t += STEP_SCALE * e;  /* Southern hemisphere observer.  */
	}
      else
	{
	  e = 1.0 - e;
	  if (southern_hemisphere == 0)
	    t += STEP_SCALE * e;
	  else
	    t -= STEP_SCALE * e;
	}
    };

  if (elevation > elevation_threshold)
    {
      /* printf ("Next inferior transit is above horizon.\n"); */
      t_set = -1.0;
      return 0;
    }
  /* Find missed set time. */
search_set:
  t_set = search_halve (t, elevation,
			t_trnsit, elevation_trnsit, func);
  f_trnsit = 1;
  return 0;
}


/* Search rise or set time by simple interval halving
   after the event has been bracketed in time.  */

double
search_halve (t1, y1, t2, y2, func)
     double t1;
     double y1;
     double t2;
     double y2;
     int (* func)();
{
  double e2, e1, em, tm, ym;

  e2 = y2 - elevation_threshold;
  e1 = y1 - elevation_threshold;
  tm = 0.5 * (t1 + t2);

while( fabs(t2 - t1) > .00001 )
  {
    /* Evaluate at middle of current interval.  */
    tm = 0.5 * (t1 + t2);
    iter_func(tm, func);
    ym = elevation;
    em = ym - elevation_threshold;
    /* Replace the interval boundary whose error has the same sign as em.  */
     if( em * e2 > 0 )
      {
	y2 = ym;
	t2 = tm;
	e2 = em;
      }
    else
      {
	y1 = ym;
	t1 = tm;
	e1 = em;
      }
  }
return (tm);
}
