/* Calculate time of transit
 * assuming RA and Dec change uniformly with time
 *
 * -- S. L. Moshier
 */

#include "kep.h"
extern double glat;

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
int f_trnsit;

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

/*  x = (J - floor(J-0.5) - 0.5) * TPI; */
 x = floor(J) - 0.5;
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
	case 0: N = COSSUN; break;

/* Moon, elevation = -34' - semidiameter + parallax
 * semidiameter = 0.272453 * parallax + 0.0799"
 */
	case 3:
		N = 1.0/(DISFAC*obpolar[2]);
		D = asin( N ); /* the parallax */
		N =  - 9.890199094634534e-3
			+ (1.0 - 0.2725076)*D
			- 3.874e-7;
		N = sin(N);
		break;

/* Other object */
	default: N = COSZEN; break;
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
  double date, date_save, date_trnsit, t0, t1;
  double rise1, set1, trnsit1, loopctr, retry;
  double TPI;
  int prtsave;

  loopctr = 0;
  prtsave = prtflg;
  TPI = PI + PI;
  JDsave = JD;
  TDTsave = TDT;
  UTsave = UT;
  date = floor(UT - 0.5) + 0.5;
  retry = 0;
  date_save = date;
  t1 = date_save;

  /* Find transit time. */
  do
    {
      date = date_save;
      t0 = t1;
      iter_func(t0, func);
      t1 = date + r_trnsit / TPI;
      if (++loopctr > 10)
	{
	  printf ("? Transit time did not converge.\n");
	  goto no_trnsit;
	}

  /* Reject transit if it is more than half a day from entered date.  */
  if (retry == 0)
    {
      if ((UTsave - t1) > 0.5)
	{
	  date_save += 1.0;
	  t1 = t_trnsit + 1.0;
	  retry = 1;
	  /*	  goto do_retry; */
	}
      if ((UTsave - t1) < -0.5)
	{
	  date_save -= 1.0;
	  t1 = t_trnsit - 1.0;
	  retry = 1;
	  /*	  goto do_retry; */
	}
    }
    }
  while (fabs(t1 - t0) > .0001);

  t_trnsit = t1;
  trnsit1 = r_trnsit;
  set1 = r_set;


  if (f_trnsit == 0)
	goto prtrnsit;

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
      /* Skip out if no event found.  Don't report rise or set.  */
      if ((f_trnsit == 0) || (++loopctr > 10))
	{
	  if (loopctr > 10)
	    {
	      printf ("? Rise time did not converge.\n");
	    }
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
      if ((f_trnsit == 0) || (++loopctr > 10))
	{
	  if (loopctr > 10)
	    {
	      printf ("? Set time did not converge.\n");
	    }
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
  if (f_trnsit != 0)
    {
      printf( "rises " );
      UT = t_rise;
      jtocal (t_rise);
      printf( "sets " );
      UT = t_set;
      jtocal (t_set);
      t0 = t_set - t_rise;
      if ((t0 > 0.0) && (t0 < 1.0))
	printf("Visible hours %.4f\n", 24.0 * t0);
      if ((fabs(JDsave - t_rise) > 0.5) && (fabs(JDsave - t_trnsit) > 0.5)
	  && (fabs(JDsave - t_set) > 0.5))
	printf("? wrong event date.\n");
    }
no_trnsit:
  JD = JDsave;
  TDT = TDTsave;
  UT = UTsave;
  /* Reset to original input date entry.  */
  prtflg = 0;
  update();
  prtflg = prtsave;
  return 0;
}
