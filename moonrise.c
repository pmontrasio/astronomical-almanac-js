/* Search program to print tables of moonrise.  */

#define STARTDATE  2448431.5
#define ENDDATE    2448831.5

#include <stdlib.h>
#include "kep.h"
/* Conversion factors between degrees and radians */
double DTR = 1.7453292519943295769e-2;
double RTD = 5.7295779513082320877e1;
double RTS = 2.0626480624709635516e5;	/* arc seconds per radian */
double STR = 4.8481368110953599359e-6;	/* radians per arc second */
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

double sqrt (), asin (), log ();

/* coordinates of object
 */
int objnum = 0;			/* I.D. number of object */

/* ecliptic polar coordinates:
 * longitude, latitude in radians
 * radius in au
 */
double FAR obpolar[3];

struct orbit objorb =
{
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
  0,				/* &ear404, */
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
struct orbit earth =
{
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

double zgetdate (), gethms ();
double search ();
static void func ();
int apparent ();

/* Results computed by domoon.c  */
/* Transit, rise, and set times in radians (2 pi = 1 day) */
extern int f_trnsit;
extern double r_trnsit;
extern double r_rise;
extern double r_set;

double t_trnsit;
double t_rise;
double t_set;

int 
main ()
{
  double u;

  kinit ();
  printf("\nTable of lunar rise, transit, and set times.\n\n");
  objnum = 0;
  u = STARTDATE;
  while (u <= ENDDATE)
    {
      prtflg = 0;
      u = search (u);
      prtflg = 1;
      if (f_trnsit)
        jtocal (t_rise);
      else
	printf("\n");
      jtocal (t_trnsit);
      if (f_trnsit)
        jtocal (t_set);
      else
	printf("\n");
      prtflg = 0;
      printf ("\n");
      u += 1.0;
    }
  exit (0);
}

/* Search for improved rise and set time estimates.  */
double
search (t)
     double t;
{
  double JDsave, TDTsave, UTsave;
  double date_save, date, t0, t1;
  double rise1, set1, trnsit1;

  JDsave = JD;
  TDTsave = TDT;
  UTsave = UT;
  date_save = floor (t - 0.5) + 0.5;

  /* Find transit time. */
  date = date_save;
  func (t);
  do
    {
      if (r_trnsit < 0.0)
	{
	  date -= 1.0;
	  r_trnsit += TPI;
	}
      if (r_trnsit > TPI)
	{
	  date += 1.0;
	  r_trnsit -= TPI;
	}
      t0 = date + r_trnsit / TPI;
      func (t0);
      t1 = date + r_trnsit / TPI;
    }
  while (fabs (t1 - t0) > .0001);

  t_trnsit = t1;
  trnsit1 = r_trnsit;
  set1 = r_set;
  if (f_trnsit == 0)
	goto sdone;

  /* Set current date to be that of the transit just found.  */
  date_save = date;
  do
    {
      if (r_rise < 0.0)
	{
	  date -= 1.0;
	  r_rise += TPI;
	}
      if (r_rise > TPI)
	{
	  date += 1.0;
	  r_rise -= TPI;
	}
      t0 = date + r_rise / TPI;
      func (t0);
      t1 = date + r_rise / TPI;
    }
  while (fabs (t1 - t0) > .0001);
  rise1 = r_rise;
  t_rise = t1;

  date = date_save;
  r_set = set1;
  do
    {
      if (r_set < 0.0)
	{
	  date -= 1.0;
	  r_set += TPI;
	}
      if (r_set > TPI)
	{
	  date += 1.0;
	  r_set -= TPI;
	}
      t0 = date + r_set / TPI;
      func (t0);
      t1 = date + r_set / TPI;
    }
  while (fabs (t1 - t0) > .0001);

  t_set = t1;
  r_trnsit = trnsit1;
  r_rise = rise1;
  /*  printf("%.15e %.15e %.15e\n", rise1, trnsit1, set1); */
sdone:
  JD = JDsave;
  TDT = TDTsave;
  UT = UTsave;
  return t_trnsit;
}

/* Compute estimate of lunar rise and set times.  */
static void
func (t)
     double t;
{
  int prtsave;

  prtsave = prtflg;
  prtflg = 0;
  objnum = 3;
  JD = t;
  update ();			/* find UT */
  kepler (TDT, &earth, rearth, eapolar);
  domoon ();
  prtflg = prtsave;
}
