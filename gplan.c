/* gplan.c
   Routines to chew through tables of perturbations.  */

#include "plantbl.h"
#if __STDC__
extern int epsiln ( double J );
double cos (double);
double sin (double);
double floor (double);
double fabs (double);
void mean_elements (double);
static int sscc (int, double, int);
double g1plan (double, struct plantbl *);
int g2plan (double, struct plantbl *, double *);
int g3plan (double, struct plantbl *, double *, int);
#else
int epsiln();
double cos (), sin (), floor (), fabs ();
static int sscc ();
#endif
#define mods3600(x) ((x) - 1.296e6 * floor ((x)/1.296e6))

static double J2000 = 2451545.0;
#if 0
/* Conversion factors between degrees and radians */
static double DTR = 1.7453292519943295769e-2;
static double RTD = 5.7295779513082320877e1;
static double RTS = 2.0626480624709635516e5;	/* arc seconds per radian */
#endif
static double STR = 4.8481368110953599359e-6;	/* radians per arc second */

/* From Simon et al (1994)  */
static double freqs[] =
{
/* Arc sec per 10000 Julian years.  */
  53810162868.8982,
  21066413643.3548,
  12959774228.3429,
  6890507749.3988,
  1092566037.7991,
  439960985.5372,
  154248119.3933,
  78655032.0744,
  52272245.1795
};

static double phases[] =
{
/* Arc sec.  */
  252.25090552 * 3600.,
  181.97980085 * 3600.,
  100.46645683 * 3600.,
  355.43299958 * 3600.,
  34.35151874 * 3600.,
  50.07744430 * 3600.,
  314.05500511 * 3600.,
  304.34866548 * 3600.,
  860492.1546,
};

static double ss[NARGS][31];
static double cc[NARGS][31];


int 
gplan (J, plan, pobj)
     double J;
     struct plantbl *plan;
     double pobj[];
{
  register double su, cu, sv, cv, T;
  double t, sl, sb, sr;
  int i, j, k, m, n, k1, ip, np, nt;
#ifdef _MSC_VER
  char FAR *p;
  double FAR *pl;
  double FAR *pb;
  double FAR *pr;
#else
  /* On some systems such as Silicon Graphics, char is unsigned
     by default.  */
#ifdef vax
  /* VAX CC rejects "signed char."  */
  char *p;
#else
#ifdef __STDC__
  signed char *p;
#else
  char *p;
#endif
#endif
  double *pl, *pb, *pr;
#endif

  T = (J - J2000) / plan->timescale;
  n = plan->maxargs;
  /* Calculate sin( i*MM ), etc. for needed multiple angles.  */
  for (i = 0; i < n; i++)
    {
      if ((j = plan->max_harmonic[i]) > 0)
	{
	  sr = (mods3600 (freqs[i] * T) + phases[i]) * STR;
	  sscc (i, sr, j);
	}
    }

  /* Point to start of table of arguments. */
  p = plan->arg_tbl;
  /* Point to tabulated cosine and sine amplitudes.  */
#ifdef _MSC_VER
  pl = (double FAR *) plan->lon_tbl;
  pb = (double FAR *) plan->lat_tbl;
  pr = (double FAR *) plan->rad_tbl;
#else
  pl = (double *) plan->lon_tbl;
  pb = (double *) plan->lat_tbl;
  pr = (double *) plan->rad_tbl;
#endif
  sl = 0.0;
  sb = 0.0;
  sr = 0.0;

  for (;;)
    {
      /* argument of sine and cosine */
      /* Number of periodic arguments. */
      np = *p++;
      if (np < 0)
	break;
      if (np == 0)
	{			/* It is a polynomial term.  */
	  nt = *p++;
	  /* Longitude polynomial. */
	  cu = *pl++;
	  for (ip = 0; ip < nt; ip++)
	    {
	      cu = cu * T + *pl++;
	    }
	  sl +=  mods3600 (cu);
	  /* Latitude polynomial. */
	  cu = *pb++;
	  for (ip = 0; ip < nt; ip++)
	    {
	      cu = cu * T + *pb++;
	    }
	  sb += cu;
	  /* Radius polynomial. */
	  cu = *pr++;
	  for (ip = 0; ip < nt; ip++)
	    {
	      cu = cu * T + *pr++;
	    }
	  sr += cu;
	  continue;
	}
      k1 = 0;
      cv = 0.0;
      sv = 0.0;
      for (ip = 0; ip < np; ip++)
	{
	  /* What harmonic.  */
	  j = *p++;
	  /* Which planet.  */
	  m = *p++ - 1;
	  if (j)
	    {
	      k = j;
	      if (j < 0)
		k = -k;
	      k -= 1;
	      su = ss[m][k];	/* sin(k*angle) */
	      if (j < 0)
		su = -su;
	      cu = cc[m][k];
	      if (k1 == 0)
		{		/* set first angle */
		  sv = su;
		  cv = cu;
		  k1 = 1;
		}
	      else
		{		/* combine angles */
		  t = su * cv + cu * sv;
		  cv = cu * cv - su * sv;
		  sv = t;
		}
	    }
	}
      /* Highest power of T.  */
      nt = *p++;
      /* Longitude. */
      cu = *pl++;
      su = *pl++;
      for (ip = 0; ip < nt; ip++)
	{
	  cu = cu * T + *pl++;
	  su = su * T + *pl++;
	}
      sl += cu * cv + su * sv;
      /* Latitiude. */
      cu = *pb++;
      su = *pb++;
      for (ip = 0; ip < nt; ip++)
	{
	  cu = cu * T + *pb++;
	  su = su * T + *pb++;
	}
      sb += cu * cv + su * sv;
      /* Radius. */
      cu = *pr++;
      su = *pr++;
      for (ip = 0; ip < nt; ip++)
	{
	  cu = cu * T + *pr++;
	  su = su * T + *pr++;
	}
      sr += cu * cv + su * sv;
    }
  pobj[0] = STR * sl;
  pobj[1] = STR * sb;
  pobj[2] = STR * plan->distance * sr + plan->distance;
  return (0);
}


/* Prepare lookup table of sin and cos ( i*Lj )
 * for required multiple angles
 */
static int 
sscc (k, arg, n)
     int k;
     double arg;
     int n;
{
  double cu, su, cv, sv, s;
  int i;

  su = sin (arg);
  cu = cos (arg);
  ss[k][0] = su;		/* sin(L) */
  cc[k][0] = cu;		/* cos(L) */
  sv = 2.0 * su * cu;
  cv = cu * cu - su * su;
  ss[k][1] = sv;		/* sin(2L) */
  cc[k][1] = cv;
  for (i = 2; i < n; i++)
    {
      s = su * cv + cu * sv;
      cv = cu * cv - su * sv;
      sv = s;
      ss[k][i] = sv;		/* sin( i+1 L ) */
      cc[k][i] = cv;
    }
  return (0);
}

/* Compute mean elements at Julian date J.  */
static double Args[NARGS];
static double LP_equinox;
static double NF_arcsec;
static double Ea_arcsec;
static double pA_precession;

void
mean_elements (J)
     double J;
{
  double x, T, T2;

  /* Time variables.  T is in Julian centuries.  */
  T = (J - 2451545.0) / 36525.0;
  T2 = T*T;

  /* Mean longitudes of planets (Simon et al, 1994)
     .047" subtracted from constant term for offset to DE403 origin. */

  /* Mercury */
  x = mods3600( 538101628.6889819 * T + 908103.213 );
  x += (6.39e-6 * T
	 - 0.0192789) * T2;
  Args[0] = STR * x;

  /* Venus */
  x = mods3600( 210664136.4335482 * T + 655127.236 );
  x += (-6.27e-6  * T
	 + 0.0059381) * T2;
  Args[1] = STR * x;

  /* Earth  */
  x = mods3600( 129597742.283429 * T + 361679.198 );
  x += (-5.23e-6 * T
	 - 2.04411e-2 ) * T2;
  Ea_arcsec = x;
  Args[2] = STR * x;

  /* Mars */
  x = mods3600(  68905077.493988 * T +  1279558.751 );
  x += (-1.043e-5 * T
	 + 0.0094264) * T2;
  Args[3] = STR * x;

  /* Jupiter */
  x = mods3600( 10925660.377991 * T + 123665.420 );
  x += ((((-3.4e-10 * T
	    + 5.91e-8) * T
	   + 4.667e-6) * T
	  + 5.706e-5) * T
         - 3.060378e-1)*T2;
  Args[4] = STR * x;

   /* Saturn */
  x = mods3600( 4399609.855372 * T + 180278.752 );
  x += (((( 8.3e-10 * T
	  - 1.452e-7) * T
	  - 1.1484e-5) * T
	   - 1.6618e-4) * T
	 + 7.561614E-1)*T2;
  Args[5] = STR * x;

  /* Uranus */
  x = mods3600( 1542481.193933 * T + 1130597.971 )
       + (0.00002156*T - 0.0175083)*T2;
  Args[6] = STR * x;

  /* Neptune */
  x = mods3600( 786550.320744 * T + 1095655.149 )
       + (-0.00000895*T + 0.0021103)*T2;
  Args[7] = STR * x;

  /* Copied from cmoon.c, DE404 version.  */
  /* Mean elongation of moon = D */
  x = mods3600( 1.6029616009939659e+09 * T + 1.0722612202445078e+06 );
  x += (((((-3.207663637426e-013 * T
	       + 2.555243317839e-011) * T
	      + 2.560078201452e-009) * T
	     - 3.702060118571e-005) * T
            + 6.9492746836058421e-03) * T /* D, t^3 */
           - 6.7352202374457519e+00) * T2; /* D, t^2 */
  Args[9] = STR * x;

  /* Mean distance of moon from its ascending node = F */
  x = mods3600( 1.7395272628437717e+09 * T + 3.3577951412884740e+05 );
  x += ((((( 4.474984866301e-013 * T
	       + 4.189032191814e-011) * T
	       - 2.790392351314e-009) * T
	      - 2.165750777942e-006) * T
	      - 7.5311878482337989e-04) * T /* F, t^3 */
	     - 1.3117809789650071e+01) * T2; /* F, t^2 */
  NF_arcsec = x;
  Args[10] = STR * x;

/* Mean anomaly of sun = l' (J. Laskar) */
  x = mods3600(1.2959658102304320e+08 * T + 1.2871027407441526e+06);
  x += ((((((((
	       1.62e-20 * T
	       - 1.0390e-17 ) * T
	      - 3.83508e-15 ) * T
	     + 4.237343e-13 ) * T
	    + 8.8555011e-11 ) * T
	   - 4.77258489e-8 ) * T
	  - 1.1297037031e-5 ) * T
	 + 8.7473717367324703e-05) * T
	- 5.5281306421783094e-01) * T2;
  Args[11] = STR * x;

  /* Mean anomaly of moon = l */
  x = mods3600( 1.7179159228846793e+09 * T + 4.8586817465825332e+05 );
  x += (((((-1.755312760154e-012 * T
		+ 3.452144225877e-011) * T
		- 2.506365935364e-008) * T
	       - 2.536291235258e-004) * T
              + 5.2099641302735818e-02) * T /* l, t^3 */
             + 3.1501359071894147e+01) * T2; /* l, t^2 */
  Args[12] = STR * x;

  /* Mean longitude of moon, re mean ecliptic and equinox of date = L  */
  x = mods3600( 1.7325643720442266e+09 * T + 7.8593980921052420e+05);
  x += ((((( 7.200592540556e-014 * T
	     + 2.235210987108e-010) * T
	    - 1.024222633731e-008) * T
	   - 6.073960534117e-005) * T
	 + 6.9017248528380490e-03) * T /* L, t^3 */
	- 5.6550460027471399e+00) * T2; /* L, t^2 */
  LP_equinox = x;
  Args[13] = STR * x;

  /* Precession of the equinox  */
 x = ((((((((( -8.66e-20*T - 4.759e-17)*T
           + 2.424e-15)*T
           + 1.3095e-12)*T
           + 1.7451e-10)*T
           - 1.8055e-8)*T
           - 0.0000235316)*T
           + 0.000076)*T
           + 1.105414)*T
           + 5028.791959)*T;
  /* Moon's longitude re fixed J2000 equinox.  */
 /*
   Args[13] -= x;
 */
   pA_precession = STR * x;

  /* Free librations.  */
  /* longitudinal libration 2.891725 years */
  x = mods3600( 4.48175409e7 * T + 8.060457e5 );
  Args[14] = STR * x;
  /* libration P, 24.2 years */
  x = mods3600(  5.36486787e6 * T - 391702.8 );
  Args[15] = STR * x;

#if 0
  Args[16] = 0.0;
#endif
  /* libration W, 74.7 years. */
 x = mods3600( 1.73573e6 * T );
 Args[17] = STR * x;
}


/* Generic program to accumulate sum of trigonometric series
   in three variables (e.g., longitude, latitude, radius)
   of the same list of arguments.  */

int
g3plan (J, plan, pobj, objnum)
     double J;
     struct plantbl *plan;
     double pobj[];
     int objnum;
{
  int i, j, k, m, n, k1, ip, np, nt;
#ifdef _MSC_VER
  char FAR *p;
  long FAR *pl;
  long FAR *pb;
  long FAR *pr;
#else
  /* On some systems such as Silicon Graphics, char is unsigned
     by default.  */
#ifdef vax
  /* VAX CC rejects "signed char."  */
  char *p;
#else
#ifdef __STDC__
  signed char *p;
#else
  char *p;
#endif
#endif
  long *pl, *pb, *pr;
#endif
  double su, cu, sv, cv;
  double T, t, sl, sb, sr;

  mean_elements (J);
#if 0
  /* For librations, moon's longitude is sidereal.  */
  if (flag)
    Args[13] -= pA_precession;
#endif
  T = (J - J2000) / plan->timescale;
  n = plan->maxargs;
  /* Calculate sin( i*MM ), etc. for needed multiple angles.  */
  for (i = 0; i < n; i++)
    {
      if ((j = plan->max_harmonic[i]) > 0)
	{
	  sscc (i, Args[i], j);
	}
    }

  /* Point to start of table of arguments. */
  p = plan->arg_tbl;
  /* Point to tabulated cosine and sine amplitudes.  */
#ifdef _MSC_VER
  pl = (long FAR *) plan->lon_tbl;
  pb = (long FAR *) plan->lat_tbl;
  pr = (long FAR *) plan->rad_tbl;
#else
  pl = (long *) plan->lon_tbl;
  pb = (long *) plan->lat_tbl;
  pr = (long *) plan->rad_tbl;
#endif
  sl = 0.0;
  sb = 0.0;
  sr = 0.0;

  for (;;)
    {
      /* argument of sine and cosine */
      /* Number of periodic arguments. */
      np = *p++;
      if (np < 0)
	break;
      if (np == 0)
	{			/* It is a polynomial term.  */
	  nt = *p++;
	  /* "Longitude" polynomial (phi). */
	  cu = *pl++;
	  for (ip = 0; ip < nt; ip++)
	    {
	      cu = cu * T + *pl++;
	    }
	  /*	  sl +=  mods3600 (cu); */
	  sl += cu;
	  /* "Latitude" polynomial (theta). */
	  cu = *pb++;
	  for (ip = 0; ip < nt; ip++)
	    {
	      cu = cu * T + *pb++;
	    }
	  sb += cu;
	  /* Radius polynomial (psi). */
	  cu = *pr++;
	  for (ip = 0; ip < nt; ip++)
	    {
	      cu = cu * T + *pr++;
	    }
	  sr += cu;
	  continue;
	}
      k1 = 0;
      cv = 0.0;
      sv = 0.0;
      for (ip = 0; ip < np; ip++)
	{
	  /* What harmonic.  */
	  j = *p++;
	  /* Which planet.  */
	  m = *p++ - 1;
	  if (j)
	    {
/*	      k = abs (j); */
		if (j < 0)
			k = -j;
		else
			k = j;
	      k -= 1;
	      su = ss[m][k];	/* sin(k*angle) */
	      if (j < 0)
		su = -su;
	      cu = cc[m][k];
	      if (k1 == 0)
		{		/* set first angle */
		  sv = su;
		  cv = cu;
		  k1 = 1;
		}
	      else
		{		/* combine angles */	
		  t = su * cv + cu * sv;
		  cv = cu * cv - su * sv;
		  sv = t;
		}
	    }
	}
      /* Highest power of T.  */
      nt = *p++;
      /* Longitude. */
      cu = *pl++;
      su = *pl++;
      for (ip = 0; ip < nt; ip++)
	{
	  cu = cu * T + *pl++;
	  su = su * T + *pl++;
	}
      sl += cu * cv + su * sv;
      /* Latitiude. */
      cu = *pb++;
      su = *pb++;
      for (ip = 0; ip < nt; ip++)
	{
	  cu = cu * T + *pb++;
	  su = su * T + *pb++;
	}
      sb += cu * cv + su * sv;
      /* Radius. */
      cu = *pr++;
      su = *pr++;
      for (ip = 0; ip < nt; ip++)
	{
	  cu = cu * T + *pr++;
	  su = su * T + *pr++;
	}
      sr += cu * cv + su * sv;
    }
  t = plan->trunclvl;
  pobj[0] = Args[objnum - 1] + STR * t * sl;
  pobj[1] = STR * t * sb;
  pobj[2] = plan->distance * (1.0 + STR * t * sr);
  return (0);
}


/* Generic program to accumulate sum of trigonometric series
   in two variables (e.g., longitude, radius)
   of the same list of arguments.  */
int 
g2plan (J, plan, pobj)
     double J;
     struct plantbl *plan;
     double pobj[];
{
  int i, j, k, m, n, k1, ip, np, nt;
#if _MSC_VER
  char FAR *p;
  long FAR *pl;
  long FAR *pr;
#else
  /* On some systems such as Silicon Graphics, char is unsigned
     by default.  */
#ifdef vax
  /* VAX CC rejects "signed char."  */
  char *p;
#else
#ifdef __STDC__
  signed char *p;
#else
  char *p;
#endif
#endif
  long *pl, *pr;
#endif
  double su, cu, sv, cv;
  double T, t, sl, sr;

  mean_elements (J);
#if 0
  /* For librations, moon's longitude is sidereal.  */
  if (flag)
    Args[13] -= pA_precession;
#endif
  T = (J - J2000) / plan->timescale;
  n = plan->maxargs;
  /* Calculate sin( i*MM ), etc. for needed multiple angles.  */
  for (i = 0; i < n; i++)
    {
      if ((j = plan->max_harmonic[i]) > 0)
	{
	  sscc (i, Args[i], j);
	}
    }

  /* Point to start of table of arguments. */
  p = plan->arg_tbl;
  /* Point to tabulated cosine and sine amplitudes.  */
#ifdef _MSC_VER
  pl = (long FAR *) plan->lon_tbl;
  pr = (long FAR *) plan->rad_tbl;
#else
  pl = (long *) plan->lon_tbl;
  pr = (long *) plan->rad_tbl;
#endif
  sl = 0.0;
  sr = 0.0;

  for (;;)
    {
      /* argument of sine and cosine */
      /* Number of periodic arguments. */
      np = *p++;
      if (np < 0)
	break;
      if (np == 0)
	{			/* It is a polynomial term.  */
	  nt = *p++;
	  /* Longitude polynomial. */
	  cu = *pl++;
	  for (ip = 0; ip < nt; ip++)
	    {
	      cu = cu * T + *pl++;
	    }
	  /*	  sl +=  mods3600 (cu); */
	  sl += cu;
	  /* Radius polynomial. */
	  cu = *pr++;
	  for (ip = 0; ip < nt; ip++)
	    {
	      cu = cu * T + *pr++;
	    }
	  sr += cu;
	  continue;
	}
      k1 = 0;
      cv = 0.0;
      sv = 0.0;
      for (ip = 0; ip < np; ip++)
	{
	  /* What harmonic.  */
	  j = *p++;
	  /* Which planet.  */
	  m = *p++ - 1;
	  if (j)
	    {
/*	      k = abs (j); */
		if (j < 0)
			k = -j;
		else
			k = j;
	      k -= 1;
	      su = ss[m][k];	/* sin(k*angle) */
	      if (j < 0)
		su = -su;
	      cu = cc[m][k];
	      if (k1 == 0)
		{		/* set first angle */
		  sv = su;
		  cv = cu;
		  k1 = 1;
		}
	      else
		{		/* combine angles */
		  t = su * cv + cu * sv;
		  cv = cu * cv - su * sv;
		  sv = t;
		}
	    }
	}
      /* Highest power of T.  */
      nt = *p++;
      /* Longitude. */
      cu = *pl++;
      su = *pl++;
      for (ip = 0; ip < nt; ip++)
	{
	  cu = cu * T + *pl++;
	  su = su * T + *pl++;
	}
      sl += cu * cv + su * sv;
      /* Radius. */
      cu = *pr++;
      su = *pr++;
      for (ip = 0; ip < nt; ip++)
	{
	  cu = cu * T + *pr++;
	  su = su * T + *pr++;
	}
      sr += cu * cv + su * sv;
    }
  t = plan->trunclvl;
  pobj[0] = t * sl;
  pobj[2] = t * sr;
  return (0);
}



/* Generic program to accumulate sum of trigonometric series
   in one variable.  */

double
g1plan (J, plan)
     double J;
     struct plantbl *plan;
{
  int i, j, k, m, k1, ip, np, nt;
#ifdef _MSC_VER
  char FAR *p;
  long FAR *pl;
#else
  /* On some systems such as Silicon Graphics, char is unsigned
     by default.  */
#ifdef vax
  /* VAX CC rejects "signed char."  */
  char *p;
#else
#ifdef __STDC__
  signed char *p;
#else
  char *p;
#endif
#endif
  long *pl;
#endif
  double su, cu, sv, cv;
  double T, t, sl;

  T = (J - J2000) / plan->timescale;
  mean_elements (J);
  /* Calculate sin( i*MM ), etc. for needed multiple angles.  */
  for (i = 0; i < NARGS; i++)
    {
      if ((j = plan->max_harmonic[i]) > 0)
	{
	  sscc (i, Args[i], j);
	}
    }

  /* Point to start of table of arguments. */
  p = plan->arg_tbl;
  /* Point to tabulated cosine and sine amplitudes.  */
#ifdef _MSC_VER
  pl = (long FAR *) plan->lon_tbl;
#else
  pl = (long *) plan->lon_tbl;
#endif
  sl = 0.0;

  for (;;)
    {
      /* argument of sine and cosine */
      /* Number of periodic arguments. */
      np = *p++;
      if (np < 0)
	break;
      if (np == 0)
	{			/* It is a polynomial term.  */
	  nt = *p++;
	  cu = *pl++;
	  for (ip = 0; ip < nt; ip++)
	    {
	      cu = cu * T + *pl++;
	    }
	  /*	  sl +=  mods3600 (cu); */
	  sl += cu;
	  continue;
	}
      k1 = 0;
      cv = 0.0;
      sv = 0.0;
      for (ip = 0; ip < np; ip++)
	{
	  /* What harmonic.  */
	  j = *p++;
	  /* Which planet.  */
	  m = *p++ - 1;
	  if (j)
	    {
/*	      k = abs (j); */
		if (j < 0)
			k = -j;
		else
			k = j;
	      k -= 1;
	      su = ss[m][k];	/* sin(k*angle) */
	      if (j < 0)
		su = -su;
	      cu = cc[m][k];
	      if (k1 == 0)
		{		/* set first angle */
		  sv = su;
		  cv = cu;
		  k1 = 1;
		}
	      else
		{		/* combine angles */
		  t = su * cv + cu * sv;
		  cv = cu * cv - su * sv;
		  sv = t;
		}
	    }
	}
      /* Highest power of T.  */
      nt = *p++;
      /* Cosine and sine coefficients.  */
      cu = *pl++;
      su = *pl++;
      for (ip = 0; ip < nt; ip++)
	{
	  cu = cu * T + *pl++;
	  su = su * T + *pl++;
	}
      sl += cu * cv + su * sv;
    }
  return (plan->trunclvl * sl);
}


/* Compute geocentric moon.  */

extern struct plantbl moonlr, moonlat;
extern double coseps, sineps;

int 
gmoon (J, rect, pol)
     double J;
     double rect[], pol[];
{
  double x, cosB, sinB, cosL, sinL;

  g2plan (J, &moonlr, pol);
  x = pol[0];
  x += LP_equinox;
  if (x < -6.48e5)
    x += 1.296e6;
  if (x > 6.48e5)
    x -= 1.296e6;
  pol[0] = STR * x;
  x = g1plan (J, &moonlat);
  pol[1] = STR * x;
  x = (1.0 + STR * pol[2]) * moonlr.distance;
  pol[2] = x;
  /* Convert ecliptic polar to equatorial rectangular coordinates.  */
  epsiln(J);
  cosB = cos(pol[1]);
  sinB = sin(pol[1]);
  cosL = cos(pol[0]);
  sinL = sin(pol[0]);
  rect[0] = cosB * cosL * x;
  rect[1] = (coseps * cosB * sinL - sineps * sinB) * x;
  rect[2] = (sineps * cosB * sinL + coseps * sinB) * x;
  return 0;
}
