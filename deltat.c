/* DeltaT = Ephemeris Time - Universal Time
 *
 * The tabulated values of deltaT, in hundredths of a second,
 * were taken from The Astronomical Almanac, page K8.  The program
 * adjusts for a value of secular tidal acceleration ndot. It is -25.8
 * arcsec per century squared for JPL's DE403 ephemeris.
 * ELP2000 and DE200 use the value -23.8946.
 *
 * The tabulated range is 1620.0 through 2003.0.  Bessel's interpolation
 * formula is implemented to obtain fourth order interpolated values at
 * intermediate times.
 *
 * Updated deltaT predictions can be obtained from this network archive:
 *    http://maia.usno.navy.mil
 * Currently (as of 2002) available series are
 *    tai-utc.dat  Changes by 1 whenever there is a leap second
 *    finals.all   EOP including UT1-UTC, always less than 1 second
 * from which deltaT = 32.184 + (tai-utc) - (UT1-UTC)
 *
 * For dates earlier than the tabulated range, the program
 * calculates approximate formulae of Stephenson and Morrison
 * or K. M. Borkowski.  These approximations have an estimated
 * error of 15 minutes at 1500 B.C.  They are not adjusted for small
 * improvements in the current estimate of ndot because the formulas
 * were derived from studies of ancient eclipses and other historical
 * information, whose interpretation depends only partly on ndot.
 *
 * A quadratic extrapolation formula, that agrees in value and slope with
 * current data, predicts future values of deltaT.
 *
 * Input Y is the Julian epoch expressed in Julian years.  Y can be
 * found from the Julian date JD by
 *     Y = 2000.0 + (JD - 2451545.0)/365.25.
 * See AA page B4.
 *
 * Output double deltat(Y) is ET-UT in seconds.
 *
 *
 * References:
 *
 * Stephenson, F. R., and L. V. Morrison, "Long-term changes
 * in the rotation of the Earth: 700 B.C. to A.D. 1980,"
 * Philosophical Transactions of the Royal Society of London
 * Series A 313, 47-70 (1984)
 *
 * Borkowski, K. M., "ELP2000-85 and the Dynamical Time
 * - Universal Time relation," Astronomy and Astrophysics
 * 205, L8-L10 (1988)
 * Borkowski's formula is derived from eclipses going back to 2137 BC
 * and uses lunar position based on tidal coefficient of -23.9 arcsec/cy^2.
 *
 * Chapront-Touze, Michelle, and Jean Chapront, _Lunar Tables
 * and Programs from 4000 B.C. to A.D. 8000_, Willmann-Bell 1991
 * Their table agrees with the one here, but the entries are
 * rounded to the nearest whole second.
 *
 * Stephenson, F. R., and M. A. Houlden, _Atlas of Historical
 * Eclipse Maps_, Cambridge U. Press (1986)
 */

#include "kep.h"

/* If the following number (read from the file aa.ini)
 * is nonzero, then the program will return it
 * and not calculate anything.
 */
double dtgiven = 0.0;
extern double dtgiven;


#define DEMO 0
#define TABSTART 1620
#define TABEND 2013
#define TABSIZ (TABEND - TABSTART + 1)

/* Note, Stephenson and Morrison's table starts at the year 1630.
 * The Chapronts' table does not agree with the Almanac prior to 1630.
 * The actual accuracy decreases rapidly prior to 1780.
 */
short dt[TABSIZ] = {
/* 1620.0 thru 1659.0 */
12400, 11900, 11500, 11000, 10600, 10200, 9800, 9500, 9100, 8800,
8500, 8200, 7900, 7700, 7400, 7200, 7000, 6700, 6500, 6300,
6200, 6000, 5800, 5700, 5500, 5400, 5300, 5100, 5000, 4900,
4800, 4700, 4600, 4500, 4400, 4300, 4200, 4100, 4000, 3800,
/* 1660.0 thru 1699.0 */
3700, 3600, 3500, 3400, 3300, 3200, 3100, 3000, 2800, 2700,
2600, 2500, 2400, 2300, 2200, 2100, 2000, 1900, 1800, 1700,
1600, 1500, 1400, 1400, 1300, 1200, 1200, 1100, 1100, 1000,
1000, 1000, 900, 900, 900, 900, 900, 900, 900, 900,
/* 1700.0 thru 1739.0 */
900, 900, 900, 900, 900, 900, 900, 900, 1000, 1000,
1000, 1000, 1000, 1000, 1000, 1000, 1000, 1100, 1100, 1100,
1100, 1100, 1100, 1100, 1100, 1100, 1100, 1100, 1100, 1100,
1100, 1100, 1100, 1100, 1200, 1200, 1200, 1200, 1200, 1200,
/* 1740.0 thru 1779.0 */
1200, 1200, 1200, 1200, 1300, 1300, 1300, 1300, 1300, 1300,
1300, 1400, 1400, 1400, 1400, 1400, 1400, 1400, 1500, 1500,
1500, 1500, 1500, 1500, 1500, 1600, 1600, 1600, 1600, 1600,
1600, 1600, 1600, 1600, 1600, 1700, 1700, 1700, 1700, 1700,
/* 1780.0 thru 1799.0 */
1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
1700, 1700, 1600, 1600, 1600, 1600, 1500, 1500, 1400, 1400,
/* 1800.0 thru 1819.0 */
1370, 1340, 1310, 1290, 1270, 1260, 1250, 1250, 1250, 1250,
1250, 1250, 1250, 1250, 1250, 1250, 1250, 1240, 1230, 1220,
/* 1820.0 thru 1859.0 */
1200, 1170, 1140, 1110, 1060, 1020, 960, 910, 860, 800,
750, 700, 660, 630, 600, 580, 570, 560, 560, 560,
570, 580, 590, 610, 620, 630, 650, 660, 680, 690,
710, 720, 730, 740, 750, 760, 770, 770, 780, 780,
/* 1860.0 thru 1899.0 */
788, 782, 754, 697, 640, 602, 541, 410, 292, 182,
161, 10, -102, -128, -269, -324, -364, -454, -471, -511,
-540, -542, -520, -546, -546, -579, -563, -564, -580, -566,
-587, -601, -619, -664, -644, -647, -609, -576, -466, -374,
/* 1900.0 thru 1939.0 */
-272, -154, -2, 124, 264, 386, 537, 614, 775, 913,
1046, 1153, 1336, 1465, 1601, 1720, 1824, 1906, 2025, 2095,
2116, 2225, 2241, 2303, 2349, 2362, 2386, 2449, 2434, 2408,
2402, 2400, 2387, 2395, 2386, 2393, 2373, 2392, 2396, 2402,
/* 1940.0 thru 1979.0 */
 2433, 2483, 2530, 2570, 2624, 2677, 2728, 2778, 2825, 2871,
 2915, 2957, 2997, 3036, 3072, 3107, 3135, 3168, 3218, 3268,
 3315, 3359, 3400, 3447, 3503, 3573, 3654, 3743, 3829, 3920,
 4018, 4117, 4223, 4337, 4449, 4548, 4646, 4752, 4853, 4959,
/* 1980.0 thru 2003.0 */
 5054, 5138, 5217, 5296, 5379, 5434, 5487, 5532, 5582, 5630,
 5686, 5757, 5831, 5912, 5998, 6078, 6163, 6230, 6297, 6347,
 6383, 6409, 6430, 6447,
/* Extrapolated values */
                         6456, 6600, 6700, 6800, 6900, 7000,
 7100, 7200, 7300, 7400,
};



double deltat(Y)
double Y;
{
double ans, p, B;
int d[6];
int i, iy, k;

if( dtgiven != 0.0 )
	return( dtgiven );

if( Y > TABEND )
	{
#if 0
/* Morrison, L. V. and F. R. Stephenson, "Sun and Planetary System"
 * vol 96,73 eds. W. Fricke, G. Teleki, Reidel, Dordrecht (1982)
 */
	B = 0.01*(Y-1800.0) - 0.1;
	ans = -15.0 + 32.5*B*B;
	return(ans);
#else
/* Extrapolate forward by a second-degree curve that agrees with
   the most recent data in value and slope, and vaguely fits
   over the past century.  This idea communicated by Paul Muller,
   who says NASA used to do something like it.  */
	B = Y - TABEND;
	/* slope */
	p = dt[TABSIZ-1] - dt[TABSIZ-2];
	/* square term */
	ans = (dt[TABSIZ - 101] - (dt[TABSIZ - 1] - 100.0 * p)) * 1e-4;
        ans = 0.01 * (dt[TABSIZ-1] + p * B + ans * B * B);
	if( prtflg )
	  printf("[extrapolated deltaT] ");
	return(ans);
#endif
	}

if( Y < TABSTART )
	{
	if( Y >= 948.0 )
		{
/* Stephenson and Morrison, stated domain is 948 to 1600:
 * 25.5(centuries from 1800)^2 - 1.9159(centuries from 1955)^2
 */
		B = 0.01*(Y - 2000.0);
		ans = (23.58 * B + 100.3)*B + 101.6;
		}
	else
		{
/* Borkowski */
		B = 0.01*(Y - 2000.0)  +  3.75;
		ans = 35.0 * B * B  +  40.;
		}
	return(ans);
	}

/* Besselian interpolation from tabulated values.
 * See AA page K11.
 */

/* Index into the table.
 */
p = floor(Y);
iy = (int) (p - TABSTART);
/* Zeroth order estimate is value at start of year
 */
ans = dt[iy];
k = iy + 1;
if( k >= TABSIZ )
	goto done; /* No data, can't go on. */

/* The fraction of tabulation interval
 */
p = Y - p;

/* First order interpolated value
 */
ans += p*(dt[k] - dt[iy]);
if( (iy-1 < 0) || (iy+2 >= TABSIZ) )
	goto done; /* can't do second differences */

/* Make table of first differences
 */
k = iy - 2;
for( i=0; i<5; i++ )
	{
	if( (k < 0) || (k+1 >= TABSIZ) )
		{
		d[i] = 0;
		}
	else
		d[i] = dt[k+1] - dt[k];
	k += 1;
	}

/* Compute second differences
 */
for( i=0; i<4; i++ )
	d[i] = d[i+1] - d[i];
B = 0.25*p*(p-1.0);
ans += B*(d[1] + d[2]);
#if DEMO
printf( "B %.4lf, ans %.4lf\n", B, ans );
#endif
if( iy+2 >= TABSIZ )
	goto done;

/* Compute third differences
 */
for( i=0; i<3; i++ )
	d[i] = d[i+1] - d[i];
B = 2.0*B/3.0;
ans += (p-0.5)*B*d[1];
#if DEMO
printf( "B %.4lf, ans %.4lf\n", B*(p-0.5), ans );
#endif
if( (iy-2 < 0) || (iy+3 > TABSIZ) )
	goto done;

/* Compute fourth differences
 */
for( i=0; i<2; i++ )
	d[i] = d[i+1] - d[i];
B = 0.125*B*(p+1.0)*(p-2.0);
ans += B*(d[0] + d[1]);
#if DEMO
printf( "B %.4lf, ans %.4lf\n", B, ans );
#endif

done:
/* Astronomical Almanac table is corrected by adding the expression
 *     -0.000091 (ndot + 26)(year-1955)^2  seconds
 * to entries prior to 1955 (AA page K8), where ndot is the secular
 * tidal term in the mean motion of the Moon.
 *
 * Entries after 1955 are referred to atomic time standards and
 * are not affected by errors in Lunar or planetary theory.
 */
ans *= 0.01;
if( Y < 1955.0 )
	{
	B = (Y - 1955.0);
#if 1
	ans += -0.000091 * (-25.8 + 26.0) * B * B;
#else
	ans += -0.000091 * (-23.8946 + 26.0) * B * B;
#endif
	}
return( ans );
}




/* Routine to fill in time values for TDT and UT.
 * These and the input date, JD are global variables.
 * jdflag is set up on reading the initialization file aa.ini.
 */

int update()
{
double T;

/* Convert Julian date to Julian years re J2000.0.
 */
T = 2000.0 + (JD - J2000)/365.25;

switch( jdflag )
	{
	case 0:
		TDT = JD;
		UT = JD;
		break;
	case 1:
		TDT = JD;
		UT = TDT - deltat(T)/86400.0;
		jtocal(UT); /* display the other date */
		break;
	case 2:
		UT = JD;
		TDT = UT + deltat(T)/86400.0;
		jtocal(TDT);
		break;
	}
jtocal(JD);
return(0);
}





/* Exercise program.
 */
#if DEMO
main()
{
char s[20];
double ans, y;

loop:
printf( "year ? " );
gets(s);
sscanf( s, "%lf", &y );
ans = deltat(y);
printf( "%.4lf\n", ans );
goto loop;
}
#endif
