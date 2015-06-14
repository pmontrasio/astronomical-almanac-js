/* Apply diurnal aberrations and calculate topocentric
 * altitude and azimuth, given the geocentric apparent
 * right ascension and declination.
 *
 * Ephemeris transit times can be obtained by modifying
 * aa.ini to declare TDT = UT.
 * 
 * This program assumes that deltaT, the difference
 * between Ephemeris Time (TDT) and Universal Time (UT),
 * has already been calculated.  The global variables
 * TDT and UT must already have been loaded with the
 * correct values of TDT and UT, respectively.  See deltat.c.
 *
 * Inputs are polar coordinates:
 * dist is the true geocentric distance in au.
 * ra and dec are in radians.
 *
 * J is the Julian date in UT measure.
 *
 * AA page B60 and D3.
 */

#include "kep.h"
extern double tlong, tlat, glat;


int altaz( pol, J )
double pol[3];
double J;
{
double dec, cosdec, sindec, lha, coslha, sinlha;
double ra, dist, last, alt, az, coslat, sinlat;
double N, D, x, y, z, TPI;

ra = pol[0];
dec = pol[1];
dist = pol[2];
TPI = 2.0*PI;

/* local apparent sidereal time, seconds converted to radians
 */
last = sidrlt( J, tlong ) * DTR/240.0;
lha = last - ra; /* local hour angle, radians */
if( prtflg )
	{
	printf( "Local apparent sidereal time " );
	hms( last );
	printf( "\n" );
	}
/* Display rate at which ra and dec are changing
 */
/*
 *if( prtflg )
 *	{
 *	x = RTS/24.0;
 *	N = x*dradt;
 *	D = x*ddecdt;
 *	if( N != 0.0 )
 *		printf( "dRA/dt %.2f\"/h, dDec/dt %.2f\"/h\n", N, D );
 *	}
 */

diurab( last, &ra, &dec );
/* Do rise, set, and transit times
   trnsit.c takes diurnal parallax into account,
   but not diurnal aberration.  */
lha = last - ra;
trnsit( J, lha, dec );

/* Diurnal parallax
 */
diurpx( last, &ra, &dec, dist );

/* Diurnal aberration
 */
/*diurab( last, &ra, &dec );*/

/* Convert ra and dec to altitude and azimuth
 */
cosdec = cos(dec);
sindec = sin(dec);
lha = last - ra;
coslha = cos(lha);
sinlha = sin(lha);

/* Use the geodetic latitude for altitude and azimuth */
x = DTR * glat;
coslat = cos(x);
sinlat = sin(x);

N = -cosdec * sinlha;
D =  sindec * coslat  -  cosdec * coslha * sinlat;
az = RTD * zatan2( D, N );
alt = sindec * sinlat  +  cosdec * coslha * coslat;
alt = RTD * asin(alt);

/* Correction for atmospheric refraction
 * unit = degrees
 */
D = refrac( alt );
alt += D;

/* Convert back to R.A. and Dec.
 */
y = sin(DTR*alt);
x = cos(DTR*alt);
z = cos(DTR*az);
sinlha = -x * sin(DTR*az);
coslha = y*coslat - x*z*sinlat;
sindec = y*sinlat + x*z*coslat;
lha = zatan2( coslha, sinlha );

y = ra; /* save previous values, before refrac() */
z = dec;
dec = asin( sindec );
ra = last - lha;
y = ra - y; /* change in ra */
while( y < -PI )
	y += TPI;
while( y > PI )
	y -= TPI;
y = RTS*y/15.0;
z = RTS*(dec - z);
if( prtflg )
	{
	printf( "atmospheric refraction %.3f deg  dRA %.3fs dDec %.2f\"\n",
		 D, y, z );
	printf( "Topocentric:  Altitude %.3f deg, ", alt );
	printf( "Azimuth %.3f deg\n", az );

	printf( "Topocentric: R.A." );
	hms( ra );
	printf( " Dec." );
	dms( dec );
	printf( "\n" );
	}
return(0);
}
