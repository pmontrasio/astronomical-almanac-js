/* Calculate and display apparent coordinates of the sun at the
 * time given by the external variables TDT, UT.
 * Before calling this routine, the geometric position of the
 * earth must be calculated and put in rearth[].
 */

#include "kep.h"
extern struct orbit earth;

int dosun()
{
double r, x, y, t;
double ecr[3], rec[3], pol[3];
int i;
double asin(), modtp(), sqrt(), cos(), sin();


/* Display ecliptic longitude and latitude.
 */
for( i=0; i<3; i++ )
	ecr[i] = -rearth[i];
r = eapolar[2];

if( prtflg )
	{
	lonlat( ecr, TDT, pol, 1 );
	}

/* Philosophical note: the light time correction really affects
 * only the Sun's barymetric position; aberration is due to
 * the speed of the Earth.  In Newtonian terms the aberration
 * is the same if the Earth is standing still and the Sun moving
 * or vice versa.  Thus the following is actually wrong, but it
 * differs from relativity only in about the 8th decimal.
 * It should be done the same way as the corresponding planetary
 * correction, however.
 */
pol[2] = r;
for( i=0; i<2; i++ )
	{
	t = pol[2]/173.1446327;
/* Find the earth at time TDT - t */
	kepler( TDT-t, &earth, ecr, pol );
	}
r = pol[2];

for( i=0; i<3; i++ )
	{
	x = -ecr[i];
	y = -rearth[i];
	ecr[i] = x;	/* position t days ago */
	rec[i] = y;	/* position now */
	pol[i] = y - x; /* change in position */
	}

if( prtflg )
	{
	printf( "light time %.4fm,  ", 1440.0*t );
	showcor( "aberration", ecr, pol );
	}

/* Estimate rate of change of RA and Dec
 * for use by altaz().
 */
deltap( ecr, rec, &dradt, &ddecdt );  /* see dms.c */
dradt /= t;
ddecdt /= t;


/* There is no light deflection effect.
 * AA page B39.
 */

/* precess to equinox of date
 */
precess( ecr, TDT, -1 );

for( i=0; i<3; i++ )
    rec[i] = ecr[i];

/* Nutation.
 */
epsiln( TDT );
nutate( TDT, ecr );

/* Display the final apparent R.A. and Dec.
 * for equinox of date.
 */
if( prtflg )
	printf ("%s.", whatconstel (ecr, TDT));
showrd( "    Apparent:", ecr, pol );

/* Show it in ecliptic coordinates */
if( prtflg )
	{
	y  =  coseps * rec[1]  +  sineps * rec[2];
	y = zatan2( rec[0], y ) + nutl;
	printf( "Apparent longitude %.3f deg\n", RTD*y );
	}

/* Report altitude and azimuth
 */

altaz( pol, UT );
return(0);
}
