/* Diurnal aberration
 * This formula is less rigorous than the method used for
 * annual aberration.  However, the correction is small.
 */

#include "kep.h"

/* distance from observer to center of Earth, in Earth radii: */
extern double trho;

/* geocentric latitude of observer, degrees: */
extern double tlat;

int diurab( last, ra, dec )
double last;	/* local apparent sidereal time, radians */
double *ra;	/* right ascension, radians */
double *dec;	/* declination, radians */
{
double lha, coslha, sinlha, cosdec, sindec;
double coslat, N, D;

lha = last - *ra;
coslha = cos(lha);
sinlha = sin(lha);
cosdec = cos(*dec);
sindec = sin(*dec);
coslat = cos( DTR*tlat );

if( cosdec != 0.0 )
	N = 1.5472e-6*trho*coslat*coslha/cosdec;
else
	N = 0.0;
*ra += N;

D = 1.5472e-6*trho*coslat*sinlha*sindec;
*dec += D;

if( prtflg )
	{
printf( "diurnal aberration dRA %.3fs dDec %.2f\"\n", RTS*N/15.0, RTS*D );
	}
return(0);
}
