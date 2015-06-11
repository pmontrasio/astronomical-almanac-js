/* Diurnal parallax, AA page D3
 */


#include "kep.h"
/* Earth radii per au. */
static double DISFAC;
extern double aearth, au;
/* distance to center of Earth, in Earth radii */
extern double trho;

/* observer's geocentric latitude, degrees */
extern double tlat;


int diurpx( last, ra, dec, dist )
double last;  /* local apparent sidereal time, radians */
double *ra;  /* right ascension, radians */
double *dec; /* declination, radians */
double dist; /* Earth - object distance, au */
{
double cosdec, sindec, coslat, sinlat;
double p[3], dp[3], x, y, z, D;

/* Don't bother with this unless the equatorial horizontal parallax
 * is at least 0.005"
 */
if( dist > 1758.8 )
	return(-1);

DISFAC = au / (0.001 * aearth);
cosdec = cos(*dec);
sindec = sin(*dec);

/* Observer's astronomical latitude
 */
x = tlat * DTR;
coslat = cos(x);
sinlat = sin(x);

/* Convert to equatorial rectangular coordinates
 * in which unit distance = earth radius
 */
D = dist * DISFAC;
p[0] = D*cosdec*cos(*ra);
p[1] = D*cosdec*sin(*ra);
p[2] = D*sindec;

dp[0] = -trho*coslat*cos(last);
dp[1] = -trho*coslat*sin(last);
dp[2] = -trho*sinlat;

x = p[0] + dp[0];
y = p[1] + dp[1];
z = p[2] + dp[2];
D = x*x + y*y + z*z;
D = sqrt(D);	/* topocentric distance */

/* recompute ra and dec */
*ra = zatan2(x,y);
*dec = asin(z/D);
showcor( "diurnal parallax", p, dp );
return(0);
}
