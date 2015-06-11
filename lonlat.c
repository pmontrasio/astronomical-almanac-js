
/* Display ecliptic longitude and latitude for
 * equinox of date.  Input is equatorial rectangular
 * coordinates for equinox J2000.
 */

#include "kep.h"

int lonlat( pp, J, polar, ofdate )
double pp[], J, polar[];
int ofdate; /* 1 means precess from J2000 to date J.  */
{
double s[3], x, y, z, yy, zz, r;
int i;

/* Make local copy of position vector
 * and calculate radius.
 */
r = 0.0;
for( i=0; i<3; i++ )
	{
	x = pp[i];
	s[i] = x;
	r += x * x;
	}
r = sqrt(r);

/* Precess to equinox of date J
 */
if( ofdate )
	precess( s, J, -1 );

/* Convert from equatorial to ecliptic coordinates
 */
epsiln(J);
yy = s[1];
zz = s[2];
x  = s[0];
y  =  coseps * yy  +  sineps * zz;
z  = -sineps * yy  +  coseps * zz;

yy = zatan2( x, y );
zz = asin( z/r );

polar[0] = yy;
polar[1] = zz;
polar[2] = r;

if( prtflg == 0 )
	return(0);

printf( "ecliptic long" );
dms( yy );
printf( " lat" );
dms( zz );
printf( " rad %.6E\n", r );
return(0);
}
