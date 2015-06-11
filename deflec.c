/* Correct for light deflection due to solar gravitation.
 * AA page B37
 *
 * SE, pq, etc. were computed earlier by angles.c.
 */

#include "kep.h"

int relativity( p, q, e )
/* unit vector from earth to object:
 */
double p[];
/* heliocentric ecliptic rectangular coordinates
 * of earth and object:
 */
double q[], e[];
{
double C;
int i;

C = 1.974e-8/(SE*(1.0+qe));
for( i=0; i<3; i++ )
	{
	dp[i] = C*(pq*e[i]/SE - ep*q[i]/SO);
	p[i] += dp[i];
	}
if( prtflg )
	printf( "elongation from sun %.2f degrees, ", acos( -ep )/DTR );
showcor( "light defl.", p, dp );
return(0);
}

