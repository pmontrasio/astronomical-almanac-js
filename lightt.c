/* Correction for light time from object to earth
 * including gravitational retardation due to the Sun.
 * AA page B36.
 */

#include "kep.h"

int lightt( elemnt, q, e )
double e[], q[];	/* rectangular position vectors */
struct orbit *elemnt;	/* orbital elements of object q */
{
double p[3], p0[3], ptemp[3];
double P, Q, E, t, x, y;
int i, k;


/* save initial q-e vector for display */
for( i=0; i<3; i++ )
	{
	p0[i] = q[i] - e[i];
	}

E = 0.0;
for( i=0; i<3; i++ )
	E += e[i]*e[i];
E = sqrt(E);

for( k=0; k<2; k++ )
	{
	P = 0.0;
	Q = 0.0;
	for( i=0; i<3; i++ )
		{
		y = q[i];
		x = y - e[i];
		p[i] = x;
		Q += y * y;
		P += x * x;
		}
	P = sqrt(P);
	Q = sqrt(Q);
/* Note the following blows up if object equals sun. */
	t = (P + 1.97e-8 * log( (E+P+Q)/(E-P+Q) ) )/173.1446327;
	kepler( TDT-t, elemnt, q, ptemp );
	}

if( prtflg )
	printf( "light time %.4fm,  ", 1440.0*t );

/* Final object-earth vector and the amount by which it changed.
 */
for( i=0; i<3; i++ )
	{
	x = q[i] - e[i];
	p[i] = x;
	dp[i] = x - p0[i];
	}
showcor( "aberration", p0, dp );

/* Calculate dRA/dt and dDec/dt.
 * The desired correction of apparent coordinates is relative
 * to the equinox of date, but the coordinates here are
 * for J2000.  This introduces a slight error.
 *
 * Estimate object-earth vector t days ago.  We have
 * p(?) = q(J-t) - e(J), and must adjust to
 * p(J-t)  =  q(J-t) - e(J-t)  =  q(J-t) - (e(J) - Vearth * t)
 *         =  p(?) + Vearth * t.
 */
velearth(TDT);
for( i=0; i<3; i++ )
	p[i] += vearth[i]*t;

deltap( p, p0, &dradt, &ddecdt );  /* see dms.c */
dradt /= t;
ddecdt /= t;
return(0);
}
