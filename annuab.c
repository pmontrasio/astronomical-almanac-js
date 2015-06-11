/* Annual aberration - AA pages B17, B37, C24
 */

#include "kep.h"
extern double Clightaud;

int annuab( p )
double p[]; /* unit vector pointing from earth to object */
{
double A, B, C;
double betai, pV;
double x[3], V[3];
int i;

/* Calculate the velocity of the earth (see vearth.c).
 */
velearth( TDT );
betai = 0.0;
pV = 0.0;
for( i=0; i<3; i++ )
	{
	A = vearth[i]/Clightaud;
	V[i] = A;
	betai += A*A;
	pV += p[i] * A;
	}
/* Make the adjustment for aberration.
 */
betai = sqrt( 1.0 - betai );
C = 1.0 + pV;
A = betai/C;
B = (1.0  +  pV/(1.0 + betai))/C;

for( i=0; i<3; i++ )
	{
	C = A * p[i]  +  B * V[i];
	x[i] = C;
	dp[i] = C - p[i];
	}

showcor( "annual aberration", p, dp );
for( i=0; i<3; i++ )
	p[i] = x[i];
return(0);
}
