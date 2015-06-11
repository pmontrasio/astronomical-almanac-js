/* Atmospheric refraction
 * Returns correction in degrees to be added to true altitude
 * to obtain apparent altitude.
 *
 * -- S. L. Moshier
 */

extern double atpress; /* millibars */
extern double attemp; /* degrees C */
extern double DTR; /* pi/180 */

#if __STDC__
double tan (double);
#else
double tan();
#endif

double refrac(alt)
double alt;	/* altitude in degrees */
{
double y, y0, D0, N, D, P, Q;
int i;

if( (alt < -2.0) || (alt >= 90.0) )
	return(0.0);

/* For high altitude angle, AA page B61
 * Accuracy "usually about 0.1' ".
 */
if( alt > 15.0 )
	{
	D = 0.00452*atpress/((273.0+attemp)*tan( DTR*alt ));
	return(D);
	}

/* Formula for low altitude is from the Almanac for Computers.
 * It gives the correction for observed altitude, so has
 * to be inverted numerically to get the observed from the true.
 * Accuracy about 0.2' for -20C < T < +40C and 970mb < P < 1050mb.
 */

/* Start iteration assuming correction = 0
 */
y = alt;
D = 0.0;
/* Invert Almanac for Computers formula numerically
 */
P = (atpress - 80.0)/930.0;
Q = 4.8e-3 * (attemp - 10.0);
y0 = y;
D0 = D;

for( i=0; i<4; i++ )
	{
	N = y + (7.31/(y+4.4));
	N = 1.0/tan(DTR*N);
	D = N*P/(60.0 + Q * (N + 39.0));
	N = y - y0;
	y0 = D - D0 - N; /* denominator of derivative */

	if( (N != 0.0) && (y0 != 0.0) )
/* Newton iteration with numerically estimated derivative */
		N = y - N*(alt + D - y)/y0;
	else
/* Can't do it on first pass */
		N = alt + D;

	y0 = y;
	D0 = D;
	y = N;
	}
return( D );
}
