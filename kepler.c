 
/* Program to solve Keplerian orbit
 * given orbital parameters and the time.
 * Returns Heliocentric equatorial rectangular coordinates of
 * the object.
 *
 * This program detects several cases of given orbital elements.
 * If a program for perturbations is pointed to, it is called
 * to calculate all the elements.
 * If there is no program, then the mean longitude is calculated
 * from the mean anomaly and daily motion.
 * If the daily motion is not given, it is calculated
 * by Kepler's law.
 * If the eccentricity is given to be 1.0, it means that
 * meandistance is really the perihelion distance, as in a comet
 * specification, and the orbit is parabolic.
 *
 * Reference: Taff, L.G., "Celestial Mechanics, A Computational
 * Guide for the Practitioner."  Wiley, 1985.
 */

#include "kep.h"
#ifdef ANSIPROT
extern int embofs ( double J, double emb[], double *r );
int gplan (double, struct plantbl *, double *);
int g3plan (double, struct plantbl *, double *, int);
int gmoon (double, double *, double *);
#else
int embofs(), g3plan(), gplan(), gmoon();
#endif
extern struct orbit earth;	/* orbital elements of the earth */
extern double eps, coseps, sineps; /* obliquity of ecliptic */

int kepler(J, e, rect, polar)
double J, rect[], polar[];
struct orbit *e;
{
double alat, E, M, W, temp;
double epoch, inclination, ascnode, argperih;
double meandistance, dailymotion, eccent, meananomaly;
double r, coso, sino, cosa;

/* Call program to compute position, if one is supplied.  */
if( e->ptable )
	{
	if( e == &earth )
		g3plan (J, e->ptable, polar, 3);
	else
		gplan (J, e->ptable, polar);
	E = polar[0]; /* longitude */
	e->L = E;
	W = polar[1]; /* latitude */
	r = polar[2]; /* radius */
	e->r = r;
	e->epoch = J;
	e->equinox = J2000;
	goto kepdon;
	}

/* Decant the parameters from the data structure
 */
epoch = e->epoch;
inclination = e->i;
ascnode = e->W * DTR;
argperih = e->w;
meandistance = e->a; /* semimajor axis */
dailymotion = e->dm;
eccent = e->ecc;
meananomaly = e->M;
/* Check for parabolic orbit. */
if( eccent == 1.0 )
	{
/* meandistance = perihelion distance, q
 * epoch = perihelion passage date
 */
	temp = meandistance * sqrt(meandistance);
	W = (J - epoch ) * 0.0364911624 / temp;
/* The constant above is 3 k / sqrt(2),
 * k = Gaussian gravitational constant = 0.01720209895 . */
	E = 0.0;
	M = 1.0;
	while( fabs(M) > 1.0e-11 )
		{
		temp = E * E;
		temp = (2.0 * E * temp + W)/( 3.0 * (1.0 + temp));
		M = temp - E;
		if( temp != 0.0 )
			M /= temp;
		E = temp;
		}
	r = meandistance * (1.0 + E * E );
	M = atan( E );
	M = 2.0 * M;
	alat = M + DTR*argperih;
	goto parabcon;
	}
if( eccent > 1.0 )
	{
/* The equation of the hyperbola in polar coordinates r, theta
 * is r = a(e^2 - 1)/(1 + e cos(theta))
 * so the perihelion distance q = a(e-1),
 * the "mean distance"  a = q/(e-1).
 */
	meandistance = meandistance/(eccent - 1.0);
	temp = meandistance * sqrt(meandistance);
	W = (J - epoch ) * 0.01720209895 / temp;
/* solve M = -E + e sinh E */
	E = W/(eccent - 1.0);
	M = 1.0;
	while( fabs(M) > 1.0e-11 )
		{
		M = -E + eccent * sinh(E) - W;
		E += M/(1.0 - eccent * cosh(E));
		}
	r = meandistance * (-1.0 + eccent * cosh(E));
	temp = (eccent + 1.0)/(eccent - 1.0);
	M = sqrt(temp) * tanh( 0.5*E );
	M = 2.0 * atan(M);	
	alat = M + DTR*argperih;
	goto parabcon;
	}
/* Calculate the daily motion, if it is not given.
 */
if( dailymotion == 0.0 )
	{
/* The constant is 180 k / pi, k = Gaussian gravitational constant.
 * Assumes object in heliocentric orbit is massless.
 */
	dailymotion = 0.9856076686/(e->a*sqrt(e->a));
	}
dailymotion *= J - epoch;
/* M is proportional to the area swept out by the radius
 * vector of a circular orbit during the time between
 * perihelion passage and Julian date J.
 * It is the mean anomaly at time J.
 */
M = DTR*( meananomaly + dailymotion );
M = modtp(M);
/* If mean longitude was calculated, adjust it also
 * for motion since epoch of elements.
 */
if( e->L )
	{
	e->L += dailymotion;
	e->L = mod360( e->L );
	}

/* By Kepler's second law, M must be equal to
 * the area swept out in the same time by an
 * elliptical orbit of same total area.
 * Integrate the ellipse expressed in polar coordinates
 *     r = a(1-e^2)/(1 + e cosW)
 * with respect to the angle W to get an expression for the
 * area swept out by the radius vector.  The area is given
 * by the mean anomaly; the angle is solved numerically.
 * 
 * The answer is obtained in two steps.  We first solve
 * Kepler's equation
 *    M = E - eccent*sin(E)
 * for the eccentric anomaly E.  Then there is a
 * closed form solution for W in terms of E.
 */

E = M; /* Initial guess is same as circular orbit. */
temp = 1.0;
do
	{
/* The approximate area swept out in the ellipse */
	temp = E - eccent * sin(E)
/* ...minus the area swept out in the circle */
		- M;
/* ...should be zero.  Use the derivative of the error
 * to converge to solution by Newton's method.
 */
	E -= temp/(1.0 - eccent*cos(E));
	}
while( fabs(temp) > 1.0e-11 );

/* The exact formula for the area in the ellipse is
 *    2.0*atan(c2*tan(0.5*W)) - c1*eccent*sin(W)/(1+e*cos(W))
 * where
 *    c1 = sqrt( 1.0 - eccent*eccent )
 *    c2 = sqrt( (1.0-eccent)/(1.0+eccent) ).
 * Substituting the following value of W
 * yields the exact solution.
 */
temp = sqrt( (1.0+eccent)/(1.0-eccent) );
W = 2.0 * atan( temp * tan(0.5*E) );

/* The true anomaly.
 */
W = modtp(W);

meananomaly *= DTR;
/* Orbital longitude measured from node
 * (argument of latitude)
 */
if( e->L )
	alat = (e->L)*DTR + W - meananomaly - ascnode;
else
	alat = W + DTR*argperih; /* mean longitude not given */

/* From the equation of the ellipse, get the
 * radius from central focus to the object.
 */
r = meandistance*(1.0-eccent*eccent)/(1.0+eccent*cos(W));

parabcon:

/* The heliocentric ecliptic longitude of the object
 * is given by
 *   tan( longitude - ascnode )  =  cos( inclination ) * tan( alat ).
 */
coso = cos( alat );
sino = sin( alat );
inclination *= DTR;
W = sino * cos( inclination );
E = zatan2( coso, W ) + ascnode;

/* The ecliptic latitude of the object
 */
W = sino * sin( inclination );
W = asin(W);

kepdon:

/* Convert to rectangular coordinates,
 * using the perturbed latitude.
 */
rect[2] = r * sin(W);
cosa = cos(W);
rect[1] = r * cosa * sin(E);
rect[0] = r * cosa * cos(E);

/* Convert from heliocentric ecliptic rectangular
 * to heliocentric equatorial rectangular coordinates
 * by rotating eps radians about the x axis.
 */
epsiln( e->equinox );
W = coseps*rect[1] - sineps*rect[2];
M = sineps*rect[1] + coseps*rect[2];
rect[1] = W;
rect[2] = M;

/* Precess the position
 * to ecliptic and equinox of J2000.0
 * if not already there.
 */
precess( rect, e->equinox, 1 );

/* If earth, adjust from earth-moon barycenter to earth
 * by AA page E2.
 */
if( e == &earth )
	{
	embofs( J, rect, &r ); /* see below */
	}

/* Rotate back into the ecliptic.  */
epsiln( J2000 );
W = coseps*rect[1] + sineps*rect[2];
M = -sineps*rect[1] + coseps*rect[2];

/* Convert to polar coordinates */
E = zatan2( rect[0], W );
W = asin( M/r );

/* Output the polar cooordinates
 */
polar[0] = E; /* longitude */
polar[1] = W; /* latitude */
polar[2] = r; /* radius */

return(0);
}


/* Adjust position from Earth-Moon barycenter to Earth
 *
 * J = Julian day number
 * emb = Equatorial rectangular coordinates of EMB.
 * pr = Earth's distance to the Sun (au)
 */
extern double emrat;

int embofs( J, ea, pr )
double J;
double ea[];
double *pr;
{
double pm[3], polm[3];
double a, b;
int i;

/* Compute the vector Moon - Earth.  */
gmoon( J, pm, polm );

/* Precess the lunar position
 * to ecliptic and equinox of J2000.0
 */
precess( pm, J, 1 );

/* Adjust the coordinates of the Earth
 */
a = 1.0 / (emrat +  1.0);
b = 0.0;
for( i=0; i<3; i++ )
	{
	ea[i] = ea[i] - a * pm[i];
	b = b + ea[i] * ea[i];
	}
/* Sun-Earth distance.  */
*pr = sqrt(b);
return(0);
}
