/* This program reduces catalogue coordinates of a star
 * to apparent place.
 *
 * - S. L. Moshier, November, 1987
 */

#include "kep.h"

int rstar( el )
struct star *el;
{
double p[3], q[3], e[3], m[3], temp[3], polar[3];
double T, vpi, epoch;
double cosdec, sindec, cosra, sinra;
int i;
double log();


/* Convert from RA and Dec to equatorial rectangular direction
 */
loop:
cosdec = cos( el->dec );
sindec = sin( el->dec );
cosra = cos( el->ra );
sinra = sin( el->ra );
q[0] = cosra * cosdec;
q[1] = sinra * cosdec;
q[2] = sindec;

/* space motion */
vpi = 21.094952663 * el->v * el->px;
m[0] =    -el->mura * cosdec * sinra
	- el->mudec * sindec * cosra
	+       vpi * q[0];

m[1] =     el->mura * cosdec * cosra
	- el->mudec * sindec * sinra
	+       vpi * q[1];

m[2] =    el->mudec * cosdec
	+       vpi * q[2];

epoch = el->epoch;

/* Convert FK4 to FK5 catalogue */

if( epoch == B1950 )
	{
	fk4fk5( q, m, el );
	goto loop;
	}


for( i=0; i<3; i++ )
	e[i] = rearth[i];

/* precess the earth to the star epoch */
precess( e, epoch, -1 );

/* Correct for proper motion and parallax
 */
T = (TDT - epoch)/36525.0;
for( i=0; i<3; i++ )
	{
	p[i] = q[i]  +  T * m[i]  -  el->px * e[i];
	}

/* precess the star to J2000 */
precess( p, epoch, 1 );
/* reset the earth to J2000 */
for( i=0; i<3; i++ )
	e[i] = rearth[i];

/* Find Euclidean vectors between earth, object, and the sun
 * angles( p, q, e );
 */
angles( p, p, e );

/* Find unit vector from earth in direction of object
 */
for( i=0; i<3; i++ )
	{
	p[i] /= EO;
	temp[i] = p[i];
	}

if( prtflg )
	{
	printf( "approx. visual magnitude %.1f\n", el->mag );
/* Report astrometric position
 */
	showrd( "Astrometric J2000.0:", p, polar );

/* Also in 1950 coordinates
 */
	precess( temp, B1950, -1 );
	showrd( "Astrometric B1950.0:", temp, polar );

/* For equinox of date: */
	for( i=0; i<3; i++ )
		temp[i] = p[i];
	precess( temp, TDT, -1 );
	showrd( "Astrometric of date:", temp, polar );
	}

/* Correct position for light deflection
 * relativity( p, q, e );
 */
relativity( p, p, e );


/* Correct for annual aberration
 */
annuab( p );

/* Precession of the equinox and ecliptic
 * from J2000.0 to ephemeris date
 */
precess( p, TDT, -1 );

/* Ajust for nutation
 * at current ecliptic.
 */
epsiln( TDT );
nutate( TDT, p );

/* Display the final apparent R.A. and Dec.
 * for equinox of date.
 */
showrd( "    Apparent:", p, polar );

/* Go do topocentric reductions.
 */
dradt = 0.0;
ddecdt = 0.0;
polar[2] = 1.0e38; /* make it ignore diurnal parallax */
altaz( polar, UT );
return(0);
}


/* Compute stellar ephemeris for ELOBJECT at time TDT.  */
 
extern struct orbit *elobject;

int dostar ()
{
  rstar( (struct star *) elobject );
  return 0;
}
