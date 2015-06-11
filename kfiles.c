
/* Disc file input routines to read initialization parameters
 * or file containing orbital elements.
 */

#if __BORLANDC__
#include <stdlib.h>
#endif

#include "kep.h"

#ifdef ANSIPROT
#include <string.h>
#endif

#ifdef _MSC_VER
#if _MSC_VER >= 1000
#include <stdlib.h>
#include <string.h>
#endif
#endif

#ifndef ANSIPROT
FILE *fincat();
#endif

extern char *intfmt, *strfmt;/* see dms.c */

static char starnam[80] = {'s','t','a','r','.','c','a','t','\0'};
static char orbnam[80] = {'o','r','b','i','t','.','c','a','t','\0'};
static int linenum = 1;

/* Read initialization file aa.ini
 * and adjust topocentric coordinates of observer.
 *
 * The following items will be read in automatically from the disc file
 * named aa.ini, if one is provided.  The file contains one ASCII
 * string number per line so is easily edited.
 *
 * Terrestrial geocentric latitude and longitude of observer
 * in degrees East of Greenwich, North of equator
 * (kfiles.c converts from geodetic latitude input.)
 */
double tlong = -71.13;	/* Cambridge, Massachusetts */
double tlat = 42.38; /* geocentric */
double glat = 42.17; /* geodetic */
/* Parameters for calculation of azimuth and elevation
 */
double attemp = 20.0;	/* atmospheric temperature, degrees Centigrade */
double atpress = 1013.0; /* atmospheric pressure, millibars */
/* Distance from observer to center of earth, in earth radii
 */
double trho = 0.9985;
double flat = 298.257222;
double height = 0.0;

/* Constants used elsewhere. These are DE403 values. */
double aearth = 6378137.;  /* Radius of the earth, in meters.  */
double au = 1.49597870691e8; /* Astronomical unit, in kilometers.  */
double emrat = 81.300585;  /* Earth/Moon mass ratio.  */
double Clight = 2.99792458e5;  /* Speed of light, km/sec  */
double Clightaud; /* C in au/day  */

extern double tlong, tlat, glat, trho, attemp, atpress, dtgiven;
extern double Rearth;

int kinit()
{
double a, b, fl, co, si, u;
FILE *f, *fopen();
char s[84];

printf( "\n\tSteve Moshier's Ephemeris Program v5.6\n\n" );
printf( "Planetary and lunar positions approximate DE404.\n" );

f = fopen( "aa.ini", "r" );
if( f )
	{
	fgets( s, 80, f );
	sscanf( s, "%lf", &tlong );
	fgets( s, 80, f );
	sscanf( s, "%lf", &glat );
	fgets( s, 80, f );
	sscanf( s, "%lf", &height );
	u = glat * DTR;

/* Reduction from geodetic latitude to geocentric latitude
 * AA page K5
 */
	co = cos(u);
	si = sin(u);
	fl = 1.0 - 1.0/flat;
	fl = fl*fl;
	si = si*si;
	u = 1.0/sqrt( co*co + fl*si );
	a = aearth*u + height;
	b = aearth*fl*u  +  height;
	trho = sqrt( a*a*co*co + b*b*si );
	tlat = RTD * acos( a*co/trho );
	if( glat < 0.0 )
		tlat = -tlat;
	trho /= aearth;

/* Reduction from geodetic latitude to geocentric latitude
 * AA page K5
 */
/*
	tlat = glat
		- 0.19242861 * sin(2.0*u)
		+ 0.00032314 * sin(4.0*u)
		- 0.00000072 * sin(6.0*u);

	trho =    0.998327073
		+ 0.001676438 * cos(2.0*u)
		- 0.000003519 * cos(4.0*u)
		+ 0.000000008 * cos(6.0*u);
	trho += height/6378160.;
*/
	printf( "Terrestrial east longitude %.4f deg\n", tlong );
	printf( "geocentric latitude %.4f deg\n", tlat );
	printf( "Earth radius %.5f\n", trho );

	fgets( s, 80, f );
	sscanf( s, "%lf", &attemp );
	printf( "temperature %.1f C\n", attemp );
	fgets( s, 80, f );
	sscanf( s, "%lf", &atpress );
	printf( "pressure %.0f mb\n", atpress );
	fgets( s, 80, f );
	sscanf( s, "%d", &jdflag );
	switch( jdflag )
		{
		case 0: printf("TDT and UT assumed equal.\n");
			break;
		case 1: printf("Input time is TDT.\n" );
			break;
		case 2: printf("Input time is UT.\n" );
			break;
		default: printf("Illegal jdflag\n" );
		exit(0);
		}
	fgets( s, 80, f );
	sscanf( s, "%lf", &dtgiven );
	if( dtgiven != 0.0 )
		printf( "Using deltaT = %.2fs.\n", dtgiven );
	fclose(f);
	}
Clightaud = 86400.0 * Clight / au;
/* Radius of the earth in au
   Thanks to Min He <Min.He@businessobjects.com> for pointing out
   this needs to be initialized early.  */
Rearth = 0.001 * aearth / au;
return(0);
}



/* Program to read in a file containing orbital parameters
 */
extern struct orbit earth;

int getorbit(el)
struct orbit *el;
{
FILE *f;
char s1[128], s2[128], *u, *v;
int i;


getnum( "Name of orbit catalogue file: ", orbnam, strfmt );
f = fincat( orbnam, 2, s1, s2 );
if( f == 0 )
	return(-1); /* failure flag */


printf( "%s\n", s1 );
printf( "%s\n", s2 );

/* Read in ASCII floating point numbers
 */
sscanf( s1, "%lf %lf %lf %lf %lf %lf",
	&el->epoch, &el->i, &el->W, &el->w, &el->a, &el->dm );

sscanf( s2, "%lf %lf %lf %lf %lf %15s", &el->ecc, &el->M,
	&el->equinox, &el->mag, &el->sdiam, &el->obname[0] );

el->obname[15] = '\0';

/* Clear out the rest of the data structure
 */
el->ptable = 0;
el->L = 0.0;
el->r = 0.0;
el->plat = 0.0;
if( strcmp( &el->obname[0], "Earth" )  )
	{
	return(0);
	}
else
	{
	u = (char *)&earth;
	v = (char *)el;
	for( i=0; i < (int) sizeof(struct orbit); i++ )
		*u++ = *v++;
	printf( "Read in earth orbit\n" );
	return(1);
	}
}



int getstar(el)
struct star *el;
{
int sign;
char s[128];
double rh, rm, rs, dd, dm, ds, x, z;
FILE *f;
char *p;
int i;

getnum( "Name of star catalogue file: ", starnam, strfmt );
f = fincat( starnam, 1, s, (char *)0 );
if( f == 0 )
	return(-1); /* failure flag */


printf( "%s\n", s );
/* Read in the ASCII string data and name of the object
 */
sscanf( s, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s",
	&el->epoch, &rh, &rm, &rs, &dd, &dm, &ds,
	&el->mura, &el->mudec, &el->v, &el->px, &el->mag, &el->obname[0] );

x = el->epoch;
if( x == 2000.0 )
	x = J2000;
else if( x == 1950.0 )
	x = B1950;
else if( x == 1900.0 )
	x = J1900;
else
	x = J2000  +  365.25 * (x - 2000.0);
el->epoch = x;

/* read the right ascension */
el->ra = 2.0 * PI * (3600.0*rh + 60.0*rm + rs)/86400.0;

/* read the declination */
sign = 1;
if( (dd < 0.0) || (dm < 0.0) || (ds < 0.0) )
	sign = -1;
z = (3600.0*fabs(dd) + 60.0*fabs(dm) + fabs(ds))/RTS;
if( dd == 0.0 )
	{
/* Scan the text for possible minus sign in front of declination 0 */
	p = s;
/* skip over 4 fields */
	for( i=0; i<4; i++ )
		{
		while( *p++ == ' ' )
			;
		while( *p++ != ' ' )
			;
		}
	while( *p++ == ' ' )
		;
	--p;
	if( *p == '-' )
		sign = -1;
	}
if( sign < 0 )
	z = -z;
el->dec = z;

#if DEBUG
printf( "%.2f\n", el->epoch );
printf( "%.0f %.0f %.3f\n", rh, rm, rs );
printf( "%.8f\n", el->ra );
printf( "%.0f %.0f %.3f\n", dd, dm, ds );
printf( "%.8f\n", el->dec );
printf( "d %.3f mua %.3f mud %.3f v %.3f\n",
	el->px, el->mura, el->mudec, el->v );
#endif

el->mura *= 15.0/RTS;	/* s/century -> "/century -> rad/century */
el->mudec /= RTS;
z = el->px;
if( z < 1.0 )
	{
	if( z <= 0.0 )
		el->px = 0.0;
	else
		el->px = STR * z;  /* assume px in arc seconds */
	}
else
	{
	el->px = 1.0/(RTS * z);	/* parsecs -> radians */
	}
return(0);
}




/* Open catalogue and find line number
 */
FILE *fincat( name, n, str1, str2 )
char *name;
int n;   /* number of lines per catalogue entry */
char *str1, *str2;
{
int i;
FILE *f, *fopen();

f = fopen( name, "r" );
if( f == 0 )
	{
	printf( "Can't find file %s\n", name );
	return(0); /* failure flag */
	}

getnum( "Line number", &linenum, intfmt );
if( linenum <= 0 )
	goto failure;

for( i=0; i<linenum; i++ )
	{

	fgets( str1, 126, f );
	if( *str1 == '-' )
		goto endf;
	if( n > 1 )
		{
		fgets( str2, 126, f );
		if( *str2 == '-' )
			goto endf;
		}
	}
fclose(f);
return( f );


endf:
printf( "End of file reached.\n" );
failure:
fclose(f);
return(0);
}

