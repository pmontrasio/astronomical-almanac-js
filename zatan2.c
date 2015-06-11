/*							atan2()
 *
 *	Quadrant correct inverse circular tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, z, atan2();
 *
 * z = atan2( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns radian angle between 0 and +2pi whose tangent
 * is y/x.
 *
 *
 *
 * ACCURACY:
 *
 * See atan.c.
 *
 */


/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
Certain routines from the Library, including this one, may
be used and distributed freely provided this notice is retained
and source code is included with all distributions.
*/

#if __STDC__
double atan (double);
#else
double atan();
#endif

extern double PI;

double zatan2( x, y )
double x, y;
{
double z, w;
short code;


code = 0;

if( x < 0.0 )
	code = 2;
if( y < 0.0 )
	code |= 1;

if( x == 0.0 )
	{
	if( code & 1 )
		return( 1.5*PI );
	if( y == 0.0 )
		return( 0.0 );
	return( 0.5*PI );
	}

if( y == 0.0 )
	{
	if( code & 2 )
		return( PI );
	return( 0.0 );
	}


switch( code )
	{
	default:
	case 0: w = 0.0; break;
	case 1: w = 2.0 * PI; break;
	case 2:
	case 3: w = PI; break;
	}

z = atan( y/x );

return( w + z );
}
