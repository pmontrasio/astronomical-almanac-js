
/* This is based on reduce (rplanet.c) and returns the distance
 * in AU of the planet from Earth
 */

#include "consts.h"
#include "kep.h"

double distance( planet, q, e )
     struct orbit *planet;	/* orbital elements of q */
     double q[], e[];	/* heliocentric coordinates */
{
  double p[3], temp[3];
  double a, b;
  int i;
  double sqrt();

  /* Save the geometric coordinates at TDT
   */
  for( i=0; i<3; i++ )
    temp[i] = q[i];

  /* Adjust for light time (planetary aberration)
   */
  lightt( planet, q, e );

  /* Find Euclidean vectors between earth, object, and the sun
   */
  for( i=0; i<3; i++ )
    p[i] = q[i] - e[i];

  angles( p, q, e );

  a = 0.0;
  for( i=0; i<3; i++ )
    {
      b = temp[i] - e[i];
      a += b * b;
    }
  a = sqrt(a);
  return a;
}


/*extern double robject[];*/

double planet_distance(struct orbit *planet)
{
  /* calculate heliocentric position of the object */
  kepler( TDT, planet, robject, obpolar );
  /* apply correction factors and print apparent place */
  return distance( planet, robject, rearth );
}
