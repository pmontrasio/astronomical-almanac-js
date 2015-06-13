#include <stdio.h>
#include "consts.h"

/* Main program starts here.
 */
int main(int argc, char* argv[])
{
  double distance;
  double julian_day;
  int planet;

  sscanf(argv[1], "%lf", &julian_day);
  sscanf(argv[2], "%d", &planet);
  distance = km_to_planet(julian_day, planet);
  printf("%f\n", distance);
}
