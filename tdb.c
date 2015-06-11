/* Find Barycentric Dynamical Time from Terrestrial Dynamical Time.
   Reference: Astronomical Almanac, page B5.  */

/* radians per arc second */
#define STR 4.8481368110953599359e-6

/* 2000 January 1.5 */
#define J2000 2451545.0

#if __STDC__
double sin (double);
double floor (double);
#else
double sin(), floor();
#endif

/* Argument JED is a Julian date, in TDT.
   Output is the corresponding date in TDB.  */

double tdb(JED)
double JED;
{
double M, T;

/* Find time T in Julian centuries from J2000.  */
T = (JED - J2000)/36525.0;

/* Mean anomaly of sun = l' (J. Laskar) */
M = 129596581.038354 * T +  1287104.76154;

/* Reduce arc seconds mod 360 degrees.  */
M = M - 1296000.0 * floor( M/1296000.0 );

M += ((((((((
  1.62e-20 * T
- 1.0390e-17 ) * T
- 3.83508e-15 ) * T
+ 4.237343e-13 ) * T
+ 8.8555011e-11 ) * T
- 4.77258489e-8 ) * T
- 1.1297037031e-5 ) * T
+ 1.4732069041e-4 ) * T
- 0.552891801772 ) * T * T;
M *= STR;
/* TDB - TDT, in seconds.  */
T = 0.001658 * sin(M) + 0.000014 * sin(M+M);

T = JED + T / 86400.0;
return(T);
}
