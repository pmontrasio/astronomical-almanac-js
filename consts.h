/* November 28, 4057 B.C.*/
/* #define STARTDATE  239935.5 */

#ifdef __EMSCRIPTEN__
#define FAR
#endif

#define STARTDATE  625379.5

/* January 31, 1982 */
#define ENDDATE    2445000.5

/* Define one of these nonzero, the others zero.  */
#define NEWMOON 1
#define FULLMOON 0
#define FIRST_QUARTER_MOON 0
#define THIRD_QUARTER_MOON 0
#define EQUINOX 0

/* Conversion factors between degrees and radians */
extern double DTR;
extern double RTD;
extern double RTS;
extern double STR;
extern double PI;
extern double TPI;

extern double PI;

/* Standard epochs.  Note Julian epochs (J) are measured in
 * years of 365.25 days.
 */
extern double J2000;
extern double B1950;
extern double J1900;

double sqrt(), asin(), log();

/* coordinates of object
 */
extern int objnum;

/* ecliptic polar coordinates:
 * longitude, latitude in radians
 * radius in au
 */
extern double FAR obpolar[3];


extern struct orbit objorb;

extern double appobj[3];


/* coordinates of Earth
 */
/* Heliocentric rectangular equatorial position
 * of the earth at time TDT re equinox J2000
 */
extern double FAR rearth[3];
extern double vearth[3];

/* Corresponding polar coordinates of earth:
 * longitude and latitude in radians, radius in au
 */
extern double FAR eapolar[3];

extern struct orbit earth;

/* Julian date of ephemeris */
extern double JD;
extern double TDT;
extern double UT;

/* flag = 0 if TDT assumed = UT,
 *      = 1 if input time is TDT,
 *      = 2 if input time is UT.
 */
extern int jdflag;

/* correction vector, saved for display  */
extern double dp[3];

/* display enable flag */
extern int prtflg;

extern double dradt, ddecdt;

double zgetdate(), gethms();
double func(), search();
int apparent();

extern struct orbit *elobject;
extern double robject[3];

extern double km_to_planet(double, int);
extern double planet_distance(struct orbit *);
