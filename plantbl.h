
#ifdef _MSC_VER
#if _MSC_VER < 1000
/* Microsoft C version 7 or earlier */
#define FAR far
#else
/* Microsoft Visual C, 32 bit compiler */
#define FAR 
#endif
#else
#ifdef __BORLANDC__
#if __BORLANDC__ < 0x0550
/* 16-bit Borland compiler */
#define FAR __far
#else
#define FAR 
#endif
#else
#define FAR 
#endif
#endif

#define NARGS 18

#ifdef _MSC_VER
#define SIGNED
#else
  /* On some systems such as Silicon Graphics, char is unsigned
     by default.  */
#ifdef vax
  /* VAX CC rejects "signed char."  */
#define SIGNED
#else
#ifdef __STDC__
#define SIGNED signed
#else
#define SIGNED
#endif
#endif
#endif

struct plantbl {
  char maxargs;
  char max_harmonic[NARGS];
  char max_power_of_t;
  char SIGNED FAR *arg_tbl;
  void FAR *lon_tbl;
  void FAR *lat_tbl;
  void FAR *rad_tbl;
  double distance;
  double timescale;
  double trunclvl;
};
