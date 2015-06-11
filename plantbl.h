
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

struct plantbl {
  char maxargs;
  char max_harmonic[NARGS];
  char max_power_of_t;
  char FAR *arg_tbl;
  void FAR *lon_tbl;
  void FAR *lat_tbl;
  void FAR *rad_tbl;
  double distance;
  double timescale;
  double trunclvl;
};
