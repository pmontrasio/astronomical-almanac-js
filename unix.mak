CC= gcc
CFLAGS=  -O2 -Wall
OBJS = altaz.o angles.o annuab.o constel.o deflec.o deltat.o diurab.o \
diurpx.o dms.o epsiln.o fk4fk5.o kepler.o kfiles.o lightt.o lonlat.o \
nutate.o precess.o refrac.o rplanet.o rstar.o sidrlt.o sun.o domoon.o \
trnsit.o vearth.o zatan2.o \
gplan.o mer404.o ven404.o ear404.o mar404.o jup404.o \
sat404.o ura404.o nep404.o plu404.o mlr404.o mlat404.o

INCS = kep.h plantbl.h

aa: aa.o $(OBJS) $(INCS)
	$(CC) -o aa aa.o $(OBJS) -lm
#	coff2exe aa

aa.o: aa.c $(INCS)

conjunct: conjunct.o $(OBJS) $(INCS)
	$(CC) -o conjunct conjunct.o $(OBJS) -lm
#	coff2exe conjunct

conjunct.o: conjunct.c $(INCS)

moonrise: moonrise.o $(OBJS) $(INCS)
	$(CC) -o moonrise moonrise.o $(OBJS) -lm

moonrise.o: moonrise.c $(INCS)

altaz.o: altaz.c $(INCS)

angles.o: angles.c $(INCS)

annuab.o: annuab.c $(INCS)

constel.o: constel.c $(INCS)

deflec.o: deflec.c $(INCS)

deltat.o: deltat.c $(INCS)

diurab.o: diurab.c $(INCS)

diurpx.o: diurpx.c $(INCS)

dms.o: dms.c $(INCS)

epsiln.o: epsiln.c $(INCS)

fk4fk5.o: fk4fk5.c $(INCS)

kepler.o: kepler.c $(INCS)

kfiles.o: kfiles.c $(INCS)

lightt.o: lightt.c $(INCS)

lonlat.o: lonlat.c $(INCS)

nutate.o: nutate.c $(INCS)

precess.o: precess.c $(INCS)

refrac.o: refrac.c $(INCS)

rplanet.o: rplanet.c $(INCS)

rstar.o: rstar.c $(INCS)

sidrlt.o: sidrlt.c $(INCS)

sun.o: sun.c $(INCS)

trnsit.o: trnsit.c $(INCS)

vearth.o: vearth.c $(INCS)

zatan2.o: zatan2.c $(INCS)

domoon.o: domoon.c $(INCS)

ear404.o: ear404.c $(INCS)

jup404.o: jup404.c $(INCS)

mar404.o: mar404.c $(INCS)

mer404.o: mer404.c $(INCS)

nep404.o: nep404.c $(INCS)

sat404.o: sat404.c $(INCS)

ura404.o: ura404.c $(INCS)

ven404.o: ven404.c $(INCS)

plu404.o: plu404.c $(INCS)

mlr404.o: mlr404.c $(INCS)

mlat404.o: mlat404.c $(INCS)
