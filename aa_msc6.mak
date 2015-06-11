# MSDOS make file
CFLAGS=/c /AL
aa.obj: aa.c kep.h plantbl.h
	cl $(CFLAGS) aa.c

altaz.obj: altaz.c kep.h
	cl $(CFLAGS) altaz.c

angles.obj: angles.c kep.h
	cl $(CFLAGS) angles.c

annuab.obj: annuab.c kep.h
	cl $(CFLAGS) annuab.c

constel.obj: constel.c kep.h
	cl $(CFLAGS) constel.c

deflec.obj: deflec.c kep.h
	cl $(CFLAGS) deflec.c

deltat.obj: deltat.c kep.h
	cl $(CFLAGS) deltat.c

diurab.obj: diurab.c kep.h
	cl $(CFLAGS) diurab.c

diurpx.obj: diurpx.c kep.h
	cl $(CFLAGS) diurpx.c

dms.obj: dms.c kep.h
	cl $(CFLAGS) dms.c

epsiln.obj: epsiln.c kep.h
	cl $(CFLAGS) epsiln.c

fk4fk5.obj: fk4fk5.c kep.h
	cl $(CFLAGS) fk4fk5.c

kepler.obj: kepler.c kep.h
	cl $(CFLAGS) kepler.c

kfiles.obj: kfiles.c kep.h
	cl $(CFLAGS) kfiles.c

lightt.obj: lightt.c kep.h
	cl $(CFLAGS) lightt.c

lonlat.obj: lonlat.c kep.h
	cl $(CFLAGS) lonlat.c

nutate.obj: nutate.c kep.h
	cl $(CFLAGS) nutate.c

precess.obj: precess.c kep.h
	cl $(CFLAGS) precess.c

refrac.obj: refrac.c kep.h
	cl $(CFLAGS) refrac.c

rplanet.obj: rplanet.c kep.h
	cl $(CFLAGS) rplanet.c

rstar.obj: rstar.c kep.h
	cl $(CFLAGS) rstar.c

sidrlt.obj: sidrlt.c kep.h
	cl $(CFLAGS) sidrlt.c

sun.obj: sun.c kep.h
	cl $(CFLAGS) sun.c

trnsit.obj: trnsit.c kep.h
	cl $(CFLAGS) trnsit.c

vearth.obj: vearth.c kep.h
	cl $(CFLAGS) vearth.c

zatan2.obj: zatan2.c kep.h
	cl $(CFLAGS) zatan2.c

mer404.obj: mer404.c kep.h plantbl.h
	cl $(CFLAGS) mer404.c

ven404.obj: ven404.c kep.h plantbl.h
	cl $(CFLAGS) ven404.c
 
ear404.obj: ear404.c kep.h plantbl.h
	cl $(CFLAGS) ear404.c

arcdot.obj: arcdot.c kep.h
	cl $(CFLAGS) arcdot.c

domoon.obj: domoon.c kep.h
	cl $(CFLAGS) domoon.c

mar404.obj: mar404.c kep.h plantbl.h
	cl $(CFLAGS) mar404.c

jup404.obj: jup404.c kep.h plantbl.h
	cl $(CFLAGS) jup404.c

sat404.obj: sat404.c plantbl.h
	cl $(CFLAGS) sat404.c

ura404.obj: ura404.c plantbl.h
	cl $(CFLAGS) ura404.c

nep404.obj: nep404.c kep.h plantbl.h
	cl $(CFLAGS) nep404.c

plu404.obj: plu404.c kep.h plantbl.h
	cl $(CFLAGS) plu404.c

mlr404.obj: mlr404.c kep.h plantbl.h
	cl $(CFLAGS) mlr404.c

mlat404.obj: mlat404.c kep.h plantbl.h
	cl $(CFLAGS) mlat404.c

gplan.obj: gplan.c kep.h plantbl.h
	cl $(CFLAGS) gplan.c

aa.exe: aa.obj altaz.obj angles.obj annuab.obj constel.obj deflec.obj \
deltat.obj diurab.obj diurpx.obj dms.obj epsiln.obj fk4fk5.obj \
kepler.obj kfiles.obj lightt.obj lonlat.obj nutate.obj precess.obj \
refrac.obj rplanet.obj rstar.obj sidrlt.obj sun.obj trnsit.obj \
arcdot.obj vearth.obj zatan2.obj domoon.obj \
mer404.obj ven404.obj ear404.obj mar404.obj \
jup404.obj sat404.obj ura404.obj nep404.obj plu404.obj gplan.obj \
mlr404.obj mlat404.obj
	link @aa.rsp
