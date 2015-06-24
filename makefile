CC= emcc

# for node
NODE_CFLAGS=  -O2 --closure 1 --emit-symbol-map

# for the browser
BROWSER_CFLAGS = -O2 --closure 1 --emit-symbol-map -s MODULARIZE=1 -s EXPORTED_FUNCTIONS="['_km_to_planet']"

OBJS = altaz.o angles.o annuab.o constel.o deflec.o deltat.o diurab.o \
diurpx.o dms.o epsiln.o fk4fk5.o kepler.o kfiles.o lightt.o lonlat.o \
nutate.o precess.o refrac.o rplanet.o rstar.o sidrlt.o sun.o domoon.o \
trnsit.o vearth.o zatan2.o \
gplan.o mer404.o ven404.o ear404.o mar404.o jup404.o \
sat404.o ura404.o nep404.o plu404.o mlr404.o mlat404.o

PLANET_OBJS = planet.o planet_distance.o distance.o consts.o \
mer404.o ven404.o ear404.o mar404.o jup404.o \
sat404.o ura404.o nep404.o plu404.o mlr404.o mlat404.o \
lightt.o kepler.o angles.o vearth.o dms.o gplan.o zatan2.o deltat.o \
precess.o epsiln.o planet_kinit.o

INCS = kep.h plantbl.h

all: aa conjunct moonrise planet

clean:
	rm -f *.o aa aa.html aa.html.mem aa.js aa.js.mem aa.js.symbols \
	conjunct conjunct.js conjunct.js.mem conjunct.js.symbols conjunct.html conjunct.html.mem \
	moonrise moonrise.js moonrise.js.mem moonrise.js.symbols moonrise.html moonrise.html.mem \
	planet planet.js planet.js planet.js.mem planet.js.symbols \
	planet-node.js planet-node.js.mem planet-node.js.symbols \
	planet.html planet.html.mem

consts.o: consts.c $(INC)

aa: aa.o $(OBJS)
	$(CC) $(NODE_CFLAGS) -o aa.js aa.o $(OBJS) -lm

aa.o: aa.c $(INCS)

conjunct: conjunct.o $(OBJS) $(INCS)
	$(CC) $(NODE_CFLAGS) -o conjunct.js conjunct.o $(OBJS) -lm

conjunct.o: conjunct.c $(INCS)

moonrise: moonrise.o $(OBJS) $(INCS)
	$(CC) $(NODE_CFLAGS) -o moonrise.js moonrise.o $(OBJS) -lm

moonrise.o: moonrise.c $(INCS)

planet: $(PLANET_OBJS) $(INCS)
	$(CC) $(NODE_CFLAGS) -o planet-node.js $(PLANET_OBJS) -lm
	$(CC) $(BROWSER_CFLAGS) -o planet.js $(PLANET_OBJS) -lm

planet.o: planet.c $(INCS)

altaz.o: altaz.c $(INCS)

angles.o: angles.c $(INCS)

annuab.o: annuab.c $(INCS)

constel.o: constel.c $(INCS)

deflec.o: deflec.c $(INCS)

deltat.o: deltat.c $(INCS)

distance.o: distance.c $(INCS)

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

planet_distance.o: planet_distance.c $(INCS)

planet_kinit.o: planet_kinit.c $(INCS)

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
