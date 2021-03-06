astronomical-almanac-js
=======================

This is the port to asm.js of Debian's [astronomical-almanac](https://packages.debian.org/source/sid/astronomical-almanac).

The original license is GPL v2.0 and is found in the ```debian/copyright``` file.
This package is licensed under the same licence.

The original README is the ```readme.404``` file in the main directory.

The changes made are:

* The ```makefile``` generates JavaScript targets using emcc and modified CFLAGS
* Added ```planet.js``` and ```planet-node.js``` targets. They wrap the algorithm that computes the distance between two planets. ```planet.js``` can be used in the browser using ```Module.cwrap``` and ```planet-node.js``` can be used from the command line.

All programs compile but the only one guaranteed to work is the new ```planet.js``` library and the ```planet-node.js``` standalone program.

Building
--------

    $ npm install
    $ make clean
    $ make -j  # or just make

The distance algorithm is embedded in ```planet.js```, ```planet.js.mem``` and ```planet.js.symbols``` for browser based projects.  You'll need ```planet-node.js```, ```planet-node.js.mem``` and ```planet-node.js.symbols``` for running it from the command line.

Examples
--------

In the browser

    <script type="text/javascript" src="planet.js"></script>

    var kmToPlanet = Module().cwrap("km_to_planet", "number", ["number", "number"]);
    var now = new Date().getTime() / 86400000 + 2440587.5; // Julian day
    console.log(kmToPlanet(now, 4)); // Earth to Mars

A live example: http://connettiva.eu/newton (source at [newton-says](https://github.com/pmontrasio/newton-says))

From the command line

    $ node planet-node.js $(node julian_day.js) 4  # Earth to Mars

Available planet numbers
------------------------

* 1 Mercury
* 2 Venus
* 4 Mars
* 5 Jupiter
* 6 Saturn
* 7 Uranus
* 8 Neptune
* 9 Pluto


