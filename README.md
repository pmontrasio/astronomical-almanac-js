astronomical-almanac-js
=======================

This is the port to asm.js of Debian's [astronomical-almanac](https://packages.debian.org/source/sid/astronomical-almanac).

The original license is GPL v2.0 and is found in the debian/ directory.
This package is licensed under the same licence.

The original README is the readme.404 file in the main directory.

The changes made are:

* The makefile generates JavaScript targets using emcc and modified CFLAGS
* Added planet.js and planet-node.js targets. They wrap the algorithm that computes the distance between two planets. planet.js can be used in the browser using Module.cwrap and planet-node.js can be used from the command line.

Examples
--------

In the browser

    <script type="text/javascript" src="planet.js"></script>

    var kmToPlanet = Module().cwrap("km_to_planet", "number", ["number", "number"]);
    console.log(kmToPlanet(now, 4)); // Earth to Mars

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


