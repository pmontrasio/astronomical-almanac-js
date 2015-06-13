var Planet = require("./planet.js");
var printf = require('printf');

var km_to_planet = Planet.cwrap('km_to_planet', 'number', ['number', 'number']);

const mercury = 1;
const venus = 2;
const mars = 4;
const jupiter = 5;
const saturn = 6;
const neptune = 7;
const uranus = 8;
const pluto = 9;

var planets = [{name: "Mercury", id: mercury}, {name: "Venus", id: venus}, {name: "Mars", id: mars},
               {name: "Jupiter", id: jupiter}, {name: "Saturn", id: saturn}, {name: "Neptune", id: neptune},
               {name: "Uranus", id: uranus}, {name: "Pluto", id: pluto}];

var julian_day = new Date().getTime()/86400000 + 2440587.5;

for (var planet of planets.values()) {
  console.log(printf("%7s: %14.3f km", planet.name, km_to_planet(julian_day, planet.id)));
}
