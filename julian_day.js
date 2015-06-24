function calculateJD(cdDate)
{
   var year = cdDate.getUTCFullYear()
   var month = cdDate.getUTCMonth() + 1
   if (month < 3)
   {
      year = year - 1
      month = month + 12
   }
   var day = cdDate.getUTCDate()
   var hour = cdDate.getUTCHours()
   var min = cdDate.getUTCMinutes()
   var sec = cdDate.getUTCSeconds()

   var A = parseInt(year / 100)
   var B = 2 - A + parseInt(A / 4)

   var E = hour / 24 + min / 1440 + sec / (8.64 * Math.pow(10, 4))

   var jd = parseInt(365.25 * (year + 4716)) + parseInt(30.6001 * (month + 1)) + day + B - 1524.5 + E

   // correct for dynamical time
   // 2011 : 66.7 seconds. 66.7 / (24 * 3600) = .000771991

   // jd = jd + .000771991

   return jd;
}

// http://www.physics.emory.edu/astronomy/events/mars/mars_calc_V4.js

console.log(calculateJD(new Date()));
