/*  Find angle between two vectors.

    For two vectors p, q, the angle A between them is given by

       p.q / (|p| |q|)  = cos A  .

    If the angle is small, an expression in sin A would be preferred.
    Set r = q - p.  Then

       p.q = p.p + p.r ,

       |p|^2 = p.p ,

       |q|^2 = p.p + 2 p.r + r.r ,

                    p.p^2 + 2 p.p p.r + p.r^2
       cos^2 A  =  ----------------------------
                      p.p (p.p + 2 p.r + r.r)

                    p.p + 2 p.r + p.r^2 / p.p
                =  --------------------------- ,
                       p.p + 2 p.r + r.r

       sin^2 A  =  1 - cos^2 A

                     r.r - p.r^2 / p.p
                =  --------------------
                    p.p + 2 p.r + r.r

                =   (r.r - p.r^2 / p.p) / q.q  .

*/

#if __STDC__
double sqrt (double);
double acos (double);
double asin (double);
double atan (double);
#else
double sqrt(), acos(), asin(), atan();
#endif

/* extern double PIO2, PI; */

double arcdot(p,q)
double p[], q[];
{
double pp, pr, qq, rr, rt, pt, qt, pq;
int i;

pq = 0.0;
qq = 0.0;
pp = 0.0;
pr = 0.0;
rr = 0.0;
for (i=0; i<3; i++)
  {
    pt = p[i];
    qt = q[i];
    pq += pt * qt;
    qq += qt * qt;
    pp += pt * pt;
    rt = qt - pt;
    pr += pt * rt;
    rr += rt * rt;
  }
if (rr == 0.0 || pp == 0.0 || qq == 0.0)
  return 0.0;
/*
if (pq == 0.0)
  return PIO2;
if (pr == 0.0)
  return (atan(sqrt(rr/pp)));
*/
rt = (rr - (pr * pr) / pp) / qq;
#if DEBUG
pt = pq / sqrt(pp*qq);
pt = acos(pt);
qt = sqrt(rt);
qt = asin(qt);
if( pq < 0.0)
  qt = PI - qt;
printf("%.16e %.16e\n", qt, pt);
#endif
if (rt <= 0.25 && pq > 0.0)
  {
    rt = sqrt(rt);
    qt = asin(rt);
  }
else
  {
    pt = pq / sqrt(pp*qq);
    qt = acos(pt);
  }
return qt;
}

#if DEBUG
main()
{
double p[3], q[3], a;
int i, j;
for(j=0; j<10; j++)
  {
    for (i=0; i<3; i++)
      {
	drand(&p[i]);
	drand(&q[i]);
	p[i] -= 1.5;
	q[i] -= 1.5;
      }
    a = arcdot(p,q);
  }
}
#endif
