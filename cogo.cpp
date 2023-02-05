/******************************************************/
/*                                                    */
/* cogo.cpp - coordinate geometry                     */
/*                                                    */
/******************************************************/
/* Copyright 2012,2015,2017,2023 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */

#include "cogo.h"
using namespace std;

#define CMPSWAP(m,n,o) if (fabs(m)>fabs(n)) {o=m;m=n;n=o;}
#define EQOPP(m,n) if (m+n==0 && m>0) {m=-m;n=-n;}
double area3(complex<double> a,complex<double> b,complex<double> c)
{
  double surface,temp,xl,xm,xh,yl,ym,yh,area[6];
  complex<double> m;
  // Translate the points near the origin for greater precision.
  xl=a.real();xm=b.real();xh=c.real();
  yl=a.imag();ym=b.imag();yh=c.imag();
  CMPSWAP(xl,xm,surface);
  CMPSWAP(yl,ym,temp);
  CMPSWAP(xm,xh,surface);
  CMPSWAP(ym,yh,temp);
  CMPSWAP(xl,xm,surface);
  CMPSWAP(yl,ym,temp);
  m=complex<double>(xm,ym);
  a-=m;
  b-=m;
  c-=m;
  // Compute the six areas.
  area[0]=a.real()*b.imag();
  area[1]=-b.real()*a.imag();
  area[2]=b.real()*c.imag();
  area[3]=-c.real()*b.imag();
  area[4]=c.real()*a.imag();
  area[5]=-a.real()*c.imag();
  // Sort the six areas into absolute value order for numerical stability.
  CMPSWAP(area[0],area[1],surface); // sorting network
  CMPSWAP(area[2],area[3],temp);
  CMPSWAP(area[4],area[5],surface);
  CMPSWAP(area[0],area[2],temp);
  CMPSWAP(area[1],area[4],surface);
  CMPSWAP(area[3],area[5],temp);
  CMPSWAP(area[0],area[1],surface);
  CMPSWAP(area[2],area[3],temp);
  CMPSWAP(area[4],area[5],surface);
  CMPSWAP(area[1],area[2],temp);
  CMPSWAP(area[3],area[4],surface);
  CMPSWAP(area[2],area[3],temp);
  // Make signs of equal-absolute-value areas alternate.
  EQOPP(area[0],area[5]);
  EQOPP(area[0],area[3]);
  EQOPP(area[4],area[1]);
  EQOPP(area[2],area[5]);
  EQOPP(area[0],area[1]);
  EQOPP(area[2],area[1]);
  EQOPP(area[2],area[3]);
  EQOPP(area[4],area[3]);
  EQOPP(area[4],area[5]);
  surface=(((((area[0]+area[1])+area[2])+area[3])+area[4])+area[5])/2;
  return surface;
}

double pldist(complex<double> a,complex<double> b,complex<double> c)
/* Signed distance from a to the line bc. */
{
  return area3(a,b,c)/abs(b-c)*2;
}
