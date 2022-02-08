/******************************************************/
/*                                                    */
/* khe.cpp - compute khe function                     */
/*                                                    */
/******************************************************/
/* Copyright 2021,2022 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */

/* The khe function խ(z) is defined on the left half-plane as follows:
 * խ(z) is asymptotic to 4*exp(z)+1 as Re(z)-> -∞.
 * That is, lim[Re(z)-> -∞] (խ(z)-1)/exp(z) =4.
 * խ(z)*խ(z+πi)=խ(2z+πi)²
 * խ(z)+խ(z+πi)=2խ(2z)
 * խ(z) is continuous on the left half-plane.
 *
 * The function խ(x+yi), as a function of y where x is fixed, is called a loop
 * (it's periodic). To compute a loop, this program starts with 36 points on
 * a circle of radius 65 ulps centered on p in [?,2), where the circle centered
 * at p=1 is the loop for x=-ln(2^54/65) (-33.25556048034141) within machine
 * precision, assuming 8-byte floats. It then applies the functional equations,
 * doubling the number of points each time.
 */

#include <iostream>
#include <cfloat>
#include <cassert>
#include "agm.h"
#include "khe.h"
#include "pairwisesum.h"
using namespace std;
#define ULPRAD 85
/* This can be 65 (5*13), 85 (5*17), or 221 (13*17).
 * These numbers have 36 integral points on their circles.
 * The normal setting is 65.
 */

#if ULPRAD==65
const signed char c65[]=
{
  65,63,60,56,52,39,33,25,16,
  0,-16,-25,-33,-39,-52,-56,-60,-63,
  -65,-63,-60,-56,-52,-39,-33,-25,-16,
  0,16,25,33,39,52,56,60,63
};

const double a65[]=
{
  atan2(0,65),atan2(16,63),atan2(25,60),atan2(33,56),atan2(39,52),
  atan2(52,39),atan2(56,33),atan2(60,25),atan2(63,16),atan2(65,0)
};
#endif

#if ULPRAD==85
const signed char c65[]=
{
  85,84,77,75,68,51,40,36,13,
  0,-13,-36,-40,-51,-68,-75,-77,-84,
  -85,-84,-77,-75,-68,-51,-40,-36,-13,
  0,13,36,40,51,68,75,77,84
};

const double a65[]=
{
  atan2(0,85),atan2(13,84),atan2(36,77),atan2(40,75),atan2(51,68),
  atan2(68,51),atan2(75,40),atan2(77,36),atan2(84,13),atan2(85,0)
};
#endif

map<double,vector<vector<complex<double> > > > loopCache;
/* The key is the circle center used to make the 65-ulp loop, which loops
 * must be divided by when fetching them from cache. The value is a sequence
 * of loops, the 0th being the 65-ulp loop with 36 points, and each successive
 * loop having twice as many points.
 */

vector<complex<double> > tinyCircle(complex<double> center)
/* Returns a circle with radius 65 ulps. center must be between 2 and -2,
 * exclusive, or the result will be inaccurate.
 */
{
  vector<complex<double> > ret;
  int i;
  for (i=0;i<36;i++)
    ret.push_back(center+complex<double>(c65[i],c65[(i+27)%36])*DBL_EPSILON);
  return ret;
}

double circleCenter(double x)
/* Returns the center to pass to tinyCircle to compute խ(x+iy).
 * The result should be at most 2-65ulp, and if called with 2x, the result
 * should be greater than 2-65ulp.
 */
{
  return exp(-x)*ULPRAD/4*DBL_EPSILON;
}

vector<complex<double> > agmExpand(vector<complex<double> > loop)
{
  vector<complex<double> > ret;
  int i,sz=loop.size();
  array<complex<double>,2> agpair;
  assert(sz%2==0);
  ret.resize(sz*2);
  for (i=0;i<sz;i++)
  {
    agpair=invAgm1(loop[i],loop[(i+sz/2)%sz]);
    //if (i && abs(agpair[1]-ret[i-1])==abs(agpair[0]-ret[i-1]))
      //cout<<"=\n"; // This starts happening when the loop is expanded 11 times.
    if (i && abs(agpair[1]-ret[i-1])+abs(agpair[0]-ret[i+sz-1])<
	     abs(agpair[0]-ret[i-1])+abs(agpair[1]-ret[i+sz-1]))
      swap(agpair[0],agpair[1]);
    ret[i]=agpair[0];
    ret[i+sz]=agpair[1];
  }
  return ret;
}

KheCachedLoop _getLoop(double x)
{
  double center=0,tryCenter=0;
  KheCachedLoop ret;
  ret.center=0;
  ret.loop=nullptr;
  int nExpand=-1,i;
  for (i=0;x<0 && tryCenter<2-ULPRAD*DBL_EPSILON;i++)
  {
    tryCenter=circleCenter(ldexp(x,i));
    if (tryCenter<2-ULPRAD*DBL_EPSILON)
    {
      nExpand=i;
      center=tryCenter;
    }
  }
  if (center)
  {
    if (!loopCache.count(center))
      loopCache[center].push_back(tinyCircle(center));
    while (loopCache[center].size()-1<nExpand)
      loopCache[center].push_back(agmExpand(loopCache[center].back()));
    ret.loop=&loopCache[center][nExpand];
    ret.center=center;
  }
  return ret;
}

vector<complex<double> > getLoop(double x)
/* Returns the loop with real part equal to x, which must be negative.
 * If x results in a circle center greater than 2, returns an empty vector;
 * խ(z) is then within an ulp of 4*exp(z)+1. If x>-1/60., it may return
 * the numbers in the loop in the wrong order. If x>=0, returns empty.
 */
{
  KheCachedLoop cloop;
  vector<complex<double> > ret;
  int i;
  cloop=_getLoop(x);
  if (cloop.center)
  {
    ret=*cloop.loop;
    for (i=0;i<ret.size();i++)
      ret[i]/=cloop.center;
  }
  return ret;
}

KheInterp getInterp(complex<double> z)
/* Returns 12 numbers from the loop, of which points[1] and points[10] come
 * from the quadrants of the original circle of radius 65 ulps, and the
 * amount by which z.imag() is along the interval from points[1] to points[10].
 */
{
  KheCachedLoop cloop;
  KheInterp ret;
  double y;
  int pointStart,i,sz;
  cloop=_getLoop(z.real());
  if (cloop.center)
  {
    y=z.imag()-(2*M_PI)*rint(z.imag()/(2*M_PI));
    sz=cloop.loop->size();
    pointStart=lrint(y*(sz/18)/M_PI)*9;
    ret.along=y*(sz/36)-(pointStart/9)*(M_PI/2);
    if (ret.along<0)
    {
      ret.along+=M_PI/2;
      pointStart-=9;
    }
    if (pointStart<2)
      pointStart+=sz;
    for (i=0;i<12;i++)
      ret.points[i]=(*cloop.loop)[(pointStart+i-1)%sz]/cloop.center;
  }
  else
    ret.along=NAN;
  return ret;
}

double avgRadius(vector<complex<double> > loop)
{
  int i,sz=loop.size();
  vector<double> diams;
  diams.resize(sz/2);
  for (i=0;i<sz/2;i++)
    diams[i]=abs(loop[i]-loop[i+sz/2]);
  return pairwisesum(diams)/sz;
}

vector<double> vecLog(vector<complex<double> > loop)
{
  vector<double> ret;
  int i;
  for (i=0;i<loop.size();i++)
    ret.push_back(log(abs(loop[i])));
  return ret;
}

vector<double> vecArg(vector<complex<double> > loop)
{
  vector<double> ret;
  double a,lasta=0;
  int i;
  for (i=0;i<loop.size();i++)
  {
    a=arg(loop[i]);
    a+=2*M_PI*rint((lasta-a)/2/M_PI);
    ret.push_back(a);
    lasta=a;
  }
  return ret;
}

double xt(int n)
{
  return a65[9]*(n/9)+a65[n%9];
}

complex<double> khe(complex<double> z)
/* Uses linear interpolation. TODO use cubic interpolation.
 * Computes the khe function of z. If z is too close to the imaginary axis,
 * may give wrong answers.
 */
{
  KheInterp interp=getInterp(z);
  complex<double> ret;
  int i,n;
  double subalong,interval[3],p,q;
  complex<double> off,pnt[4],ctrl[2],slp[2];
  if (z.real()>=0)
    ret=complex<double>(NAN,NAN);
  else if (isnan(interp.along))
    ret=4.*exp(z)+1.;
  else
  {
    for (i=0;i<9;i++)
      if (a65[i]<=interp.along)
      {
	n=i;
	subalong=interp.along-a65[i];
	interval[1]=a65[i+1]-a65[i];
      }
    if (abs(interp.points[n+1])<abs(interp.points[n+2]))
      off=interp.points[n+1];
    else
      off=interp.points[n+2];
    for (i=0;i<4;i++)
      pnt[i]=interp.points[n+i]-off;
    assert(subalong<interval[1]);
    if (n==0)
      interval[0]=a65[1];
    else
      interval[0]=a65[n]-a65[n-1];
    if (n==8)
      interval[2]=a65[1];
    else
      interval[2]=a65[n+2]-a65[n+1];
    slp[0]=((pnt[2]-pnt[1])*interval[0]+
	    (pnt[1]-pnt[0])/interval[0]*interval[1]*interval[1])/
	    (interval[0]+interval[1]);
    slp[1]=((pnt[2]-pnt[1])*interval[2]+
	    (pnt[3]-pnt[2])/interval[2]*interval[1]*interval[1])/
	    (interval[1]+interval[2]);
    ctrl[0]=pnt[1]+slp[0]/3.;
    ctrl[1]=pnt[2]-slp[1]/3.;
    p=subalong/interval[1];
    q=1-p;
    p=1-q;
    ret=(pnt[1]*q*q*q+3.*ctrl[0]*p*q*q+3.*ctrl[1]*p*p*q+pnt[2]*p*p*p)+off;
  }
  return ret;
}

void outMaxMag(vector<complex<double> > &loop)
/* Outputs all local maxima of the absolute value of the loop.
 * loop[0] is the global maximum.
 */
{
  int i,sz=loop.size();
  for (i=0;i<sz-1;i++)
    if (i==0 || (abs(loop[i])>abs(loop[i-1]) && abs(loop[i])>abs(loop[i+1])))
      cout<<xt(i)*1./xt(sz)<<' '<<loop[i]/loop[0]*2.<<endl;
}
