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
  return exp(-x)*65/4*DBL_EPSILON;
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
  for (i=0;x<0 && tryCenter<2-65*DBL_EPSILON;i++)
  {
    tryCenter=circleCenter(ldexp(x,i));
    if (tryCenter<2-65*DBL_EPSILON)
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
  int nExpand=-1,i;
  cloop=_getLoop(x);
  if (cloop.center)
  {
    ret=*cloop.loop;
    for (i=0;i<ret.size();i++)
      ret[i]/=cloop.center;
  }
  return ret;
}

//KheInterp getInterp(complex<double> z)

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
