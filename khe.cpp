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

const char a65[]={0,7,11,15,18,26,29,33,37,44};

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
  return exp(x)/65*4/DBL_EPSILON;
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

int xt(int n)
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
