/******************************************************/
/*                                                    */
/* khe.cpp - compute khe function                     */
/*                                                    */
/******************************************************/
/* Copyright 2021-2023 Pierre Abbat
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
 *
 * The computation is accurate within 1.5e-8 times the maximum absolute value
 * of the loop as long as x<-0.133, and probably as long as x<-0.0664. At the
 * -0.0663 seam, something goes wrong. Between -33 (-0.016113) and -32/2048
 * (-0.015625), the loop goes out of order.
 */

#include <iostream>
#include <cfloat>
#include <cassert>
#include "agm.h"
#include "khe.h"
#include "pairwisesum.h"
using namespace std;

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

#if ULPRAD==221
const short c65[]=
{
  221,220,204,195,171,140,104,85,21,
  0,-21,-85,-104,-140,-171,-195,-204,-220,
  -221,-220,-204,-195,-171,-140,-104,-85,-21,
  0,21,85,104,140,171,195,204,220
};

const double a65[]=
{
  atan2(0,221),atan2(21,220),atan2(85,204),atan2(104,195),atan2(140,171),
  atan2(171,140),atan2(195,104),atan2(204,85),atan2(220,21),atan2(221,0)
};
#endif

int mostInnings=0;

map<double,vector<vector<complex<double> > > > loopCache;
/* The key is the circle center used to make the 65-ulp loop, which loops
 * must be divided by when fetching them from cache. The value is a sequence
 * of loops, the 0th being the 65-ulp loop with 36 points, and each successive
 * loop having twice as many points.
 */

unsigned gcd(unsigned a,unsigned b)
{
  while (a&&b)
  {
    if (a>b)
    {
      b^=a;
      a^=b;
      b^=a;
    }
    b%=a;
  }
  return a+b;
}

KheSwapStep::KheSwapStep(int n,int d,vector<complex<double> > &loop)
/* loop[a] and loop[b] are both real, in which case loop[a] should be near
 * loop[0]/k where k is in [1,-3,5,-7,...], or loop[a] and loop[b] are
 * both imaginary, in which case loop[a] should have positive imaginary part.
 */
{
  double realness,imagness;
  a=n;
  b=(n+loop.size()/2)%loop.size();
  dir=(d<0)?-1:1;
  realness=fabs(real(loop[a]))+fabs(real(loop[b]));
  imagness=fabs(imag(loop[a]))+fabs(imag(loop[b]));
  assert(realness>imagness);
  if ((realness>imagness && abs(loop[b])>abs(loop[a])) ||
      (imagness>realness && imag(loop[b])>imag(loop[a])))
    ::swap(loop[a],loop[b]);
  dist=abs(loop[a]-loop[b]);
  lastDiff=loop[a]-loop[b];
}

void KheSwapStep::step(vector<complex<double> > &loop)
{
  a+=dir;
  b+=dir;
  if (a<0)
    a+=loop.size();
  if (a>=loop.size())
    a-=loop.size();
  if (b<0)
    b+=loop.size();
  if (b>=loop.size())
    b-=loop.size();
  dist=abs(loop[a]-loop[b]);
}

void KheSwapStep::swap(vector<complex<double> > &loop)
{
  if (real((loop[a]-loop[b])/lastDiff)<0)
    ::swap(loop[a],loop[b]);
  if (loop[a]!=loop[b])
    lastDiff=loop[a]-loop[b];
}

bool operator<(const KheSwapStep &a,const KheSwapStep &b)
{
  return a.dist<b.dist;
}

bool operator>(const KheSwapStep &a,const KheSwapStep &b)
{
  return a.dist>b.dist;
}

bool meet(KheSwapStep &n,KheSwapStep &s)
{
  return n.dir!=s.dir && (n.a==s.a || n.a==s.b);
}

void findMeeters(vector<KheSwapStep *> &swapStep)
/* Arrange the steppers in order of a or b, whichever is less, except that
 * the stepper with a=0 going backward is placed at the end. This puts
 * steppers that will meet in pairs. A quadratic sort is fast enough, as
 * the loop is much bigger than the set of steppers;
 */
{
  int i,j,a0,a1,sz=swapStep.size();
  for (i=0;i<sz-1;i++)
    for (j=0;j<sz-1;j++)
    {
      a0=swapStep[j]->a;
      if (a0>swapStep[j]->b || (a0==0 && swapStep[j]->dir<0))
	a0=swapStep[j]->b;
      a1=swapStep[j+1]->a;
      if (a1>swapStep[j+1]->b || (a1==0 && swapStep[j+1]->dir<0))
	a1=swapStep[j+1]->b;
      if (a0>a1 || (a0==a1 && swapStep[j]->dir>swapStep[j+1]->dir))
	swap(swapStep[j],swapStep[j+1]);
    }
  for (i=0;i<sz;i++)
    swapStep[i]->partner=swapStep[i^1];
}

void sortSteppers(vector<KheSwapStep *> &swapStep)
/* Arrange the steppers in order of increasing distance.
 */
{
  int i,j,sz=swapStep.size();
  for (i=0;i<sz-1;i++)
    for (j=0;j<sz-1;j++)
    {
      if (*swapStep[j]>*swapStep[j+1])
	swap(swapStep[j],swapStep[j+1]);
    }
}

vector<complex<double> > agmExpand(vector<complex<double>> loop,double center)
/* Starting angles and where they end up:
 * 0°		(1,0)
 * 15°		(0,1/12)
 * 18°		(0,1/10)
 * 22.5°	(0,1/8)
 * 30°		(0,1/6)
 * 40°		(1/9,0)
 * 45°		(0,1/4)
 * 54°		(0,-1/10)
 * 60°		(0,0)
 * 72°		(1/5,0)
 * 75°		(0,1/12)
 * 80°		(1/9,0)
 * 90°		(0,1/2)
 * 105°		(0,-1/12)
 * 120°		(-1/3,0)
 * 126°		(0,-1/10)
 * 135°		(0,-1/4)
 * 144° 	(1/5,0)
 * 150°		(0,1/6)
 * 157.5°	(0,-1/8)
 * 160°		(1/9,0)
 * 162°		(0,1/10)
 * 165°		(0,-1/12)
 * 180°		(0,0)
 *
 * General rule for 2πa/b, where a and b are relatively prime:
 * b=4n+1	(1/b,0)
 * b=4n-1	(-1/b,0)
 * b=4n+2	(0,0)
 * b=4n,a=4n+1	(0,2/b)
 * b=4n,a=4n-1	(0,-2/b)
 *
 * Expansion, where the first is the AM, the second is the GM, and the
 * third and fourth are the results of expansion:
 * 0° (1,0), 180° (0,0) -> 0° (2,0), 180° (0,0)
 * 60° (0,0), 240° (-1/3,0) -> 30° (0,1/3), 210° (0,-1/3)
 * 90° (0,1/2), 270° (0,-1/2) -> 45° (0,1/2), 225° (0,1/2)
 * 120° (-1/3,0), 300° (0,0) -> 60° (0,0), 240° (-2/3,0)
 * 180° (0,0), 0° (1,0) -> 90° (0,1), 270° (0,-1)
 * 240° (-1/3,0), 60° (0,0) -> 120° (-2/3,0), 300° (0,0)
 * 270° (0,-1/2), 90° (0,1/2) -> 135° (0,-1/2), 315° (0,-1/2)
 * 300° (0,0), 120° (-1/3,0) -> 150° (0,1/3), 330° (0,-1/3)
 * Swapping should start at 0°/180° and 90°/270° and proceed in both directions,
 * ending at 45°/225° and 135°/315°, where the numbers being swapped
 * end up equal.
 */
{
  vector<complex<double> > ret;
  vector<int> fractInx;
  int i,j,sz=loop.size(),innings=0;
  array<complex<double>,2> agpair;
  vector<KheSwapStep *> swapStep;
  assert(sz%2==0);
  ret.resize(sz*2);
  for (i=0;i<sz;i++)
  {
    if (abs(loop[i])<center && abs(loop[i-1])>=center)
      innings++;
    agpair=invAgm1(loop[i],loop[(i+sz/2)%sz]);
    ret[i]=agpair[0];
    ret[i+sz]=agpair[1];
  }
  if (innings>mostInnings)
  {
    //cout<<innings<<" innings, "<<sz<<" loop size\n";
    mostInnings=innings;
  }
  /* The number of innings is 1, 3, 7, 13, 19, 29, ... (A099957).
   * The real arcs appear in the order 1/1, -1/3, 1/5, ..., each traced
   * as many times as the totient of the denominator. For each real arc,
   * there is an inning opposite to it. So place steppers on each tracing
   * of a real arc, and place as many pairs of steppers as there are innings.
   */
  for (i=1;;i+=2)
  {
    fractInx.clear();
    for (j=0;j<i;j++)
      if (gcd(i,j)==1)
	fractInx.push_back(lrint(2.*j*sz/i));
    if (swapStep.size()/2+fractInx.size()<=innings)
      for (j=0;j<fractInx.size();j++)
      {
	swapStep.push_back(new KheSwapStep(fractInx[j],1,ret));
	swapStep.push_back(new KheSwapStep(fractInx[j],-1,ret));
      }
    else
      break;
  }
  findMeeters(swapStep);
  sortSteppers(swapStep);
  while (swapStep.size())
  {
    swapStep.back()->step(ret);
    swapStep.back()->swap(ret);
    sz=swapStep.size();
    if (meet(*swapStep[sz-1],*swapStep[sz-1]->partner))
    {
      for (i=0;i<sz-2;i++)
	if (swapStep[i]==swapStep[sz-1]->partner)
	  swap(swapStep[i],swapStep[i+1]);
      delete swapStep[sz-1];
      delete swapStep[sz-2];
      swapStep.resize(sz-2);
    }
    sz=swapStep.size();
    i=sz-1;
    while (i>0 && *swapStep[i]<*swapStep[i-1])
    {
      swap(swapStep[i],swapStep[i-1]);
      i--;
    }
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
/* last2a is an attempt to follow the loop near -1/128, in which some steps
 * are bigger than 180°. It didn't work.
 */
{
  vector<double> ret;
  double a,lasta=0,last2a=0;
  int i;
  for (i=0;i<loop.size();i++)
  {
    a=arg(loop[i]);
    a+=2*M_PI*rint((lasta-a)/2/M_PI);
    ret.push_back(a);
    last2a=lasta;
    lasta=a;
  }
  return ret;
}

void Khe::init(int circleSize)
{
  int64_t i=abs(circleSize),j=0,n=0,sq=i*i;
  radius=i;
  while (i>j)
    if (i*i+j*j==sq)
    {
      arcTan[n]=atan2(j,i);
      arcTan[9-n]=atan2(i,j);
      cirCoord[n]=i--;
      cirCoord[27+n]=cirCoord[9-n]=j++;
      cirCoord[18+n]=cirCoord[18-n]=-cirCoord[n];
      cirCoord[9+n]=cirCoord[27-n]=-cirCoord[9-n];
      if (n)
	cirCoord[36-n]=cirCoord[n];
      n++;
    }
    else if (i*i+j*j>radius*radius)
      i--;
    else
      j++;
  assert(n==5);
  /* circleSize must be such that there are 36 integral points on a circle 
   * of radius circleSize. Such numbers are 65, 85, 145, 185, 205, 221, etc.
   * See http://oeis.org/A131574 .
   */
}

Khe::Khe()
{
  init(65);
}

Khe::Khe(int circleSize)
{
  init(circleSize);
}

vector<complex<double> > Khe::tinyCircle(complex<double> center)
/* Returns a circle of size radius ulps. center must be between 2 and -2,
 * exclusive, or the result will be inaccurate.
 */
{
  vector<complex<double> > ret;
  int i;
  for (i=0;i<36;i++)
    ret.push_back(center+complex<double>(cirCoord[i],cirCoord[(i+27)%36])*DBL_EPSILON);
  return ret;
}

double Khe::circleCenter(double x)
/* Returns the center to pass to tinyCircle to compute խ(x+iy).
 * The result should be at most 2-radius*ulp, and if called with 2x, the result
 * should be greater than 2-radius*ulp.
 */
{
  return exp(-x)*radius/4*DBL_EPSILON;
}

KheCachedLoop Khe::_getLoop(double x)
{
  double center=0,tryCenter=0;
  KheCachedLoop ret;
  ret.center=0;
  ret.loop=nullptr;
  int nExpand=-1,i;
  for (i=0;x<0 && tryCenter<2-radius*DBL_EPSILON;i++)
  {
    tryCenter=circleCenter(ldexp(x,i));
    if (tryCenter<2-radius*DBL_EPSILON)
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
      loopCache[center].push_back(agmExpand(loopCache[center].back(),center));
    ret.loop=&loopCache[center][nExpand];
    ret.center=center;
  }
  return ret;
}

vector<complex<double> > Khe::getLoop(double x)
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

KheInterp Khe::getInterp(complex<double> z)
/* Returns 12 numbers from the loop, of which points[1] and points[10] come
 * from the quadrants of the original circle of size radius ulps, and the
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

double Khe::xt(int n)
{
  return arcTan[9]*(n/9)+arcTan[n%9];
}

complex<double> Khe::operator()(complex<double> z)
/* uses cubic interpolation.
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
      if (arcTan[i]<=interp.along)
      {
	n=i;
	subalong=interp.along-arcTan[i];
	interval[1]=arcTan[i+1]-arcTan[i];
      }
    if (abs(interp.points[n+1])<abs(interp.points[n+2]))
      off=interp.points[n+1];
    else
      off=interp.points[n+2];
    for (i=0;i<4;i++)
      pnt[i]=interp.points[n+i]-off;
    assert(subalong<interval[1]);
    if (n==0)
      interval[0]=arcTan[1];
    else
      interval[0]=arcTan[n]-arcTan[n-1];
    if (n==8)
      interval[2]=arcTan[1];
    else
      interval[2]=arcTan[n+2]-arcTan[n+1];
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

void Khe::outMaxMag(vector<complex<double> > &loop)
/* Outputs all local maxima of the absolute value of the loop.
 * loop[0] is the global maximum.
 */
{
  int i,sz=loop.size();
  for (i=0;i<sz-1;i++)
    if (i==0 || (abs(loop[i])>abs(loop[i-1]) && abs(loop[i])>abs(loop[i+1])))
      cout<<xt(i)*1./xt(sz)<<' '<<loop[i]/loop[0]*2.<<endl;
}
