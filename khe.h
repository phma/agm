/******************************************************/
/*                                                    */
/* khe.h - compute khe function                       */
/*                                                    */
/******************************************************/
/* Copyright 2021-2023 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */

#ifndef KHE_H
#define KHE_H

#include <cmath>
#include <complex>
#include <vector>
#include <array>
#include <map>

unsigned gcd(unsigned a,unsigned b);

struct KheCachedLoop
{
  double center;
  std::vector<std::complex<double> > *loop;
};

struct KheInterp
{
  std::array<std::complex<double>,12> points;
  double along;
};

class KheSwapStep
{
public:
  KheSwapStep()=default;
  KheSwapStep(int n,int d,std::vector<std::complex<double> > &loop);
  void step(std::vector<std::complex<double> > &loop);
  void swap(std::vector<std::complex<double> > &loop);
  std::complex<double> lastDiff;
  double dist;
  int a,b;
  int dir;
  friend bool operator<(const KheSwapStep &a,const KheSwapStep &b);
  friend bool operator>(const KheSwapStep &a,const KheSwapStep &b);
  friend bool meet(KheSwapStep &n,KheSwapStep &s);
};

std::vector<std::complex<double> > agmExpand(std::vector<std::complex<double> > loop);
std::vector<double> vecLog(std::vector<std::complex<double> > loop);
std::vector<double> vecArg(std::vector<std::complex<double> > loop);
double avgRadius(std::vector<std::complex<double> > loop);
double xt(int n);

class Khe
{
public:
  Khe();
  Khe(int circleSize); // Must be a member of http://oeis.org/A131574
  std::vector<std::complex<double> > getLoop(double x);
  double xt(int n);
  std::complex<double> operator()(std::complex<double> z);
  void outMaxMag(std::vector<std::complex<double> > &loop);
private:
  int cirCoord[36];
  double arcTan[10];
  int radius; // Radius of circle returned by tinyCircle
  std::map<double,std::vector<std::vector<std::complex<double> > > > loopCache;
  /* The key is the circle center used to make the circular loop, which loops
  * must be divided by when fetching them from cache. The value is a sequence
  * of loops, the 0th being the circular loop with 36 points, and each
  * successive loop having twice as many points.
  */
  void init(int circleSize);
  std::vector<std::complex<double> > tinyCircle(std::complex<double> center);
  double circleCenter(double x);
  KheCachedLoop _getLoop(double x);
  KheInterp getInterp(std::complex<double> z);
};

#endif
