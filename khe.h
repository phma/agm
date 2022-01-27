/******************************************************/
/*                                                    */
/* khe.h - compute khe function                       */
/*                                                    */
/******************************************************/
/* Copyright 2021,2022 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */

#ifndef KHE_H
#define KHE_H

#include <cmath>
#include <complex>
#include <vector>
#include <array>

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

double circleCenter(double x);
std::vector<std::complex<double> > tinyCircle(std::complex<double> center);
std::vector<std::complex<double> > agmExpand(std::vector<std::complex<double> > loop);
std::vector<double> vecLog(std::vector<std::complex<double> > loop);
std::vector<double> vecArg(std::vector<std::complex<double> > loop);
std::vector<std::complex<double> > getLoop(double x);
//KheInterp getInterp(complex<double> z);
double avgRadius(std::vector<std::complex<double> > loop);
int xt(int n);
void outMaxMag(std::vector<std::complex<double> > &loop);

#endif
