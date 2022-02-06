/******************************************************/
/*                                                    */
/* deriv4.cpp - derivatives from four points          */
/*                                                    */
/******************************************************/
/* Copyright 2014,2022 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */

#include "deriv4.h"
using namespace std;

// Computes derivatives of f(x) from f(-1.5), f(-0.5), f(0.5), and f(1.5).

complex<double> deriv0(array<complex<double>,4> xsect)
{
  return (-xsect[3]+9.*xsect[2]+9.*xsect[1]-xsect[0])/16.;
}

complex<double> deriv1(array<complex<double>,4> xsect)
{
  return (xsect[0]-27.*xsect[1]+27.*xsect[2]-xsect[3])/24.;
}

complex<double> deriv2(array<complex<double>,4> xsect)
{
  return (xsect[3]-xsect[2]-xsect[1]+xsect[0])/2.;
}

complex<double> deriv3(array<complex<double>,4> xsect)
{
  return xsect[3]-3.*xsect[2]+3.*xsect[1]-xsect[0];
}
