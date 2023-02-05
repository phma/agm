/******************************************************/
/*                                                    */
/* relprime.cpp - relatively prime numbers            */
/*                                                    */
/******************************************************/
/* Copyright 2015,2021 Pierre Abbat.
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */
#include <map>
#include <cmath>
#include "relprime.h"
#include "khe.h"

using namespace std;

map<unsigned,unsigned> relprimes;

unsigned relprime(unsigned n)
// Returns the integer closest to n/φ of those relatively prime to n.
{
  unsigned ret,twice;
  double phin;
  ret=relprimes[n];
  if (!ret)
  {
    phin=n*M_1PHI;
    ret=rint(phin);
    twice=2*ret-(ret>phin);
    while (gcd(ret,n)!=1)
      ret=twice-ret+(ret<=phin);
    relprimes[n]=ret;
  }
  return ret;
}
