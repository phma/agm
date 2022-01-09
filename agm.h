/******************************************************/
/*                                                    */
/* agm.h - compute arithmetic-geometric mean          */
/*                                                    */
/******************************************************/
/* Copyright 2021,2022 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */

#include <cmath>
#include <complex>
#include <string>
#include <vector>
#include <array>

struct AgmRec
{
  std::complex<double> a,g;
  int32_t d;
};

struct AgmResult
{
  std::complex<double> m;
  std::string branch;
};

AgmRec agm1(AgmRec ag);
AgmResult agm(std::complex<double> a,std::complex<double> g=1,std::string branch="");
std::vector<std::complex<double> > agmLattice(std::complex<double> a,std::complex<double> g=1,unsigned depth=0,unsigned level=0,std::string branch="");
std::complex<double> pvAgm(std::complex<double> a,std::complex<double> g);
std::array<std::complex<double>,2> invAgm1(std::complex<double> a,std::complex<double> g);
