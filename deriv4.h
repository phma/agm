/******************************************************/
/*                                                    */
/* deriv4.h - derivatives from four points            */
/*                                                    */
/******************************************************/
/* Copyright 2014,2022 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */

#include <complex>
#include <array>

std::complex<double> deriv0(std::array<std::complex<double>,4> xsect);
std::complex<double> deriv1(std::array<std::complex<double>,4> xsect);
std::complex<double> deriv2(std::array<std::complex<double>,4> xsect);
std::complex<double> deriv3(std::array<std::complex<double>,4> xsect);
