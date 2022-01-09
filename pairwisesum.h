/******************************************************/
/*                                                    */
/* pairwisesum.h - add many numbers                   */
/*                                                    */
/******************************************************/
/* Copyright 2015,2016,2018,2021 Pierre Abbat.
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */
#ifndef MANYSUM_H
#define MANYSUM_H
#include <map>
#include <vector>
#include <cmath>
/* Adds together many numbers (like millions) accurately.
 * pairwisesum takes an array or vector with the numbers already computed.
 */

double pairwisesum(double *a,unsigned n);
double pairwisesum(std::vector<double> &a);
long double pairwisesum(long double *a,unsigned n);
long double pairwisesum(std::vector<long double> &a);

#endif
