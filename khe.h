/******************************************************/
/*                                                    */
/* khe.h - compute khe function                       */
/*                                                    */
/******************************************************/

#ifndef KHE_H
#define KHE_H

#include <cmath>
#include <complex>
#include <vector>
#include <array>

std::vector<std::complex<double> > tinyCircle(std::complex<double> center);
std::vector<std::complex<double> > agmExpand(std::vector<std::complex<double> > loop);
std::vector<double> vecLog(std::vector<std::complex<double> > loop);
std::vector<double> vecArg(std::vector<std::complex<double> > loop);
double avgRadius(std::vector<std::complex<double> > loop);
int xt(int n);
void outMaxMag(std::vector<std::complex<double> > &loop);

#endif
