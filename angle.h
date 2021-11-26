/******************************************************/
/*                                                    */
/* angle.h - angles as binary fractions of rotation   */
/*                                                    */
/******************************************************/

// Angles are represented as integers, with 2147483648 representing 360°.
// This file overloads the functions sin and cos.

#ifndef ANGLE_H
#define ANGLE_H
#include <cmath>
#include <complex>

double sin(int angle);
double cos(int angle);
std::complex<double> cossin(int angle);
#endif
