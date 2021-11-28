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
int argi(std::complex<double> z);

#define SEC1 1657
#define FURMAN1 32768
#define MIN1 99421
#define DEG1 5965232
#define AT0512 0x80ae90e
// AT0512 is arctangent of 5/12, 22.619865°
#define DEG30 0x0aaaaaab
#define AT34 0x0d1bfae2
// AT34 is arctangent of 3/4, 36.8698976°
#define DEG40 0x0e38e38e
#define DEG45 0x10000000
#define DEG50 0x11c71c72
#define DEG60 0x15555555
#define DEG72 0x1999999a
#define DEG90 0x20000000
#define DEG120 0x2aaaaaab
#define DEG144 0x33333333
#define DEG150 0x35555555
#define DEG180 0x40000000
#define DEG270 0x60000000
#define DEG360 0x80000000
#define PHITURN 0xcf1bbcdd

#endif
