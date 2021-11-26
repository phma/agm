/******************************************************/
/*                                                    */
/* angle.cpp - angles as binary fractions of rotation */
/*                                                    */
/******************************************************/

#include "angle.h"
using namespace std;

double sin(int angle)
{
  return sinl(angle*M_PIl/1073741824.);
}

double cos(int angle)
{
  return cosl(angle*M_PIl/1073741824.);
}

complex<double> cossin(int angle)
{
  return complex<double>(cos(angle),sin(angle));
}
