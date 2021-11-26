/******************************************************/
/*                                                    */
/* agm.cpp - compute arithmetic-geometric mean        */
/*                                                    */
/******************************************************/

#include "agm.h"
#include "angle.h"
using namespace std;

AgmRec agm1(AgmRec ag)
{
  AgmRec ret;
  complex<double> c;
  ret.a=(ag.a+ag.g)/2.;
  ret.g=sqrt(ag.a*ag.g);
  c=ret.a*cossin(ag.d);
  if (abs(c-ret.g)>abs(c+ret.g))
    ret.g=-ret.g;
  ret.d=argi(ret.g/ret.a);
  return ret;
}
