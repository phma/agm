/******************************************************/
/*                                                    */
/* agm.cpp - compute arithmetic-geometric mean        */
/*                                                    */
/******************************************************/

#include "agm.h"

AgmRec agm1(AgmRec ag)
{
  AgmRec ret;
  ret.a=(ag.a+ag.g)/2.;
  ret.g=sqrt(ag.a*ag.g);
  return ret;
}
