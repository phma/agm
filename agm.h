/******************************************************/
/*                                                    */
/* agm.h - compute arithmetic-geometric mean          */
/*                                                    */
/******************************************************/

#include <cmath>
#include <complex>

struct AgmRec
{
  std::complex<double> a,g;
  int32_t d;
};

AgmRec agm1(AgmRec ag);
