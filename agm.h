/******************************************************/
/*                                                    */
/* agm.h - compute arithmetic-geometric mean          */
/*                                                    */
/******************************************************/

#include <cmath>
#include <complex>
#include <string>

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
