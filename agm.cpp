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

AgmResult agm(complex<double> a,complex<double> g,string branch)
{
  AgmResult ret;
  AgmRec in,out;
  size_t pos;
  int i=0;
  in.a=a;
  in.g=g;
  pos=branch.find_last_not_of('\0');
  if (pos<branch.length())
    branch.erase(pos+1);
  else
    branch.clear();
  while (true)
  {
    if (i<branch.length())
      in.d=branch[i]<<24;
    else
      in.d=0;
    out=agm1(in);
    ret.branch+=(char)((out.d+0x800000)>>24);
    ++i;
    if (i>=branch.length() && (out.a==out.g || out.g==0.
			       || (out.a==in.a && out.g==in.g)))
      break;
    in=out;
  }
  ret.m=out.g;
  return ret;
}
