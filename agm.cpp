/******************************************************/
/*                                                    */
/* agm.cpp - compute arithmetic-geometric mean        */
/*                                                    */
/******************************************************/

#include "agm.h"
#include "angle.h"
using namespace std;

/* The arithmetic-geometric mean is a transcendental function of two arguments
 * with three branch points and countably infinitely many branches. Since
 * agm(ax,ay)=a*agm(x,y), it can be curried into a function of one argument
 * by holding the other fixed. With the argument fixed at 1, the branch points
 * are 0, -1, and 1.
 *
 * The branch is passed as the third argument to the agm function and is
 * represented as a string of bytes. It could be represented as a bit string,
 * but a byte string allows one to pass smoothly from one branch to another.
 * If it's a bit string, only finitely many of the bits can be 1, as for every
 * 1 bit, abs(a)+abs(g) decreases by at least 14%. The bits or bytes represent
 * the angle between a and g at each iteration.
 */

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
      in.d=branch[i]<<23;
    else
      in.d=0;
    out=agm1(in);
    ret.branch+=(char)((out.d+0x400000)>>23);
    ++i;
    if (i>=branch.length() && (out.a==out.g || out.g==0.
			       || (out.a==in.a && out.g==in.g)))
      break;
    in=out;
  }
  ret.m=out.g;
  return ret;
}
