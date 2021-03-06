/******************************************************/
/*                                                    */
/* agm.cpp - compute arithmetic-geometric mean        */
/*                                                    */
/******************************************************/
/* Copyright 2021,2022 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */

#include <cassert>
#include <iostream>
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

vector<complex<double> > agmLattice(complex<double> a,complex<double> g,unsigned depth,unsigned level,string branch)
{
  vector<complex<double> > ret,otherSide;
  AgmResult result;
  int i;
  if (level<depth)
  {
    ret=agmLattice(a,g,depth,level+1,branch);
    while (branch.length()<level)
      branch+='\0';
    branch[level]^=0x80;
    result=agm(a,g,branch);
    otherSide=agmLattice(a,g,depth,level+1,result.branch);
    for (i=0;i<otherSide.size();i++)
      ret.push_back(otherSide[i]);
  }
  else
  {
    result=agm(a,g,branch);
    ret.push_back(result.m);
  }
  return ret;
}

complex<double> pvAgm(complex<double> a,complex<double> g)
{
  AgmResult res=agm(a,g,"");
  return res.m;
}

array<complex<double>,2> invAgm1(complex<double> a,complex<double> g)
/* Returns two numbers whose arithmetic mean is a and geometric mean is g.
 * You are responsible for swapping them if necessary.
 */
{
  complex<double> diffsq=(a+g)*(a-g); // a*a-g*g is imprecise if a is close to g
  complex<double> rt=sqrt(diffsq);
  array<complex<double>,2> ret;
  if (abs(a-rt)>abs(a+rt))
    rt=-rt;
  ret[0]=a+rt;
  if (abs(a+rt)>2*abs(a-rt))
    ret[1]=g*g/ret[0];
  else
    ret[1]=a-rt;
  return ret;
}
