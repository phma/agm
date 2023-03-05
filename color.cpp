/******************************************************/
/*                                                    */
/* color.cpp - colors for domain coloring             */
/*                                                    */
/******************************************************/
/* Copyright 2020,2021,2023 Pierre Abbat.
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */

#include <cfloat>
#include "color.h"
using namespace std;

double compand(double mag)
/* խ(z) has two types of singularities on the imaginary axis. At 0, it
 * approaches ∞ like -π/z, but at πi, it approaches 0 somewhat like exp(1/(z-πi)).
 * If not companded, this would make the region near πi all black, obscuring
 * the detail. The companding function is the reciprocal of the Lambert W
 * function of the reciprocal.
 */
{
  double ret,pexp,corr=1;
  int i=0;
  mag=1/mag;
  if (mag>3)
    ret=log(mag)-log(log(mag));
  else
    ret=mag;
  while (i<100 && fabs(corr)>1e-15)
  {
    pexp=ret*exp(ret);
    corr=(pexp-mag)/(ret+1)/exp(ret);
    ret-=corr;
    i++;
  }
  return 1/ret;
}

Color::Color(double red,double green,double blue)
{
  r=min(1.,max(0.,red));
  g=min(1.,max(0.,green));
  b=min(1.,max(0.,blue));
}

Color::Color()
{
  r=g=b=0;
}

string Color::ppm()
{
  string ret("rgb");
  double r8,g8,b8;
  r8=trunc(r*256);
  g8=trunc(g*256);
  b8=trunc(b*256);
  if (r8>255)
    r8=255;
  if (g8>255)
    g8=255;
  if (b8>255)
    b8=255;
  ret[0]=r8;
  ret[1]=g8;
  ret[2]=b8;
  return ret;
}

void Color::mix(const Color &diluent,double part)
{
  r=(r*(1-part)+diluent.r*part);
  g=(g*(1-part)+diluent.g*part);
  b=(b*(1-part)+diluent.b*part);
}

const Color white(1,1,1);

Colorize::Colorize()
{
  high=low=NAN;
  ori=0;
  scheme=CS_HV;
}

void Colorize::setLimits(double l,double h)
/* l and h are the minimum and maximum absolute values of the khe function
 * in the area being plotted, truncated to the rightmost line of dots. If l
 * is 0, it is replaced with the smallest positive number.
 */
{
  if (l<DBL_MIN)
    l=DBL_MIN;
  low=log(compand(l));
  high=log(compand(h));
}

void Colorize::setOrientation(int o)
{
  ori=o;
}

void Colorize::setScheme(int s)
{
  scheme=s;
}

int Colorize::getScheme()
{
  return scheme;
}

Color Colorize::operator()(complex<double> z)
{
  double a=arg(z);
  double r=abs(z);
  double cp,bright,chroma,octave,modOctave;
  if (r==0)
    cp=low;
  else if (isnan(r))
    cp=(low+high)/2;
  else
    cp=log(compand(r));
  bright=(cp-low)/(high-low);
  octave=log(r)/log(2);
  modOctave=octave-floor(octave);
  if (isnan(r))
    chroma=0;
  else
    chroma=bright*(1-bright)*(1-modOctave/3);
  return Color(bright+chroma*cos(a),
	       bright+chroma*cos(a+2*M_PI/3),
	       bright+chroma*cos(a-2*M_PI/3));
}
