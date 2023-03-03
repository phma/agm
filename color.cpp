/******************************************************/
/*                                                    */
/* color.cpp - colors for domain coloring             */
/*                                                    */
/******************************************************/
/* Copyright 2020,2021,2023 Pierre Abbat.
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */

#include "color.h"
using namespace std;

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
{
  low=l;
  high=h;
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
