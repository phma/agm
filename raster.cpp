/******************************************************/
/*                                                    */
/* raster.cpp - raster image output                   */
/*                                                    */
/******************************************************/
/* Copyright 2023 Pierre Abbat.
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */
#include <cstdio>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "raster.h"
#include "color.h"

using namespace std;
fstream rfile;

complex<double> FordCircle::pole(complex<double> z)
{
  return residue/(z-y);
}

bool FordCircle::in(complex<double> z)
{
  return abs(z-complex<double>(-radius,y))<=radius;
}

double FordCircle::farIn(complex<double> z)
{
  return radius-abs(z-complex<double>(-radius,y));
}

void ropen(string fname)
{
  if (fname=="")
    fname="/dev/stdout";
  rfile.open(fname.c_str(),ios_base::out|ios_base::binary);
}

void rclose()
{
  rfile.close();
}

void ppmheader(int width,int height)
{
  rfile<<"P6\n"<<width<<" "<<height<<endl<<255<<endl;
}

void rasterplot(Khe &khe,int width,int height,string filename)
/* The image extends from -πi to πi and from 0 as far left as determined
 * by width and height, the pixels being square.
 */
{
  int i,j,k;
  Color pixel;
  Colorize col;
  complex<double> pnt,z;
  double scale;
  char letter;
  ropen(filename);
  if (width<0 || height<=0)
    throw(range_error("rasterdraw: paper size must be nonnegative"));
  scale=2*M_PI/height;
  col.setLimits(abs(khe(complex<double>(-scale/2,M_PI))),abs(khe(-scale/2)));
  ppmheader(width,height);
  for (i=0;i<height;i++)
    for (j=0;j<width;j++)
    {
      pnt=complex<double>((j-width+0.5)*scale,(height/2.-i-0.5)*scale);
      z=khe(pnt);
      pixel=col(z);
      rfile<<pixel.ppm();
    }
  rclose();
}
