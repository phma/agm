/******************************************************/
/*                                                    */
/* color.h - colors for domain coloring               */
/*                                                    */
/******************************************************/
/* Copyright 2020,2021,2023 Pierre Abbat.
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */
#ifndef COLOR_H
#define COLOR_H
class ContourInterval;
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>

#define CS_HV 0

class Color
{
public:
  Color(double red,double green,double blue);
  Color();
  double fr()
  {
    return r;
  }
  double fg()
  {
    return g;
  }
  double fb()
  {
    return b;
  }
  int br()
  {
    return std::min(255,(int)floor(255*r));
  }
  int bg()
  {
    return std::min(255,(int)floor(255*g));
  }
  int bb()
  {
    return std::min(255,(int)floor(255*b));
  }
  std::string ppm();
  void mix(const Color &diluent,double part);
private:
  double r,g,b;
};

extern const Color white;

class Colorize
{
public:
  Colorize();
  void setLimits(double l,double h);
  void setOrientation(int o);
  void setScheme(int s);
  int getScheme();
  Color operator()(std::complex<double> z);
private:
  double low;
  double high;
  int ori;
  int scheme;
};

#endif
