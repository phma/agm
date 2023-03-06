/******************************************************/
/*                                                    */
/* raster.h - raster image output                     */
/*                                                    */
/******************************************************/
/* Copyright 2023 Pierre Abbat.
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */
#include "khe.h"

class FordCircle
{
public:
  double y,radius;
  std::complex<double> residue; // set to 0 for essential singularity
  std::complex<double> pole(std::complex<double> z);
  bool in(std::complex<double> z);
  double farIn(std::complex<double> z);
};

void rasterplot(Khe &khe,int width,int height,std::string filename);
