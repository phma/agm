/******************************************************/
/*                                                    */
/* pairwisesum.cpp - add many numbers                 */
/*                                                    */
/******************************************************/
/* Copyright 2015,2016,2019,2021 Pierre Abbat.
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */
#include <cfloat>
#include <cstring>
#include <iostream>
#include <cmath>
#include <vector>
#include "pairwisesum.h"
using namespace std;

double pairwisesum(double *a,unsigned n)
/* Adds up n numbers at a pairwise, for less roundoff error.
 * i^(i+1) is 1 if i is even, 3 if i%4==1, 7 if i%8==3, etc.
 * i^(i+8) is 8 if i%16==0, 24 if i%32==16, etc.
 */
{
  unsigned i,j,b;
  double sums[32],sum=0;
  for (i=0;i+7<n;i+=8) // Add up groups of 8, leaving 0 to 7 numbers at the end
  {
    b=i^(i+8);
    if (b==8)
      sums[3]=(((a[i]+a[i+1])+(a[i+2]+a[i+3]))+((a[i+4]+a[i+5])+(a[i+6]+a[i+7])));
    else
    {
      sums[3]+=(((a[i]+a[i+1])+(a[i+2]+a[i+3]))+((a[i+4]+a[i+5])+(a[i+6]+a[i+7])));
      for (j=4;b>>(j+1);j++)
	sums[j]+=sums[j-1];
      sums[j]=sums[j-1];
    }
  }
  for (;i<n;i++) // Add up the last 0 to 7 numbers
  {
    b=i^(i+1);
    if (b==1)
      sums[0]=a[i];
    else
    {
      sums[0]+=a[i];
      for (j=1;b>>(j+1);j++)
	sums[j]+=sums[j-1];
      sums[j]=sums[j-1];
    }
  }
  for (i=0;i<32;i++)
    if ((n>>i)&1)
      sum+=sums[i];
  return sum;
}

long double pairwisesum(long double *a,unsigned n)
{
  unsigned i,j,b;
  long double sums[32],sum=0;
  for (i=0;i+7<n;i+=8)
  {
    b=i^(i+8);
    if (b==8)
      sums[3]=(((a[i]+a[i+1])+(a[i+2]+a[i+3]))+((a[i+4]+a[i+5])+(a[i+6]+a[i+7])));
    else
    {
      sums[3]+=(((a[i]+a[i+1])+(a[i+2]+a[i+3]))+((a[i+4]+a[i+5])+(a[i+6]+a[i+7])));
      for (j=4;b>>(j+1);j++)
	sums[j]+=sums[j-1];
      sums[j]=sums[j-1];
    }
  }
  for (;i<n;i++)
  {
    b=i^(i+1);
    if (b==1)
      sums[0]=a[i];
    else
    {
      sums[0]+=a[i];
      for (j=1;b>>(j+1);j++)
	sums[j]+=sums[j-1];
      sums[j]=sums[j-1];
    }
  }
  for (i=0;i<32;i++)
    if ((n>>i)&1)
      sum+=sums[i];
  return sum;
}

/* This is the original version of pairwisesum.
 * This code is left here to show what the optimized version is trying
 * to accomplish and as a reference for unit tests.
 */
#if 0
double pairwisesum(double *a,unsigned n)
// a is clobbered.
{
  unsigned i,j;
  if (n)
  {
    for (i=1;i<n;i*=2)
      for (j=0;j+i<n;j+=2*i)
	a[j]+=a[j+i];
    return a[0];
  }
  else
    return 0;
}
#endif

double pairwisesum(vector<double> &a)
{
  if (a.size())
    return pairwisesum(&a[0],a.size());
  else
    return 0;
}

long double pairwisesum(vector<long double> &a)
{
  if (a.size())
    return pairwisesum(&a[0],a.size());
  else
    return 0;
}
