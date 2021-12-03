/******************************************************/
/*                                                    */
/* angle.cpp - angles as binary fractions of rotation */
/*                                                    */
/******************************************************/

#include "angle.h"
using namespace std;

double sin(int angle)
{
  return sinl(angle*M_PIl/1073741824.);
}

double cos(int angle)
{
  return cosl(angle*M_PIl/1073741824.);
}

complex<double> cossin(int angle)
{
  return complex<double>(cos(angle),sin(angle));
}

int argi(complex<double> z)
{
  return rint(arg(z)/M_PIl*1073741824.);
}

double bintorot(int angle)
{
  return angle/2147483648.;
}

double bintogon(int angle)
{
  return bintorot(angle)*400;
}

double bintodeg(int angle)
{
  return bintorot(angle)*360;
}

double bintomin(int angle)
{
  return bintorot(angle)*21600;
}

double bintosec(int angle)
{
  return bintorot(angle)*1296000;
}

double bintorad(int angle)
{
  return bintorot(angle)*M_PIl*2;
}

int rottobin(double angle)
{
  double iprt=0,fprt;
  fprt=2*modf(angle/2,&iprt);
  if (fprt>=1)
    fprt-=2;
  if (fprt<-1)
    fprt+=2;
  return (int)llrint(2147483648.*fprt);
}

int degtobin(double angle)
{
  return rottobin(angle/360);
}

int mintobin(double angle)
{
  return rottobin(angle/21600);
}

int sectobin(double angle)
{
  return rottobin(angle/1296000);
}

int gontobin(double angle)
{
  return rottobin(angle/400);
}

int radtobin(double angle)
{
  return rottobin(angle/M_PIl/2);
}

double radtodeg(double angle)
{
  return angle*180/M_PIl;
}

double degtorad(double angle)
{
  return angle/180*M_PIl;
}

double radtomin(double angle)
{
  return angle*10800/M_PIl;
}

double mintorad(double angle)
{
  return angle/10800*M_PIl;
}

double radtosec(double angle)
{
  return angle*648000/M_PIl;
}

double sectorad(double angle)
{
  return angle/648000*M_PIl;
}

double radtogon(double angle)
{
  return angle*200/M_PIl;
}

double gontorad(double angle)
{
  return angle/200*M_PIl;
}
