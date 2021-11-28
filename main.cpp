/******************************************************/
/*                                                    */
/* main.cpp - main program                            */
/*                                                    */
/******************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "agm.h"
#include "angle.h"
using namespace std;

int main(int argc,char **argv)
{
  AgmResult ag;
  ag=agm(cossin(0x3fffffff)*1.000000001,1);
  cout<<ag.m<<endl;
  return 0;
}
  
