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
  int i;
  ag=agm(cossin(0x3fffffff)*0.999999999,1);
  cout<<ag.m;
  for (i=0;i<ag.branch.length();i++)
    cout<<' '<<hex<<setw(2)<<(int)(ag.branch[i]&0xff)<<dec;
  cout<<endl;
  return 0;
}
  
