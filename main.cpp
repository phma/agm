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

void printag(AgmResult ag)
{
  int i;
  cout<<ag.m;
  for (i=0;i<ag.branch.length();i++)
    cout<<' '<<hex<<setw(2)<<(int)(ag.branch[i]&0xff)<<dec;
  cout<<endl;
}

int main(int argc,char **argv)
{
  AgmResult ag;
  int i;
  for (i=0;i<32;i++)
  {
    ag=agm(cossin(DEG45*i)*0.5,1,ag.branch);
    printag(ag);
  }
  return 0;
}
  
