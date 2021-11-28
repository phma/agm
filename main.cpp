/******************************************************/
/*                                                    */
/* main.cpp - main program                            */
/*                                                    */
/******************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "agm.h"
using namespace std;

int main(int argc,char **argv)
{
  AgmResult ag;
  ag=agm(sqrt(2),1);
  cout<<ag.m<<endl;
  return 0;
}
  
