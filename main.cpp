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
#include "ps.h"
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
  PostScript ps;
  int i;
  ps.open("agm.ps");
  ps.setpaper(papersizes["A4 landscape"],0);
  ps.prolog();
  ps.startpage();
  ps.setscale(-2,-1.5,2,1.5,0);
  ps.setcolor(0,0,1);
  ps.startline();
  for (i=0;i<32;i++)
  {
    ag=agm(cossin(DEG45*i)*0.5,1,ag.branch);
    ps.lineto(ag.m);
    printag(ag);
  }
  ps.endline();
  return 0;
}
  
