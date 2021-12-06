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
  vector<complex<double> > lattice;
  int i;
  ps.open("agm.ps");
  ps.setpaper(papersizes["A4 landscape"],0);
  ps.prolog();
  ps.startpage();
  ps.setscale(-2,-1.5,2,1.5,0);
  ps.setcolor(0,0,1);
  ps.startline();
  for (i=0;i<1440;i+=1)
  {
    ag=agm(cossin(degtobin(i))*0.9,1,ag.branch);
    ps.lineto(ag.m);
    //printag(ag);
  }
  ps.endline();
  ps.endpage();
  ps.startpage();
  ps.setcolor(0,0,1);
  ps.startline();
  ag.branch="";
  ag.m=1;
  for (i=0;i<1440;i+=1)
  {
    ag=agm(cossin(degtobin(i))*1.1,1,ag.branch);
    ps.lineto(ag.m);
    //printag(ag);
  }
  ps.endline();
  ps.endpage();
  ps.startpage();
  ps.setcolor(0,0,1);
  ps.startline();
  ag.branch="";
  ag.m=1;
  for (i=0;i>-180;i--)
  {
    ag=agm(cossin(degtobin(i)),1,ag.branch);
  }
  for (i=-180;i<=180;i+=1)
  {
    ag=agm(cossin(degtobin(i)),1,ag.branch);
    ps.lineto(ag.m);
    //printag(ag);
  }
  ps.endline();
  ps.endpage();
  ps.startpage();
  ps.setcolor(0,0,1);
  lattice=agmLattice(sqrt(2),1,8);
  for (i=0;i<lattice.size();i++)
  {
    ps.dot(lattice[i]);
  }
  ps.endpage();
  return 0;
}
  
