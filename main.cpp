/******************************************************/
/*                                                    */
/* main.cpp - main program                            */
/*                                                    */
/******************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cfloat>
#include "agm.h"
#include "angle.h"
#include "pairwisesum.h"
#include "ldecimal.h"
#include "ps.h"
using namespace std;

const signed char c65[]=
{
  65,63,60,56,52,39,33,25,16,
  0,-16,-25,-33,-39,-52,-56,-60,-63,
  -65,-63,-60,-56,-52,-39,-33,-25,-16,
  0,16,25,33,39,52,56,60,63
};
const bool inverted=false;
const char a65[]={0,7,11,15,18,26,29,33,37,44};

vector<complex<double> > tinyCircle(complex<double> center)
/* Returns a circle with radius 65 ulps. center must be between 2 and -2,
 * exclusive, or the result will be inaccurate.
 */
{
  vector<complex<double> > ret;
  int i;
  for (i=0;i<36;i++)
    ret.push_back(center+complex<double>(c65[i],c65[(i+27)%36])*DBL_EPSILON);
  return ret;
}

double avgRadius(vector<complex<double> > loop)
{
  int i,sz=loop.size();
  vector<double> diams;
  diams.resize(sz/2);
  for (i=0;i<sz/2;i++)
    diams[i]=abs(loop[i]-loop[i+sz/2]);
  return pairwisesum(diams)/sz;
}

void printag(AgmResult ag)
{
  int i;
  cout<<ag.m;
  for (i=0;i<ag.branch.length();i++)
    cout<<' '<<hex<<setw(2)<<(int)(ag.branch[i]&0xff)<<dec;
  cout<<endl;
}

vector<double> vecLog(vector<complex<double> > loop)
{
  vector<double> ret;
  int i;
  for (i=0;i<loop.size();i++)
    ret.push_back(log(abs(loop[i])));
  return ret;
}

vector<double> vecArg(vector<complex<double> > loop)
{
  vector<double> ret;
  double a,lasta=0;
  int i;
  for (i=0;i<loop.size();i++)
  {
    a=arg(loop[i]);
    a+=2*M_PI*rint((lasta-a)/2/M_PI);
    ret.push_back(a);
    lasta=a;
  }
  return ret;
}

int xt(int n)
{
  return a65[9]*(n/9)+a65[n%9];
}

void plotLogArg(PostScript &ps,vector<double> &logloop,vector<double> &argloop)
{
  int i,sz=argloop.size();
  vector<int> mismatch,mismatchse;
  double minlog=0,maxlog=0,minarg=0,maxarg=0,maxx;
  assert(logloop.size()==argloop.size());
  for (i=0;i<sz;i++)
  {
    if (logloop[i]>maxlog)
      maxlog=logloop[i];
    if (argloop[i]>maxarg)
      maxarg=argloop[i];
    if (logloop[i]<minlog)
      minlog=logloop[i];
    if (argloop[i]<minarg)
      minarg=argloop[i];
    if (i && (logloop[i]!=logloop[sz-i] || argloop[i]+argloop[sz-i]!=0))
      mismatch.push_back(i);
  }
  cout<<mismatch.size()<<" mismatches\n";
  for (i=0;i<mismatch.size();i++)
  {
    if (i==0)
      mismatchse.push_back(mismatch[i]);
    else if (mismatch[i]>mismatch[i-1]+1)
    {
      mismatchse.push_back(mismatch[i-1]+1);
      mismatchse.push_back(mismatch[i]);
    }
  }
  if (i)
    mismatchse.push_back(mismatch[i-1]+1);
  maxx=xt(sz);
  ps.startpage();
  ps.setscale(0,-1,3,1,0);
  ps.setcolor(1,1,0);
  for (i=0;i<mismatchse.size();i+=2)
  {
    ps.startline();
    ps.lineto(complex<double>((xt(mismatchse[i])-2)/maxx*3,-1));
    ps.lineto(complex<double>((xt(mismatchse[i+1])+2)/maxx*3,-1));
    ps.lineto(complex<double>((xt(mismatchse[i+1])+2)/maxx*3,1));
    ps.lineto(complex<double>((xt(mismatchse[i])-2)/maxx*3,1));
    ps.endline(true);
  }
  ps.setcolor(0,0,1);
  ps.startline();
  for (i=0;i<sz;i++)
    ps.lineto(complex<double>(xt(i)/maxx*3,(logloop[i]-minlog)/(maxlog-minlog)));
  ps.lineto(complex<double>(3,(logloop[0]-minlog)/(maxlog-minlog)));
  ps.endline();
  ps.setcolor(1,0,0);
  ps.startline();
  for (i=0;i<sz;i++)
    ps.lineto(complex<double>(xt(i)/maxx*3,(argloop[i]-maxarg)/(maxarg-minarg)));
  ps.lineto(complex<double>(3,(argloop[0]-maxarg)/(maxarg-minarg)));
  ps.endline();
  ps.endpage();
}

int main(int argc,char **argv)
{
  AgmResult ag;
  PostScript ps;
  vector<complex<double> > lattice;
  vector<vector<complex<double> > > loops;
  vector<double> logloop,argloop;
  int i,j;
  double minreal=1,maxreal=1,maximag=0,diam[3];
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
  loops.resize(12);
  for (i=0;i<11;i++) //12
  {
    if (!i)
      loops[i]=tinyCircle(1);
    else
      loops[i]=agmExpand(loops[i-1]);
    ps.startpage();
    ps.setcolor(0,0,1);
    for (j=0;j<loops[i].size();j++)
    {
      if (inverted)
      {
	if ((1./loops[i][j]).real()>maxreal)
	  maxreal=(1./loops[i][j]).real();
	if ((1./loops[i][j]).real()<minreal)
	  minreal=(1./loops[i][j]).real();
	if ((1./loops[i][j]).imag()>maximag)
	  maximag=(1./loops[i][j]).imag();
      }
      else
      {
	if (loops[i][j].real()>maxreal)
	  maxreal=loops[i][j].real();
	if (loops[i][j].real()<minreal)
	  minreal=loops[i][j].real();
	if (loops[i][j].imag()>maximag)
	  maximag=loops[i][j].imag();
      }
    }
    logloop=vecLog(loops[i]);
    argloop=vecArg(loops[i]);
    ps.setscale(minreal,-maximag,maxreal,maximag);
    cout<<"Iter "<<i<<" Bounds "<<minreal<<' '<<maxreal<<' '<<maximag<<endl;
    for (j=0;j<loops[i].size();j++)
    {
      ps.dot(inverted?1./loops[i][j]:loops[i][j]);
    }
    ps.endpage();
    plotLogArg(ps,logloop,argloop);
  }
  for (i=0;i<3;i++)
  {
    diam[i]=avgRadius(loops[i]);
    cout<<i<<' '<<ldecimal(diam[i])<<' '<<ldecimal(log(diam[i]))<<endl;
  }
  cout<<ldecimal(log(diam[0]/diam[1])/log(diam[1]/diam[2]))<<" should be 2\n";
  cout<<"y for 65-ulp loop around 1 is "<<ldecimal(2*(log(diam[1]/diam[0])))<<endl;
  cout<<pvAgm(pvAgm(2,3),pvAgm(5,7))<<' ';
  cout<<pvAgm(pvAgm(2,5),pvAgm(3,7))<<' ';
  cout<<pvAgm(pvAgm(2,7),pvAgm(5,3))<<endl;
  return 0;
}
  
