/******************************************************/
/*                                                    */
/* main.cpp - main program                            */
/*                                                    */
/******************************************************/
/* Copyright 2021,2022 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */
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
#include "khe.h"
using namespace std;

const bool inverted=false;

void printag(AgmResult ag)
{
  int i;
  cout<<ag.m;
  for (i=0;i<ag.branch.length();i++)
    cout<<' '<<hex<<setw(2)<<(int)(ag.branch[i]&0xff)<<dec;
  cout<<endl;
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

string complexStr(complex<double> z)
{
  return ldecimal(z.real())+','+ldecimal(z.imag());
}

void outMismatch(const vector<complex<double> > &prevLoop,const vector<complex<double> > &thisLoop)
{
  int i,sz=thisLoop.size();
  ofstream file("mismatch.html");
  file<<"<html><head><title>Mismatches</title></head>\n";
  file<<"<body>Max "<<complexStr(thisLoop[0]);
  file<<"\n<table border>\n";
  assert(sz==73728);
  for (i=9200;i<9233;i++)
  {
    file<<"<tr><th>"<<i<<"</th>";
    file<<"<td>"<<complexStr(thisLoop[i])<<"</td>";
    file<<"<td>"<<complexStr(thisLoop[sz/2-i])<<"</td>";
    file<<"<td>"<<complexStr(thisLoop[sz/2+i])<<"</td>";
    file<<"<td>"<<complexStr(thisLoop[sz-i])<<"</td>";
    if (i%2==0)
    {
      file<<"<td>"<<complexStr(prevLoop[i/2]*2.)<<"</td>";
      file<<"<td>"<<complexStr(prevLoop[sz/4-i/2]*2.)<<"</td>";
      file<<"<td>"<<complexStr(prevLoop[sz/4+i/2]*2.)<<"</td>";
      file<<"<td>"<<complexStr(prevLoop[sz/2-i/2]*2.)<<"</td>";
    }
    else
      file<<"<td></td><td></td><td></td><td></td>";
    file<<"</tr>\n";
  }
  file<<"</table></body></html>\n";
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
  double x65;
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
  for (i=0;i<12;i++)
  {
    //if (!i)
      //loops[i]=tinyCircle(1);
    //else
      //loops[i]=agmExpand(loops[i-1]);
    loops[i]=getLoop(-pow(0.5,i/12.));
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
    outMaxMag(loops[i]);
  }
  ps.startpage();
  ps.setscale(-2/3.,-1,2,1,0);
  ps.setcolor(0,0,1);
  for (i=1;i<64;i++)
    if (i&1)
      if (i&2)
	ps.circle(-1./i,1./i);
      else
	ps.circle(1./i,1./i);
    else
    {
      ps.circle(complex<double>(0.,1./i),1./i);
      ps.circle(complex<double>(0.,-1./i),1./i);
    }
  ps.endpage();
  for (i=0;i<3;i++)
  {
    diam[i]=avgRadius(loops[i]);
    cout<<i<<' '<<ldecimal(diam[i])<<' '<<ldecimal(log(diam[i]))<<endl;
  }
  cout<<ldecimal(log(diam[0]/diam[1])/log(diam[1]/diam[2]))<<" should be 2\n";
  x65=-2*(log(diam[1]/diam[0]));
  cout<<"x for 65-ulp loop around 1 is "<<ldecimal(x65)<<endl;
  cout<<"Center for x="<<x65<<" is "<<ldecimal(circleCenter(x65))<<endl;
  cout<<"Center for x="<<x65/2<<" is "<<ldecimal(circleCenter(x65/2))<<endl;
  cout<<pvAgm(pvAgm(2,3),pvAgm(5,7))<<' ';
  cout<<pvAgm(pvAgm(2,5),pvAgm(3,7))<<' ';
  cout<<pvAgm(pvAgm(2,7),pvAgm(5,3))<<endl;
  return 0;
}
  
