/******************************************************/
/*                                                    */
/* main.cpp - main program                            */
/*                                                    */
/******************************************************/
/* Copyright 2021-2023 Pierre Abbat
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
#include "deriv4.h"
#include "pairwisesum.h"
#include "ldecimal.h"
#include "ps.h"
#include "khe.h"
using namespace std;

const bool inverted=false;
Khe khe,khe85(85),khe221(221);
// the following should succeed
//Khe kheneg(-221); // same as 221
//Khe khe255(255); // same as 85 but thrice as big
// the following should abort
//Khe khebig(13*17*29*37*41*53);
//Khe khe1105(1105);
//Khe khemin(2147483648);

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
  maxx=khe.xt(sz);
  ps.startpage();
  ps.setscale(0,-1,3,1,0);
  ps.setcolor(1,1,0);
  for (i=0;i<mismatchse.size();i+=2)
  {
    ps.startline();
    ps.lineto(complex<double>((khe.xt(mismatchse[i])-2)/maxx*3,-1));
    ps.lineto(complex<double>((khe.xt(mismatchse[i+1])+2)/maxx*3,-1));
    ps.lineto(complex<double>((khe.xt(mismatchse[i+1])+2)/maxx*3,1));
    ps.lineto(complex<double>((khe.xt(mismatchse[i])-2)/maxx*3,1));
    ps.endline(true);
  }
  ps.setcolor(0,0,1);
  ps.startline();
  for (i=0;i<sz;i++)
    ps.lineto(complex<double>(khe.xt(i)/maxx*3,(logloop[i]-minlog)/(maxlog-minlog)));
  ps.lineto(complex<double>(3,(logloop[0]-minlog)/(maxlog-minlog)));
  ps.endline();
  ps.setcolor(1,0,0);
  ps.startline();
  for (i=0;i<sz;i++)
    ps.lineto(complex<double>(khe.xt(i)/maxx*3,(argloop[i]-maxarg)/(maxarg-minarg)));
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

void plotSquare(PostScript &ps,Khe f,complex<double> cen,complex<double> h)
{
  complex<double> fcen;
  vector<complex<double> > values;
  int i,j;
  double max=0,rad,r,g,b;
  const int sz=3;
  for (i=-sz;i<=sz;i++)
    for (j=-sz;j<=sz;j++)
      values.push_back(f(cen+complex<double>(i,j)*h));
  fcen=values[values.size()/2];
  for (i=0;i<values.size();i++)
  {
    values[i]-=fcen;
    if (abs(values[i])>max)
      max=abs(values[i]);
  }
  ps.startpage();
  ps.setscale(-max,-max,max,max,0);
  rad=max/sz/10;
  for (i=0;i<values.size();i++)
  {
    r=((i/(2*sz+1))+0.5)/(2*sz+1);
    b=((i%(2*sz+1))+0.5)/(2*sz+1);
    g=(2-r-b)/2;
    ps.setcolor(r,g,b);
    ps.circle(values[i],rad,true);
  }
  ps.endpage();
}

double relativeError(double x)
/* Computes the RMS error of the khe function. x should be a seam.
 */
{
  array<complex<double>,4> xsect;
  double h=x/1e7; // small enough that the third derivative should be less than an ulp
  double maxabs;
  complex<double> ctr;
  vector<double> errors;
  int i,j;
  for (i=0;i<355;i++)
  {
    ctr=khe(complex<double>(x,i));
    if (i==0)
      maxabs=ctr.real();
    for (j=0;j<4;j++)
      xsect[j]=khe(complex<double>(x+(j-1.5)*h,i));
    errors.push_back(norm(deriv3(xsect)/4./maxabs));
  }
  return sqrt(pairwisesum(errors)/errors.size());
}

array<double,3> zoomBounds(double x)
/* [0] is right bound, [1] is imag part of top and bottom, [2] is left.
 * x is negative.
 */
{
  array<double,3> ret;
  double infAsymp=4*exp(x); // approaches 1
  double zeroAsymp=-M_PI/x; // measured from 0
  double khex=khe(x).real();
  double weight=1-exp(-x*x/20);
  double inc;
  ret[0]=weight*(infAsymp+1)+(1-weight)*zeroAsymp;
  ret[1]=weight*(infAsymp)+(1-weight)*zeroAsymp/2;
  ret[2]=weight*(1-infAsymp)-(1-weight)*zeroAsymp/3;
  ret[1]*=4/3.; // enlarge bounds
  inc=(ret[0]-ret[2])/6.;
  ret[2]-=inc;
  ret[0]+=inc;
  return ret;
}

void drawGrid(PostScript &ps,const array<double,3> &bounds)
{
  int i,right,top,left;
  double gridOrigin;
  double gridSpacing=1;
  double mantissa;
  while (gridSpacing<bounds[1])
    gridSpacing*=10;
  while (gridSpacing>bounds[1]/3)
    gridSpacing/=10;
  if (bounds[2]>0)
    gridOrigin=1;
  else
    gridOrigin=0;
  mantissa=log10(bounds[1]/3/gridSpacing);
  left=lrint((bounds[2]-gridOrigin)/gridSpacing);
  top=lrint(bounds[1]/gridSpacing);
  right=lrint((bounds[0]-gridOrigin)/gridSpacing);
  for (i=left;i<=right;i++)
  {
    if (i%10==0)
      ps.setcolor(1,0,1);
    else
      ps.setcolor(1,mantissa,1);
    ps.line2p(complex<double>(gridOrigin+i*gridSpacing,-bounds[1]),
	      complex<double>(gridOrigin+i*gridSpacing,bounds[1]));
  }
  for (i=-top;i<=top;i++)
  {
    if (i%10==0)
      ps.setcolor(1,0,1);
    else
      ps.setcolor(1,mantissa,1);
    ps.line2p(complex<double>(bounds[2],i*gridSpacing),
	      complex<double>(bounds[0],i*gridSpacing));
  }
}

void zoomOut()
/* Draw graphs of loops from x=-32 to -1/15 in geometric progression.
 * This is almost 9 octaves; at 30 frames per second, it should take 60 frames
 * per octave for 18 seconds.
 */
{
  PostScript ps;
  const int framesPerOctave=60;
  vector<double> xcoord;
  array<double,3> bounds;
  double x;
  int i,j,n;
  for (i=0;i<framesPerOctave;i++)
    xcoord.push_back(-32*pow(0.5,(double)i/framesPerOctave));
  ps.open("zoom.ps");
  ps.setpaper(papersizes["A4 landscape"],0);
  ps.prolog();
  for (i=0;;i++)
  {
    x=xcoord[i%framesPerOctave]/(1<<(i/framesPerOctave));
    if (x>-1/15.)
      break;
    cout<<x<<endl;
    bounds=zoomBounds(x);
    ps.startpage();
    ps.setscale(bounds[2],-bounds[1],bounds[0],bounds[1]);
    drawGrid(ps,bounds);
    ps.setcolor(0,0,0);
    n=lrint(-1024/x);
    ps.startline();
    for (j=0;j<n;j++)
      ps.lineto(khe(complex<double>(x,j*2*M_PI/n)));
    ps.endline(true);
    ps.endpage();
  }
}

void modform1(complex<double> z)
{
  complex<double> invz=4.*M_PI*M_PI/z;
  complex<double> ratio=khe(invz)/khe(z);
  cout<<"At scaled z "<<-z/2./M_PI<<", ratio "<<ratio<<endl;
}

void modform()
/* Check whether խ(z), with suitable scaling, is a modular form.
 * Modular forms are defined in the upper half plane and have period 1.
 * խ(z) is defined in the left half plane and has period 2πi.
 */
{
  cout<<"խ(-2π)="<<khe(-2*M_PI)<<endl;
  modform1(-2*M_PI);
  modform1(-3*M_PI);
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
  double hi,lo,mid;
  complex<double> pnt0,pnt1,avg;
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
  hi=-12;
  lo=-24;
  mid=-16;
  while (hi>mid && mid>lo)
  {
    loops[0]=khe.getLoop(mid);
    if (loops[0].size()>64)
      hi=mid;
    else
      lo=mid;
    mid=(hi+lo)/2;
  }
  cout<<"Loop size doubles at "<<ldecimal(mid)<<endl;
  /* Draw a 355-pointed star with each point two radians after the previous one.
   * The last line is not exact, but as 355/113 is a good approximation to π,
   * it is really close.
   */
  ps.startpage();
  ps.setcolor(0,0,1);
  ps.setscale(1-2e-7,-2e-7,1+2e-7,2e-7);
  ps.startline();
  for (i=0;i<355;i++)
    ps.lineto(khe(complex<double>(-17,i*2)));
  ps.endline(true);
  ps.endpage();
  /* Draw just before and just after the loop doubles in size, with the distance
   * between them exaggerated. Most of the difference is from interpolating
   * between points that are closer together. The difference between the loop
   * just before, which is 36 points of an exact circle, and the loop just after,
   * in which [0] and [180] have slightly greater real part and [90] and [270]
   * slightly less, is negligible.
   */
  ps.startpage();
  ps.setscale(1-2e-7,-2e-7,1+2e-7,2e-7);
  loops[0].clear();
  loops[1].clear();
  for (i=0;i<360;i++)
  {
    pnt0=khe(complex<double>(mid-1e-8,degtorad(i)));
    pnt1=khe(complex<double>(mid+1e-8,degtorad(i)));
    avg=(pnt0+pnt1)/2.;
    loops[0].push_back(avg+2e2*(pnt0-avg));
    loops[1].push_back(avg+2e2*(pnt1-avg));
  }
  ps.setcolor(0,0,1);
  ps.startline();
  for (i=0;i<360;i++)
    ps.lineto(loops[0][i]);
  ps.endline(true);
  ps.setcolor(1,0,0);
  ps.startline();
  for (i=0;i<360;i++)
    ps.lineto(loops[1][i]);
  ps.endline(true);
  ps.endpage();
  // Compute the error of the khe function.
  for (i=0;i<12;i++)
  {
    cout<<"Relative error at "<<mid/(1<<i)<<" is "<<
	  ldecimal(relativeError(mid/(1<<i)))<<endl;
    cout<<"Relative error at "<<mid/1.5/(1<<i)<<" is "<<
	  ldecimal(relativeError(mid/1.5/(1<<i)))<<endl;
  }
  for (i=0;i<12;i++)
  {
    loops[i]=khe.getLoop(-33./(1<<i));
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
    cout<<"Iter "<<i<<" Size "<<logloop.size()<<" Bounds "<<minreal<<' '<<maxreal<<' '<<maximag<<endl;
    for (j=0;j<loops[i].size();j++)
    {
      ps.dot(inverted?1./loops[i][j]:loops[i][j]);
    }
    ps.endpage();
    plotLogArg(ps,logloop,argloop);
    khe.outMaxMag(loops[i]);
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
  plotSquare(ps,khe,complex<double>(-16,0.),0.00002);
  plotSquare(ps,khe,complex<double>(-16,0.),0.002);
  plotSquare(ps,khe,complex<double>(-16.974354,0.),0.002);
  plotSquare(ps,khe,complex<double>(-16.974354,0.1),0.002);
  plotSquare(ps,khe,complex<double>(-16.974354,0.),0.000001);
  for (i=0;i<3;i++)
  {
    diam[i]=avgRadius(loops[i]);
    cout<<i<<' '<<ldecimal(diam[i])<<' '<<ldecimal(log(diam[i]))<<endl;
  }
  cout<<pvAgm(pvAgm(2,3),pvAgm(5,7))<<' ';
  cout<<pvAgm(pvAgm(2,5),pvAgm(3,7))<<' ';
  cout<<pvAgm(pvAgm(2,7),pvAgm(5,3))<<endl;
  /* The loop at x, as x approaches 0 from below, fits in a box {u+vi |
   * -t/3<u<t, -t/2<v<t/2}, where t*x=-π. This does not hold for x much less
   * than -1.
   */
  cout<<"Scale factor for loop box "<<ldecimal(khe(-1/64.).real()/64)<<endl;
  cout<<"65:  "<<khe(complex<double>(-1,-1))<<endl;
  cout<<"85:  "<<khe85(complex<double>(-1,-1))<<endl;
  cout<<"221: "<<khe221(complex<double>(-1,-1))<<endl;
  zoomOut();
  modform();
  return 0;
}
