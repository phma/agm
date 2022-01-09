/******************************************************/
/*                                                    */
/* ps.cpp - PostScript output                         */
/*                                                    */
/******************************************************/
/* Copyright 2012-2019,2021,2022 Pierre Abbat
 * Licensed under the Apache License, Version 2.0.
 * This file is part of AGM.
 */
#include <cstdio>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cassert>
#include <iomanip>
#include "ldecimal.h"
#include "config.h"
#include "ps.h"
#include "angle.h"
using namespace std;

#define PAPERRES 0.004

char rscales[]={10,12,15,20,25,30,40,50,60,80};
const double PSPoint=25.4/72;
map<string,papersize> papersizes=
/* These mean the physical orientation of the paper in the printer. If you
 * want to print in landscape, but the paper is portrait in the printer,
 * set pageorientation to 1.
 */
{
  {"A4 portrait",{210000,297000}},
  {"A4 landscape",{297000,210000}},
  {"US Letter portrait",{215900,279400}},
  {"US Letter landscape",{279400,215900}},
  {"US Legal portrait",{215900,355600}},
  {"US Legal landscape",{355600,215900}}
};

int fibmod3(int n)
{
  int i,a,b;
  for (i=a=0,b=1;a<n;i++)
  {
    b+=a;
    a=b-a;
  }
  return (a==n)?(i%3):-1;
}

complex<double> turn(complex<double> pnt,int ori)
{
  return pnt*cossin(ori);
}

PostScript::PostScript()
{
  oldr=oldg=oldb=NAN;
  paper=complex<double>(210,297);
  scale=1;
  orientation=pages=0;
  indocument=inpage=inlin=false;
  psfile=nullptr;
}

PostScript::~PostScript()
{
  if (psfile)
    close();
}

void PostScript::setpaper(papersize pap,int ori)
/* ori is 0 for no rotation, 1 for 90° rotation, making portrait
 * into landscape and vice versa. Do this before each page,
 * or before calling prolog if all pages are the same size.
 */
{
  paper=complex<double>(pap.width/1e3,pap.height/1e3);
  pageorientation=ori;
}

double PostScript::aspectRatio()
// Returns >1 for landscape, <1 for portrait.
{
  if (pageorientation&1)
    return paper.imag()/paper.real();
  else
    return paper.real()/paper.imag();
}

void PostScript::open(string psfname)
{
  if (psfile)
    close();
  psfile=new ofstream(psfname);
}

bool PostScript::isOpen()
{
  return psfile!=nullptr;
}

void PostScript::prolog()
{
  if (psfile && !indocument)
  {
    *psfile<<"%!PS-Adobe-3.0\n%%BeginProlog\n%%%%Pages: (atend)"<<endl;
    *psfile<<"%%BoundingBox: 0 0 "<<rint(paper.real()/PSPoint)<<' '<<rint(paper.imag()/PSPoint)<<endl;
    *psfile<<"\n/. % ( x y )\n{ newpath 0.1 0 360 arc fill } bind def\n\n";
    *psfile<<"/- % ( x1 y1 x2 y2 )\n\n{ newpath moveto lineto stroke } bind def\n\n";
    *psfile<<"/c. % ( str )\n{ dup stringwidth -2 div exch -2 div exch\n"<<
            "3 2 roll 2 index 2 index rmoveto show rmoveto } bind def\n\n";
    *psfile<<"/mmscale { 720 254 div dup scale } bind def\n";
    *psfile<<"/col { setrgbcolor } def\n";
    *psfile<<"/n { newpath } def\n";
    *psfile<<"/m { moveto } def\n";
    *psfile<<"/l { lineto } def\n";
    *psfile<<"/c { curveto } def\n";
    *psfile<<"/s { stroke } def\n";
    *psfile<<"/af { arc fill } def\n";
    *psfile<<"/as { arc stroke } def\n";
    *psfile<<"%%EndProlog"<<endl;
    indocument=true;
    pages=0;
  }
}

void PostScript::startpage()
{
  if (psfile && indocument && !inpage)
  {
    ++pages;
    *psfile<<"%%Page: "<<pages<<' '<<pages<<"\n<< /PageSize [";
    *psfile<<paper.real()*36e1/127<<' '<<paper.imag()*36e1/127<<"] >> setpagedevice\n";
    *psfile<<"gsave mmscale 0.1 setlinewidth\n";
    *psfile<<paper.real()/2<<' '<<paper.imag()/2<<' '<<"translate ";
    *psfile<<(pageorientation&3)*90<<" rotate ";
    *psfile<<paper.real()/-2<<' '<<paper.imag()/-2<<' '<<"translate"<<endl;
    *psfile<<"/Helvetica findfont 3 scalefont setfont"<<endl;
    oldr=oldg=oldb=NAN;
    inpage=true;
  }
}

void PostScript::endpage()
{
  if (psfile && indocument && inpage)
  {
    *psfile<<"grestore showpage"<<endl;
    inpage=false;
  }
}

void PostScript::trailer()
{
  if (inpage)
    endpage();
  if (psfile && indocument)
  {
    *psfile<<"%%BeginTrailer\n%%Pages: "<<pages<<"\n%%EndTrailer"<<endl;
    indocument=false;
  }
}

void PostScript::close()
{
  if (indocument)
    trailer();
  delete(psfile);
  psfile=nullptr;
}

int PostScript::getPages()
{
  return pages;
}

double PostScript::xscale(double x)
{
  return scale*(x-modelcenter.real())+paper.real()/2;
}

double PostScript::yscale(double y)
{
  return scale*(y-modelcenter.imag())+paper.imag()/2;
}

string PostScript::escape(string text)
{
  string ret;
  int ch;
  while (text.length())
  {
    ch=text[0];
    if (ch=='(' || ch==')')
      ret+='\\';
    ret+=ch;
    text.erase(0,1);
  }
  return ret;
}

void PostScript::setcolor(double r,double g,double b)
{
  if (r!=oldr || g!=oldg || b!=oldb)
  {
    *psfile<<fixed<<setprecision(3)<<r<<' '<<g<<' '<<b<<" col"<<endl;
    oldr=r;
    oldg=g;
    oldb=b;
  }
}

void PostScript::setscale(double minx,double miny,double maxx,double maxy,int ori)
/* To compute minx etc. using dirbound on e.g. a pointlist pl:
 * minx=pl.dirbound(-ori);
 * miny=pl.dirbound(DEG90-ori);
 * maxx=-pl.dirbound(DEG180-ori);
 * maxy=-pl.dirbound(DEG270-ori);
 */
{
  double xsize,ysize,papx,papy;
  int i;
  orientation=ori;
  modelcenter=complex<double>(minx+maxx,miny+maxy)/2.;
  xsize=fabs(minx-maxx);
  ysize=fabs(miny-maxy);
  papx=paper.real();
  papy=paper.imag();
  if (pageorientation&1)
    swap(papx,papy);
  for (scale=1;scale*xsize/10<papx && scale*ysize/10<papy;scale*=10);
  for (;scale*xsize/80>papx*0.9 || scale*ysize/80>papy*0.9;scale/=10);
  for (i=0;i<9 && (scale*xsize/rscales[i]>papx*0.9 || scale*ysize/rscales[i]>papy*0.9);i++);
  scale/=rscales[i];
  *psfile<<"% minx="<<minx<<" miny="<<miny<<" maxx="<<maxx<<" maxy="<<maxy<<" scale="<<scale<<endl;
}

double PostScript::getscale()
{
  return scale;
}

void PostScript::dot(complex<double> pnt,string comment)
{
  assert(psfile);
  pnt=turn(pnt,orientation);
  if (isfinite(pnt.real()) && isfinite(pnt.imag()))
  {
    *psfile<<ldecimal(xscale(pnt.real()),PAPERRES)<<' '<<ldecimal(yscale(pnt.imag()),PAPERRES)<<" .";
    if (comment.length())
      *psfile<<" %"<<comment;
    *psfile<<endl;
  }
}

void PostScript::circle(complex<double> pnt,double radius)
{
  assert(psfile);
  pnt=turn(pnt,orientation);
  if (isfinite(pnt.real()) && isfinite(pnt.imag()))
    *psfile<<ldecimal(xscale(pnt.real()),PAPERRES)<<' '<<ldecimal(yscale(pnt.imag()),PAPERRES)
    <<" n "<<ldecimal(scale*radius,PAPERRES)<<" 0 360 as %"
    <<ldecimal(radius*radius,radius*radius/1000)<<endl;
}

void PostScript::line2p(complex<double> pnt1,complex<double> pnt2)
{
  pnt1=turn(pnt1,orientation);
  pnt2=turn(pnt2,orientation);
  if (isfinite(pnt1.real()) && isfinite(pnt1.imag()) && isfinite(pnt2.real()) && isfinite(pnt2.imag()))
    *psfile<<ldecimal(xscale(pnt1.real()),PAPERRES)<<' '<<ldecimal(yscale(pnt1.imag()),PAPERRES)
    <<' '<<ldecimal(xscale(pnt2.real()),PAPERRES)<<' '<<ldecimal(yscale(pnt2.imag()),PAPERRES)<<" -"<<endl;
}

void PostScript::startline()
{
  assert(psfile);
  *psfile<<"n"<<endl;
}

void PostScript::lineto(complex<double> pnt)
{
  assert(psfile);
  pnt=turn(pnt,orientation);
  *psfile<<ldecimal(xscale(pnt.real()),PAPERRES)<<' '<<ldecimal(yscale(pnt.imag()),PAPERRES)<<(inlin?" l":" m");
  *psfile<<endl;
  inlin=true;
}

void PostScript::endline(bool closed)
{
  assert(psfile);
  if (closed)
    *psfile<<"closepath ";
  *psfile<<"s"<<endl;
  inlin=false;
}

void PostScript::widen(double factor)
{
  *psfile<<"currentlinewidth "<<ldecimal(factor)<<" mul setlinewidth"<<endl;
}

void PostScript::write(complex<double> pnt,string text)
{
  pnt=turn(pnt,orientation);
  *psfile<<ldecimal(xscale(pnt.real()),PAPERRES)<<' '<<ldecimal(yscale(pnt.imag()),PAPERRES)
  <<" m ("<<escape(text)<<") show"<<endl;
}

void PostScript::centerWrite(complex<double> pnt,string text)
{
  pnt=turn(pnt,orientation);
  *psfile<<ldecimal(xscale(pnt.real()),PAPERRES)<<' '<<ldecimal(yscale(pnt.imag()),PAPERRES)
  <<" m ("<<escape(text)<<") c."<<endl;
}

void PostScript::comment(string text)
{
  *psfile<<'%'<<text<<endl;
}
