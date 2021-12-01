/******************************************************/
/*                                                    */
/* ps.h - PostScript output                           */
/*                                                    */
/******************************************************/
#ifndef PS_H
#define PS_H
#include <string>
#include <iostream>
#include <map>
#include <complex>

struct papersize
{
  int width,height; // in micrometers
};
extern std::map<std::string,papersize> papersizes;

class PostScript
{
protected:
  std::ostream *psfile;
  int pages;
  bool indocument,inpage,inlin;
  double scale; // paper size is in millimeters, but model space is in meters
  int orientation,pageorientation;
  double oldr,oldg,oldb;
  std::complex<double> paper,modelcenter;
public:
  PostScript();
  ~PostScript();
  void setpaper(papersize pap,int ori);
  double aspectRatio();
  void open(std::string psfname);
  bool isOpen();
  void prolog();
  void startpage();
  void endpage();
  void trailer();
  void close();
  int getPages();
  double xscale(double x);
  double yscale(double y);
  std::string escape(std::string text);
  void setcolor(double r,double g,double b);
  void setscale(double minx,double miny,double maxx,double maxy,int ori=0);
  double getscale();
  void dot(std::complex<double> pnt,std::string comment="");
  void circle(std::complex<double> pnt,double radius);
  void line2p(std::complex<double> pnt1,std::complex<double> pnt2);
  void startline();
  void lineto(std::complex<double> pnt);
  void endline(bool closed=false);
  void widen(double factor);
  void write(std::complex<double> pnt,std::string text);
  void centerWrite(std::complex<double> pnt,std::string text);
  void comment(std::string text);
};

#endif
