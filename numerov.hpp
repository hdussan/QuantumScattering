/* 


In case the function in the right hand side
needs another argument, that is xx.
   Numerov algorithm needs
         y at r0    : y0
         y at r0+dr : y1
 */
 #ifndef numerov_hpp
 #define numerov_hpp
 #include "righthandside.hpp"
 #include <cmath>
 #include <iostream>
 #include <cstdlib>
 #include <vector>

 using namespace std;
 class diffEquation
 {
  private:
    int points;
    double r0;
    double y0;
    double dr;
    double (*function)(int l,double r);
    double (*Function)(int l,double xx,double r);
  public:
   diffEquation(int numberPoints,
                double rinitial,
                double wf_bc,double rstep,
                double (*f)(int,double));

   diffEquation(int numberPoints,
                double rinitial,
                double wf_bc,double rstep,
                double (*f)(int,double,double));

   void numerovSolves(int l,double y1,
                      vector<double> &rn,vector<double> &yn);

   void numerovSolves(int l,double y1,double xx,
                      vector<double> &rn,vector<double> &yn);
 };
 #endif
