 #include "numerov.hpp"
 using namespace std;
/* */
diffEquation::diffEquation(int numberPoints,
                           double rinitial,double wf_bc,double step,
			   double (*f)(int,double))
 {
   points = numberPoints;
   r0     = rinitial;
   y0     = wf_bc;
   dr     = step;
   function = f;
 }

/*************************
   Constructor for right hand side function that
   needs 3 arguments
 *************************/
diffEquation::diffEquation(int numberPoints,
                           double rinitial,double wf_bc,double step,
			   double (*f)(int,double,double))
 {
   points = numberPoints;
   r0     = rinitial;
   y0     = wf_bc;
   dr     = step;
   Function = f;
 }

/* 
   Numerov algorithm needs
         y at r0    : y0
         y at r0+dr : y1
 */
void diffEquation::numerovSolves(int l,
                                 double y1,
                                 vector<double> &rn,
                                 vector<double> &yn)
 {
   int i;
   double ri,yi,h2;
   h2 = dr*dr;
   ri = r0;
   rn.push_back(r0);
   yn.push_back(y0);
   yn.push_back(y1);
   for(i=1;i<points;i++)
     { 
       ri+=dr;
       yi=(2.-5.*h2*function(l,ri)/6.)*yn[i] 
     	      - ( 1.+ h2*function(l,rn[i-1])/12. )*yn[i-1];
       yi/=(1. + h2*function(l,ri+dr)/12.);
       rn.push_back(ri);
       yn.push_back(yi);
     }
 }

/* 
In case the function in the right hand side
needs another argument, that is xx.
   Numerov algorithm needs
         y at r0    : y0
         y at r0+dr : y1
 */
void diffEquation::numerovSolves(int l,
                                 double y1,
				 double xx,
                                 vector<double> &rn,
                                 vector<double> &yn)
 {
   int i;
   double ri,yi,h2;
   h2 = dr*dr;
   ri = r0;
   rn.push_back(r0);
   yn.push_back(y0);
   yn.push_back(y1);
   for(i=1;i<points;i++)
     { 
       ri+=dr;
       yi=(2.-5.*h2*Function(l,xx,ri)/6.)*yn[i] 
	   - ( 1.+ h2*Function(l,xx,rn[i-1])/12. )*yn[i-1];
       yi/=(1. + h2*Function(l,xx,ri+dr)/12.);
       rn.push_back(ri);
       yn.push_back(yi);
     }
 }
