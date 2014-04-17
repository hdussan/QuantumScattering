/***************************************************************
 *                      Gauss_q.hpp                              *
 *HEADER FILE TO SOLVE INTEGRALS USING GAUSSIAN QUADRATURE     *
 *GIVEN THE LIMITS OF THE INTEGRAL X1 AND X2, THE ROUTINE      *
 *CALCULATES THE ABSCISAS X_J AND THE WEIGHTS W_J TO DO THE    * 
 *NEEDED INTEGRAL (n-POINTS QUADRATURE FORMULA).               * 
 *n IS THE SIZE OF THE VECTOR X.                               * 
 ***************************************************************/
#ifndef Gauss_q_hpp
#define Gauss_q_hpp
#include "constants.hpp"
using namespace std;
/***************************************************************
 *             GAUSSIAN QUADRATURE                             *
 ***************************************************************/
template<class T>
void GausLeg(const T x1,const T x2,vector<T> &x,vector<T> &w)
{const double epsilon=1.0e-15;
 int m,j,i;
 T z1,z,xm,xl,pp,p3,p2,p1;
 int n=x.size();
 m=(n+1)/2;
 xm=0.5*(x2+x1);
 xl=0.5*(x2-x1);
 for(i=0;i<m;i++)                     /*************************/
   {z=cos(pi*(i+0.75)/(n+0.5));       /*Approx. to the ith root*/
   do{p1=1.0;                         /*************************/
      p2=0.0;
      for(j=0;j<n;j++)                /*************************/
	{p3=p2;                       /*Legendre Rec. Rel.     */ 
    	 p2=p1;                       /*************************/
	 p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
        }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
      }
   while((fabs(z-z1)) >epsilon);
   x[i]=xm-xl*z;
   x[n-1-i]=xm+xl*z;
   w[i]=2.0*xl/((1.0-z*z)*pp*pp);
   w[n-1-i]=w[i];
   }
}

/************************************
 *      GAUSSIAN QUADRATURE         *
 * WITH TANGENT MAPPING TO INTEGRATE*
 * FROM  'A'   TO INFINITY          *
 ************************************/
 template <class T>
 void GaussTang(const T A,const T C,vector<T> &X,vector<T> &W)
 {
  const T pi_2=0.5*pi;  
  int n=X.size();
  vector <T> x(n),w(n);
  GausLeg(0.0,1.0,x,w);
  for(int i=0; i<n; i++)
     {
       T cos2x=cos(pi_2*x[i])*cos(pi_2*x[i]);
       X[i]=A+C*tan(pi_2*x[i]);
       W[i]=pi_2*C*w[i]/cos2x;
     }
 }

/*
Generation of tangential mapping with specified intervals:
    a:     lower limit
    infty: Good estimation of infinity (depends on the 
           function to integrate)
 */
 template <class T>
 void GaussMapping(T a,T infty,vector<T> &X,vector<T> &W)
 {
   const T pi_2 = pi/2.0;
   int n = X.size();
   vector<T> x(n),w(n);
   T cinf,cos_pix2;
   GausLeg(0.0,1.0,x,w);
   cinf = (infty-a)/tan( pi_2*x[n-1] );
   for(int i=0; i<n; i++)
      {
	cos_pix2 = cos(pi_2*x[i]);
        X[i] =  cinf*tan( pi_2*x[i] ) + a;
	W[i] =  cinf*w[i]*pi_2/( cos_pix2*cos_pix2 );
      }
 }
#endif
