/***
     righthandside
    right hand side of the differential equation
    to use numerov.
 ***/ 
#include "righthandside.hpp"
 using namespace std;

 double wellPotential(double r)
 {
   double const r0 = 6.0;
   return (r<r0)? -50.:0.0;
 }

 double f(int l,double Ecm,double r)
 {
   double const Mu = M; //reduced mass in MeV
   //Well potential
   double V = wellPotential(r);
   double centrifugal = l*(l+1.)/(r*r);
   double RHS = 2.*Mu*(Ecm -  V )/hbc2 - centrifugal; 
   return RHS;
 }

/*****************************
  Woods-Saxon Potential
 *****************************/
 double wsPotential(double r)
 {
   double const r0=1.21;//fm
   double const a = 0.5;
   int const A  = 40; //Calcium
   double const V0 = -50.0;
  
   double R = r0*pow(double(A),1./3.);
   return V0/(1. + exp((r-R)/a) );
 }

 double g(int l,double Ecm,double r)
 {
   double const Mu = M; //reduced mass in MeV
   //Well potential
   double V = wsPotential(r);
   double centrifugal = l*(l+1.)/(r*r);
   double RHS = 2.*Mu*( Ecm -  V )/hbc2 - centrifugal; 
   return RHS;
 }
 
