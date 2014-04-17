 #ifndef righthandside_hpp
 #define righthandside_hpp
 #include "constants.hpp"
 using namespace std;

 double const hbc2 = hbarc*hbarc;
// well potential and right hand side (for Numerov) 
// for well potential
 double wellPotential(double r);
 double f(int l,double Ecm,double r);

// Woods-Saxon potential and right hand side (for Numerov)
// for this potential
 double wsPotential(double r);
 double g(int l,double Ecm,double r);
 #endif
