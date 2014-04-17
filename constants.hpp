/********************************************************
 *        constants.hpp
 * DECLARATION OF MATH CONSTANTS
 * AND SOME PHYSICS CONSTANTS
 *
 ********************************************************/
 #ifndef constants_hpp
 #define constants_hpp
 #include<cmath> 
 #include<vector>

 using namespace std;
/*********************************************
 *      Numerical CONSTANTS    
 *********************************************/
  double const pi=M_PI;
  double const e_=2.718281828459;
  double const ln_2=0.6931471805;

/*********************************************
 *        Physical Constants           
 *********************************************/
  double const M  = 938.95;          // Nucleon Mass [MeV]
  double const Mp = 938.27231;       // Proton Mass  [MeV]
  double const me=0.511;             // Electron Mass [MeV]
  double const m_mu=105.6583;        // Muon mass     [MeV]
  double const mu_tau=1776.84;       // Tau mass         [MeV]
  double const c_=2.99792e+10;      //speed of light    [cm/s]
  double const G_=6.6726e-8;           // Gravitation const [cm3/gs2]
  double const hbarc=197.327053;    //  (hbar*c)       [MeV*fm]
  double const alpha =0.0072973531062704997;
//double const alpha=0.091701236; //(e^2/(4pi*hbarc)) Structure fine constant
  double const e2=1.440292;             //e^2 [MeV*fm]
  double const m0=931.5;           //Other Nucleon Mass?
  double const kconstant=0.048192;   //  2/hbar^2
#endif
