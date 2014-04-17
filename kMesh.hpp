/*
 FILE :   kMesh.hpp
 PURPOSE: 
     -CREATES A GRID OF GAUSS
      POINTS FOR GAUSS INTEGRATIONS
      AVOIDING POLE FOR ScattEnergy>0.
     -ONCE THE REDUCIBLE SELF ENERGY 
      IS FOUND,FINDS THE CORRECTED
      PROPAGATOR IN k-SPACE.
     -ADDING CONSTRUCTION OF ANGULAR 
      MESH FOR CALCUALTION IN A VECTOR
      BASIS REPRESENTATION.
 */
 #ifndef kMesh_hpp
 #define kMesh_hpp
 #include <cmath>
 #include <vector>
 #include <iostream>
 #include <complex>
 #include <gsl/gsl_sf_bessel.h>
 #include <Eigen/LU>
 #include <Eigen/Dense>
 #include "constants.hpp"
 #include "Gauss_q.hpp"


 using namespace std;
 using Eigen::MatrixXd;
 using Eigen::MatrixXcd;

 class kgrid
 {
  private:
    double ScattEnergy;
    double Mass;
  public:
      //DEFAULT CONSTRUCTOR NUCLEON MASS
    kgrid(double Ecm) {ScattEnergy=Ecm; Mass=M;}

    kgrid(double Ecm,double mu){ScattEnergy=Ecm; Mass=mu;}
    double getPole();
    void makeKmesh(vector<double> &kmesh,vector<double> &wkmesh);
    void makeKmesh2(double kmax, 
                    vector<double> &kmesh,vector<double> &wkmesh);
    void makeKmesh2(double kmax, double kcutoff,
                    vector<double> &kmesh,vector<double> &wkmesh,
                    vector<double> &kSub,vector<double> &wkSub);

    void makeKmeshNegE(vector<double> &kmesh,vector<double> &wkmesh);

    void cosTheta(vector<double> &Xk,vector<double> &WXk);
    void makeThetaMesh(vector<double> &thetak,vector<double> &wtheta);

    void getPropagator(vector<double> &kmesh,vector<double> &Go_k);
    void G_ink_space(vector<double> &kmesh,
		     MatrixXcd &Sigma,
                     MatrixXcd &Gij);
   
    void get_G_r1r2(int l,double r1,double r2,
		    vector<double> &kmesh,
		    vector<double> &wmesh,
		    MatrixXcd &G_k,
                    complex<double> &G_r1r2);
 }; 
 #endif
