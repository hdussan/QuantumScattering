#ifndef scattering_hpp
#define scattering_hpp
#include "constants.hpp"
#include "numerov.hpp"
#include "kMesh.hpp"
#include <boost/math/special_functions.hpp>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <cmath>
#include <vector>
using namespace std;
using namespace boost::math;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;

class scattering
{
 private:
    int l;
    double Ecm; 
    double Mu;
 public: 
    scattering(int angularL,double EcentreMass, double reducedMass);
    vector<double> scatteredWave(double r0,double dr);
    double phaseShift(double r0,double dr,int numberPoints);
    MatrixXd SpectralFunction(double r0,double dr);
    double phaseShift(double dr,vector<double> &u_lr);
    
    void scattered_wf_in_k(vector<double> &k,
                           vector<double> &dk,
                           vector<double> &u_lk);
};

#endif
