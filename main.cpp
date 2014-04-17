#include "DifferentialEq.hpp"
using namespace std;

int main()
{
  int l,nr;
  double Ecm;
  
  cin>>Ecm;

  l=2;

  vector<double> rn,un;
  scattering woodssaxon(l,Ecm,M);
  double r0=0.1;
  double dr=0.1;
 
  vector<double> ur;
  ur =  woodssaxon.scatteredWave(r0,dr);
  double r=r0;
  int rpoints =ur.size();
  cout<<"\t    E_cm = "<<Ecm<<" MeV \t";
  cout<<"   l = "<<l<<"\n";
  cout<<"delta_l =\t";
  cout<<woodssaxon.phaseShift(dr,ur)<<"\n";

  MatrixXd Srr = woodssaxon.SpectralFunction(r0,dr);
  r=r0;
  for(int i=0;i<rpoints;i++)
     {
       cout<<r<<"\t"<<Srr(i,i)<<"\n";
       r+=dr;
      }
  /**********************************
    1. Use  ur := Scattered Wave Function
    2. Now calculate the bound wave function
       (solve eigen-state eq.)

    3. Product with Bound State
       Orthogonality tested
   
  **********************************/
  /**
  nr = 0;
  double SpectralDensity=0.0;
  vector<double> ubound;
  double j =l+0.5;
  bound_wf_data boundwf(nr,l,j,64);
  ubound = boundwf.wf_in_rspace(rpoints,r0,dr);
  r=r0;
  for(int i=0;i<rpoints;i++)
    {
      cout<<r<<"\t"<<ubound[i]<<"\t"<<ur[i]<<"\n";
      SpectralDensity +=dr*r*r*ubound[i]*ur[i]; 
      r+=dr;
    }
  cout<<"  Spectral Density ("<<Ecm<<" MeV ) = ";
  cout<<SpectralDensity*SpectralDensity<<"\n";
  **/
}
