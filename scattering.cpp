 #include "scattering.hpp"
 using namespace std;
 scattering::scattering(int angularL,
                        double EcentreMass, double reducedMass)
 {
  l   = angularL;
  Ecm = EcentreMass;
  Mu  = reducedMass;
 }
/**
   r0 := minimum value of r to start integrating 
         the differential equation
   dr := step in disctrete position space  //fm

**/
 vector<double> scattering::scatteredWave(double r0,double dr)
 {
  vector<double> u_l;
  int numberPoints = 200;
  double rmax = numberPoints*dr;
  double wf_bc =0.; //boundary condition at r -> 0 fm

  // Solve Schroedinger Eq.
  diffEquation schroedinger(numberPoints,r0,wf_bc,dr,g);//f
  double u_1=0.01;
  vector<double> rn;
  schroedinger.numerovSolves(l,u_1,Ecm,rn,u_l);

  //Normalisation
  double delta_l,rinfty,normalise;
  double k0 =sqrt(2.*Mu*Ecm)/hbarc;
  delta_l = phaseShift(dr,u_l);
  rinfty = rn[numberPoints-3];
  normalise = sin(k0*rinfty - pi*l/2. + delta_l)/u_l[numberPoints-3];
 
  double r =r0;
  for(int i=0;i<numberPoints;i++) 
    {
      u_l[i]*=normalise/(k0*r);
      r+=dr;
    }
  return u_l;
 } 


/**
 Using Equation 98.21 in Davydov
 (Only good in a domain where the result is less or equal to 1 )

  k sin(\delta) = -2 \mu \int V(r) R_l(r) g_l(r) dr
                ~ -2 \mu \int V(r) j^2_l(kr) dr

****/
double scattering::phaseShift(double r0,double dr,int numberPoints)
 {
   double delta =0.0;
   double const k0 = sqrt(2.*Mu*Ecm)/hbarc;
   double r=r0;
   double j_l_kr,V_r;
   double sum =0.0;

   for(int i=0;i<numberPoints;i++)
      {
	j_l_kr = sph_bessel(l,k0*r);
	V_r =  wsPotential(r);
	sum   += V_r*j_l_kr*j_l_kr*r*r*dr;

	if( fabs(V_r) < 1.e-6 ) break;

	r+=dr;
      }
 
   sum *= -2.*Mu*k0/(hbc2*2.*pi);

   delta = asin(sum);
   return delta;
 }
/*************************************
      Spectral Function 
       in    1/(MeV*fm3)
 *************************************/
 MatrixXd scattering::SpectralFunction(double r0,double dr)
 {
   double const k0 = sqrt(2.*Mu*Ecm)/hbarc;
   double rho = 2.*Mu*k0/(hbarc*hbarc);
   vector<double> u_l = scatteredWave(r0,dr);
   int r_size = u_l.size();
   MatrixXd Srr(r_size,r_size);
   for(int ir=0;ir<r_size;ir++)
     {
      for(int jr=0;jr<r_size;jr++)
	Srr(ir,jr)=u_l[ir]*u_l[jr];
     }

   Srr*=rho;
   return Srr;
 }


  /****************************************

     u_l[numberPoints-5] = u_l1
     u_l[numberPoints-3] = u_l2
     u_l1 =A sinl1 + B cosl1
     u_l2 =A sinl2 + B cosl2

     denom =sinl1*cosl2 - sinl2*cosl1
     A = (u_l1*cosl2 - u_l2*cosl1)/denom
     B = (u_l2*sinl1 - u_l1*sinl2)/denom
    Phase shift satisfies:
      tan d_l = A/B 
    See for example
    F. Scheck, Quantum Physics,Eq. 2.13, p 141.

    Note: the Wave function u_lr does not have
          to be normalised since the resulting 
          phase shift is a ratio, any normali- 
          sation constant vanishes 

  ****************************************/

double scattering::phaseShift(double dr,vector<double> &u_lr)
{
  double d_l = 0.0;
  int const numberPoints = u_lr.size();
  double const k0 = sqrt(2.*Mu*Ecm)/hbarc;
  double rinf1,rinf2,sinl1,sinl2,cosl1,cosl2,u_l1,u_l2,denom;
  double Al,Bl;
  rinf1 = (numberPoints-5)*dr;
  rinf2 = (numberPoints-3)*dr;
  sinl1 = sin(k0*rinf1 - l*pi/2.);
  cosl1 = cos(k0*rinf1 - l*pi/2.);
  sinl2 = sin(k0*rinf2 - l*pi/2.);
  cosl2 = cos(k0*rinf2 - l*pi/2.);
 
   u_l1 = u_lr[numberPoints - 5];
   u_l2 = u_lr[numberPoints - 3];

   denom = sinl1*cosl2 - sinl2*cosl1;
   Al = (u_l1*cosl2 - u_l2*cosl1)/denom;
   Bl = (u_l2*sinl1 - u_l1*sinl2)/denom;
   d_l = atan(Al/Bl);
  return d_l;
}
/*********************************
  Pass scattering wave function
  to k-space
 *********************************/
  void  scattering::scattered_wf_in_k(vector<double> &k,
                                      vector<double> &dk,
                                      vector<double> &u_lk)
  {
    vector<double> rn,un_r;
    double r0 = 0.05;
    double dr = 0.05;
    un_r =scatteredWave(r0,dr);
    double ri=r0;
    /**
     for(int i=0;i<un_r.size();i++)
        {
         cout<<ri<<"\t"<<un_r[i]<<"\n";
         ri+=dr;
        }
    **/
     double utemp = 0.0;
     double j_l_kr = 0.0;
     kgrid kmesh(Ecm,Mu);
     kmesh.makeKmesh(k,dk);

     cout<<" k size = "<<k.size()<<"\n";
     /**/
     int ik=0;
     ri = r0;
     utemp = 0.0;

     for(int ir=0;ir<un_r.size();ir++)
       {
	 j_l_kr = sph_bessel(l,k[ik]*ri);
         utemp+= dr*ri*ri*ri*j_l_kr*un_r[ir];
         ri+=dr;
	 cout<<ri<<"\t"<<j_l_kr<<"\n";
       }
     u_lk.push_back(utemp);
     cout<<" Integral "<<utemp<<"\n";
     cout<<" u_lk size = "<<u_lk.size()<<"\n";
     /**/
  }
