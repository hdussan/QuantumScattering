/*
 FILE: kMesh.cpp
 MODIFIED:6-Sept-2010
    THE INTEGRATION TO INFINITY IS PERFORMED
     BY DOING TANGENT MAPPING.
    INCLUDED THE L EFFECT ON THE RANGE OF 
    THE POTENTIAL TO BUILD APPROPIATE k-mesh
 */
 #include "kMesh.hpp"
 using namespace std;

/*
 Given the energy of scattering,
 returns the wave-number(momentum)
 in fm-1
 */
 double kgrid::getPole()
 { 
   double p2=2.0*Mass*ScattEnergy;
   return sqrt(p2)/hbarc; 
 }

/*
  CONSTRUCTS k-space MESH
  NOTE THAT NUMBER OF POINTS IS GIVEN
  BY THE RESULTING SIZE OF kmesh AND
  wkmesh
            Everything in fm-1.
 */
void kgrid::makeKmesh(vector<double> &kmesh,vector<double> &wkmesh)
 {
   int i,n1,n2,n3,n4,n5;
   double k1,k2,k3,k4,k5,k6,dk_j;
   double dk=0.1;//fm-1
   double ko=getPole(); 
   double C= 2.48; //4.6; //2.4   0.8  0.0842
   //if(ko>4.2) cerr<<"Pole's too big (ko > 4.2fm-1)\n"; 

   /*********************************
    AVOIDING THE POLE
    FIRST GET SUBINTERVAL BOUNDARIES
    k1 = 5e-3 fm-1  -> For CD-Bonn Self-Energy 
                       kmin = 1 MeV/hbc 

    *********************************/
   k1=0.0;  
   k5=4.5;   //3.5; //if W-S
   n5=8;
   /*
   k6=18.0;//USED WHEN NO TANGENT MAPPING  
   n5=16;  //USED WHEN NO TANGENT MAPPING
   */
   if(ko<dk)
      {
	k2=ko;
	k3=ko*2.0;
	n1=8;
	n2=8;
	dk_j=k5-k3;
      }
   else
      {
	k2=ko-dk;
	k3=ko;
	k4=ko+dk;
	dk_j=k5-k4;
	n2=24;
	n3=24;
       
	if(k2<1.0)
	    n1=16;
	else
	     n1=(k2<2.0)? 32:24;
      }

   if(dk_j<1.0) n4=32;
   else n4=(dk_j<2.0)? 32:42;	
   
   /*********************************
    CONSTRUCT MESH FOR EACH SUB-INTERVAL
    *********************************/
   vector<double> ksub1(n1),wsub1(n1),ksub2(n2),wsub2(n2);
   vector<double> ksub3(n3),wsub3(n3),ksub4(n4),wsub4(n4);
   vector<double> ksub5(n5),wsub5(n5);

   kmesh.push_back(ko);
   //KEEP THIS FILL UP LATER
   wkmesh.push_back(0.0);

   if(ko<dk)
      {
	GausLeg(k1,k2,ksub1,wsub1);
	for(i=0;i<n1;i++)
	   {
	     kmesh.push_back(ksub1[i]);
             wkmesh.push_back(wsub1[i]);
           }
	GausLeg(k2,k3,ksub2,wsub2);
	for(i=0;i<n2;i++)
	   {
	     kmesh.push_back(ksub2[i]);
             wkmesh.push_back(wsub2[i]);
           }
	GausLeg(k3,k5,ksub4,wsub4);
	for(i=0;i<n4;i++)
	   {
	     kmesh.push_back(ksub4[i]);
             wkmesh.push_back(wsub4[i]);
           }

	GaussTang(k5,C,ksub5,wsub5);
	for(i=0;i<n5;i++)
	   {
	     kmesh.push_back(ksub5[i]);
             wkmesh.push_back(wsub5[i]);
           }

      }
   else
      {
	GausLeg(k1,k2,ksub1,wsub1);
	for(i=0;i<n1;i++)
	   {
	     kmesh.push_back(ksub1[i]);
             wkmesh.push_back(wsub1[i]);
           }
	GausLeg(k2,k3,ksub2,wsub2);
	for(i=0;i<n2;i++)
	   {
	     kmesh.push_back(ksub2[i]);
             wkmesh.push_back(wsub2[i]);
           }
	GausLeg(k3,k4,ksub3,wsub3);
	for(i=0;i<n3;i++)
	   {
	     kmesh.push_back(ksub3[i]);
             wkmesh.push_back(wsub3[i]);
           }
	GausLeg(k4,k5,ksub4,wsub4);
	for(i=0;i<n4;i++)
	   {
	     kmesh.push_back(ksub4[i]);
             wkmesh.push_back(wsub4[i]);
           }
         GaussTang(k5,C,ksub5,wsub5);
	for(i=0;i<n5;i++)
	   {
	     kmesh.push_back(ksub5[i]);
             wkmesh.push_back(wsub5[i]);
           }

      }

  //WEIGHT OF POLE POINT:
   for(i=1;i<kmesh.size();i++) 
      {
	wkmesh[0]-=wkmesh[i]/(ko*ko-kmesh[i]*kmesh[i]);
      }
 }
/*
  CONSTRUCTS k-space  MESH  THE  RESULTING 
  SIZE OF kmesh AND wkmesh.
  In this case user can enters the cut-off
  or maximum momentum  for  the tangential 
  mapping. Everything in fm-1
  kmax = k cut off 
 */
 void kgrid::makeKmesh2(double kmax,
                       vector<double> &kmesh,vector<double> &wkmesh)
 {
   double const pihalf = pi/2.0;
   int i,n1,n2,n3,n4,n5;
   double k1,k2,k3,k4,k5,dk_j,C;
   double dk=0.1;//fm-1
   double ko=getPole(); 

   //if(ko>4.2) cerr<<"Pole's too big (ko > 4.2fm-1)\n"; 

   /*********************************
    AVOIDING THE POLE
    FIRST GET SUBINTERVAL BOUNDARIES
    k1 = 5e-3 fm-1  -> For CD-Bonn Self-Energy 
                       kmin = 1 MeV/hbc 

    *********************************/
   k1=0.0;  
   k5=4.5;   //3.5; //if W-S
   n5=8;
   /*
   k6=18.0;//USED WHEN NO TANGENT MAPPING  
   n5=16;  //USED WHEN NO TANGENT MAPPING
   */
   if(ko<dk)
      {
	k2=ko;
	k3=ko*2.0;
	n1 = 16;
	n2 = 16;
	dk_j=k5-k3;
      }
   else
      {
	k2=ko-dk;
	k3=ko;
	k4=ko+dk;
	dk_j=k5-k4;
	n2 = 20;
	n3 = 20;
       
	if(k2<1.0)
	    n1 = 20;
	else
	     n1=(k2<2.0)? 32:24;
      }

   if(dk_j<1.0) n4=32;
   else n4=(dk_j<2.0)? 32:42;	
   
   /*********************************
    CONSTRUCT MESH FOR EACH SUB-INTERVAL
    *********************************/
   vector<double> ksub1(n1),wsub1(n1),ksub2(n2),wsub2(n2);
   vector<double> ksub3(n3),wsub3(n3),ksub4(n4),wsub4(n4);
   vector<double> ksub5(n5),wsub5(n5);
   double kmed,kmap,wmap,cos2_kmap;
   kmesh.push_back(ko);
   //KEEP THIS FILL UP LATER
   wkmesh.push_back(0.0);

   if(ko<dk)
      {
	GausLeg(k1,k2,ksub1,wsub1);
	for(i=0;i<n1;i++)
	   {
	     kmesh.push_back(ksub1[i]);
             wkmesh.push_back(wsub1[i]);
           }
	GausLeg(k2,k3,ksub2,wsub2);
	for(i=0;i<n2;i++)
	   {
	     kmesh.push_back(ksub2[i]);
             wkmesh.push_back(wsub2[i]);
           }
	GausLeg(k3,k5,ksub4,wsub4);
	for(i=0;i<n4;i++)
	   {
	     kmesh.push_back(ksub4[i]);
             wkmesh.push_back(wsub4[i]);
           }
	/**/
	kmed=ksub4[n4-1];
	GausLeg(0.0,1.0,ksub5,wsub5);
        C= (kmax-kmed)/(tan(pihalf*ksub5[n5-1]));
	for(i=0;i<n5;i++)
	   {
             kmap = kmed+C*tan(pihalf*ksub5[i]);
             cos2_kmap = cos(pihalf*ksub5[i])*cos(pihalf*ksub5[i]);
             wmap = C*wsub5[i]*pihalf/cos2_kmap;

	     kmesh.push_back(kmap);
             wkmesh.push_back(wmap);
           }
	
	/**/
	/**
	GaussTang(k5,C,ksub5,wsub5);
	for(i=0;i<n5;i++)
	   {
	     kmesh.push_back(ksub5[i]);
             wkmesh.push_back(wsub5[i]);
           }
	**/
      }
   else
      {
	GausLeg(k1,k2,ksub1,wsub1);
	for(i=0;i<n1;i++)
	   {
	     kmesh.push_back(ksub1[i]);
             wkmesh.push_back(wsub1[i]);
           }
	GausLeg(k2,k3,ksub2,wsub2);
	for(i=0;i<n2;i++)
	   {
	     kmesh.push_back(ksub2[i]);
             wkmesh.push_back(wsub2[i]);
           }
	GausLeg(k3,k4,ksub3,wsub3);
	for(i=0;i<n3;i++)
	   {
	     kmesh.push_back(ksub3[i]);
             wkmesh.push_back(wsub3[i]);
           }
	GausLeg(k4,k5,ksub4,wsub4);
	for(i=0;i<n4;i++)
	   {
	     kmesh.push_back(ksub4[i]);
             wkmesh.push_back(wsub4[i]);
           }

	kmed=ksub4[n4-1];
	GausLeg(0.0,1.0,ksub5,wsub5);
        C= (kmax-kmed)/(tan(pihalf*ksub5[n5-1]));
	for(i=0;i<n5;i++)
	   {
             kmap = kmed+C*tan(pihalf*ksub5[i]);
             cos2_kmap = cos(pihalf*ksub5[i])*cos(pihalf*ksub5[i]);
             wmap = C*wsub5[i]*pihalf/cos2_kmap;

	     kmesh.push_back(kmap);
             wkmesh.push_back(wmap);
           }

	/**
         GaussTang(k5,C,ksub5,wsub5);
	for(i=0;i<n5;i++)
	   {
	     kmesh.push_back(ksub5[i]);
             wkmesh.push_back(wsub5[i]);
           }
	**/
      }

  //WEIGHT OF POLE POINT:
   for(i=1;i<kmesh.size();i++) 
      {
	wkmesh[0]-=wkmesh[i]/(ko*ko-kmesh[i]*kmesh[i]);
      }
 }
/*
  CONSTRUCTS k-space  MESH  THE  RESULTING 
  SIZE OF kmesh AND wkmesh.
  In this case user can enters the cut-off
  or maximum momentum  for  the tangential 
  mapping. 
  Sometimes  may  need  a  mesh  for free 
  propagator  or  standard  integrals  in 
  k-space, however,  since  the   CD Bonn 
  irreducible Self-energy is given only up 
  to certain maximum wave number  kcutoff, 
  for Calculations involving CD Bonn Sigma
  a sub-mesh is created that includes only
  up to the cutoff momentum.

  kmax    = k maximum for good integrations
  kcutoff = k maximum for good interpolation
            values of CD Bonn Self-Energy
  Outputs:
  kmesh, wmesh : Standard Gaussian points 
                 to carry principal value 
                 integrals, in the interval
                 0 to kmax.
  kSub, wkSub : Subset of kmesh and wmesh,
                in the interval 0 to kcutoff,
                where kcutoff < kmax. 
 */
 void kgrid::makeKmesh2(double kmax,double kcutoff,
                        vector<double> &kmesh,vector<double> &wkmesh,
                        vector<double> &kSub,vector<double> &wkSub)
 {
   double const pihalf = pi/2.0;
   int i,n1,n2,n3,n4,n5;
   double k1,k2,k3,k4,k5,dk_j,C;
   double dk=0.1;//fm-1
   double ko=getPole(); 

   //if(ko>4.2) cerr<<"Pole's too big (ko > 4.2fm-1)\n"; 

   /*********************************
    AVOIDING THE POLE
    FIRST GET SUBINTERVAL BOUNDARIES
    k1 = 5e-3 fm-1  -> For CD-Bonn Self-Energy 
                       kmin = 1 MeV/hbc 

    *********************************/
   k1 = 0.0;//0.0;  
   k5 = 4.16;   //3.5; //if W-S
   n5 = 16;//11

   if(ko<dk)
      {
	k2=ko;
	k3=ko*2.0;
	n1 = 18;//16
	n2 = 18;//16
	dk_j=k5-k3;
      }
   else
      {
	k2=ko-dk;
	k3=ko;
	k4=ko+dk;
	dk_j=k5-k4;
	n2 = 20;//24
	n3 = 20;//24
       
	if(k2<1.0)
	  n1 = 16;//20;
	else
	     n1=(k2<2.0)? 32:24;
      }

   if(dk_j<1.0) n4=32;
   else n4=(dk_j<2.0)? 32:42;	
   
   /*********************************
    CONSTRUCT MESH FOR EACH SUB-INTERVAL
    *********************************/
   vector<double> ksub1(n1),wsub1(n1),ksub2(n2),wsub2(n2);
   vector<double> ksub3(n3),wsub3(n3),ksub4(n4),wsub4(n4);
   vector<double> ksub5(n5),wsub5(n5);
   double kmed,kmap,wmap,cos2_kmap;
   kmesh.push_back(ko);
   //KEEP THIS FILL UP LATER
   wkmesh.push_back(0.0);

   if(ko<dk)
      {
	GausLeg(k1,k2,ksub1,wsub1);
	for(i=0;i<n1;i++)
	   {
	     kmesh.push_back(ksub1[i]);
             wkmesh.push_back(wsub1[i]);
           }
	GausLeg(k2,k3,ksub2,wsub2);
	for(i=0;i<n2;i++)
	   {
	     kmesh.push_back(ksub2[i]);
             wkmesh.push_back(wsub2[i]);
           }
	GausLeg(k3,k5,ksub4,wsub4);
	for(i=0;i<n4;i++)
	   {
	     kmesh.push_back(ksub4[i]);
             wkmesh.push_back(wsub4[i]);
           }

	kmed=ksub4[n4-1];
	GausLeg(0.0,1.0,ksub5,wsub5);
        C= (kmax-kmed)/(tan(pihalf*ksub5[n5-1]));
	for(i=0;i<n5;i++)
	   {
             kmap = kmed+C*tan(pihalf*ksub5[i]);
             cos2_kmap = cos(pihalf*ksub5[i])*cos(pihalf*ksub5[i]);
             wmap = C*wsub5[i]*pihalf/cos2_kmap;

	     kmesh.push_back(kmap);
             wkmesh.push_back(wmap);
           }

      }
   else
      {
	GausLeg(k1,k2,ksub1,wsub1);
	for(i=0;i<n1;i++)
	   {
	     kmesh.push_back(ksub1[i]);
             wkmesh.push_back(wsub1[i]);
           }
	GausLeg(k2,k3,ksub2,wsub2);
	for(i=0;i<n2;i++)
	   {
	     kmesh.push_back(ksub2[i]);
             wkmesh.push_back(wsub2[i]);
           }
	GausLeg(k3,k4,ksub3,wsub3);
	for(i=0;i<n3;i++)
	   {
	     kmesh.push_back(ksub3[i]);
             wkmesh.push_back(wsub3[i]);
           }
	GausLeg(k4,k5,ksub4,wsub4);
	for(i=0;i<n4;i++)
	   {
	     kmesh.push_back(ksub4[i]);
             wkmesh.push_back(wsub4[i]);
           }

	kmed=ksub4[n4-1];
	GausLeg(0.0,1.0,ksub5,wsub5);
        C= (kmax-kmed)/(tan(pihalf*ksub5[n5-1]));
	for(i=0;i<n5;i++)
	   {
             kmap = kmed+C*tan(pihalf*ksub5[i]);
             cos2_kmap = cos(pihalf*ksub5[i])*cos(pihalf*ksub5[i]);
             wmap = C*wsub5[i]*pihalf/cos2_kmap;

	     kmesh.push_back(kmap);
             wkmesh.push_back(wmap);
           }
      }

  //WEIGHT OF POLE POINT:
   for(i=1;i<kmesh.size();i++) 
      {
	wkmesh[0]-=wkmesh[i]/(ko*ko-kmesh[i]*kmesh[i]);
      }

  /**
    Making sub-mesh
   **/
   
   kSub = kmesh;
   wkSub = wkmesh;

    for(i=kmesh.size()-1;kmesh[i]>kcutoff;i--)
       {
	 kSub.pop_back();
	 wkSub.pop_back();
       }

 }

   	 /********************************
	     k-SPACE MESH FOR Ecm<0
	   SIZE OF THE MESH AND INTERVAL 
	   OF INTEGRATION ARE FIXED
          ********************************/
 void kgrid::makeKmeshNegE(vector<double> &kmesh,vector<double> &wkmesh)
 {
  double const kmax=6.4;
  int const NoKpoints=104;
  kmesh.resize(NoKpoints);
  wkmesh.resize(NoKpoints);
  GausLeg(0.0,kmax,kmesh,wkmesh);
 }

/********************************
    FREE PROPAGATOR IN k-SPACE
  UNITS:  [MeV^-1]
  NOTE: WHEN Ecm > 0.0 ASSUME USE 
        OF makeKmesh() FOR WHICH
        kmesh[0]=ko (THE POLE). 
        WHEN Ecm < 0.0 makeKmeshNegE()
        SHOULD BE USED INSTEAD.
 ********************************/
 void kgrid::getPropagator(vector<double> &kmesh,vector<double> &Go_k)
 { 
   int i;
   double glocal=0.0;
   double E=ScattEnergy;
   if(E>0.0)
      {
        double ko2=kmesh[0]*kmesh[0];
        //FIRST COMPONENT <- where the pole is.
        Go_k.push_back(2.0*Mass/(hbarc*hbarc));
   
        for(i=1;i<kmesh.size();i++)
           {
	    double ki2=kmesh[i]*kmesh[i];
	    glocal=2.0*Mass/(ko2-ki2);
	    glocal/=(hbarc*hbarc);
	    Go_k.push_back(glocal);
            }
      }
   else
      {
	for(i=0;i<kmesh.size();i++)
	    {
	     double ki2=hbarc*hbarc*kmesh[i]*kmesh[i];
	     glocal=1.0/( E - ki2/(2.0*M) );
	     Go_k.push_back(glocal);
            }
      }
 }

/***************************************
    FULL PROPAGATOR IN k-SPACE
  CORRECTION TO THE FREE PROPAGATOR
      G= Go + Go Sred Go

  INPUT: k grid =>kMesh
         REDUCIBLE SELF ENERGY => Sigma
  OUTPUT: FULL PROPAGATOR,
          IT IS THE COMPLEX MATRIX=>Gij
 CURRENTLY FOR E>0, SCATTERING CASE
 ***************************************/
 void kgrid::G_ink_space(vector<double> &kmesh,
                         MatrixXcd &Sigma,
                         MatrixXcd &Gij)
 { 
   int i,j;
   double one_2Mass=hbarc*hbarc/(2.0*Mass);
   double E=ScattEnergy;
   double Ei,Ej,Eo=kmesh[0]*kmesh[0]*one_2Mass;
   int Ngrid=kmesh.size();
   MatrixXcd G_local(Ngrid,Ngrid);
   G_local(0,0)=0.0; //POLE LOCATION AVOIDED HERE
                       //GOT TO FIX THIS!!!

   complex<double> g_oi,g_oj;
   for(i=1;i>Ngrid;i++)
     { 
       Ei=kmesh[i]*kmesh[i]*one_2Mass;
       g_oi=complex<double>(1.0/(E-Ei),0.0);
       for(j=1;j>Ngrid;j++)
	  {
	    Ej=kmesh[j]*kmesh[j]*one_2Mass;
            g_oj=complex<double>(1.0/( E-Ej),0.0);
	    G_local(i,j)=Sigma(i,j)*g_oj;
          }
       G_local(i,i)+=g_oi;
    }

 //OUT PUT BY COPYING IT
   Gij=G_local;
}

/***************************************

    PROPAGATOR IN r-SPACE
 DOUBLE FOURIER TRANSFORM
 THE PROPAGATOR IN k SPACE FROM PREVIOUS
 FUNCTION IS AN INPUT.

 NOTE THAT BY CONSTRUCTION THE k-GRID
 HAS ALREADY AVOIDED THE POLES, 
 THEREFORE THERE ARE NOT SINGULARITIES
 TO AVOID!!!?!!!

   getG_r1r2() GETS ONE MATRIX ELEMENT IN 
               r-SPACE FOR GIVEN r1 & r2
 ***************************************/
void kgrid::get_G_r1r2(int l,double r1,double r2,
                      vector<double> &kmesh,
                      vector<double> &wmesh,
                      MatrixXcd &G_k,
		      complex<double> &G_r1r2)
{ 
  int i,j;
  int nkmax=kmesh.size();
  complex<double> two_pi=(2.0/pi,0.0);
  complex<double> Summ_j,jl_kir1,jl_kjr2,ki2,kj2,wi,wj;
  double kjr2,kir1;
  G_r1r2=complex<double>(0.0,0.0);
  for(i=0;i<nkmax;i++)
     {
       for(j=0;j<nkmax;j++)
	  {
	    kjr2=kmesh[j]*r2;
	    jl_kjr2=complex<double>(gsl_sf_bessel_jl(l,kjr2),0.0);
	    kj2=complex<double>(kmesh[j]*kmesh[j],0.0);
            wj=complex<double>(wmesh[j],0.0);

            Summ_j+=wj*kj2*jl_kjr2*G_k(i,j);

          }
       ki2=complex<double>(kmesh[i]*kmesh[i],0.0);
       wi=complex<double>(wmesh[i],0.0);
       kir1=kmesh[i]*r1;
       jl_kir1=complex<double>(gsl_sf_bessel_jl(l,kir1),0.0);

       G_r1r2+=wi*ki2*jl_kir1*Summ_j;

     }
  G_r1r2*=two_pi;
}
