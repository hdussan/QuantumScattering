 #include "readBound_wf.hpp"
 using namespace std;
bound_wf_data::bound_wf_data(int radialQuantum_n, 
                             int orbitalAngularl,
                             double totalAngularj,
                             int numberOfLines)
{
  n_quantum =  radialQuantum_n;
  l         =  orbitalAngularl;
  j         =    totalAngularj;
  lines     =    numberOfLines;
}

 string bound_wf_data::namefile()
 {
   string name;
   string nprinc, angularl, angularj;
   string prefix ="wf_";
   string sufix = "half.dat";
   ostringstream nconv1,nconv2;

   nconv1 << n_quantum;
   nprinc = nconv1.str();

   switch(l)
     {
       case 1:
          angularl = "p";
          break;
       case 2:
          angularl = "d";
          break;
       case 3:
          angularl = "f";
          break;
       case 4:
          angularl = "g";
          break;
      default:
          angularl = "s";
     }

   nconv2<< int(2.*j);
   angularj = nconv2.str();
   name = prefix + nprinc + angularl + angularj + sufix;   

   return name;
 }

vector<double> bound_wf_data::wf_in_rspace(int rpoints,
                                           double r0, 
                                           double dr)
 {
   vector<double> wf_r;
   /*
     1. read file
     2. loop over r wf in r-space
       For each r need to carry an integral 
       over the momentum
    */
   double tempdk,tempk,tempwf,tempwfn,ignore1,ignore2;
   vector<double> dk,k,wf_lj;
   string wf_name =namefile();
   ifstream wf_file;
   cout<<"\t  reading "<<wf_name<<"\n";
   wf_file.open(wf_name.c_str(),ios::in);
   for(int i=0;i<lines;i++)
      {
        wf_file>>tempdk>>tempk>>tempwfn>>ignore1>>tempwf>>ignore2;
	dk.push_back(tempdk);
         k.push_back(tempk);
	 wf_lj.push_back(tempwfn);
      }
   wf_file.close();

   /* Transforming to r-space*/
   double const fact =sqrt(2./pi);
   int ir,ik;
   double r,dk_k2,j_l_kr,sum;

   r = r0;
   for(ir=0;ir<rpoints;ir++)
     {
      sum = 0.0;
      for(ik =0; ik<lines;ik++)
	 { 
	   dk_k2  = dk[ik]*k[ik]*k[ik];
	   j_l_kr = sph_bessel(l,k[ik]*r);
	   sum   += dk_k2*j_l_kr*wf_lj[ik];
	 }
     
      sum*=fact;
      wf_r.push_back(sum);
      r += dr;
     }

   return wf_r;
 }
