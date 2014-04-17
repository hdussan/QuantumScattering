 #ifndef readBound_wf_hpp
 #define readBound_wf_hpp
 #include "constants.hpp"
 #include <iostream>
 #include <fstream>
 #include <cstring>
 #include <vector>
 #include <cmath>
 #include <boost/math/special_functions.hpp>
 using namespace std;
 using namespace boost::math;
 class bound_wf_data
 {
  private:
     int n_quantum;
     int l;
     double j;
     int lines;
  public:
     bound_wf_data(int radialQuantum_n, int orbitalAngularl,
                   double totalAngularj, int numberOfLines);
     string namefile();
     vector<double> wf_in_rspace(int rpoints,double r0,double dr);
 };
 #endif
