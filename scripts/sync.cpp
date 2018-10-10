#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(y_group);
  DATA_MATRIX(y_ind);
  
  int N = y_group.size(); // get number of data points to loop over
  int G = y_ind.cols(); // get number of columns to loop over

  // parameters:
  PARAMETER(log_group_sigma);
  PARAMETER_VECTOR(log_ind_sigma);
  PARAMETER(group_mean);
  PARAMETER_VECTOR(ind_mean); 

  Type group_sigma = exp(log_group_sigma);
  vector<Type> ind_sigma = exp(log_ind_sigma);

  Type nll = 0.0; // initialize negative log likelihood
  for(int i = 0; i < G; i++)  {
    nll -= dnorm(y_group(i), group_mean, group_sigma, true);
  }
    for(int i = 0; i < G; i++) {
      for(int j = 0; j < N; j++) {
        nll -= dnorm(y_ind(j,i), ind_mean(i), ind_sigma(i), true);
      }
    }
    
 // Derived values:
 
 Type phi;
 vector<Type> cv_s_ind;
 vector<Type> cv_s_ind_w;
 Type cv_s;
 Type cv_c;
 
 phi = pow(group_sigma,2) / pow(sum(ind_sigma),2);
 
// for(int i = 0; i < G; i++)  {
   cv_s_ind = ind_sigma / ind_mean;
   cv_s_ind_w = cv_s_ind * (ind_mean / group_mean);
// }
   cv_s = sum(cv_s_ind_w);
   cv_c = group_sigma / group_mean;
 
 Type logit_phi = logit(phi);
 REPORT(phi);
 REPORT(group_sigma);
 REPORT(ind_sigma);
 REPORT(logit_phi);
 REPORT(cv_s);
 REPORT(cv_c);
 
 ADREPORT(phi);
 ADREPORT(logit_phi);
 ADREPORT(cv_s);
 ADREPORT(cv_c);

  return nll;
}
