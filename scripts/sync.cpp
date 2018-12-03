#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(y_time);
  DATA_MATRIX(y_ind);
  
  int N = y_time.size(); // get number of data points to loop over
  int G = y_ind.cols(); // get number of columns to loop over

  // parameters:
  PARAMETER(log_numerator_sigma);
  PARAMETER_VECTOR(log_denominator_sigma);
  PARAMETER(numerator_mean);
  PARAMETER_VECTOR(denominator_mean);

  Type numerator_sigma = exp(log_numerator_sigma);
  vector<Type> denominator_sigma = exp(log_denominator_sigma);

  Type nll = 0.0; // initialize negative log likelihood
  for(int i = 0; i < N; i++)  {
    nll -= dnorm(y_time(i), numerator_mean, numerator_sigma, true);
  }
    for(int i = 0; i < G; i++) {
      for(int j = 0; j < N; j++) {
        nll -= dnorm(y_ind(j,i), denominator_mean(i), denominator_sigma(i), true);
      }
    }
    
 // Derived values:
 
 Type phi;
 vector<Type> cv_s_ind;
 vector<Type> cv_s_ind_w;
 Type cv_s;
 Type cv_c;
 Type log_cv_s;
 Type log_cv_c;
 
 phi = pow(numerator_sigma, 2.) / pow(sum(denominator_sigma), 2.);
//
// for(int i = 0; i < G; i++)  {
 cv_s_ind = denominator_sigma / denominator_mean;
 cv_s_ind_w = cv_s_ind * (denominator_mean / numerator_mean);
// }
 cv_s = sum(cv_s_ind_w);
 cv_c = numerator_sigma / numerator_mean;

 log_cv_s = log(cv_s);
 log_cv_c = log(cv_c);
 
 Type logit_phi = logit(phi);
 REPORT(phi);
 REPORT(logit_phi);
 REPORT(log_cv_s);
 REPORT(log_cv_c);
 
 ADREPORT(phi);
 ADREPORT(logit_phi);
 ADREPORT(log_cv_s);
 ADREPORT(log_cv_c);

  return nll;
}
