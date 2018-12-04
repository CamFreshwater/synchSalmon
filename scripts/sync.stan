data {
  int<lower=1> N; // rows of data
  int<lower=1> G; // groups
  vector[N] y_time; // vector to hold observations of summed abundance each year
  row_vector[N] y_ind[G]; // vectors to hold observations
}
parameters {
  real<lower=0> denominator_mean[G];
  real<lower=0> denominator_sd[G];
  real<lower=0> numerator_mean;
  real<lower=0> numerator_sd;
}
model {
  y_time ~ normal(numerator_mean, numerator_sd);
  for (j in 1:G)
    y_ind[j] ~ normal(denominator_mean[j], denominator_sd[j]);

}
generated quantities {
  real phi_logit;
  real p;
  p = pow(numerator_sd, 2.0) / pow(sum(denominator_sd), 2.0);


  // vector[G] cv_s_ind;
  // vector[G] cv_s_ind_w;
  // real cv_s;
  // real cv_c;
  // for (i in 1:G) {
  //   cv_s_ind[i] = (ind_sigma[i] / (ind_mean[i]));
  //   cv_s_ind_w[i] = cv_s_ind[i] * (ind_mean[i] / group_mean);
  // }
  // cv_s = sum(cv_s_ind_w);
  // cv_c = group_sigma / group_mean;
}

