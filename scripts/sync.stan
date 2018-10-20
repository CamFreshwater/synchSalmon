data {
  int<lower=1> N; // rows of data
  int<lower=1> G; // groups
  vector[N] y_group; // vector to hold observations
  vector[N] y_time; // vector to hold observations of summed abundance each year
  row_vector[N] y_ind[G]; // vectors to hold observations
}
parameters {
  real<lower=0> group_sigma;
  real<lower=0> ind_sigma[G];
  real<lower=0> time_sigma[G];
  real<lower=0> group_mean;
  real<lower=0> ind_mean[G];
  real<lower=0> time_mean;
}
model {
  // very weak priors:
  // group_mean ~ normal(0, 100);
  // ind_mean ~ normal(0, 100);
  // group_sigma ~ normal(0, 100);
  // ind_sigma ~ normal(0, 100);

  // likelihood model:
  y_group ~ normal(group_mean, group_sigma);
  y_time ~ normal(time_mean, time_sigma);
  for (i in 1:G) {
    y_ind[i] ~ normal(ind_mean[i], ind_sigma[i]);
  }
}
generated quantities {
  real phi;
  vector[G] cv_s_ind;
  vector[G] cv_s_ind_w;
  real cv_s;
  real cv_c;
  phi = time_sigma^2 / sum(ind_sigma)^2;
  for (i in 1:G) {
    cv_s_ind[i] = (ind_sigma[i] / (ind_mean[i]));
    cv_s_ind_w[i] = cv_s_ind[i] * (ind_mean[i] / group_mean);
  }
  cv_s = sum(cv_s_ind_w);
  cv_c = group_sigma / group_mean;
}

