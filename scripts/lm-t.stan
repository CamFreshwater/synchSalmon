data {
  int<lower=0> N; // rows of data
  int<lower=1> J; // n fixed covariates on mean
  matrix[N,J] X_ij; // covariate model matrix
  vector[N] y_i; // vector to hold observations
}
parameters {
  vector[J] b_j;
  real<lower=0> sigma;
  real<lower=2> nu;
}
model {
  nu ~ exponential(0.1);
  b_j ~ normal(0, 5);
  sigma ~ student_t(3, 0, 3);
  y_i ~ student_t(nu, X_ij * b_j, sigma);
}
