functions {
  real skew_student_t_log(real y, real nu, real mu, real sigma, real skew) {
  real lp;
  if (skew <= 0)
    reject("Skew has to be positive. Found skew=", skew);
  if (sigma <= 0)
    reject("Scale has to be positive.  Found sigma=", sigma);
  lp = log(skew) - log1p(square(skew));
  if (y < mu)
    return lp + student_t_lpdf(y * skew | nu, mu * skew, sigma);
  else
    return lp + student_t_lpdf(y / skew | nu, mu / skew, sigma);
  }
}
data {
  int<lower=0> N; // rows of data
  int<lower=1> J; // n fixed covariates on mean
  matrix[N,J] X_ij; // covariate model matrix
  vector[N] y_i; // vector to hold observations
}
parameters {
  vector[J] b_j;
  real<lower=0> sigma;
  real log_skew;
  real<lower=2> nu;
}
model {
  vector[N] mu_i;
  nu ~ exponential(0.1);
  b_j ~ normal(0, 5);
  sigma ~ student_t(3, 0, 3);
  log_skew ~ student_t(3, 0, 3);
  mu_i = X_ij * b_j;
  for (i in 1:N)
    y_i[i] ~ skew_student_t(nu, mu_i[i], sigma, exp(log_skew));
}
