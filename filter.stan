data {
  int<lower=1> N;
  int<lower=1> Nfrog;
  int<lower=0, upper=1> genus[N];
  int<lower=1, upper=Nfrog> id[N];
  int week[N]; //days since exposure
  real<lower=0> lzsp[N]; // nat log of zsp count
}
parameters {
  vector[Nfrog] alpha;
  vector[Nfrog] beta;
  real mu_alpha;
  real mu_beta;
  real alpha_g; // genus effect
  real beta_g;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[N] y_hat;

  for( i in 1:N){
    y_hat[i] = alpha_g*genus[i] + beta_g*genus[i]*week[i] +
      alpha[id[i]] + beta[id[i]]*week[i];   
      } 
}
model {
  #priors:
  mu_alpha ~ normal(0,5);
  mu_beta ~ normal(0,5);
  alpha ~ normal(mu_alpha, sigma_alpha);
  beta ~ normal(mu_beta, sigma_beta);
  alpha_g ~ normal(0,5);
  beta_g ~ normal(0,5);
  lzsp ~ normal(y_hat, sigma_y);
 }
