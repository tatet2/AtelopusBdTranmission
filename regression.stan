data {
  int<lower=1> N;
  int<lower=1> Nfrog;
  int<lower=0, upper=1> genus[N];
  int<lower=0, upper=1> mixed[N];
  int<lower=0, upper=1> cons[N];    
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
  real alpha_m; // mixed effect
  real beta_m;
  real alpha_c; // conspecific effect
  real beta_c;
  real alpha_gc; // conspecific genus interaction
  real beta_gc;
  real alpha_gm; // mixed genus interaction
  real beta_gm;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[N] y_hat;

  for( i in 1:N){
    y_hat[i] = alpha_c*cons[i] + alpha_m*mixed[i] + alpha_g*genus[i] +
      alpha_gc*genus[i]*cons[i] + alpha_gm*genus[i]*mixed[i] + 
      beta_c*cons[i]*week[i] +  beta_m*mixed[i]*week[i] + beta_g*genus[i]*week[i] +
      beta_gc*genus[i]*cons[i]*week[i] + beta_gm*genus[i]*mixed[i]*week[i] +
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
  alpha_c ~ normal(0,5);
  beta_c ~ normal(0,5);
  alpha_m ~ normal(0,5);
  beta_m ~ normal(0,5);
  alpha_gc ~ normal(0,5);
  beta_gc ~ normal(0,5);
  alpha_gm ~ normal(0,5);
  beta_gm ~ normal(0,5);
  lzsp ~ normal(y_hat, sigma_y);
 }
