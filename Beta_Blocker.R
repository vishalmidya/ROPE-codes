library("rstan")

stanString = "

data{
int<lower=0> N;
int<lower=0> n0[N];
int<lower=0> y0[N];
int<lower=0> n1[N];
int<lower=0> y1[N];
real<lower=0> delta;
real<lower=0> tau;
real eta;
}

parameters{
real logit_p0[N];
real <lower=0> sigma2;
real theta;
real theta_i[N];
}

model {

 real a;
 real sigma;
 real p0[N];
 real p1[N];


 if(fabs(theta) < delta) a = uniform_lpdf(theta|-delta, delta) + log(1.0 - eta);
   else a = normal_lpdf(theta|0, tau) - log(2.0*normal_cdf(-delta, 0, tau)) + log(eta);

  target += a;

  sigma2 ~ inv_gamma_lpdf(0.01,0.01);
  sigma = sqrt(sigma2);


for(i in 1:N) {

theta_i[i] ~ normal( theta, sigma);

logit_p0[i] ~ normal(0, 10);
p0[i] = exp(logit_p0[i]);
p0[i] = p0[i]/(1.0+p0[i]);

p1[i] = logit_p0[i] + theta_i[i];//logit_p1
p1[i] = exp(p1[i]);
p1[i] = p1[i]/(1.0+p1[i]);

}
 y0 ~ binomial(n0, p0);
 y1 ~ binomial(n1, p1);

}
"

d.truncated.normal = function(y, delta, tau)
{
  if(abs(y) < delta) 0
  else dnorm(y, 0, tau)/(2*pnorm(-delta, 0, tau))
}

beta_blocker_example = function(delta, tau){

  eta = dunif(delta,-delta,delta)/(dunif(delta,-delta,delta) + d.truncated.normal(delta, delta, tau))
  cat("tau, eta =", tau, eta, "\n")

  blockerDat<-list(y1 = c(3,  7,  5,  102,  28, 4,  98,  60, 25, 138, 64, 45,  9, 57, 25, 33, 28, 8, 6, 32, 27, 22 ),
                 n1 = c(38, 114, 69, 1533, 355, 59, 945, 632, 278,1916, 873, 263, 291, 858, 154, 207, 251, 151, 174, 209, 391, 680),
                 y0 = c(3, 14, 11, 127, 27, 6, 152, 48, 37, 188, 52, 47, 16, 45, 31, 38, 12, 6, 3, 40, 43, 39),
                 n0 = c(39, 116, 93, 1520, 365, 52, 939, 471, 282, 1921, 583, 266, 293, 883, 147, 213, 122, 154, 134, 218, 364, 674),
                 N = 22, delta = delta, tau=tau)

  stan.model = stan(model_code = stanString, data = blockerDat, iter = 5e5, pars="theta", warmup = 500, chains = 4, thin = 1, refresh = 0)

  xx = extract(stan.model, pars="theta")
  xx = xx[[1]]

  prob1 = mean(abs(xx)> blockerDat$delta)
  Bayes_factor = ((1-eta)/eta)* (prob1/(1-prob1))
  Posterior_odds = prob1/(1-prob1)
  cat("delta,eta,Posterior_Prob,log_Posterior_odds,log_Bayes_factor",
              delta,eta,prob1,log(Posterior_odds),log(Bayes_factor),"\n")
  c(delta,eta,prob1,log(Posterior_odds),log(Bayes_factor))

}

delta.saved = seq(0.1,0.4,0.05)
result_beta = array(0, c(length(delta.saved), 5))
for(i in 1:length(delta.saved)) {
  result_beta[i,] = beta_blocker_example(delta.saved[i],2.5*delta.saved[i])
}

result_beta = as.data.frame(result_beta)
names(result_beta) = c("delta","eta","Posterior_Prob","log_Posterior_odds","log_Bayes_factor")


