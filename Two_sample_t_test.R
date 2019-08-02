library("rstan")

stanString = "  

functions
{ 
   real b1_lpdf(real b1, real delta, real tau, real eta)
  {
   if(fabs(b1) <= delta) return uniform_lpdf(b1|-delta, delta) + log(1.0 - eta);
   else return normal_lpdf(b1|0, tau) - log(2.0*normal_cdf(-delta, 0, tau)) + log(eta); //truncated normal
  }
}

 
data {
 int n;
 vector[n] y;
 vector[n] x;
 real delta;
 real <lower=0> tau;
 real eta;
}

 parameters {
 real b0;
 real b1;
 real <lower=0> sigma2;
 }

 
 model {
 vector[n] mu;
 
  b0 ~ normal(0,100);

  target += b1_lpdf(b1| delta, tau, eta);

   sigma2 ~ inv_gamma(0.01,0.01);
 
   mu = b0 + sqrt(sigma2)*b1*x; //Cohen's effect size
 
   y ~ normal(mu, sqrt(sigma2));
 }
 
 "

d.truncated.normal = function(y, delta, tau)
{
  if(abs(y) < delta) 0
  else dnorm(y, 0, tau)/(2*pnorm(-delta, 0, tau))
}

two_sample_t_interval_hypothesis = function(y1, y2, delta, tau)
{
  y = c(y1, y2)
  x1 = rep(0, length(y1))
  x2 = rep(1, length(y2))
  x = c(x1, x2)
  
  eta = dunif(delta,-delta,delta)/(dunif(delta,-delta,delta) + d.truncated.normal(delta, delta, tau))
  
  cat("tau, eta =", tau, eta, "\n")
  
  data.list <- list(y = y, x = x, n = length(y), delta = delta, tau=tau, eta=eta)
  
  library(rstan)
  stan.model = stan(model_code = stanString, data = data.list, iter = 5e5, pars="b1",
                    warmup = 500, chains = 4, thin = 1, refresh = 0)
  
  xx = extract(stan.model, pars="b1")
  xx = xx[[1]]
  
  prob1 = mean(abs(xx)> delta)
  Bayes_factor = ((1-eta)/eta)* (prob1/(1-prob1))
  Posterior_odds = prob1/(1-prob1)
  cat("tau, Bayes_factor","Posterior_odds",tau, Bayes_factor, Posterior_odds,"\n")
  c(tau, Bayes_factor,Posterior_odds)
}

delta = 0.1
calcium=c(7,-4,18,17,-3,-5,1,10,11,-2)
placebo=c(-1,12,-1,-3,3,-5,5,2,-11,-1,-3)

tau.saved = seq(10,50,2)
result1 = array(0, c(length(tau.saved), 3))
for(i in 1:length(tau.saved)) result1[i,] = two_sample_t_interval_hypothesis(calcium, placebo, delta, tau.saved[i]*delta)

result1 = as.data.frame(result1)
names(result1) = c("tau", "Bayes_factor","Posterior_odds")

print(t.test(calcium, placebo))

