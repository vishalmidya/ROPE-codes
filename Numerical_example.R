# Numerical Example Comparing ROPE vs. Bayes Factor

delta = 1/2
n.simu = 10000000


T_saved = c(50, 200, 800, 3200)

liao = function(y_bar, n, TT)
{
 
  sigma = 1/sqrt(n)
  
  
  x1 = seq(delta, TT, length=n.simu)
  x2 = seq(-TT, -delta, length=n.simu)
  x = c(x1, x2)
  top1 = mean(dnorm(x, y_bar, sigma))
  
  x = seq(-delta, delta, length=n.simu)
  bottom1 = mean(dnorm(x, y_bar, sigma))
  
  Bayes_factor = top1/bottom1
  posterior_odds = Bayes_factor*(TT-delta)/delta
  posterior_prob =  posterior_odds/(1+ posterior_odds)
  
  c(TT, posterior_prob,  posterior_odds, Bayes_factor)
  
}


n = 50
y_bar = 1

for(TT in T_saved)
{
  print(liao(y_bar, n, TT))
}


posterior_prob_saved = numeric(1000)
posterior_odds_saved = numeric(1000)
Bayes_factor_saved = numeric(1000)


#####################

n = 100
TT = 200
for(i in 1:1000)
{
  print(i)
  y_bar = rnorm(1,delta, 1/sqrt(n))
  result1 = liao(y_bar, n, TT)
  posterior_prob_saved[i] = result1[2]
  posterior_odds_saved[i] = result1[3]
  Bayes_factor_saved[i] = result1[4]
}


