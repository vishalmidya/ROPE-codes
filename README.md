# ROPE-codes
## Codes for comparisons between log Posterior Odds vs. log Bayes Factor in Interval Null Hypotheses

These codes are used in my paper "Connecting and contrasting Bayes Factor and Region of Practical Equivalence Procedure for testing interval null hypothesis" to illustrate and highlight 

* The stability of log Posterior Odds with respect to log Bayes Factors in interval null hypothesis testing.
* Demonstration of a simple and effective algorithm for computing Bayes factor using draws from posterior distributions generated by standard Bayesian programs such as Stan.

## Abstract

There has been strong recent interest in testing interval null  hypothesis  for improved scientific inference. For example, Lakens et al (2018) and Lakens and Harms (2017) use this approach to study if there is a pre-specified meaningful treatment effect in gerontology and clinical trials, instead of  a point null hypothesis of any effect. Two popular Bayesian approaches are available for interval null hypothesis testing. One is the standard Bayes factor and the other is the Region of Practical Equivalence (ROPE) procedure championed by Kruschke and others over many years. This paper connects key quantities in the two approaches, which in turn allow us to contrast two key differences between the approaches with huge practical implications. The first is that Bayes factor depends heavily on the prior specification while ROPE procedure is very robust. The second difference is concerned with the statistical property when data is generated under a neutral parameter value on the common boundary of competing hypotheses. Bayes factor can be severely biased while ROPE approach gives reasonable result. Finally, the connection leads to a simple and effective algorithm for computing Bayes factor using draws from posterior distributions generated by standard Bayesian programs such as BUGS, JAGS and Stan.  

