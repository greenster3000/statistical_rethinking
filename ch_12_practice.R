# 12H1. 
# In 2014, a paper was published that was entitled “Female hurricanes are deadlier than male
# hurricanes. As the title suggests, the paper claimed that hurricanes with female names have caused
# greater loss of life, and the explanation given is that people unconsciously rate female hurricanes
# as less dangerous and so are less likely to evacuate. Statisticians severely criticized the paper after
# publication. Here, you’ll explore the complete data used in the paper and consider the hypothesis
# that hurricanes with female names are deadlier.
# Acquaint yourself with the columns by inspecting the help ?Hurricanes. In this problem, you’ll focus on 
# predicting deaths using femininity of each hurricane’s name. Fit and interpret the simplest
# possible model, a Poisson model of deaths using femininity as a predictor. You can use quap or
# ulam. Compare the model to an intercept-only Poisson model of deaths. How strong is the association 
# between femininity of name and deaths? Which storms does the model fit (retrodict) well?
# Which storms does it fit poorly?
library(rethinking)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
d <- Hurricanes
dat_list <- list(fem=d$femininity, deaths=d$deaths)
p12.1 <- ulam(
alist(
deaths ~ dpois(lambda),
log(lambda) <- a + b*fem,
a ~ dnorm(3, 0.5),
b ~ dnorm(0, 0.2)
), data=dat_list, chains=4, log_lik = TRUE, iter = 2000
)
precis(p12.1)
precis_plot(precis(p12.1))

f_seq = seq(from = 1 , to = 11, length.out = 20)
lambda <- link(p12.1, data=data.frame(fem = f_seq))
plot(d$femininity , d$deaths, xlab="Femininity", ylab="Deaths", col=rangi2)
lmu <- apply(lambda , 2, mean)
lci <- apply(lambda, 2, PI, prob=0.97)
lines(f_seq, lmu, lty=2, lwd=1.5)
shade(lci , f_seq)
postcheck(p12.1, window = 46)
# Looks like a relatively weak association between feminity and deaths
# Doesn't do well at predicting Hurricanes with very high number of deaths, overdispersion probably
# Does better at predicting lower deaths, but not much

p12.1.int <- ulam(
alist(
deaths ~ dpois(lambda),
log(lambda) <- a,
a ~ dnorm(3, 0.5)
), data=list(deaths=d$deaths), chains=4, log_lik = TRUE
)
compare(p12.1, p12.1.int)
lambda.int <- link(p12.1.int, data=data.frame(fem = f_seq))
lmu.int <- apply(lambda.int , 2, mean)
lci.int <- apply(lambda.int, 2, PI, prob=0.97)
lines(f_seq, lmu.int, lty=1, lwd=1.5)
shade(lci.int , f_seq)
# gives all the weight to the model with femininity, but SE is a lot bigger than dWAIC.
# very few extra deaths predicted even for the most feminine of hurricanes

# 12H2. Counts are nearly always over-dispersed relative to Poisson. So fit a gamma-Poisson (aka 
# negative-binomial) model to predict deaths using femininity. Show that the over-dispersed model
# no longer shows as precise a positive association between femininity and deaths, with an 89% interval
# that overlaps zero. Can you explain why the association diminished in strength?

p12.2 <- ulam(
alist(
deaths ~ dgampois( lambda, phi ),
log(lambda) <- a + b*fem,
a ~ dnorm(3, 0.5),
b ~ dnorm(0, 0.2),
phi ~ dexp(1)
), data=dat_list, chains=4, log_lik = TRUE
)
precis(p12.2)
lambda.2 <- link(p12.2, data=data.frame(fem = f_seq))
lmu.2 <- apply(lambda.2 , 2, mean)
lci.2 <- apply(lambda.2, 2, PI, prob=0.97)
lines(f_seq, lmu.2, lty=3, lwd=1.5)
shade(lci.2 , f_seq)
# now b overlaps 0
# gamma poisson allows more variance in the lambda parameter across rows
# Can see this in the fact that the predictions are less certain with the 97% PI


# 12H3. In the data, there are two measures of a hurricane’s potential to cause death: damage_norm
# and min_pressure. Consult ?Hurricanes for their meanings. It makes some sense to imagine that
# femininity of a name matters more when the hurricane is itself deadly. This implies an interaction
# between femininity and either or both of damage_norm and min_pressure. Fit a series of models
# evaluating these interactions. Interpret and compare the models. In interpreting the estimates, it
# may help to generate counterfactual predictions contrasting hurricanes with masculine and feminine
# names. Are the effect sizes plausible?

# So we want femininity * damage_norm model, as well as femininity * min_pressure and 
# fem * d_norm * min_pres




