library(rethinking)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# 13E1. Which of the following priors will produce more shrinkage in the 
# estimates? 
#   (a) α tank ∼ Normal(0, 1); 
#   (b) α tank ∼ Normal(0, 2).

# (a) will do, as it's got a smaller sigma, hence allows less variation in each
# tanks individual intercept parameters

# 13E2. Rewrite the following model as a multilevel model.
# y i ∼ Binomial(1, p i )
# logit(p i ) = α group[i] + βx i
# α group ∼ Normal(0, 1.5)
# β ∼ Normal(0, 0.5)

y ~ dbinom(1 , p)
logit(p) <- a[group] + b*x
b ~ dnorm(0, 1.5)
a[group] ~ dnorm(a_bar, sigma_a)
a_bar ~ dnorm(0, 1.5)
sigma_a ~ dexp(1)


# 13E3. Rewrite the following model as a multilevel model.
# y i ∼ Normal(μ i , σ)
# μ i = α group[i] + βx i
# α group ∼ Normal(0, 5)
# β ∼ Normal(0, 1)
# σ ∼ Exponential(1)

y ~ dnorm(mu, sigma)
mu <- a[group] + b*x
a[group] ~ dnorm(a_bar, sigma_a)
b ~ dnorm(0,1)
a_bar ~ dnorm(0, 1.5)
sigma_a ~ dexp(1)
sigma ~ dexp(1)

# 13E4. Write a mathematical model formula for a Poisson regression with varying 
# intercepts.

y ~ dpois(lambda)
log(lambda)  <- a[group] + b*x
a[group] ~ dnorm(a_bar, sigma_a)
b ~ dnorm(0, 0.5)
a_bar ~ dnorm(0, 1.5)
sigma_a ~ dexp(1)


# 13E5. Write a mathematical model formula for a Poisson regression with two 
# different kinds of varying intercepts, a cross-classified model.

y ~ dpois(lambda)
log(lambda) <- a[group] + b[block] + c*x
a[group] ~ dnorm(a_bar, sigma_a)
b[block] ~ dnorm(0, sigma_b)
c ~ dnorm(0, 0.5)
a_bar ~ dnorm(0, 1.5)
sigma_a ~ dexp(1)
sigma_b ~ dexp(1)

# 13M1. Revisit the Reed frog survival data, data(reedfrogs) , and add the
# predation & size treatment variables to the varying intercepts model. Consider
# models with either main effect alone, both main effects, as well as a model 
# including both and their interaction. Instead of focusing on inferences about 
# these two predictor variables, focus on the inferred variation across tanks. 
# Explain why it changes as it does across models.

data("reedfrogs")
d <- reedfrogs
str(d)

d$tank <- 1:nrow(d)

dat <- list(
  S = as.integer(d$surv),
  pred = as.integer(d$pred),
  big=as.integer(d$size),
  N=as.integer(d$density),
  tank=as.integer(d$tank)
)
m13m1.size <- ulam(
  alist(
    S ~ dbinom( N , p),
    logit(p) <- a[tank] + b_s*big,
    b_s ~ dnorm(0, 0.5),
    a[tank] ~ dnorm(a_bar, sigma_a),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1)
  ), data=dat, chains=4, log_lik=TRUE, iter=2000
)

m13m1.pred <- ulam(
  alist(
    S ~ dbinom( N , p),
    logit(p) <- a[tank] + b_p*pred,
    b_p ~ dnorm(0, 0.5),
    a[tank] ~ dnorm(a_bar, sigma_a),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1)
  ), data=dat, chains=4, log_lik=TRUE, iter=4000
)

m13m1.both <- ulam(
  alist(
    S ~ dbinom( N , p),
    logit(p) <- a[tank] + b_p*pred + b_s*big,
    c(b_p, b_s) ~ dnorm(0, 0.5),
    a[tank] ~ dnorm(a_bar, sigma_a),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1)
  ), data=dat, chains=4, log_lik=TRUE, iter=4000
)

m13m1.interact <- ulam(
  alist(
    S ~ dbinom( N , p),
    logit(p) <- a[tank] + b_p*pred + b_s*big + b_ps * big * pred,
    c(b_p, b_s, b_ps) ~ dnorm(0, 0.5),
    a[tank] ~ dnorm(a_bar, sigma_a),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1)
  ), data=dat, chains=4, log_lik=TRUE, iter=4000
)

# basically all have equally good out of sample predictive power
precis_plot(precis(m13m1.pred))
precis_plot(precis(m13m1.size))
precis_plot(precis(m13m1.both))
precis_plot(precis(m13m1.interact))

# Looks like it's not too certain about whether size has an effect on its own
# Looks like when we include the interaction term, predator and size has 
# an effect together

# let's plot sigma_a and a_bar for all the models

post.pred <- extract.samples(m13m1.pred)
post.size <- extract.samples(m13m1.size)
post.both <- extract.samples(m13m1.both)
post.interact <- extract.samples(m13m1.interact)

dens(post.pred$sigma_a, lty=1, xlim=c(0,2.5))
dens(post.size$sigma_a, add = TRUE, lty=2)
dens(post.both$sigma_a, add = TRUE, lty=3)
dens(post.interact$sigma_a, add = TRUE, lty=4)
title("Sigma_A distributions for each model")
legend(x=0, y = 2.5, legend=c("Predator", "Size", "Both", "Interaction"), 
       lty=c(1,2,3,4))

dens(post.pred$a_bar, lty=1, xlim=c(-1,6))
dens(post.size$a_bar, add = TRUE, lty=2)
dens(post.both$a_bar, add = TRUE, lty=3)
dens(post.interact$a_bar, add = TRUE, lty=4)
title("A_bar distributions for each model")
legend(x=-1, y = 0.8, legend=c("Predator", "Size", "Both", "Interaction"), 
       lty=c(1,2,3,4))

# Can clearly see that sigma a for the size model is larger, would expect less
# shrinkage here

# Also where they are shrinking too should change. Would expect Size only model
# to be pulled only a bit more survival-y, the rest a lot more, with predator
# only model being pulled the most of the three

# See if my prediction is correct! We should see model m13m1.size give less
# shrinkage

d$propsurv.est.pred <- inv_logit(apply(post.pred$a, 2, mean))
d$propsurv.est.size <- inv_logit(apply(post.size$a, 2, mean))
d$propsurv.est.both <- inv_logit(apply(post.both$a, 2, mean))
d$propsurv.est.interact <- inv_logit(apply(post.interact$a, 2, mean))

# Predator only prediction
plot(d$propsurv, ylim=c(0,1), pch=16, xaxt="n", xlab="tank", 
     ylab="proportion survived", col=rangi2)
axis(1, at = c(1,16,32,48), labels=c(1,16,32,48))
points(d$propsurv.est.pred)
abline(h=mean(inv_logit(post.pred$a_bar)), lty=2)
abline(v=16.5, lwd=2)
abline(v=32.5, lwd=2)
text(8, 0, "small tanks")
text(16+8, 0, "medium tanks")
text(32+8, 0, "large tanks")
title("Only predator as a predictor")

# Size only prediction
plot(d$propsurv, ylim=c(0,1), pch=16, xaxt="n", xlab="tank", 
     ylab="proportion survived", col=rangi2)
axis(1, at = c(1,16,32,48), labels=c(1,16,32,48))
points(d$propsurv.est.size)
abline(h=mean(inv_logit(post.size$a_bar)), lty=2)
abline(v=16.5, lwd=2)
abline(v=32.5, lwd=2)
text(8, 0, "small tanks")
text(16+8, 0, "medium tanks")
text(32+8, 0, "large tanks")
title("Only size as a predictor")

# Both prediction
plot(d$propsurv, ylim=c(0,1), pch=16, xaxt="n", xlab="tank", 
     ylab="proportion survived", col=rangi2)
axis(1, at = c(1,16,32,48), labels=c(1,16,32,48))
points(d$propsurv.est.both)
abline(h=mean(inv_logit(post.both$a_bar)), lty=2)
abline(v=16.5, lwd=2)
abline(v=32.5, lwd=2)
text(8, 0, "small tanks")
text(16+8, 0, "medium tanks")
text(32+8, 0, "large tanks")
title("Size and predator as predictors")

# Interaction prediction
plot(d$propsurv, ylim=c(0,1), pch=16, xaxt="n", xlab="tank", 
     ylab="proportion survived", col=rangi2)
axis(1, at = c(1,16,32,48), labels=c(1,16,32,48))
points(d$propsurv.est.interact)
abline(h=mean(inv_logit(post.interact$a_bar)), lty=2)
abline(v=16.5, lwd=2)
abline(v=32.5, lwd=2)
text(8, 0, "small tanks")
text(16+8, 0, "medium tanks")
text(32+8, 0, "large tanks")
title("Size and predator as predictors with interaction term")

# Can really see the shirnkage differences due to difference in a_bar and sigma
# In size only, lots of points barely shrunk, even some "grown" in the large 
# tanks

# 13M2. Compare the models you fit just above, using WAIC. Can you reconcile 
# the differences in WAIC with the posterior distributions of the models?

compare(m13m1.pred, m13m1.size, m13m1.both, m13m1.interact)

# err, basically no difference? the dWAIC is less than the dSE. 

# 13M3. Re-estimate the basic Reed frog varying intercept model, but now using a
# Cauchy distribution in place of the Gaussian distribution for the varying 
# intercepts. That is, fit this model:

s_i ~ Binomial(n_i , p_i)
logit(p_i) = α_tank[i]
α_tank ~ Cauchy(α, σ)
α ~ Normal(0, 1)
σ ~ Exponential(1)

# (You are likely to see many divergent transitions for this model. Can you 
# figure out why? Can you fix them?) Compare the posterior means of the 
# intercepts, α tank , to the posterior means produced in the chapter, using the 
# customary Gaussian prior. Can you explain the pattern of differences? Take 
# note of any change in the mean α as well.

data("reedfrogs")
d <- reedfrogs
str(d)

d$tank <- 1:nrow(d)

dat <- list(
  S = as.integer(d$surv),
  pred = as.integer(d$pred),
  big=as.integer(d$size),
  N=as.integer(d$density),
  tank=as.integer(d$tank)
)

m13m3 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank],
    a[tank] ~ dcauchy(alpha, sigma),
    alpha ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ), data=dat, chains=4, log_lik = TRUE
)
# from the chapter

m13.2 <- ulam(
  alist(
    S ~ dbinom( N , p ),
    logit(p) <- a[tank],
    a[tank] ~ dnorm(a_bar, sigma),
    a_bar ~ dnorm(0 , 1.5),
    sigma ~ dexp(1)
  ), data=dat, chains=4, log_lik = TRUE
)

precis_plot(precis(m13.2, depth=2), xlim = c(-2,20))
precis_plot(precis(m13m3, depth=2), xlim=c(-2,20))
precis(m13.2)
precis(m13m3)
# basically those a_tanks that were already quite high are now allowed to grow
# crazily big. Unsure as to what the cauchy distirbution is. It's learnt a
# slightly higher alpha bar now, and also oddly enough a lower sigma, despite
# this being far bigger variation in those crazy alphas

# to be honest though, the large individual alphas are all basically "always"
# on the log-odds scale already before cauchy.

# Apparently Cauchy has an undefined mean cos of tose crazy tails! That's where
# we get those large a values from

m13m3nc <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a_bar + z[tank]*sigma_a,
    z[tank] ~ dcauchy(0, 1),
    a_bar ~ dnorm(0, 1),
    sigma_a ~ dexp(1),
    gq > vector [tank]:a <<- a_bar + z*sigma_a
  ), data=dat, chains=4, log_lik = TRUE
)
compare(m13m3m, m13m3nc)
precis(m13m3, m13m3nc)

# re-parameterised to have no divergent transitions!
