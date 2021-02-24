library(rethinking)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# 13E1. Which of the following priors will produce more shrinkage in the 
# estimates? 
#   (a) α tank ~ Normal(0, 1); 
#   (b) α tank ~ Normal(0, 2).

# (a) will do, as it's got a smaller sigma, hence allows less variation in each
# tanks individual intercept parameters

# 13E2. Rewrite the following model as a multilevel model.
# y i ~ Binomial(1, p i )
# logit(p i ) = α group[i] + βx i
# α group ~ Normal(0, 1.5)
# β ~ Normal(0, 0.5)

y ~ dbinom(1 , p)
logit(p) <- a[group] + b*x
b ~ dnorm(0, 1.5)
a[group] ~ dnorm(a_bar, sigma_a)
a_bar ~ dnorm(0, 1.5)
sigma_a ~ dexp(1)


# 13E3. Rewrite the following model as a multilevel model.
# y i ~ Normal(μ i , σ)
# μ i = α group[i] + βx i
# α group ~ Normal(0, 5)
# β ~ Normal(0, 1)
# σ ~ Exponential(1)

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
# Specifically, these were tanks with 100% survival that had crazy values, 
# since cauchy's tails are infinite, the posterior can extend forever...

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
compare(m13m3, m13m3nc)
precis(m13m3)
precis(m13m3nc)

# re-parameterised to have no divergent transitions! orrr... not. hmmm.

# 13M4. Now use a Student-t distribution with ν = 2 for the intercepts:
  α_tank ~ Student(2, α, σ)
# Refer back to the Student-t example in Chapter 7 (page 234), if necessary. 
# Compare the resulting posterior to both the original model and the Cauchy 
# model in 13M3. Can you explain the differences and similarities in shrinkage 
# in terms of the properties of these distributions?

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

m13m4 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank],
    a[tank] ~ dstudent(2, a_bar, sigma_a),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1)
  ), data=dat, chains=4, log_lik = TRUE
)

# comparing shrinkages between the three models:

post.norm <- extract.samples(m13.2)
post.cauchy <- extract.samples(m13m3)
post.student <- extract.samples(m13m4)

d$propsurv.est.norm <- inv_logit(apply(post.norm$a, 2, mean))
d$propsurv.est.cauchy <- inv_logit(apply(post.cauchy$a, 2, mean))
d$propsurv.est.student <- inv_logit(apply(post.student$a, 2, mean))

# normal adaptive prior
plot(d$propsurv, ylim=c(0,1), pch=16, xaxt="n", xlab="tank", 
     ylab="proportion survived", col=rangi2)
axis(1, at = c(1,16,32,48), labels=c(1,16,32,48))
points(d$propsurv.est.norm)
abline(h=mean(inv_logit(post.norm$a_bar)), lty=2)
abline(v=16.5, lwd=2)
abline(v=32.5, lwd=2)
text(8, 0, "small tanks")
text(16+8, 0, "medium tanks")
text(32+8, 0, "large tanks")
title("Adaptive Prior with Normal Distribution")

# Cauchy adaptive prior
plot(d$propsurv, ylim=c(0,1), pch=16, xaxt="n", xlab="tank", 
     ylab="proportion survived", col=rangi2)
axis(1, at = c(1,16,32,48), labels=c(1,16,32,48))
points(d$propsurv.est.cauchy)
abline(h=mean(inv_logit(post.cauchy$alpha)), lty=2)
abline(v=16.5, lwd=2)
abline(v=32.5, lwd=2)
text(8, 0, "small tanks")
text(16+8, 0, "medium tanks")
text(32+8, 0, "large tanks")
title("Adaptive Prior with Cauchy Distribution")

# Student adaptive prior
plot(d$propsurv, ylim=c(0,1), pch=16, xaxt="n", xlab="tank", 
     ylab="proportion survived", col=rangi2)
axis(1, at = c(1,16,32,48), labels=c(1,16,32,48))
points(d$propsurv.est.student)
abline(h=mean(inv_logit(post.student$a_bar)), lty=2)
abline(v=16.5, lwd=2)
abline(v=32.5, lwd=2)
text(8, 0, "small tanks")
text(16+8, 0, "medium tanks")
text(32+8, 0, "large tanks")
title("Adaptive Prior with Student-t Distribution")

# With normal distribution, more shrinkage in smaller tanks, and more shrinkage
# farther away from a_bar

# With Cauchy dist, do get more shrinkage with smaller tanks, but seems like the
# inverse is true, in that we have more shrinkage FARTHER AWAY from alpha.
# Probably due to the fat tails of the Cauchy dist. 

# Student-t similar to cauchy


# 13M5. Modify the cross-classified chimpanzees model m13.4 so that the adaptive 
# prior for blocks contains a parameter γ̄ for its mean:
  γ_j  ~ Normal(γ̄, σ γ )
  γ̄ ~ Normal(0, 1.5)
# Compare this model to m13.4 . What has including γ̄ done?

data("chimpanzees")
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition

dat_list <- list(
  pulled_left = d$pulled_left,
  actor=d$actor,
  block_id=d$block,
  treatment=as.integer(d$treatment)
)

m13.4 <- ulam(
  alist(
    pulled_left ~ dbinom( 1, p ),
    logit(p) <- a[actor] + g[block_id] + b[treatment],
    b[treatment] ~ dnorm( 0 , 0.5 ),
    
    ## adaptive priors
    a[actor] ~ dnorm(a_bar, sigma_a),
    g[block_id] ~ dnorm(0, sigma_g),
    
    ## hyper-priors
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1)
  )
  , data = dat_list, chains = 4 , cores = 4 , iter=2000, log_lik = TRUE)

m13m5 <- ulam(
  alist(
    pulled_left ~ dbinom( 1, p),
    logit(p) <- a[actor] + g[block_id] + b[treatment],
    b[treatment] ~ dnorm(0, 0.5),
    a[actor] ~ dnorm(a_bar, sigma_a),
    g[block_id] ~ dnorm(g_bar, sigma_g),
    a_bar ~ dnorm(0, 1.5),
    g_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1)
  ), dat=dat_list, chains = 4, iter=2000, log_lik = TRUE
)
compare(m13.4, m13m5)
# treatment posteriors are in the same place, as are all of the actor
# and block intercepts.

# both give the same out-of-sample predictive power

# However, it's like the left leg/right leg example from before.
# variation around a_bar and g_bar is a lot more, gives us massive uncertainty
# around each of these posteriors. LESSON LEARNED - ONLY ONE ADAPTIVE MEAN


# 13M6. Sometimes the prior and the data (through the likelihood) are in
# conflict, because they concentrate around different regions of parameter
# space. What happens in these cases depends a lot upon the shape of the tails 
# of the distributions. Likewise, the tails of distributions strongly influence
# whether outliers are shrunk or not towards the mean. I want you to consider 
# four different models to fit to one observation at y = 0. The models differ 
# only in the distributions assigned to the likelihood and prior. Here are the 
# four models:
  
# Model NN:
  y ~ Normal(μ, 1)
  μ ~ Normal(10, 1)
# Model NT: 
  y ~ Normal(μ, 1)
  μ ~ Student(2, 10, 1)
# Model TN: 
  y ~ Student(2, μ, 1)
  μ ~ Normal(10, 1)
# Model TT: 
  y ~ Student(2, μ, 1)
  μ ~ Student(2, 10, 1)

# Estimate the posterior distributions for these models and compare them. Can 
# you explain the results, using the properties of the distributions
m13nn <- ulam(
  alist(
    y ~ dnorm(mu, 1),
    mu ~ dnorm(10, 1)
  ), data=list(y=1), chains=4, log_lik = TRUE
)

m13nt <- ulam(
  alist(
    y ~ dnorm(mu, 1),
    mu ~ dstudent(2, 10, 1)
  ), data=list(y=1), chains=4, log_lik = TRUE
)

m13tn <- ulam(
  alist(
    y ~ dstudent(2, mu, 1),
    mu ~ dnorm(10, 1)
  ), data=list(y=1), chains=4, log_lik = TRUE
)

m13tt <- ulam(
  alist(
    y ~ dstudent(2, mu, 1),
    mu ~ dstudent(2, 10, 1)
  ), data=list(y=1), chains=4, log_lik = TRUE
)

# 13H1. In 1980, a typical Bengali woman could have 5 or more children in her 
# lifetime. By the year 2000, a typical Bengali woman had only 2 or 3. You’re 
# going to look at a historical set of data, when contraception was widely 
# available but many families chose not to use it. These data reside in 
# data(bangladesh) and come from the 1988 Bangladesh Fertility Survey. Each row 
# is one of 1934 women. There are six variables, but you can focus on two of 
# them for this practice problem:
#   (1) district : ID number of administrative district each woman resided in
#   (2) use.contraception : An indicator (0/1) of whether the woman was using 
# contraception

# The first thing to do is ensure that the cluster variable, district , is a 
# contiguous set of integers. Recall that these values will be index values 
# inside the model. If there are gaps, you’ll have parameters for which there is 
# no data to inform them. Worse, the model probably won’t run. Look at the 
# unique values of the district variable:
  sort(unique(d$district))
# District 54 is absent. So district isn’t yet a good index variable, because 
# it’s not contiguous. This is easy to fix. Just make a new variable that is 
# contiguous. This is enough to do it:
  d$district_id <- as.integer(as.factor(d$district))
  sort(unique(d$district_id))

# Now there are 60 values, contiguous integers 1 to 60. Now, focus on predicting 
# use.contraception, clustered by district_id. Fit both (1) a traditional 
# fixed-effects model that uses an index variable for district and (2) a 
# multilevel model with varying intercepts for district. Plot the predicted 
# proportions of women in each district using contraception, for both the 
# fixed-effects model and the varying-effects model. That is, make a plot in 
# which district ID is on the horizontal axis and expected proportion using
# contraception is on the vertical. Make one plot for each model, or layer them 
# on the same plot, as you prefer. How do the models disagree? Can you explain
# the pattern of disagreement? In particular, can you explain the most extreme 
# cases of disagreement, both why they happen where they do and why the models
# reach different inferences?
library(rethinking)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

data("bangladesh")
d <- bangladesh
?bangladesh
d$district_id <- as.integer(as.factor(d$district))


dat <- list(
  N = aggregate(d$use.contraception, list(d$district_id), length)$x,
  C = aggregate(d$use.contraception, list(d$district_id), sum)$x,
  did = 1:60
)
m13h1.fix <- ulam(
  alist(
    C ~ dbinom(N, p),
    logit(p) <- a[did],
    a[did] ~ dnorm(0, 1.5)
  ), data=dat, chains=4, log_lik = TRUE
)

m13h1.var <- ulam(
  alist(
    C ~ dbinom(N, p),
    logit(p) <- a[did],
    a[did] ~ dnorm(a_bar, sigma_a),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1)
), data=dat, chains=4, log_lik = TRUE
)

post.fix <- extract.samples(m13h1.fix)
post.var <- extract.samples(m13h1.var)

dat$fix.pconc <- inv_logit(apply(post.fix$a, 2, mean))
dat$var.pconc <- inv_logit(apply(post.var$a, 2, mean))
dat$pconc <- dat$C / dat$N
plot(dat$pconc, ylim=c(0,1), pch=16, xlab="District ID", 
     ylab="Proportion using Contraception", col=rangi2)
points(dat$fix.pconc)
points(dat$var.pconc, pch=4)
abline(h=sum(dat$C) / sum(dat$N))
abline(h=mean(inv_logit(post.var$a_bar)), lty=2)
# crosses (varying effects model) have shrunk more than the fixed ones
# they've shrunk more when farther away from the mean

# I think if we were to plot by size of district we'd see an effect as well

# 13H2. Return to data(Trolley) from Chapter 12. Define and fit a varying 
# intercepts model for these data. Cluster intercepts on individual participants 
# as indicated by the unique values in the id variable. Include action, 
# intention, and contact as ordinary terms. Compare the varying intercepts model 
# and a model that ignores individuals, using both WAIC and posterior 
# predictions. What is the impact of individual variation in these data?

data("Trolley")
d <- Trolley

dat <- list(
  R=d$response,
  C=d$contact,
  I=d$intention,
  A=d$action,
  story=as.integer(d$story),
  indiv=as.integer(d$id)
)

m13h2.base <- ulam(
  alist(
    R ~ dordlogit(phi, cutpoints),
    phi <- bC*C + bI*I + bA*A,
    c(bC, bI, bA) ~ dnorm(0, 0.5),
    cutpoints ~ dnorm(0, 1.5)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)

m13h2.indiv <- ulam(
  alist(
    R ~ dordlogit(phi, cutpoints),
    phi <- bC*C + bI*I + bA*A + a[indiv],
    c(bC, bI, bA) ~ dnorm(0, 0.5),
    cutpoints ~ dnorm(0, 1.5),
    a[indiv] ~ dnorm(a_bar, sigma_a),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)

compare(m13h2.base, m13h2.indiv)
precis(m13h2.base)
precis(m13h2.indiv)

# 13H3. The Trolley data are also clustered by story, which indicates a unique 
# narrative for each vignette. Define and fit a cross-classified varying 
# intercepts model with both id and story. Use the same ordinary terms as in the 
# previous problem. Compare this model to the previous models. What do you infer 
# about the impact of different stories on responses?

m13h2.stor <- ulam(
  alist(
    R ~ dordlogit(phi, cutpoints),
    phi <- bC*C + bI*I + bA*A + a[indiv] + g[story],
    c(bC, bI, bA) ~ dnorm(0, 0.5),
    cutpoints ~ dnorm(0, 1.5),
    a[indiv] ~ dnorm(a_bar, sigma_a),
    g[story] ~ dnorm(0, sigma_g),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)

compare(m13h2.base, m13h2.indiv, m13h2.stor)
precis(m13h2.base)
precis(m13h2.indiv)
precis(m13h2.stor)
