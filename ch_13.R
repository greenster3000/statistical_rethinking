library(rethinking)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
data(reedfrogs)
d <- reedfrogs
str(d)
?reedfrogs

# Varying intercepts model

# Don't have just one intercept (lose info by averaging survivals across tanks)

# Although we have one intercept per tank, also estimate variation among tanks
# Different to dummy variable, in which tanks can't learn from each other

# Instead have adaptive pooling 

# First model is just a repeat of earlier models with an intercept per tank

d$tank <- 1:nrow(d)

dat <- list(
  S = d$surv,
  N = d$density,
  tank = d$tank
)

m13.1 <- ulam(
  alist(
    S ~ dbinom( N, p),
    logit(p) <- a[tank],
    a[tank] ~ dnorm(0, 1.5)
  ), data=dat, chains=4, log_lik = TRUE
)

precis(m13.1, depth=2)
precis_plot(precis(m13.1, depth=2))
postcheck(m13.1, window=24)

# Now for multilevel model with adaptive pooling
# To do this, make the prior for A into a function of some new parameters
# Above, the prior for each of the 48 a[tank]s is dnorm(0, 1.5)
# Now, the prior for each tank is a normal dist with alpha-bar and sigma.
# These priors of alpha-bar and sigma have their own priors
# Called an adaptive prior

m13.2 <- ulam(
  alist(
    S ~ dbinom( N , p ),
    logit(p) <- a[tank],
    a[tank] ~ dnorm(a_bar, sigma),
    a_bar ~ dnorm(0 , 1.5),
    sigma ~ dexp(1)
  ), data=dat, chains=4, log_lik = TRUE
)

compare(m13.1, m13.2)

precis(m13.2)

# Ordinary model has 48 parameters, but 25 effective parameters
# Multilevel model has 50 parameters, but 22 effective ones
# This is because it can adaptively regularize with this multilevel feature
# Results in a less flexible posterior, fewer effective params

# Has learnt a similar regularizing sigma

post <- extract.samples(m13.2)
d$propsurv.est <- inv_logit(apply(post$a, 2, mean))
plot(d$propsurv, ylim=c(0,1), pch=16, xaxt="n", xlab="tank", 
     ylab="proportion survived", col=rangi2)
axis(1, at = c(1,16,32,48), labels=c(1,16,32,48))
points(d$propsurv.est)
abline(h=mean(inv_logit(post$a_bar)), lty=2) 
# estimated median survival
# NOT mean survival
abline(v=16.5, lwd=2)
abline(v=32.5, lwd=2)
text(8, 0, "small tanks")
text(16+8, 0, "medium tanks")
text(32+8, 0, "large tanks")

# Three things due to adaptive pooling:
# All the estimates of survival (open points) are closer to the line than actual
# Small tanks have shrunk to the median estimate more (less data in small tanks)
# Points farther away from median shrink more

plot(NULL , xlim = c(-3,4), ylim=c(0,0.35),
     xlab="log-odds survive", ylab="density")

for ( i in 1:100)
  curve(dnorm(x,post$a_bar[i], post$sigma[i]), add = TRUE,
        col=col.alpha("black", 0.2))

sim_tanks <- rnorm(8000, post$a_bar, post$sigma)
dens(inv_logit(sim_tanks), lwd=2, adj=0.1)

# Here there's simulation of 100 log-odds survival, and the implied probability

# Now on to simulating some data. Simulate so we know the no-pooling and partial
# pooling estimate of survival, and compare to our original values that we sim

# Partial pooling (varying effects/adaptive priors) is better than complete
# pooling (one intercept for all ponds) as complete pooling tends to underfit
# Also better than no pooling (one intercept per pond) as this model doesn't 
# share any information between ponds, and tends to overfit to each pond.
# Especially if one pond has very little data. 

# These are the values we're going to try and recover later
a_bar <- 1.5
sigma <- 1.5

# Simulating number of survivors for each pond
nponds <- 60
Ni <- as.integer(rep(c(5,10,25,35), each=15))
set.seed(5005)
# Creating a true a_pond we're going to try and get back
a_pond <- rnorm(nponds, mean=a_bar, sd=sigma)
dsim <- data.frame(pond=1:nponds, Ni=Ni, true_a=a_pond)

# And simulating number of survivors using this true a_pond
dsim$Si <- rbinom( nponds, prob=logistic(dsim$true_a), size = dsim$Ni)

# and the probability of just taking the direct survival ratio per pond
# This is the no-pooling estimate of probability of survival
dsim$p_nopool <- dsim$Si / dsim$Ni

dat <- list(Si=dsim$Si, Ni=dsim$Ni, pond=dsim$pond)
m13.3 <- ulam(
  alist(
    Si ~ dbinom(Ni, p),
    logit(p) <- a_pond[pond],
    a_pond[pond] ~ dnorm(a_bar, sigma),
    a_bar ~ dnorm(0, 1.5),
    sigma ~ dexp(1)
  ), data=dat, chains=4
)

precis(m13.3, depth=2)

post <- extract.samples(m13.3)
dsim$p_partpool <- apply(inv_logit(post$a_pond), 2, mean)
dsim$p_true <- inv_logit(dsim$true_a)
nopool_error <- abs(dsim$p_nopool - dsim$p_true)
partpool_error <- abs(dsim$p_partpool - dsim$p_true)

mean(nopool_error)
mean(partpool_error)

plot(1: 60, nopool_error, xlab="pond", ylab="absolute error",
     col=rangi2, pch=16)
points(1:60, partpool_error)
# This function aggregate applies the given function (here mean) to the data,
# and aggregates the mean over the different values passed as the second
# argument in a nice table (here 4 different values of Ni for pond size)
# Length of first two arguments has to be the same size

nopool_avg <- aggregate(nopool_error, list(dsim$Ni),mean)
partpool_avg <- aggregate(partpool_error, list(dsim$Ni),mean)


for ( i in 1:4 ) {
  xval = (i-1) * 15
  abline(v=15.5 + xval)
  y_val_np = nopool_avg[i,2]
  y_val_pp = partpool_avg[i, 2]
  segments(x0 = xval , y0=y_val_np, x1=15.5 + xval, y1=y_val_np, lty=2)
  segments(x0 = xval , y0=y_val_pp, x1=15.5 + xval, y1=y_val_pp, lty=1)
}
  
# Same deal as before, error of partial pooling is, on average, lower than that
# for no pooling. Get better results for small ponds on left than bigger ponds
# on right. No harm in doing partial pooling at right, same result, model just
# learns to pay less attention to the grand mean and more attention to the 
# data in each pond

#########

# CHIMPS
# Back to chimp data where we're looking at more than one kind of cluster,
# both actor (seven chimps) and experimental block (6 days of the experiment)

# This is cross-classified data, not all chimps nested in the same experimental
# block. If it were, would be called hierarchical. However, can also make 
# adaptive priors for adaptive priors, and these are also called hierarchical

# Just repeat the recipe as before for the block varying intercept the same as 
# the actor (or pond before) 

# The reason we only add an alpha bar for the actor and not the block is because
# each row of the data gets its own parameters here. All the shift is
# accommodated in the one alpha bar. 

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
precis(m13.4, depth=2)
precis_plot(precis(m13.4, depth=2))

# sigma_g spends time around 0 (bounded here) hence inefficient and low n_eff

dens(extract.samples(m13.4)$sigma_a, xlim=c(0, 4), ylim=c(0, 3.5))
dens(extract.samples(m13.4)$sigma_g, add=TRUE)

# sigma_g a lot closer to zero and a lot more concentrated. sigma_a larger and 
# more spread out. Can see in the precis_plot

# if we ignore the block, do we see any evidence of overfitting cf with block?

m13.5 <- ulam(
  alist(
    pulled_left ~ dbinom( 1, p),
    logit(p) <-  a[actor] + b[treatment],
    b[treatment] ~ dnorm(0, 0.5),
    a[actor] ~ dnorm(a_bar, sigma_a),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1)
  ), data=dat_list, chains=4, cores=4, iter=2000, log_lik = TRUE
)

compare(m13.4, m13.5)

# pWAIC only 2 smaller despite having 7 fewer actual parameters.
# WAIC almost exactly the same, less than the dSE 
# No need to "select" the model. Do that from the experimental design

m13.6 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + g[block_id] + b[treatment] ,
    b[treatment] ~ dnorm( 0 , sigma_b ),
    a[actor] ~ dnorm( a_bar , sigma_a ),
    g[block_id] ~ dnorm( 0 , sigma_g ),
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1),
    sigma_b ~ dexp(1)
  ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )
coeftab( m13.4 , m13.6 )
precis_plot(precis(m13.6, depth=2))
compare(m13.4, m13.5, m13.6)
# we can do this! Can use partial pooling on "treatment" effects, despite that
# being the thing we're controlling for. Doesn't matter how the clusters arise,
# we can use partial pooling to get better inferences. Here, it barely makes a 
# difference to the treatment parameters at all, as there's lots of data. 

########
# DIVERGENT TRANSITIONS
# In the models above, we get divergent transitions. This is due to difficult to
# sample parts of the posterior. Inefficient sampling. Can either increase
# Stan's adapt_delta param to give smaller step-sizes (but won't always work) or
# re-parameterise our model so that it's mathematically identical, but 
# numerically easier to sample from. Just like standardising data and giving it
# mean 0 and sd 1

m13.7 <- ulam(
  alist(
    v ~ normal(0 ,3),
    x ~ normal(0,exp(v))
  ), data=list(N=1), chains=4
)
divergent(m13.7)
precis(m13.7)
traceplot(m13.7)

# 170 divergent transitions, far too high r-hat, low n_eff and horrible trace

m13.7nc <- ulam(
  alist(
    v ~ normal(0, 3),
    z ~ normal(0, 1),
    real[1]:x <- z*exp(v)
  ), data=list(N=1), chains=4
)
divergent(m13.7nc)
precis(m13.7nc)
traceplot(m13.7nc)

# all looks far better now. have to add in real[1]: at the start to make stan
# work, and remember it's not a random variable, so assign it with <-.
# Also, have to add in gq> < around real[1]:x to make stan calcualte it
# not too sure why normal and not rnorm

# this is just like dividing by the s.d. to normalise. to reverse the process, 
# have to multiply by s.d.

# now for the same but with the chimps
# treatment stays the same cos it's not multi-level. 
# for the other ones with adaptive priors (and hence multilevel) have to 
# "subtract" out the s.d. part and include it in the original model definition
# line
m13.4nc <- ulam(
  alist(
    pulled_left ~ dbinom( 1, p),
    logit(p) <- a_bar + z[actor]*sigma_a + x[block_id]*sigma_g + b[treatment],
    b[treatment] ~ dnorm(0, 0.5),
    z[actor] ~ dnorm(0 , 1),
    x[block_id] ~ dnorm(0 , 1),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1),
    gq> vector[actor]:a <<-a_bar * z*sigma_a,
    gq> vector[block_id]:g <<- x*sigma_g
  ), data=dat_list, chains=4, cores=4, log_lik = TRUE
)

compare(m13.4nc, m13.4)
divergent(m13.4)
divergent(m13.4nc)

# 13.5 - Multilevel Posterior Predictions

# Multilvel models are complex, we can't just look at posterior means and
# intervals

# Tools so far: 
# 1. Posterior prediction checks, want to do OK at retrodicting, but not perfect
# 2. Counterfactual plots a la chapter 5, varying inputs and seeing predictions
# 3. WAIC and AIC give us good estimates of out-of-sample model accuracy

# In addition with multi-level models, want to be aware of some nuance
# 1. Adaptive regularizations benefit is that we trade off poorer fit in sample
# for hopefully better fit out of sample. Shrinkage and learning a bit of the
# grand mean

# 2. Predicting either previously seen clusters, or new clusters. These are
# different and require some care. 

# 13.5.1 Posterior Predictions for same clusters
# Chimps - posterior predictions (retrodictions) for each of the 7 clusters 
# (actors)

# Sim data
chimp <- 2
d_pred <- list(
  actor = rep(chimp, 4),
  treatment = 1:4,
  block_id = rep(1, 4)
)

# Link and plot data
p <- link(m13.4, data=d_pred)
p_mu <- apply(p, 2, mean)
p_ci <- apply(p, 2, PI)
plot(NULL, xlab="treatment", ylab="proportion pulled left", ylim=c(0,1), xlim=c(1,4))
lines(1:4, p_mu)
shade(p_ci, 1:4)

# Using sim and remembering the model
post <- extract.samples(m13.4)
str(post)

p_link <- function(treatment , actor=1, block_id = 1) {
  logodds <- with(post, a[,actor] + g[,block_id] + b[,treatment])
  return (inv_logit(logodds))
}
p_raw <- sapply(1:4, function(i) p_link( i, actor=2, block_id = 1))
p_mu <- apply(p_raw, 2, mean)
p_ci <- apply(p_raw, 2, PI)

# Produce the same results. But second one allows more control

# 13.5.2
# Posterior predictions of new clusters

# No unique procedure for generalzing, have to think about your problem and 
# parameterize a simulation that embodies the target generalization

# Example, new chimpanzee (actor). Use a_bar and sigma_a!

p_link_abar <- function(treatment) {
  log_odds <- with(post, a_bar + b[,treatment])
  return(inv_logit(log_odds))
}

# Ignores block, cos we're extrapolating to new blocks
# Ignores actor, because we're not looking at any one actor

# same as before, apply this function
p_raw <- sapply(1:4, function(i) p_link_abar(i))
p_mu <- apply(p_raw, 2, mean)
p_ci <- apply(p_raw, 2, PI)

plot(NULL, xlab="Treatment", ylab="Prop pulled left",
     ylim=c(0,1), xaxt="n", xlim=c(1,4))
axis(1, at=1:4, labels=c("R/N", "L/N", "R/P", "L/P"))
lines(1:4, p_mu)
shade(p_ci, 1:4)


# Can see the effect of prosocial on the left, but doesn't show the variation 
# between actors. This just shows us the average actor and the uncertainty 
# around them.

# To do this use sigma_a and rnorm with a_bar to simulate different actors
# Simulate the random chimps, not the ones in the posterior. Also, simulate them
# first, then "link" (or our own function here)

a_sim <- with(post, rnorm(length(post$a_bar), a_bar, sigma_a))
p_link_asim <- function( treatment ) {
  logodds <- with(post, a_sim + b[,treatment])
  return(inv_logit(logodds))
}
p_raw_asim <- sapply(1:4, function(i) p_link_asim(i))
p_mu_asim <- apply(p_raw_asim, 2, mean)
p_ci_asim <- apply(p_raw_asim, 2, PI)

plot(NULL, xlab="Treatment", ylab="Prop pulled left",
     ylim=c(0,1), xaxt="n", xlim=c(1,4))
axis(1, at=1:4, labels=c("R/N", "L/N", "R/P", "L/P"))
lines(1:4, p_mu_asim)
shade(p_ci_asim, 1:4)

# These are marginal of actor, which means they average over uncertainty among
# actors

# previous ones set actor to the average, ignoring uncertainty among actors

# This shows different simulated actors, accounting for the uncertainty between
# them, as well as showing the effect of prosoc_left, and the effect of the
# interaction between treatment and actor. Those at the top and bottom show less
# variation than those in the middle, treatment has little effect and they show
# strong handedness variation. In the middle, see more effect of treatment

plot(NULL, xlab="Treatment", ylab="Prop pulled left",
     ylim=c(0,1), xaxt="n", xlim=c(1,4))
axis(1, at=1:4, labels=c("R/N", "L/N", "R/P", "L/P"))
for (i in 1:100) lines(1:4, p_raw_asim[i,], col=grau(0.25), lwd=2)