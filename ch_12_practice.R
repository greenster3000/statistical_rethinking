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

# So we want 
# 1. femininity + damage_norm
# 2. femininity + min_pressure
# 3. femininity + min_pressure + damage_norm

# 12H5. One hypothesis from developmental psychology, usually attributed to Carol Gilligan, pro-
# poses that women and men have different average tendencies in moral reasoning. Like most hypothe-
# ses in social psychology, it is descriptive, not causal. The notion is that women are more concerned
# with care (avoiding harm), while men are more concerned with justice and rights. Evaluate this hy-
# pothesis, using the Trolley data, supposing that contact provides a proxy for physical harm. Are
# women more or less bothered by contact than are men, in these data? Figure out the model(s) that is
# needed to address this question.
library(rethinking)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
data(Trolley)
d <- Trolley
dat <- list(
  R=d$response,
  C=d$contact,
  M=d$male  
)

m12h5gender <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints),
    phi <- bM * M,
    bM ~ dnorm( 0 , 0.5),
    cutpoints ~ dnorm( 0 , 1.5)
  ), data=dat, chains = 4, cores = 4
) 

m12h5gender.contact <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints),
    phi <- bM*M + bC*C,
    c(bM, bC) ~ dnorm( 0 , 0.5),
    cutpoints ~ dnorm( 0 , 1.5)
  ), data=dat, chains = 4, cores = 4
) 

m12h5gender.interact <- ulam(
  alist(
    R ~ dordlogit( phi , cutpoints),
    phi <- bM*M + bC*C + bCM*C*M,
    c(bM, bC, bCM) ~ dnorm( 0 , 0.5),
    cutpoints ~ dnorm( 0 , 1.5)
  ), data=dat, chains = 4, cores = 4
)

precis_plot(precis(m12h5gender))
precis_plot(precis(m12h5gender.contact))
precis_plot(precis(m12h5gender.interact))

kM <- 0:1
pdat <- data.frame(M=kM)
s <- sim(m12h5gender, data=pdat)
simplehist(s, xlab="response")

kC <- 0
kM <- 0:1
pdat <- data.frame(M=kM, C=kC)
s <- sim(m12h5gender.contact, data=pdat)
simplehist(s, xlab="response")

kC <- 1
kM <- 0:1
pdat <- data.frame(M=kM, C=kC)
s <- sim(m12h5gender.contact, data=pdat)
simplehist(s, xlab="response")

# Gender only model, average response from men and women 
kM <- 0:1
pdat <- data.frame(M=kM)
s <- sim(m12h5gender, data=pdat)
mean(s[,1])
mean(s[,2])
mean(s[,2]) - mean(s[,1])
# Gender and Intent model, average response from men and women when I = 0 or I = 1
kC <- 0
kM <- 0:1
pdat <- data.frame(M=kM, C=kC)
s <- sim(m12h5gender.contact, data=pdat)
mean(s[,1])
mean(s[,2])
mean(s[,2]) - mean(s[,1])

kC <- 1
kM <- 0:1
pdat <- data.frame(M=kM, C=kC)
p <- sim(m12h5gender.contact, data=pdat)
mean(p[,1])
mean(p[,2])
mean(p[,2]) - mean(p[,1])
simplehist(s[,1] - s[,2])
simplehist(p[,1] - p[,2])

# if anything, less of a difference when Contact is involved
# also, posterior for interaction term is around zero
kC <- 0:1
kM <- 0
pdat <- data.frame(C=kC, M=kM)
plot( NULL, type="n", xlab="Male", ylab="probability", 
      xlim=c(0,1), ylim=c(0,1), xaxp=c(0,1,1), yaxp=c(0,1,2))
phi <- link(m12h5gender, data=pdat)
post <- extract.samples(m12h5gender.intent)
for ( s in 1:50 ) {
  pk <- pordlogit(1:6, phi[s,], post$cutpoints[s, ])
  for (i in 1:6) lines(kC, pk[,i], col=col.alpha("blue", alpha=0.1))
}

kC <- 0:1
kM <- 1
pdat <- data.frame(C=kC, M=kM)
phi <- link(m12h5gender, data=pdat)
post <- extract.samples(m12h5gender.intent)
plot( NULL, type="n", xlab="Male", ylab="probability", 
      xlim=c(0,1), ylim=c(0,1), xaxp=c(0,1,1), yaxp=c(0,1,2))
for ( s in 1:50 ) {
  pk <- pordlogit(1:6, phi[s,], post$cutpoints[s, ])
  for (i in 1:6) lines(kC, pk[,i], col=col.alpha("red", alpha=0.1))
}

# try making triptych plots

# 12H6. The data in data(Fish) are records of visits to a national park. See ?Fish for details. The
# question of interest is how many fish an average visitor takes per hour, when fishing. The problem is
# that not everyone tried to fish, so the fish_caught numbers are zero-inflated. As with the monks
# example in the chapter, there is a process that determines who is fishing (working) and another 
# process that determines fish per hour (manuscripts per day), conditional on fishing (working). We want
# to model both. Otherwise we’ll end up with an underestimate of rate of fish extraction from the park.
# You will model these data using zero-inflated Poisson GLMs. Predict fish_caught as a function
# of any of the other variables you think are relevant. One thing you must do, however, is use a proper
# Poisson offset/exposure in the Poisson portion of the zero-inflated model. Then use the hours vari-
# able to construct the offset. This will adjust the model for the differing amount of time individuals
# spent in the park.

data(Fish)
d <- Fish

# try and simulate data and get back to our original values to check our model makes sense

prob_fish <- 0.8 
rate_fish <- 10 # per hour
N <- 1000 # 250 people fishing
set.seed(1999)
go_fish <- rnbinom( N , 1 , prob_fish) # fish or not
hours <- rlnorm(N) # hours in park
fish <- go_fish*rpois( N , rate_fish*hours)
data.test = list(fish=fish, log_hours = log(hours))

m12h6test <- ulam(
  alist(
    fish ~ dzipois( p , lambda ),
    logit(p) <- ap,
    log(lambda) <- log_hours + al,
    ap ~ dnorm(0, 1),
    al ~ dnorm(1, 0.5)
  ), data=data.test, chains=4, cores=4
  
)
post <- extract.samples(m12h6test)
mean(inv_logit(post$ap))
mean(exp(post$al))

# doesn't quite capture the rate of fish per hour, but unsure of what else to do here
# gets back the prob pretty well
d$log_hours <- log(d$hours)
dat = list(
  fish = d$fish_caught,
  persons = d$persons,
  livebait = d$livebait,
  log_hours = d$log_hours,
  kids = d$kids,
  camper=d$camper
)
m12h6 <- ulam(
  alist(
    fish ~ dzipois( p , lambda ),
    logit(p) <- ap,
    log(lambda) <- al + log_hours,
    ap ~ dnorm(0, 1),
    al ~ dnorm(3, 0.5)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)
post <- extract.samples(m12h6)
mean(inv_logit(post$ap))
mean(exp(post$al))
# prob of fishing is 0.33, fish per hour is 0.88
precis(m12h6)

m12h6.persons <- ulam(
  alist(
    fish ~ dzipois( p , lambda ),
    logit(p) <- ap,
    log(lambda) <- al + log_hours + bP*persons,
    ap ~ dnorm(0, 1),
    al ~ dnorm(3, 0.5),
    bP ~ dnorm(0, 0.2)
  ), data=dat, chains=4, cores=4, log_lik = TRUE
)
compare(m12h6, m12h6.persons)

# 12H7. In the trolley data— data(Trolley) —we saw how education level (modeled as an ordered category)
# is associated with responses. But is this association causal? One plausible confound is that
# education is also associated with age, through a causal process: People are older when they finish
# school than when they begin it. Reconsider the Trolley data in this light. Draw a DAG that repre-
# sents hypothetical causal relationships among response, education, and age. Which statical model or
# models do you need to evaluate the causal influence of education on responses? Fit these models to
# the trolley data. What do you conclude about the causal relationships among these three variables?
library(dagitty)
dag12h7 <- dagitty( "dag{ A -> E; E -> R; A -> R }" )
coordinates(dag12h7) <- list( x=c(A=0,R=1,E=2) , y=c(A=0,R=1,E=0) )
drawdag( dag12h7 )

data("Trolley")
d <- Trolley
edu_levels = c(6,1,8,4,7,2,5,3)
dat <- list(
  age=d$age,
  response=d$response,
  education=as.integer(edu_levels[d$edu]),
  alpha=rep(2, 7),
  A=d$action,
  C=d$contact,
  I=d$intention,
  M=d$male
)

m12h7age <- ulam(
  alist(
    response ~ dordlogit( phi , kappa),
    phi <- bA * age,
    kappa ~ dnorm(0, 1.5),
    bA ~ dnorm(0, 0.5)
    ), data=dat, chains=4, cores = 4
)

m12h7edu <- ulam(
  alist(
    response ~ dordlogit( phi , kappa),
    phi <- bE * sum(delta_j[1:education]),
    kappa ~ dnorm(0, 1.5),
    bE ~ dnorm(0, 0.5),
    vector[8]: delta_j <<- append_row(0, delta),
    simplex[7]: delta ~ dirichlet(alpha)
  ), data=dat, chains=4, cores = 4
)
m12h7both <- ulam(
  alist(
    response ~ dordlogit( phi , kappa),
    phi <- bE * sum(delta_j[1:education]) + bA*age,
    kappa ~ dnorm(0, 1.5),
    c(bA,bE) ~ dnorm(0, 0.5),
    vector[8]: delta_j <<- append_row(0, delta),
    simplex[7]: delta ~ dirichlet(alpha)
  ), data=dat, chains=4, cores = 4
)
m12h7all <- ulam(
  alist(
    response ~ dordlogit( phi , kappa),
    phi <- bE * sum(delta_j[1:education]) + bage*age + bI*I + bA*A + bC*C,
    kappa ~ dnorm(0, 1.5),
    c(bage,bE, bI, bC, bA) ~ dnorm(0, 0.5),
    vector[8]: delta_j <<- append_row(0, delta),
    simplex[7]: delta ~ dirichlet(alpha)
  ), data=dat, chains=4, cores = 4
)
m12h8 <- ulam(
  alist(
    response ~ dordlogit( phi , kappa),
    phi <- bE * sum(delta_j[1:education]) + bage*age + bI*I + bA*A + bC*C + bM*M,
    kappa ~ dnorm(0, 1.5),
    c(bage,bE, bI, bC, bA, bM) ~ dnorm(0, 0.5),
    vector[8]: delta_j <<- append_row(0, delta),
    simplex[7]: delta ~ dirichlet(alpha)
  ), data=dat, chains=4, cores = 4
)

# Consider one more variable in the trolley data: Gender. Suppose that gender might influence
# education as well as response directly. Draw the DAG now that includes response, education, age, and
# gender. Using only the DAG, is it possible that the inferences from 12H7 above are confounded by
# gender? If so, define any additional models you need to infer the causal influence of education on
# response. What do you conclude?
dag_12h8 <- dagitty( "dag {
 E -> R
 E <- A -> R
 E <- G -> R}")
drawdag(dag_12h8)
adjustmentSets(dag_12h8 , exposure="E", outcome="R")
# need to also condition on gender
m12h8 <- ulam(
  alist(
    response ~ dordlogit( phi , kappa),
    phi <- bE * sum(delta_j[1:education]) + bage*age + bI*I + bA*A + bC*C + bM*M,
    kappa ~ dnorm(0, 1.5),
    c(bage,bE, bI, bC, bA, bM) ~ dnorm(0, 0.5),
    vector[8]: delta_j <<- append_row(0, delta),
    simplex[7]: delta ~ dirichlet(alpha)
  ), data=dat, chains=4, cores = 4
)
precis_plot(precis(m12h7all))
precis_plot(precis(m12h8))
# education really doesn't have much of an effect, it's all down to gender
postcheck(m12h8)

