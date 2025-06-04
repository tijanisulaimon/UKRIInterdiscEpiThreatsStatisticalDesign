# script to estimate power to compare seroprevalence between predicted high and low risk areas
start.time <- Sys.time()

# load packages
library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
library(lme4)
library(parallel)
library(binom)
library(plyr)
library(ggplot2)

# no of data sets to simulate (1000 takes ~ 4 min)
nsim <- 1000

# sampling / design choices
n.village <- 50    # number of villages
n.hh <- 10         # number of households per village
n <- 2             # number of people per household   

# if evenly distribute high and low risk sites across
# the nine sentinel sites (depending on whether sentinel sites)
# themselves differ in risk at a higher scale e.g. Western lower than Eastern -
# can also estimate seroprevalence at sentinel-site level as well as more
# nuanced - comparing average with amount of heterogeneity at a smaller scale

# effect size: OR representing difference in seroprevalence 
# between high and low risk villages
OR <- 3

# load serology data to get parameter (intercept and variance) estimates

# load parameter estimates
par.tab <- read.csv("parameter.estimates.csv", row.names = 1)
virus <- colnames(par.tab)

#*******Example simulate serology data*******************

dat <- expand.grid(hh = 1:n.hh, village = 1:n.village, n = n)
print(sum(dat$n))

# allocate villages to high and low prevalence in 1:1 ratio 
dat$risk.level <- dat$village %% 2 - 0.5
# simulate seropositives
simdat <-
  sim.glmm(
    design.data = dat, 
    fixed.eff = 
      list(
        intercept = par.tab["(Intercept)", 1],
        risk.level = log(OR)),
    distribution = "binomial",
    rand.V = c(hh = par.tab["barcode_hh", 1], 
               village = par.tab["village", 1]))


simdat$risk.level <- as.factor(simdat$risk.level)


summary <- ddply(simdat,.(risk.level),summarise,num=sum(response),denom=sum(n))
meanPrev <- binom.confint(summary$num,summary$denom,methods="exact")
meanPrev$risk.level <- c("Low","High")

ggplot(meanPrev) +
  geom_point(aes(x=risk.level,y=mean*100)) +
  geom_errorbar(aes(x=risk.level,ymin=lower*100,ymax=upper*100))



#*************Simulate, model fit and calculate p-value************
# function to simulate data and estimate p-value for null hypothesis 
# that high and low risk areas have the same seroprevalence
res.tab.fn <- function(...) {
    # create template data set
    dat <- expand.grid(hh = 1:n.hh, village = 1:n.village, n = n)
    print(sum(dat$n))
    # allocate villages to high and low prevalence in 1:1 ratio 
    dat$risk.level <- dat$village %% 2 - 0.5
    # simulate seropositives
    simdat <-
      sim.glmm(
        design.data = dat, 
        fixed.eff = 
          list(
            intercept = par.tab["(Intercept)", sp],
            risk.level = log(OR)),
        distribution = "binomial",
        rand.V = c(hh = par.tab["barcode_hh", sp], 
                   village = par.tab["village", sp]))
    
    fit <- glmer(cbind(response, n - response) ~ risk.level + (1 | hh) +(1 | village), family = binomial, data = simdat)
    fit0 <- update(fit, ~ . - risk.level)
    #coef(summary(fit))["risk.level", "Pr(>|z|)"] # Wald P not reliable - gives inflated type 1 error
    anova(fit, fit0)[2, "Pr(>Chisq)"]
  
}


# repeat simulations many times and calculate p-value
cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterExport(cl, varlist=c("res.tab.fn"),
              envir=environment())               


start <- Sys.time()
sim.res <- parLapply(cl,1:nsim,res.tab.fn) 
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start


# estimate power
apply(do.call("rbind", sim.res) < 0.05, 2, mean)


