# # script to estimate power to compare seroprevalence and force of infection between predicted high and low risk areas
# 
# # assume a time constant and age constant FOI
# 
# 
# start.time <- Sys.time()
# 
# # load packages
# library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
# library(lme4)
# library(parallel)
# library(binom)
# library(plyr)
# library(ggplot2)
# 
# # no of data sets to simulate (1000 takes ~ 4 min)
# nsim <- 1000
# 
# # sampling / design choices
# n.village <- 50    # number of villages
# n.hh <- 10         # number of households per village
# n <- 2             # number of people per household   
# 
# 
# # if evenly distribute high and low risk sites across
# # the nine sentinel sites (depending on whether sentinel sites)
# # themselves differ in risk at a higher scale e.g. Western lower than Eastern -
# # can also estimate seroprevalence at sentinel-site level as well as more
# # nuanced - comparing average with amount of heterogeneity at a smaller scale
# 
# # effect size: OR representing difference in seroprevalence 
# # between high and low risk villages
# OR <- 3
# 
# # load serology data to get parameter (intercept and variance) estimates
# 
# # load parameter estimates
# par.tab <- read.csv("parameter.estimates.csv", row.names = 1)
# virus <- colnames(par.tab)
# 
# #*******Example simulate serology data*******************
# 
# dat <- expand.grid(hh = 1:n.hh, village = 1:n.village, n = n)
# print(sum(dat$n))
# 
# # allocate villages to high and low prevalence in 1:1 ratio 
# dat$risk.level <- dat$village %% 2 - 0.5
# # simulate seropositives
# simdat <-
#   sim.glmm(
#     design.data = dat, 
#     fixed.eff = 
#       list(
#         intercept = par.tab["(Intercept)", 1],
#         risk.level = log(OR)),
#     distribution = "binomial",
#     rand.V = c(hh = par.tab["barcode_hh", 1], 
#                village = par.tab["village", 1]))
# 
# 
# 
# simdat$risk.level <- as.factor(simdat$risk.level)
# 
# summary <- ddply(simdat,.(risk.level),summarise,num=sum(response),denom=sum(n))
# meanPrev <- binom.confint(summary$num,summary$denom,methods="exact")
# meanPrev$risk.level <- c("Low","High")
# 
# ggplot(meanPrev) +
#   geom_point(aes(x=risk.level,y=mean*100)) +
#   geom_errorbar(aes(x=risk.level,ymin=lower*100,ymax=upper*100))
# 
# #*************Simulate, model fit and calculate p-value************
# # function to simulate data and estimate p-value for null hypothesis 
# # that high and low risk areas have the same seroprevalence
# res.tab.fn <- function(...) {
#     # create template data set
#     dat <- expand.grid(hh = 1:n.hh, village = 1:n.village, n = n)
#     print(sum(dat$n))
#     # allocate villages to high and low prevalence in 1:1 ratio 
#     dat$risk.level <- dat$village %% 2 - 0.5
#     # simulate seropositives
#     simdat <-
#       sim.glmm(
#         design.data = dat, 
#         fixed.eff = 
#           list(
#             intercept = par.tab["(Intercept)", sp],
#             risk.level = log(OR)),
#         distribution = "binomial",
#         rand.V = c(hh = par.tab["barcode_hh", sp], 
#                    village = par.tab["village", sp]))
#     
#     fit <- glmer(cbind(response, n - response) ~ risk.level + (1 | hh) + (1 | village), family = binomial, data = simdat)
#     fit0 <- update(fit, ~ . - risk.level)
#     #coef(summary(fit))["risk.level", "Pr(>|z|)"] # Wald P not reliable - gives inflated type 1 error
#     anova(fit, fit0)[2, "Pr(>Chisq)"]
#   
# }
# 
# 
# # repeat simulations many times and calculate p-value
# cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
# clusterExport(cl, varlist=c("res.tab.fn"),
#               envir=environment())               
# 
# 
# start <- Sys.time()
# sim.res <- parLapply(cl,1:nsim,res.tab.fn) # lapply(1:nsim,res.tab.fn) 
# stopCluster(cl)   # stop the clusters
# end <- Sys.time()
# end - start
# 
# 
# # estimate power
# apply(do.call("rbind", sim.res) < 0.05, 2, mean)
# 
# 



# script to estimate power to compare seroprevalence and force of infection between predicted high and low risk areas

# assume a time constant and age constant FOI


start.time <- Sys.time()

# load packages 
library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
library(lme4)
library(parallel)
library(binom)
library(plyr)
library(ggplot2)
library(dplyr)

# no of data sets to simulate (1000 takes ~ 4 min)
nsim <- 1000

# sampling / design choices
n.village <- 50    # number of villages
n.hh <- 10         # number of households per village
n <- 2             # number of people per household   

# Set age groups (5-year bins) and calc. midpoints
age_bins <- seq(0, 95, by = 5)
age_midpoints <- age_bins + 2.5  # midpoint of each 5-year bin?
age_groups <- paste0(age_bins, "-", age_bins + 4)

# between high and low risk villages
OR <- 3
age_OR  <- 1.01    # per‐year odds ratio (0% increase, or set to 1?)

# age distribution (set equal for now?)
age_probs <- rep(1 / length(age_bins), length(age_bins))  

# load serology data to get parameter (intercept and variance) estimates

# load parameter estimates
par.tab <- read.csv("parameter.estimates.csv", row.names = 1)
virus <- colnames(par.tab)

sp <- 1
beta0 <- par.tab["(Intercept)",   sp]    # baseline intercept on log‐odds scale
sigma_hh   <- par.tab["barcode_hh", sp]  # household esimate
sigma_vil  <- par.tab["village",    sp]  # village estimate

#*******Example simulate serology data*******************

dat <- expand.grid(hh = 1:n.hh, village = 1:n.village, n = n)
print(sum(dat$n))

# allocate villages to high and low prevalence in 1:1 ratio 
dat$risk.level <- dat$village %% 2 - 0.5

# Assign age groups and midpoints
dat$age_group <- sample(
  age_groups, 
  size = nrow(dat), 
  replace = TRUE,
  prob = age_probs
)

dat$age_midpoint <- age_midpoints[match(dat$age_group, age_groups)]

head(dat)

# center the age_midpoint 
mean_age <- mean(dat$age_midpoint)
dat$age_mid_c <- dat$age_midpoint - mean_age

# simulate seropositives
simdat <-
  sim.glmm(
    design.data = dat, 
    fixed.eff = 
      list(
        intercept = beta0, # do we need to adjust intercept?
        risk.level = log(OR),
        age_mid_c   = log(age_OR)
        ),
    distribution = "binomial",
    rand.V = c(hh = par.tab["barcode_hh", 1], 
               village = par.tab["village", 1]
               ))


simdat$risk.level <- factor(
  simdat$risk.level,
  levels = c(-0.5, 0.5),
  labels = c("Low", "High")
)

simdat$age_group <- factor(simdat$age_group, levels = age_groups)


risk_summary <- ddply(simdat,.(risk.level),summarise,num=sum(response),denom=sum(n))
meanPrev <- binom.confint(risk_summary$num, risk_summary$denom,methods="exact")
meanPrev$risk.level <- risk_summary$risk.level
meanPrev

ggplot(meanPrev) +
  geom_point(aes(x=risk.level,y=mean*100)) +
  geom_errorbar(aes(x=risk.level,ymin=lower*100,ymax=upper*100))


# Serop by age group (across risk levels, no?)
age_summary <- ddply(
  simdat,
  .(age_group),
  summarise,
  num   = sum(response),
  denom = length(response)
)

age_ci <- binom.confint( x = age_summary$num, n = age_summary$denom, methods = "exact")

age_prev <- cbind(age_summary,mean  = age_ci$mean, lower = age_ci$lower, upper = age_ci$upper)

# Add the numeric midpoint 
age_prev$age_mid <- age_midpoints[match(age_prev$age_group, age_groups)]

# Seroprevalence vs age_midpoint
ggplot(age_prev, aes(x = age_mid, y = mean * 100)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), width = 2) +
  theme_minimal()

# by risk and age
age_risk_summary <- ddply(
  simdat,
  .(risk.level, age_group),
  summarize,
  num   = sum(response),
  denom = sum(n)
)

ci_mat <- binom.confint(x = age_risk_summary$num, n = age_risk_summary$denom, methods = "exact")

age_risk_summary <- cbind(
  age_risk_summary,
  mean  = ci_mat$mean,
  lower = ci_mat$lower,
  upper = ci_mat$upper
)

age_risk_summary$age_mid <- age_midpoints[match(age_risk_summary$age_group, age_groups)]

head(age_risk_summary)


ggplot(age_risk_summary, aes(x = age_mid, y = mean * 100, color = risk.level)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100), width = 4, alpha = 0.6) +
  facet_wrap(~risk.level, ncol = 1) +
  theme_minimal()

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
  
  fit <- glmer(cbind(response, n - response) ~ risk.level + (1 | hh) + (1 | village), family = binomial, data = simdat)
  fit0 <- update(fit, ~ . - risk.level)
  #coef(summary(fit))["risk.level", "Pr(>|z|)"] # Wald P not reliable - gives inflated type 1 error
  anova(fit, fit0)[2, "Pr(>Chisq)"]
  
}


# repeat simulations many times and calculate p-value
cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterExport(cl, varlist=c("res.tab.fn"),
              envir=environment())               


start <- Sys.time()
sim.res <- parLapply(cl,1:nsim,res.tab.fn) # lapply(1:nsim,res.tab.fn) 
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start


# estimate power
apply(do.call("rbind", sim.res) < 0.05, 2, mean)



