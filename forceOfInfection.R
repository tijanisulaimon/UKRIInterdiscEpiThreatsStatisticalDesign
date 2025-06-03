
library(serofoi)

tc <- fit_seromodel(serosurvey=chagas2012)


seromodel <- fit_seromodel(
  serosurvey = chagas2012
)

plot_seromodel(seromodel=seromodel,serosurvey=chagas2012)
plot_foi_estimates(seromodel=seromodel,serosurvey=chagas2012)

seromodel <- fit_seromodel(
  serosurvey = chagas2012,
  model_type = "time",
  foi_index = data.frame(
    year = 1935:2011,
    foi_index = c(rep(1, 46), rep(2, 31))
  ),
  iter = 1000
)

plot_foi_estimates(seromodel=seromodel,serosurvey=chagas2012)


chagas2012


#**************************************
survey_year <- 2025
n_sample_mean <- 20
