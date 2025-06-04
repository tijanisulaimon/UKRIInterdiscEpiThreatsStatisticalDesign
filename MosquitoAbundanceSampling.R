# Simulations to assess numbers of traps required to detect differences 
# in mosquito abundance across ecological strata/ land use types

library("ggplot2")
library("GLMMmisc") # available via devtools::install_github("pcdjohnson/GLMMmisc")
library("parallel")

# From the outbreak entomology report, HLC caught between 0 and 20, but often <5 of each
# species, difficult to find data on numbers in areas that have artisanal mining in other places

# Mosquito surveys - quarterly
# Gabon estimates
# https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0011501
# anthropic 482/1340 hours, forest 110/1645 hours

1340/24 # 55 days
482/55  # 8 per day in anthropic

1645/24  # 68 days
110/68   # 1.62 per day in forest 
# approximately 5x more

#*****************************Simulate data example******************
# create dataset - study design first, assuming three teams
days <- 1:4           # number of days at each area per quarterly survey
traps <- 1:2         # traps set per day 
habitat <- c("forest","cocoa1","cocoa2","palmoil","mining","urban")
area <- LETTERS[1:9]   # district

moz.data <- expand.grid(traps=traps,habitat=habitat,days=days,area=area)

teamAlloc <-  c(rep(1,length(habitat)*length(traps)*length(days)*length(area)/3)
                ,rep(2,length(habitat)*length(traps)*length(days)*length(area)/3)
                ,rep(3,length(habitat)*length(traps)*length(days)*length(area)/3)    
                 )

moz.data$teamAlloc <- teamAlloc

moz.data$daysRep <- NA
moz.data$daysRep <- rep(c(rep(1,length(habitat)*length(traps))
                     ,rep(2,length(habitat)*length(traps))
                     ,rep(3,length(habitat)*length(traps))
                     ,rep(4,length(habitat)*length(traps))
                     ,rep(5,length(habitat)*length(traps))
                     ,rep(6,length(habitat)*length(traps))
                     ,rep(7,length(habitat)*length(traps))
                     ,rep(8,length(habitat)*length(traps))
                     ,rep(9,length(habitat)*length(traps))
                     ,rep(10,length(habitat)*length(traps))
                     ,rep(11,length(habitat)*length(traps))
                     ,rep(12,length(habitat)*length(traps))
                    ),3)

                  
# simulate abundance data, assuming poisson distributed
moz.data <-
  sim.glmm(
    design.data = moz.data,
    fixed.eff =
      list(
        intercept = log(2),      
        habitat=log(
          c(cocoa1=5, 
            cocoa2=5,
            palmoil=5,
            mining=5, 
            urban=1,
            forest=1)
        )),
    rand.V =            # add in some random effects - potentially try and infer from
      inv.mor(          # outbreak report
        c(daysRep = 1.2,   # across time
          area = 1.2)), # across districts             
    distribution = "poisson")   

# Plot

ggplot(moz.data) +
  geom_boxplot(aes(x=habitat,y=response))


#********************Simulations to determine number of adult traps************#
#********************In this one, focusing on larval habitat areas********
#* per year quarter
#* assume three teams take on three sentinel sites each - 3 or 4 days at each sentinel site


pwrFunc <- function(...){
  library("GLMMmisc")
  # create dataset - study design first
  days <- 1:4            # number of days at each area per quarterly survey
  traps <- 1:2         # traps set per day 
  habitat <- c("forest","cocoa1","cocoa2","palmoil","mining","urban")
  area <- LETTERS[1:9]   # district
  
  moz.data <- expand.grid(traps=traps,habitat=habitat,days=days,area=area)
  
  teamAlloc <-  c(rep(1,length(habitat)*length(traps)*length(days)*length(area)/3)
                  ,rep(2,length(habitat)*length(traps)*length(days)*length(area)/3)
                  ,rep(3,length(habitat)*length(traps)*length(days)*length(area)/3)    
  )
  
  moz.data$teamAlloc <- teamAlloc
  
  moz.data$daysRep <- NA
  moz.data$daysRep <- rep(c(rep(1,length(habitat)*length(traps))
                            ,rep(2,length(habitat)*length(traps))
                            ,rep(3,length(habitat)*length(traps))
                            ,rep(4,length(habitat)*length(traps))
                            ,rep(5,length(habitat)*length(traps))
                            ,rep(6,length(habitat)*length(traps))
                            ,rep(7,length(habitat)*length(traps))
                            ,rep(8,length(habitat)*length(traps))
                            ,rep(9,length(habitat)*length(traps))
                            ,rep(10,length(habitat)*length(traps))
                            ,rep(11,length(habitat)*length(traps))
                            ,rep(12,length(habitat)*length(traps))
  ),3)
  
  
  # simulate abundance data, assuming poisson distributed
  moz.data<-
    sim.glmm(
      design.data = moz.data,
      fixed.eff =
        list(
          intercept = log(2),      
          habitat=log(
            c(cocoa1=1.5, 
              cocoa2=2,
              palmoil=3,
              mining=1.5, 
              urban=1,
              forest=1)
          )),
      rand.V =            # add in some random effects - potentially try and infer from
        inv.mor(          # outbreak report
          c(days = 1.2,   # across time
            area = 1.2)), # across districts             
      distribution = "poisson")   
  
  #fit glmm
  moz.pois <-
    lme4::glmer(response ~ habitat + (1 | area) + (1| days) ,
                family = "poisson", data = moz.data)
  output <- summary(moz.pois)
  cfs<-output$coefficients
  forest <-cfs[1,4]
  cocoa1 <-cfs[2,4]
  cocoa2 <-cfs[3,4]
  palmoil <-cfs[4,4]
  mining <-cfs[5,4]
  urban <- cfs[6,4]
  
  pvalforest  <- 0
  pvalcocoa1  <- 0
  pvalcocoa2  <- 0
  pvalpalmoil <- 0
  pvalurban <- 0
  pvalmining <- 0
  
  if(length(output$optinfo$conv$lme4$messages)==0){
    if(forest<=0.05){pvalforest<-1}
    if(cocoa1<=0.05){pvalcocoa1<-1}
    if(cocoa2<=0.05){pvalcocoa2<-1}
    if(palmoil<=0.05){pvalpalmoil<-1}
    if(mining<=0.05){pvalmining<-1}
    if(urban<=0.05){pvalurban<-1}
  }
  return(c(pvalforest,pvalcocoa1,pvalcocoa2,pvalpalmoil,pvalmining,pvalurban))
}


#*******************************Run sims*********************************
cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterExport(cl, varlist=c("pwrFunc"),
envir=environment())               


start <- Sys.time()
sim.res <- parLapply(cl,1:1000,pwrFunc) 
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start


l <- do.call(rbind,sim.res)
l2 <- rowSums(l)
length(l2[l2==5]) 

#*****************************************************************************
# 93% power to detect a 5x difference between urban/ forest and the other categories
# 91% power to detect a 2x difference between urban/ forest and other categories
# 87% power to detect differences of 1x in urban, 1.5x in mining, x3 palmoil, 1.5x cocoa1, 2x cocoa2
