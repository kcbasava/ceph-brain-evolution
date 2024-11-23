#load packages----
library(brms)
library(ape)
library(mice)
library(tidyverse)
library(cmdstanr)

options(scipen=999) #turn off scientific notation
getwd()

setwd("/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution")

#read in data and phylogenies
cephtreeMCC <- read.nexus("cephtreeMCC.tree") #consensus (maximum clade credibility) tree
cormat <- vcv(cephtreeMCC, corr=TRUE) #correlation matrix for consensus phylogeny

#the following code can be run to reproduce the standardized version of the dataset with additional CNS and ML measures----
#or can be read in as a csv below

# cephdat2 <- read.csv("/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/cephdat2.csv")

# #standardize function from rethinking (https://rdrr.io/github/rmcelreath/rethinking/src/R/utilities.r)
# standardize <- function(x) {
#   x <- scale(x)
#   z <- as.numeric(x)
#   attr(z,"scaled:center") <- attr(x,"scaled:center")
#   attr(z,"scaled:scale") <- attr(x,"scaled:scale")
#   return(z)
# }
# 
# #log and standardize continuous vars----
# st.cephdat2 <- cephdat2
# View(st.cephdat2[c(7:34)])
# st.cephdat2[c(7:34)] <- lapply(st.cephdat2[c(7:34)], function(x) round(standardize(log(x)), digits=3))
# st.cephdat2 <- st.cephdat2[,-c(12,13)] #remove original latitude columns with negative values
# nrow(cephdat2[!is.na(cephdat2$CNS.2),])
# 
# # 
# #nov 18 all 3 measures----
# #wide to long
# View(stcd2long)
# #pivot_longer by CNS and ML 
# stcd2long <- pivot_longer(st.cephdat2, cols = starts_with("CNS"), names_to="CNS_measure", values_to="CNS_value") %>%
#   pivot_longer(cols = starts_with("ML"), names_to = "ML_measure", values_to="ML_value")
# 
# #this returns all possible combinations of CNS and ML, drop all by matching pairs
# stcd2long <- stcd2long %>%
#   filter(substring(CNS_measure, 5) == substring(ML_measure, 4))
# 
# getwd()
# write.csv(stcd2long, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/stcd2long.csv")

#Nov 15, 2024: should use all 3 measures for this analysis, not just CNS2 and 3
#read in long data with additional CNS/ML measures----
stcd2long <- read.csv("sensitivity-checks/stcd2long.csv")
stcd2long$benthic <- factor(stcd2long$benthic)
stcd2long$sociality.bin <- factor(stcd2long$sociality.bin)
stcd2long$sociality.3 <- factor(stcd2long$sociality.3)
stcd2long$habitat3 <- factor(stcd2long$habitat3)
stcd2long$depth_cat <- factor(stcd2long$depth_cat)

#trying to drop all CNS NAs to see if mi() on ML works now
stcd2long <- stcd2long[!is.na(stcd2long$CNS_value),]

#maturity age----
## maximum age of sexual maturity ----
m2.smmax <- brm(bf(CNS_value ~ mi(ML_value) + mi(lifespan.max) + mi(matage.max) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                 bf(matage.max | mi() ~ 1 + mi(ML_value) + mi(lifespan.max) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                 bf(lifespan.max | mi() ~ 1 + mi(ML_value) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                 bf(ML_value | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
               data=stcd2long, 
               data2=list(cormat=cormat), 
               prior = c(prior(normal(0,0.5), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               backend="cmdstanr", cores=4)
summary(m2.smmax)
saveRDS(m2.smmax, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.smmax.rds")
# somewhat positive 

## mean age---
m2.smmean <- brm(bf(CNS_value ~ mi(ML_value) + mi(lifespan.mean) + mi(matage.mean) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(matage.mean | mi() ~ 1 + mi(ML_value) + mi(lifespan.mean) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(lifespan.mean | mi() ~ 1 + mi(ML_value) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                  bf(ML_value | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
                data=stcd2long, 
                data2=list(cormat=cormat), 
                prior = c(prior(normal(0,0.5), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                backend="cmdstanr", cores=4)
summary(m2.smmean)
save(m2.smmean, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.smmean.rda")

## min age----
m2.smmin <- brm(bf(CNS_value ~ mi(ML_value) + mi(lifespan.min) + mi(matage.min) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                   bf(matage.min | mi() ~ 1 + mi(ML_value) + mi(lifespan.min) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                   bf(lifespan.min | mi() ~ 1 + mi(ML_value) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                   bf(ML_value | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
                 data=stcd2long, 
                 data2=list(cormat=cormat), 
                 prior = c(prior(normal(0,0.5), class = Intercept),
                           prior(normal(0,0.5), class = b)),
                 backend="cmdstanr", cores=4)
summary(m2.smmin)
save(m2.smmin, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.smmin.rda")

#sociality----
## sociality binary----
m2.soc <- brm(bf(CNS_value ~ mi(ML_value) + sociality.bin + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
               bf(ML_value | mi() ~ 1 + benthic + depth.mean + (1|gr(phy.species, cov=cormat))), 
             data=stcd2long, 
             data2=list(cormat=cormat), 
             prior = c(prior(normal(0,1), class = Intercept),
                       prior(normal(0,0.5), class = b)),
             backend="cmdstanr", cores=4)
summary(m2.soc)
saveRDS(m2.soc, 
file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.soc.rds")
readRDS("/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.soc.rds")
# still mostly negative

## scoaility 3 category----
m2.soc3 <- brm(bf(CNS_value ~ mi(ML_value) + sociality.3 + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat)))+
                bf(ML_value | mi() ~ 1 + benthic + depth.mean + (1|gr(phy.species, cov=cormat))), 
              data=stcd2long, 
              data2=list(cormat=cormat), 
              family=gaussian("identity"),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(normal(0, 1), class = b)),
              backend="cmdstanr", cores=4)
summary(m2.soc3)
save(m2.soc3, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.soc3.rda")

# habitat----
## benthic-pelagic----
m2.benth <- brm(bf(CNS_value ~ mi(ML_value) + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML_value | mi() ~ 1 + WoS + (1 | gr(phy.species, cov = cormat))),
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               data = stcd2long,
               data2 = list(cormat = cormat),
               backend = "cmdstanr",
               cores = 4)
summary(m2.benth)
#less but still mostly positive
saveRDS(m2.benth, 
        file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.benth.rds")

## habitat 3 category----
m2.hab <- brm(bf(CNS_value ~ mi(ML_value) + habitat3 + WoS + (1 | gr(phy.species, cov = cormat))) +
               bf(ML_value | mi() ~ 1 + WoS + (1 | gr(phy.species, cov = cormat))),
             prior = c(prior(normal(0, 1), class = Intercept),
                       prior(normal(0,0.5), class = b)),
             data = stcd2long,
             data2 = list(cormat=cormat),
             backend = "cmdstanr",
             cores = 4)
summary(m2.hab)
save(m2.hab, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.hab.rda")

#depth----
## maximum depth----
m2.maxdepth <- brm(bf(CNS_value ~ mi(ML_value) + depth.max + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                    bf(ML_value | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                  prior = c(prior(normal(0,1), class = Intercept),
                            prior(normal(0,0.5), class = b)),
                  data = stcd2long,
                  data2 = list(cormat = cormat),
                  backend = "cmdstanr",
                  cores = 4)
summary(m2.maxdepth)
saveRDS(m2.maxdepth, 
        file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.maxdepth.rds")
#still pretty negative

## minimum depth----
m2.mindepth <- brm(bf(CNS_value ~ mi(ML_value) + depth.min + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML_value | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = stcd2long,
                   data2 = list(cormat = cormat),
                   backend = "cmdstanr",
                   cores = 4)
summary(m2.mindepth)
save(m2.mindepth, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.mindepth.rda")

## depth interaction with benthic-pelagic, maximum----
m2.maxdb <- brm(bf(CNS_value ~ mi(ML_value) + depth.max*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML_value | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = stcd2long,
                data2 = list(cormat = cormat),
                backend = "cmdstanr",
                cores = 4)
summary(m2.maxdb)
save(m2.maxdb, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.maxdb.rda")

## depth interaction with benthic-pelagic, minimum----
m2.mindb <- brm(bf(CNS_value ~ mi(ML_value) + depth.min*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML_value | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = stcd2long,
                data2 = list(cormat = cormat),
                backend = "cmdstanr",
                cores = 4)
summary(m2.mindb)
save(m2.mindb, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.mindb.rda")

## depth categories----
m2.depthcat <- brm(bf(CNS_value ~ mi(ML_value) + depth_cat + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                    bf(ML_value | mi() ~ 1 + depth_cat + benthic + (1 | gr(phy.species, cov = cormat))),
                  prior = c(prior(normal(0,1), class = Intercept),
                            prior(normal(0,0.5), class = b)),
                  data = stcd2long,
                  data2 = list(cormat=cormat),
                  backend = "cmdstanr",
                  cores = 4)
summary(m2.depthcat)
save(m2.depthcat, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.depthcat.rda")

#latitude range----
m2.latrange <- brm(bf(CNS_value ~ mi(ML_value) + lat.range + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                    bf(ML_value | mi() ~ 1 + lat.range + benthic + (1 | gr(phy.species, cov = cormat))),
                  prior = c(prior(normal(0,1), class = Intercept),
                            prior(normal(0,0.5), class = b)),
                  data = stcd2long,
                  data2 = list(cormat=cormat),
                  backend = "cmdstanr",
                  cores = 4)
summary(m2.latrange)
save(m2.latrange, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.latrange.rda")

#equator distance----
m2.eqdist <- brm(bf(CNS_value ~ mi(ML_value) + eq.dist + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML_value | mi() ~ 1 + pos.latmean + benthic + (1 | gr(phy.species, cov = cormat))),
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = stcd2long,
                data2 = list(cormat=cormat),
                backend = "cmdstanr",
                iter=4000,
                cores = 4)
summary(m2.eqdist)
save(m2.eqdist, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.eqdist.rda")

#TB's constrained imputation function----
constrained_imputation <- function(model, vars) {
  
  # Extract and edit Stan code
  scode <- capture.output(stancode(model))
  
  # Initialize imp_code as empty
  imp_code <- scode
  
  for (var in vars) {
    # Match brms var name
    brms_var <- gsub("[\\._]", "", var)
    stan_var <- paste0("Ymi_", brms_var)
    
    # Find lower and upper bounds from the data
    lower_bound <- min(model$data[[var]], na.rm = TRUE)
    upper_bound <- max(model$data[[var]], na.rm = TRUE)
    
    # Amend Stan code and set lower and upper bounds on the imputed variable
    imp_code <- gsub(paste0("vector[Nmi_", brms_var, "] ", stan_var, ";"),
                     paste0("vector<lower=", lower_bound, ", upper=", upper_bound, ">[Nmi_", brms_var, "] ", stan_var, ";"),
                     imp_code, fixed = TRUE)
  }
  
  # Replace and compile the model object with the amended Stan code
  attributes(model$fit)$CmdStanModel <- cmdstan_model(write_stan_file(imp_code))
  
  # Fit to data with modified model
  model_constrain <- update(model,
                            cores = 4, chains = 4, iter = 8000,
                            recompile = FALSE, 
                            control = list(adapt_delta = 0.95))
  
  return(model_constrain)
}

#diet----
## diet basic
m2.diet.empty <- brm(bf(CNS_value ~ mi(ML_value) + mi(diet.breadth) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                      bf(diet.breadth | mi() ~ 1 + mi(ML_value) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                      bf(ML_value | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
                    family=gaussian,
                    data=stcd2long,
                    data2=list(cormat=cormat),
                    prior = c(prior(normal(0,1), class = Intercept),
                              prior(normal(0,0.5), class = b)),
                    chains = 0,
                    backend="cmdstanr")
m2.diet.constrained <- constrained_imputation(model = m2.diet.empty,
                                             vars = c("diet.breadth"))
summary(m2.diet.constrained)
save(m2.diet.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.diet.rda")

## diet-benthic interaction----
m2.dhab.empty <- brm(bf(CNS_value ~ mi(ML_value) + mi(diet.breadth)*benthic + depth.mean + pos.latmean + articles.read + (1|gr(phy.species,cov=cormat))) +
                      bf(diet.breadth | mi() ~ 1 + mi(ML_value) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                      bf(ML_value | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
                    family=gaussian,
                    data=stcd2long,
                    data2=list(cormat=cormat),
                    prior = c(prior(normal(0,1), class = Intercept),
                              prior(normal(0,0.5), class = b)),
                    chains=0,
                    backend="cmdstanr", cores=4)
m2.dhab.constrained <- constrained_imputation(model = m2.dhab.empty, vars = "diet.breadth")
summary(m2.dhab.constrained)
save(m2.dhab.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.dhab.rda")

#predators----
## number of predator groups----
m2.preds.empty <- brm(bf(CNS_value ~ mi(ML_value) + mi(predator.breadth) + depth.mean + benthic + articles.read + (1 | gr(phy.species, cov = cormat))) +
                       bf(predator.breadth | mi() ~ 1 + mi(ML_value) + depth.mean + benthic + articles.read + (1 | gr(phy.species, cov = cormat))) +
                       bf(ML_value | mi() ~ 1 + benthic + depth.mean + (1 | gr(phy.species, cov = cormat))),
                     family = gaussian,
                     prior = c(prior(normal(0,1), class = Intercept),
                               prior(normal(0,0.5), class = b)),
                     data = stcd2long,
                     data2 = list(cormat = cormat),
                     chains=0,
                     backend = "cmdstanr", cores = 4)
m2.preds.constrained <- constrained_imputation(model = m2.preds.empty, vars = "predator.breadth")
summary(m2.preds.constrained)
save(m2.preds.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.preds.rda")

## predators-benthic interaction----
m2.phab.empty <- brm(bf(CNS_value ~ mi(ML_value) + mi(predator.breadth)*benthic + articles.read + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))) +
                      bf(predator.breadth | mi() ~ 1 + mi(ML_value) + articles.read + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))) +
                      bf(ML_value | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))),
                    family = gaussian,
                    prior = c(prior(normal(0,1), class = Intercept),
                              prior(normal(0,0.5), class = b)),
                    data = stcd2long,
                    data2 = list(cormat=cormat),
                    chains=0,
                    backend = "cmdstanr", cores = 4)
m2.phab.constrained <- constrained_imputation(model = m2.phab.empty, vars = "predator.breadth")
summary(m2.phab.constrained)
save(m2.phab.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/sensitivity-checks/m2.phab.rda")

#dropping S. lessoniana from some analyses for sensitivity check----
st.cephdat <- read.csv("st.cephdat.csv") #logged and standardized data
st.cephdat$benthic <- factor(st.cephdat$benthic)
st.cephdat$sociality.bin <- factor(st.cephdat$sociality.bin)
st.cephdat$sociality.3 <- factor(st.cephdat$sociality.3)
st.cephdat$habitat3 <- factor(st.cephdat$habitat3)
st.cephdat_noless$depth_cat <- factor(st.cephdat$depth_cat)
View(st.cephdat[-c("Sepioteuthis_lessoniana"),])
st.cephdat_noless <- st.cephdat[-70,]

#noless sociality
m.soc_noless <- brm(bf(CNS.1 ~ mi(ML.1) + sociality.bin + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                bf(ML.1 | mi() ~ 1 + benthic + depth.mean + (1|gr(phy.species, cov=cormat))), 
              data=st.cephdat_noless, 
              data2=list(cormat=cormat), 
              prior = c(prior(normal(0,1), class = Intercept),
                        prior(normal(0,0.5), class = b)),
              backend="cmdstanr", cores=4)
summary(m.soc_noless)

#benthic noless 
m.benth_noless <- brm(bf(CNS.1 ~ mi(ML.1) + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML.1 | mi() ~ 1 + WoS + (1 | gr(phy.species, cov = cormat))),
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               data = st.cephdat_noless,
               data2 = list(cormat = cormat),
               backend = "cmdstanr",
               cores = 4)
summary(m.benth_noless)

## benthic*depth noless
m.maxdb_noless <- brm(bf(CNS.1 ~ mi(ML.1) + depth.max*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               data = st.cephdat_noless,
               data2 = list(cormat = cormat),
               backend = "cmdstanr",
               cores = 4)
summary(m.maxdb_noless)

m.mindb_noless <- brm(bf(CNS.1 ~ mi(ML.1) + depth.min*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               data = st.cephdat_noless,
               data2 = list(cormat = cormat),
               backend = "cmdstanr",
               cores = 4)
summary(m.mindb_noless)

#diet-benthic noless
m.dhab.empty_noless <- brm(bf(CNS.1 ~ mi(ML.1) + mi(diet.breadth)*benthic + depth.mean + pos.latmean + articles.read + (1|gr(phy.species,cov=cormat))) +
                      bf(diet.breadth | mi() ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                      bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
                    family=gaussian,
                    data=st.cephdat_noless,
                    data2=list(cormat=cormat),
                    prior = c(prior(normal(0,1), class = Intercept),
                              prior(normal(0,0.5), class = b)),
                    chains=0,
                    backend="cmdstanr", cores=4)

m.dhab_noless <- constrained_imputation(model = m.dhab.empty_noless, vars = "diet.breadth")
summary(m.dhab_noless)

#predators-benthic noless----
m.phab.empty_noless <- brm(bf(CNS.1 ~ mi(ML.1) + mi(predator.breadth)*benthic + articles.read + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))) +
                      bf(predator.breadth | mi() ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))) +
                      bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))),
                    family = gaussian,
                    prior = c(prior(normal(0,1), class = Intercept),
                              prior(normal(0,0.5), class = b)),
                    data = st.cephdat_noless,
                    data2 = list(cormat=cormat),
                    chains=0,
                    backend = "cmdstanr", cores = 4)
m.phab_noless <- constrained_imputation(model = m.phab.empty_noless, vars = "predator.breadth")
summary(m.phab_noless)

#mean age mat noless----
m.smmean_noless <- brm(bf(CNS.1 ~ mi(ML.1) + mi(lifespan.mean) + mi(matage.mean) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(matage.mean | mi() ~ 1 + mi(ML.1) + mi(lifespan.mean) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(lifespan.mean | mi() ~ 1 + mi(ML.1) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
                data=st.cephdat_noless, 
                data2=list(cormat=cormat), 
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                backend="cmdstanr", cores=4)
summary(m.smmean_noless)