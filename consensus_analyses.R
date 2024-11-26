#using consensus phylogeny for all models, with the updated dataset April 10, 2024 excluding I. paradoxus
#UPDATING november 19, 2024 with corrected depth data----

#load packages----
library(brms)
library(ape)
library(mice)
library(tidyverse)
library(cmdstanr)
library(dagitty)
options(scipen=999) #turn off scientific notation

#set working directory----
setwd("/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution")
getwd()

#read in data and phylogenies
cephtreeMCC <- read.nexus("cephtreeMCC.tree") #consensus (maximum clade credibility) tree
cormat <- vcv(cephtreeMCC, corr=TRUE) #correlation matrix for consensus phylogeny
st.cephdat <- read.csv("st.cephdat.csv") #logged and standardized data
st.cephdat$benthic <- factor(st.cephdat$benthic)
st.cephdat$sociality.bin <- factor(st.cephdat$sociality.bin)
st.cephdat$sociality.3 <- factor(st.cephdat$sociality.3)
st.cephdat$habitat3 <- factor(st.cephdat$habitat3)
st.cephdat$depth_cat <- factor(st.cephdat$depth_cat)

#ML baseline----
m.MLst <- brm(bf(CNS.1 ~ mi(ML.1) + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = st.cephdat,
                   data2 = list(cormat = cormat),
                   iter = 8000, chains = 4,
                   control = list(adapt_delta = 0.89),
                   backend = "cmdstanr",
                   cores = 4)
summary(m.MLst)
save(m.MLst, file="m.MLst.rda")

#age of sexual maturity----
## with maximum age of sexual maturity and lifespan----
m.smmax <- brm(bf(CNS.1 ~ mi(ML.1) + mi(lifespan.max) + mi(matage.max) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                 bf(matage.max | mi() ~ 1 + mi(ML.1) + mi(lifespan.max) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                 bf(lifespan.max | mi() ~ 1 + mi(ML.1) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
               data=st.cephdat, 
               data2=list(cormat=cormat), 
               prior = c(prior(normal(0,0.5), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               iter = 8000, chains = 4,
               control = list(adapt_delta = 0.89),
               backend="cmdstanr", cores=4)
summary(m.smmax)
save(m.smmax, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.smmax.rda")
0.33/0.75 #0.33 instead of 0.3 mean 

## minimum age----
m.smmin <- brm(bf(CNS.1 ~ mi(ML.1) + mi(lifespan.min) + mi(matage.min) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(matage.min | mi() ~ 1 + mi(ML.1) + mi(lifespan.min) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(lifespan.min | mi() ~ 1 + mi(ML.1) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
                data=st.cephdat, 
                data2=list(cormat=cormat), 
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                iter = 8000, chains = 4,
                control = list(adapt_delta = 0.89),
                backend="cmdstanr", cores=4)
summary(m.smmin)
save(m.smmin, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.smmin.rda")
-0.11/0.75 

m.smmean <- brm(bf(CNS.1 ~ mi(ML.1) + mi(lifespan.mean) + mi(matage.mean) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                 bf(matage.mean | mi() ~ 1 + mi(ML.1) + mi(lifespan.mean) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                 bf(lifespan.mean | mi() ~ 1 + mi(ML.1) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
               data=st.cephdat, 
               data2=list(cormat=cormat), 
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               iter = 8000, chains = 4,
               control = list(adapt_delta = 0.89),
               backend="cmdstanr", cores=4)
summary(m.smmean)
save(m.smmean, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.smmean.rda")
0.26 /0.75 #instead of 0.18

#sociality----
## binary----
m.soc <- brm(bf(CNS.1 ~ mi(ML.1) + sociality.bin + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
               bf(ML.1 | mi() ~ 1 + benthic + depth.mean + (1|gr(phy.species, cov=cormat))), 
             data=st.cephdat, 
             data2=list(cormat=cormat), 
             prior = c(prior(normal(0,1), class = Intercept),
                       prior(normal(0,0.5), class = b)),
             iter = 8000, chains = 4,
             control = list(adapt_delta = 0.89),
             backend="cmdstanr", cores=4)
summary(m.soc)
save(m.soc, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.soc.rda")
-0.26/0.75 #-0.26 instead of -0.21

## 3 category----
m.soc3 <- brm(bf(CNS.1 ~ mi(ML.1) + sociality.3 + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat)))+
                bf(ML.1 | mi() ~ 1 + benthic + depth.mean + (1|gr(phy.species, cov=cormat))), 
              data=st.cephdat, 
              data2=list(cormat=cormat), 
              family=gaussian("identity"),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(normal(0, 1), class = b)),
              iter = 8000, chains = 4,
              control = list(adapt_delta = 0.89),
              backend="cmdstanr", cores=4)
summary(m.soc3)
save(m.soc3, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.soc3.rda") 

## sociality decapodiformes only----
decadat <- read.csv("decadat.csv")
m.socdec <- brm(bf(CNS.1 ~ mi(ML.1) + sociality.bin + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + (1|gr(phy.species, cov=cormat))), 
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data=decadat, 
                data2=list(cormat=cormat), 
                backend="cmdstanr", cores=4)
summary(m.socdec)
save(m.socdec, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/nov2024_fits/m.socdec.rda") 

# behavioral complexity----
#TB's constrained imputation function
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

## combined cognition----
m.cog.empty <- brm(bf(CNS.1 ~ mi(ML.1) +  mi(cog2) + benthic + articles.read + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                     bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat)))+
                     bf(cog2 | mi() ~ 1 + benthic + articles.read + (1|gr(phy.species,cov=cormat))),
                   family=gaussian,
                   data=st.cephdat,
                   data2=list(cormat=cormat),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   backend="cmdstanr", 
                   chains = 0)

m.cog.constrained <- constrained_imputation(model = m.cog.empty,
                                            vars = c("ML.1", "cog2"))
summary(m.cog.constrained)
save(m.cog.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.cog.rda")

## defense repertoire----
m.def.empty <- brm(bf(CNS.1 ~ mi(ML.1) +  mi(defense.repertoire) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                      bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat)))+
                      bf(defense.repertoire | mi() ~ 1 + benthic + articles.read + (1|gr(phy.species,cov=cormat))),
                    family=gaussian,
                    data=st.cephdat,
                    data2=list(cormat=cormat),
                    prior = c(prior(normal(0,1), class = Intercept),
                              prior(normal(0,0.5), class = b)),
                    chains = 0,
                    backend="cmdstanr")

m.def.constrained <- constrained_imputation(model = m.def.empty, vars = "defense.repertoire")
summary(m.def.constrained)
save(m.def.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.def.rda")

## foraging repertoire----
m.hunt.empty <- brm(bf(CNS.1 ~ mi(ML.1) +  mi(foraging.repertoire) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                       bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat)))+
                       bf(foraging.repertoire | mi() ~ 1 + benthic + articles.read + (1|gr(phy.species,cov=cormat))),
                     family=gaussian,
                     data=st.cephdat,
                     data2=list(cormat=cormat),
                     prior = c(prior(normal(0,1), class = Intercept),
                               prior(normal(0,0.5), class = b)),
                     chains=0,
                     backend="cmdstanr")

m.hunt.constrained <- constrained_imputation(model = m.hunt.empty,
                                              vars = c("ML.1", "foraging.repertoire"))
summary(m.hunt.constrained)
save(m.hunt.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.hunt.rda")

#ecological richness----
## dietary breadth----
m.diet.empty <- brm(bf(CNS.1 ~ mi(ML.1) + mi(diet.breadth) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                      bf(diet.breadth | mi() ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                      bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
                    family=gaussian,
                   data=st.cephdat,
                   data2=list(cormat=cormat),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   chains = 0,
                   backend="cmdstanr")

m.diet.constrained <- constrained_imputation(model = m.diet.empty,
                                            vars = c("diet.breadth"))

summary(m.diet.constrained)
save(m.diet.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/nov2024_fits/m.diet.rda")

## diet-benthic interaction----
m.dhab.empty <- brm(bf(CNS.1 ~ mi(ML.1) + mi(diet.breadth)*benthic + depth.mean + pos.latmean + articles.read + (1|gr(phy.species,cov=cormat))) +
                 bf(diet.breadth | mi() ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
               family=gaussian,
               data=st.cephdat,
               data2=list(cormat=cormat),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               chains=0,
               backend="cmdstanr", cores=4)

m.dhab.constrained <- constrained_imputation(model = m.dhab.empty, vars = "diet.breadth")
summary(m.dhab.constrained)
save(m.dhab.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.dhab.rda")

## number of predator groups----
m.preds.empty <- brm(bf(CNS.1 ~ mi(ML.1) + mi(predator.breadth) + depth.mean + benthic + articles.read + (1 | gr(phy.species, cov = cormat))) +
                  bf(predator.breadth | mi() ~ 1 + mi(ML.1) + depth.mean + benthic + articles.read + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + (1 | gr(phy.species, cov = cormat))),
                family = gaussian,
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = st.cephdat,
                data2 = list(cormat = cormat),
                chains=0,
                backend = "cmdstanr", cores = 4)
m.preds.constrained <- constrained_imputation(model = m.preds.empty, vars = "predator.breadth")
summary(m.preds.constrained)
0.08/0.75
save(m.preds.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.preds.rda")

## predators-benthic interaction----
m.phab.empty <- brm(bf(CNS.1 ~ mi(ML.1) + mi(predator.breadth)*benthic + articles.read + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))) +
                 bf(predator.breadth | mi() ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))),
               family = gaussian,
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               data = st.cephdat,
               data2 = list(cormat=cormat),
               chains=0,
               backend = "cmdstanr", cores = 4)
m.phab.constrained <- constrained_imputation(model = m.phab.empty, vars = "predator.breadth")
summary(m.phab.constrained)
save(m.phab.constrained, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.phab.rda")

## habitat binary----
m.benth <- brm(bf(CNS.1 ~ mi(ML.1) + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML.1 | mi() ~ 1 + WoS + (1 | gr(phy.species, cov = cormat))),
                prior = c(prior(normal(0, 1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = st.cephdat,
                data2 = list(cormat = cormat),
               iter = 8000, chains = 4,
               control = list(adapt_delta = 0.89),
                backend = "cmdstanr",
                cores = 4)
summary(m.benth)
save(m.benth, file="m.benth.rda")
0.58/0.75
## habitat 3 category----
m.hab <- brm(bf(CNS.1 ~ mi(ML.1) + habitat3 + WoS + (1 | gr(phy.species, cov = cormat))) +
                bf(ML.1 | mi() ~ 1 + WoS + (1 | gr(phy.species, cov = cormat))),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(normal(0,0.5), class = b)),
              data = st.cephdat,
              data2 = list(cormat=cormat),
             iter = 8000, chains = 4,
             control = list(adapt_delta = 0.89),
              backend = "cmdstanr",
              cores = 4)
summary(m.hab)
save(m.hab, file="m.hab.rda")

## latitude range----
m.latrange <- brm(bf(CNS.1 ~ mi(ML.1) + lat.range + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                    bf(ML.1 | mi() ~ 1 + lat.range + benthic + (1 | gr(phy.species, cov = cormat))),
                  prior = c(prior(normal(0,1), class = Intercept),
                            prior(normal(0,0.5), class = b)),
                  data = st.cephdat,
                  data2 = list(cormat=cormat),
                  backend = "cmdstanr",
                  iter=4000,
                  cores = 4)
summary(m.latrange)

## distance from equator----
# rerunning August 2024 for distance from equator not mean 
m.eqdist <- brm(bf(CNS.1 ~ mi(ML.1) + eq.dist + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                    bf(ML.1 | mi() ~ 1 + pos.latmean + benthic + (1 | gr(phy.species, cov = cormat))),
                  prior = c(prior(normal(0,1), class = Intercept),
                            prior(normal(0,0.5), class = b)),
                  data = st.cephdat,
                  data2 = list(cormat=cormat),
                  backend = "cmdstanr",
                  iter=4000,
                  cores = 4)
summary(m.eqdist)
#and that's still nothing
save(m.eqdist, file="m.eqdist.rda")

## mean depth----
m.meandepth <- brm(bf(CNS.1 ~ mi(ML.1) + depth.mean + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                      bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                    prior = c(prior(normal(0,1), class = Intercept),
                              prior(normal(0,0.5), class = b)),
                    data = st.cephdat,
                    data2 = list(cormat = cormat),
                   iter = 8000, chains = 4,
                   control = list(adapt_delta = 0.89),
                    backend = "cmdstanr",
                    cores = 4)
summary(m.meandepth)
save(m.meandepth, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.meandepth.rda")
load(file="m.meandepth.rda")

## maximum depth----
m.maxdepth <- brm(bf(CNS.1 ~ mi(ML.1) + depth.max + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = st.cephdat,
                   data2 = list(cormat = cormat),
                   iter = 8000, chains = 4,
                   control = list(adapt_delta = 0.89),
                   backend = "cmdstanr",
                   cores = 4)
summary(m.maxdepth)
save(m.maxdepth, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.maxdepth.rda")
10/75 #instead of -0.22

## minimum depth----
m.mindepth <- brm(bf(CNS.1 ~ mi(ML.1) + depth.min + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                    bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                  prior = c(prior(normal(0,1), class = Intercept),
                            prior(normal(0,0.5), class = b)),
                  data = st.cephdat,
                  data2 = list(cormat = cormat),
                  iter = 8000, chains = 4,
                  control = list(adapt_delta = 0.89),
                  backend = "cmdstanr",
                  cores = 4)
summary(m.mindepth)
save(m.mindepth, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.mindepth.rda")
-0.17/0.75 #instead of -0.14

## benthic*depth----
m.maxdb <- brm(bf(CNS.1 ~ mi(ML.1) + depth.max*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = st.cephdat,
                data2 = list(cormat = cormat),
               iter = 8000, chains = 4,
               control = list(adapt_delta = 0.89),
               backend = "cmdstanr",
                cores = 4)
summary(m.maxdb)
save(m.maxdb, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.maxdb.rda")

m.mindb <- brm(bf(CNS.1 ~ mi(ML.1) + depth.min*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               data = st.cephdat,
               data2 = list(cormat = cormat),
               iter = 8000, chains = 4,
               control = list(adapt_delta = 0.89),
               backend = "cmdstanr",
               cores = 4)
summary(m.mindb)
save(m.mindb, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.mindb.rda")

m.meandb <- brm(bf(CNS.1 ~ mi(ML.1) + depth.mean*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               data = st.cephdat,
               data2 = list(cormat = cormat),
               iter = 8000, chains = 4,
               control = list(adapt_delta = 0.89),
               backend = "cmdstanr",
               cores = 4)
summary(m.meandb)
save(m.meandb, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/nov2024_fits/m.meandb.rda")

## depth categories----
m.depthcat <- brm(bf(CNS.1 ~ mi(ML.1) + depth_cat + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML.1 | mi() ~ 1 + depth_cat + benthic + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = st.cephdat,
                   data2 = list(cormat=cormat),
                  iter = 8000, chains = 4,
                  control = list(adapt_delta = 0.89),
                   backend = "cmdstanr",
                   cores = 4)
summary(m.depthcat)
save(m.depthcat, file="m.depthcat.rda")

#ASR and signal----
library(phytools)
library(mice)
getwd()
#decomposition of phylogenetic distance matrix into orthogonal vectors (PVRs)
phylo.vectors = PVR::PVRdecomp(cephtreeMCC)
cephdat <- read.csv("cephdat.csv")
ML.dat <- cephdat[c("phy.species", "CNS.1", "ML.1")]
ML.dat$CNS.1 <- log(ML.dat$CNS.1)
ML.dat$ML.1 <- log(ML.dat$ML.1)
ML.dat <- complete(mice(ML.dat)) #imputation

#calculate EQ
ML.dat$EQ <- ML.dat$CNS.1/ML.dat$ML.1
logEQ <- as.vector(ML.dat$EQ)
names(logEQ) <- ML.dat$phy.species

#ancestral state reconstruction and plot
ASR <- contMap(cephtreeMCC, logEQ, plot=FALSE)

length(ASR$cols)

ASR$cols[1:1001]<-colorRampPalette(c("#feba2c","#d6556d","#2a0593"))(1001)
plot(ASR)

# Plot the mapped characters with the new colors
plot(ASR, type="fan", outline=FALSE, legend = 0.7*max(nodeHeights(cephtreeMCC)),
  fsize = c(0.5, 0.7))

save(ASR, file="ASR.rda")
load(file="ASR.rda")
plot(ASR)
