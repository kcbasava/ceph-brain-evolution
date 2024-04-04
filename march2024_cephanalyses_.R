#load packages----
library(brms)
library(ape)
library(mice)
library(tidyverse)
library(cmdstanr)
library(ggdag)
library(dagitty)

#set working directory----
setwd("~/.../cephalopod_analyses")
#read in data and phylogenies
cephtreeMCC <- read.nexus("cephtreeMCC.tree") #consensus (maximum clade credibility) tree
cephtrees100 <- read.nexus("cephtrees100.trees") #100 trees sampled from posterior distribution
cormat <- vcv(cephtreeMCC, corr=TRUE) #correlation matrix for consensus phylogeny
load("ceph100_cormat.rda") #list of 100 correlation matrices for posterior trees
st.cephdat <- read.csv("st.cephdat.csv") #logged and standardized data
st.cephdat$benthic <- factor(st.cephdat$benthic)
st.cephdat$benthic_alt <- factor(st.cephdat$benthic_alt)
st.cephdat$sociality.bin <- factor(st.cephdat$sociality.bin)
st.cephdat$sociality.3 <- factor(st.cephdat$sociality.3)
st.cephdat$habitat3 = factor(st.cephdat$habitat3)
load("cephdag.rda")

#Prediction 1: Age at sexual maturity----
adjustmentSets(cephdag, outcome="CNS", exposure="maturity")
# { ML, WoS, benthic, depth, latitude, lifespan, phylogeny }

## main model: maximum age of sexual maturity and lifespan----
m1.smmax0 <- brm(bf(CNS.1 ~ mi(ML.1) + mi(lifespan.max) + mi(matage.max) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(matage.max | mi() ~ 1 + mi(ML.1) + mi(lifespan.max) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(lifespan.max | mi() ~ 1 + mi(ML.1) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
                data=st.cephdat, 
                data2=list(cormat=ceph100cormat[[1]]), 
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                iter = 4000, chains = 4,
                backend="cmdstanr", cores=4)

#loop model over list of matrices
m1.smmax_loops <- vector("list", 100) 
for (i in seq_along(m1.smmax_loops)) {
  m1.smmax_loops[[i]] <- update(m1.smmax0,
                                data2 = list(cormat = ceph100cormat[[i]],
                                             backend = "cmdstanr",
                                             cores = 4)
  )
}

#combine models
m1.smmax_comb <- combine_models(m1.smmax_loops[[i]])

summary(m1.smmax_comb)
save(m1.smmax_comb, file="m1smmax100.rda")
plot(conditional_effects(m1.smmax_comb, effects = "matage.max", resp = "CNS1"), points=T)

## sensitivity check with minimum age on consensus phylogeny----
m1.smmin <- brm(bf(CNS.1 ~ mi(ML.1) + mi(lifespan.min) + mi(matage.min) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(matage.min | mi() ~ 1 + mi(ML.1) + mi(lifespan.min) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(lifespan.min | mi() ~ 1 + mi(ML.1) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
                data=st.cephdat, 
                data2=list(cormat=cormat), 
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                iter = 4000, chains = 4,
                backend="cmdstanr", cores=4)
summary(m1.smmin)
save(m1.smmin, file="m1smmin.rda")

#Prediction # 2: Sociality-----
adjustmentSets(cephdag, outcome="CNS", exposure="sociality")
# { WoS, benthic, depth, phylogeny }

## main model: binary sociality----
m2.soc0 <- brm(bf(CNS.1 ~ mi(ML.1) + sociality.bin + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                bf(ML.1 | mi() ~ 1 + benthic + depth.mean + (1|gr(phy.species, cov=cormat))), 
              data=st.cephdat, 
              data2=list(cormat=ceph100cormat[[1]]), 
              prior = c(prior(normal(0,1), class = Intercept),
                        prior(normal(0,0.5), class = b)),
              iter = 4000, chains = 4,
              backend="cmdstanr", cores=4)

#loop model over list of matrices
m2.soc_loops <- vector("list", 100) 
for (i in seq_along(m2.soc_loops)) {
  m2.soc_loops[[i]] <- update(m2.soc0,
                              data2 = list(cormat = ceph100cormat[[i]],
                                           backend = "cmdstanr",
                                           cores = 4)
  )
}

#combine models
m2.soc_comb <- combine_models(m2.soc_loops[[i]])
plot(conditional_effects(m2.soc_comb, effects = "sociality.bin", resp = "CNS1"), points = TRUE)
summary(m2.soc_comb)
save(m2.soc_comb, file="m2.soc100.rda")

## sensitivity check w/ 3 category sociality on consensus phylogeny----
class(st.cephdat$sociality.3)
m2.soc3 <- brm(bf(CNS.1 ~ mi(ML.1) + sociality.3 + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat)))+
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + (1|gr(phy.species, cov=cormat))), 
               data=st.cephdat, 
               data2=list(cormat=cormat), 
               family=gaussian("identity"),
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(normal(0, 1), class = b)),
               iter = 4000, chains = 4,
               backend="cmdstanr", cores=4)
summary(m2.soc3)
save(m2.soc3, file="m2soc3.rda") 

hypothesis(m2.soc3, hypothesis="CNS1_sociality.3.L < 0")

#Prediction 3: Behavioral repertoire----

#function to constrain imputation
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

#these use articles read as the proxy for research effort because this is count-based data of behavior more heavily confounded by # of observations
adjustmentSets(cephdag, outcome="CNS", exposure="behavior")
# { WoS, depth, latitude, phylogeny }

##main model: summed cognition and social behaviors----
min(na.omit(st.cephdat$cog2)) #lower bound for observed cog2 variable to constrain imputed values >1
m3.cog2.empty <- brm(bf(CNS.1 ~ mi(ML.1) +  mi(cog2) + benthic + articles.read + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat)))+
                 bf(cog2 | mi() ~ 1 + benthic + articles.read + (1|gr(phy.species,cov=cormat))),
               family=gaussian,
               data=st.cephdat,
               data2=list(cormat=cormat),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               backend="cmdstanr",chains=0)

m3.cog_cons <- constrained_imputation(model = m3.cog2.empty,
                                              vars = c("ML.1", "cog2"))
summary(m3.cog_cons)
save(m3.cog_cons, file="m3.cog_cons.rda")

## main model: defense repertoire----
m3.def.empty <- brm(bf(CNS.1 ~ mi(ML.1) +  mi(defense.repertoire) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                      bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat)))+
                      bf(defense.repertoire | mi() ~ 1 + benthic + articles.read + (1|gr(phy.species,cov=cormat))),
                    family=gaussian,
                    data=st.cephdat,
                    data2=list(cormat=cormat),
                    prior = c(prior(normal(0,1), class = Intercept),
                              prior(normal(0,0.5), class = b)),
                    backend="cmdstanr", chains = 0)

m3.def_cons <- constrained_imputation(model = m3.def.empty,
                                             vars = c("ML.1", "defense.repertoire"))

summary(m3.def_cons)
save(m3.def_cons, file="m3.def_cons.rda")

## main model: hunting repertoire----
m3.hunt.empty <- brm(bf(CNS.1 ~ mi(ML.1) +  mi(foraging.repertoire) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                       bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat)))+
                       bf(foraging.repertoire | mi() ~ 1 + benthic + articles.read + (1|gr(phy.species,cov=cormat))),
                     family=gaussian,
                     data=st.cephdat,
                     data2=list(cormat=cormat),
                     prior = c(prior(normal(0,1), class = Intercept),
                               prior(normal(0,0.5), class = b)),
                     backend="cmdstanr", chains=0)

m3.hunt_cons <- constrained_imputation(model = m3.hunt.empty,
                                              vars = c("ML.1", "foraging.repertoire"))

summary(m3.hunt_cons)
save(m3.hunt_cons, file="m3.hunt_cons.rda")

#Prediction 4: Ecological complexity----
adjustmentSets(cephdag, outcome="CNS", exposure="diet")
# { ML, WoS, benthic, depth, latitude, phylogeny }
## main model: diet breadth ----
min(na.omit(st.cephdat$diet.breadth))
m4.diet0 <- brm(bf(CNS.1 ~ mi(ML.1) + mi(diet.breadth) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                 bf(diet.breadth | mi() + resp_trunc(lb=-2.02) ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
               family=gaussian,
               data=st.cephdat,
               data2=list(cormat=ceph100cormat[[1]]),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               backend="cmdstanr", cores=4)

#loop model over list of matrices
m4.diet_loops <- vector("list", 100) 
for (i in seq_along(m4.diet_loops)) {
  m4.diet_loops[[i]] <- update(m4.diet0,
                               data2 = list(cormat = ceph100cormat[[i]],
                                            backend = "cmdstanr",
                                            cores = 4)
  )
}

#combine models
m4.diet_comb <- combine_models(m4.diet_loops[[i]])
save(m4.diet_comb, file="m4diet100.rda")
summary(m4.diet_comb)

## diet-benthic interaction consensus phylogeny----
m4.dhab0 <- brm(bf(CNS.1 ~ mi(ML.1) + mi(diet.breadth)*benthic + depth.mean + pos.latmean + articles.read + (1|gr(phy.species,cov=cormat))) +
                 bf(diet.breadth | mi() + resp_trunc(lb=-2.1) ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
               family=gaussian,
               data=st.cephdat,
               data2=list(cormat=ceph100cormat[[1]]),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               backend="cmdstanr", cores=4)

#loop model over list of matrices
m4.dhab_loops <- vector("list", 100) 
for (i in seq_along(m4.dhab_loops)) {
  m4.dhab_loops[[i]] <- update(m4.dhab0,
                               data2 = list(cormat = ceph100cormat[[i]],
                                            backend = "cmdstanr",
                                            cores = 4)
  )
}

#combine models
m4.dhab_comb <- combine_models(m4.dhab_loops[[i]])
save(m4.dhab_comb, file="m4.dhab100.rda")
summary(m4.dhab_comb)


## main model: predator breadth ----
adjustmentSets(cephdag, outcome="CNS", exposure="predators")
# { ML, WoS, benthic, depth, latitude }
min(na.omit(st.cephdat$predator.breadth))
m5.preds0 <- brm(bf(CNS.1 ~ mi(ML.1) + mi(predator.breadth) + depth.mean + pos.latmean + benthic + articles.read + (1 | gr(phy.species, cov = cormat))) +
                  bf(predator.breadth | mi() + resp_trunc(lb=-1.55) ~ 1 + mi(ML.1) + depth.mean + pos.latmean + benthic + articles.read + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))),
                family = gaussian,
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = st.cephdat,
                data2 = list(cormat = ceph100cormat[[1]]),
                backend = "cmdstanr", cores = 4)

#loop model over list of matrices
m5.preds_loops <- vector("list", 100) 
for (i in seq_along(m5.preds_loops)) {
  m5.preds_loops[[i]] <- update(m5.preds0,
                                data2 = list(cormat = ceph100cormat[[i]],
                                             backend = "cmdstanr",
                                             cores = 4)
  )
}

#combine models
m5.preds_comb <- combine_models(m5.preds_loops[[i]])
summary(m5.preds_comb)
plot(conditional_effects(m5.preds, effects = "predator.breadth", resp = "CNS1"), points = TRUE)

## main model: predator breadth habitat interaction----
m5.phab <- brm(bf(CNS.1 ~ mi(ML.1) + mi(predator.breadth)*benthic + articles.read + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))) +
                 bf(predator.breadth | mi() + resp_trunc(lb=-1.6) ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))),
               family = gaussian,
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               data = st.cephdat,
               data2 = list(cormat=ceph100cormat[[1]]),
               backend = "cmdstanr", cores = 4)

#loop model over list of matrices
m5.phab_loops <- vector("list", 100) 
for (i in seq_along(m5.phab_loops)) {
  m5.phab_loops[[i]] <- update(m5.phab,
                               data2 = list(cormat = ceph100cormat[[i]],
                                            backend = "cmdstanr",
                                            cores = 4)
  )
}

#combine models
m5.phab_comb <- combine_models(m5.phab_loops[[i]])
save(m5.phab_comb, file="m5.phab100.rda")

## main model: habitat----
m5.hab0 <- brm(bf(CNS.1 ~ mi(ML.1) + habitat3 + WoS + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML.1 | mi() ~ 1 + WoS + (1 | gr(phy.species, cov = cormat))),
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               data = st.cephdat,
               data2 = list(cormat=ceph100cormat[[1]]),
               backend = "cmdstanr",
               cores = 4)

#loop model over list of matrices
m5.hab_loops <- vector("list", 100) 
for (i in seq_along(m5.hab_loops)) {
  m5.hab_loops[[i]] <- update(m5.hab0,
                               data2 = list(cormat = ceph100cormat[[i]],
                                            backend = "cmdstanr",
                                            cores = 4)
  )
}

#combine models
m5.hab_comb <- combine_models(m5.hab_loops[[i]])
summary(m5.hab_comb)
conditional_effects(m5.hab_comb, effects = "habitat3", resp = "CNS1")
save(m5.hab_comb, file="m5.hab100.rda")

## sensitivity check: benthic vs. pelagic----
m5.benth <- brm(bf(CNS.1 ~ mi(ML.1) + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML.1 | mi() ~ 1 + WoS + (1 | gr(phy.species, cov = cormat))),
                prior = c(prior(normal(0, 1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = st.cephdat,
                data2 = list(cormat = cormat),
                iter=4000,
                backend = "cmdstanr",
                cores = 4)
summary(m5.benth)
save(m5.benth, file="m5benth.rda")

## latitude range on consensus phylogeny----
adjustmentSets(cephdag, exposure="latitude", outcome="CNS") 
# { WoS, benthic, phylogeny }
m5.latrange <- brm(bf(CNS.1 ~ mi(ML.1) + lat.range + WoS + benthic + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML.1 | mi() ~ 1 + lat.range + benthic + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = st.cephdat,
                   data2 = list(cormat=cormat),
                   backend = "cmdstanr",
                   cores = 4)
summary(m5.latrange)
save(m5.latrange, file="m5.latrange.rda")

load(file="m5.latrange.rda")

# mean latitude on consensus phylogeny----
m5.latmean <- brm(bf(CNS.1 ~ mi(ML.1) + pos.latmean + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                    bf(ML.1 | mi() ~ 1 + pos.latmean + benthic + (1 | gr(phy.species, cov = cormat))),
                  prior = c(prior(normal(0,1), class = Intercept),
                            prior(normal(0,0.5), class = b)),
                  data = st.cephdat,
                  data2 = list(cormat=cormat),
                  backend = "cmdstanr",
                  cores = 4)
summary(m5.latmean)
save(m5.latmean, file="m5.latmean.rda")

## minimum depth on consensus phylogeny----
m5.mindepth <- brm(bf(CNS.1 ~ mi(ML.1) + depth.min + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML.1 | mi() ~ 1 + depth.min  + benthic + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = st.cephdat,
                   data2 = list(cormat=cormat),
                   backend = "cmdstanr",
                   iter=4000, chains=4,
                   cores = 4)
summary(m5.mindepth) 
save(m5.mindepth, file="m5mindepth.rda")

#max depth on consensus phylogeny----
m5.maxdepth <- brm(bf(CNS.1 ~ mi(ML.1) + depth.max + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = st.cephdat,
                   data2 = list(cormat = cormat),
                   backend = "cmdstanr",
                   cores = 4)
plot(conditional_effects(m5.maxdepth, effects = "depth.max", resp = "CNS1"), points = TRUE)
summary(m5.maxdepth)
save(m5.maxdepth, file="m5.maxdepth.rda")

## main model: minimum depth interaction w/ habitat----
m5.mindb0 <- brm(bf(CNS.1 ~ mi(ML.1) + depth.min*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = st.cephdat,
                data2 = list(cormat=ceph100cormat[[1]]),
                backend = "cmdstanr",
                cores = 4)

#loop model over list of matrices
m5.mindb_loops <- vector("list", 100) 
for (i in seq_along(m5.mindb_loops)) {
  m5.mindb_loops[[i]] <- update(m5.mindb0,
                                data2 = list(cormat = ceph100cormat[[i]],
                                             backend = "cmdstanr",
                                             cores = 4)
  )
}

#combine models
m5mindb_comb <- combine_models(m5.mindb_loops[[i]])
save(m5mindb_comb, file="m5mindb100.rda")
summary(m5mindb_comb)

load(file="m5mindb100.rda")

## main model: maximum depth interaction w/ habitat----
m5.maxdb0 <- brm(bf(CNS.1 ~ mi(ML.1) + depth.max*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = st.cephdat,
                   data2 = list(cormat = ceph100cormat[[1]]),
                   backend = "cmdstanr",
                   cores = 4)

#loop model over list of matrices
m5.maxdb_loops <- vector("list", 100) 
for (i in seq_along(m5.maxdb_loops)) {
  m5.maxdb_loops[[i]] <- update(m5.maxdb0,
                                data2 = list(cormat = ceph100cormat[[i]],
                                             backend = "cmdstanr",
                                             cores = 4)
  )
}

#combine models
m5maxdb_comb <- combine_models(m5.maxdb_loops[[i]])
save(m5maxdb_comb, file="m5maxdb100.rda")
summary(m5maxdb_comb)


