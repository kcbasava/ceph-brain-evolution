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
m1.smmax <- brm(bf(CNS.1 ~ mi(ML.1) + mi(lifespan.max) + mi(matage.max) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(matage.max | mi() ~ 1 + mi(ML.1) + mi(lifespan.max) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                  bf(lifespan.max | mi() ~ 1 + mi(ML.1) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
                data=st.cephdat, 
                data2=list(cormat=ceph100cormat[[1]]), 
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                iter = 2000, warmup = 1000, chains = 4,
                backend="cmdstanr", cores=4)

#loop model over list of matrices
m1.smmax_loops <- vector("list", 100) 
for (i in seq_along(m1.smmax_loops)) {
  m1.smmax_loops[[i]] <- update(m1.smmax,
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

## sensitivity checks with consensus phylogeny----
m1.sm <- brm(bf(CNS.1 ~ mi(ML.1) + mi(lifespan.mean) + mi(matage.mean) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
               bf(matage.mean | mi() ~ 1 + mi(ML.1) + mi(lifespan.mean) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
               bf(lifespan.mean | mi() ~ 1 + mi(ML.1) + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
               bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species, cov=cormat))), 
             data=st.cephdat, 
             data2=list(cormat=cormat), 
             prior = c(prior(normal(0,1), class = Intercept),
                       prior(normal(0,0.5), class = b)),
             iter = 2000, warmup = 1000, chains = 4,
             backend="cmdstanr", cores=4)
summary(m1.sm)

#Prediction # 2: Sociality-----
adjustmentSets(cephdag, outcome="CNS", exposure="sociality")
# { WoS, benthic, depth, phylogeny }

## main model: binary sociality----
m2.soc <- brm(bf(CNS.1 ~ mi(ML.1) + sociality.bin + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                bf(ML.1 | mi() ~ 1 + benthic + depth.mean + (1|gr(phy.species, cov=cormat))), 
              data=st.cephdat, 
              data2=list(cormat=ceph100_cormat[[1]]), 
              prior = c(prior(normal(0,1), class = Intercept),
                        prior(normal(0,0.5), class = b)),
              iter = 2000, warmup = 1000, chains = 4,
              backend="cmdstanr", cores=4)

#loop model over list of matrices
m2.soc_loops <- vector("list", 100) 
for (i in seq_along(m2.soc_loops)) {
  m2.soc_loops[[i]] <- update(m2.soc,
                              data2 = list(cormat = ceph100_cormat[[i]],
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
               iter = 2000, warmup = 1000, chains = 4,
               backend="cmdstanr", cores=4)
summary(m2.soc3)
hypothesis(m2.soc3, hypothesis="CNS1_sociality.3.L < 0")

#Prediction 3: Behavioral repertoire----
#these use articles read as the proxy for research effort because this is count-based data of behavior more heavily confounded by # of observations
adjustmentSets(cephdag, outcome="CNS", exposure="behavior")
# { WoS, depth, latitude, phylogeny }

##main model: summed cognition and social behaviors----
min(na.omit(st.cephdat$cog2)) #lower bound for observed cog2 variable to constrain imputed values >1
m3.cog2 <- brm(bf(CNS.1 ~ mi(ML.1) +  mi(cog2) + benthic + articles.read + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat)))+
                 bf(cog2 | mi() + resp_trunc(lb=-0.963) ~ 1 + benthic + articles.read + (1|gr(phy.species,cov=cormat))),
               family=gaussian,
               data=st.cephdat,
               data2=list(cormat=ceph100_cormat[[1]]),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               backend="cmdstanr", cores=4)

#loop model over list of matrices
m3.cog2_loops <- vector("list", 100) 
for (i in seq_along(m3.cog2_loops)) {
  m3.cog2_loops[[i]] <- update(m3.cog2,
                                 data2 = list(cormat = ceph100_cormat[[i]],
                                              backend = "cmdstanr",
                                              cores = 4)
  )
}

#combine models
m3.cog2_comb <- combine_models(m3.cog2_loops[[i]])
summary(m3.cog2_comb)
save(m3.cog2_comb, file="m3.antipr100.rda")

## main model: defense repertoire----
min(na.omit(st.cephdat$defense.repertoire))
m3.def <- brm(bf(CNS.1 ~ mi(ML.1) +  mi(defense.repertoire) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                   bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat)))+
                   bf(defense.repertoire | mi() ~ 1 + benthic + articles.read + (1|gr(phy.species,cov=cormat))),
                 family=gaussian,
                 data=st.cephdat,
                 data2=list(cormat=ceph100_cormat[[1]]),
                 prior = c(prior(normal(0,1), class = Intercept),
                           prior(normal(0,0.5), class = b)),
                 backend="cmdstanr", cores=4)

#loop model over list of matrices
m3.def_loops <- vector("list", 100) 
for (i in seq_along(m3.def_loops)) {
  m3.def_loops[[i]] <- update(m3.def,
                                 data2 = list(cormat = ceph100_cormat[[i]],
                                              backend = "cmdstanr",
                                              cores = 4)
  )
}

#combine models
m3.def_comb <- combine_models(m3.def_loops[[i]])
summary(m3.def_comb)
save(m3.def_comb, file="m3.def100.rda")

## main model: hunting repertoire----
min(na.omit(st.cephdat$foraging.repertoire))
m3.hunt <- brm(bf(CNS.1 ~ mi(ML.1) +  mi(foraging.repertoire) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat)))+
                 bf(foraging.repertoire | mi() + resp_trunc(lb=-1.3) ~ 1 + benthic + articles.read + (1|gr(phy.species,cov=cormat))),
               family=gaussian,
               data=st.cephdat,
               data2=list(cormat=ceph100_cormat[[1]]),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               backend="cmdstanr", cores=4)
summary(m3.hunt)

#loop model over list of matrices
m3.hunt_loops <- vector("list", 100) 
for (i in seq_along(m3.hunt_loops)) {
  m3.hunt_loops[[i]] <- update(m3.hunt,
                               data2 = list(cormat = ceph100_cormat[[i]],
                                            backend = "cmdstanr",
                                            cores = 4)
  )
}

#combine models
m3.hunt_comb <- combine_models(m3.hunt_loops[[i]])
plot(conditional_effects(m3.hunt_comb, effects = "foraging.repertoire", resp = "CNS1"), points = TRUE)
save(m3.hunt_comb, file="m3.hunt100.rda")

#Prediction 4: Ecological complexity----
adjustmentSets(cephdag, outcome="CNS", exposure="diet")
# { ML, WoS, benthic, depth, latitude, phylogeny }
## main model: diet breadth ----
min(na.omit(st.cephdat$diet.breadth))
m4.diet <- brm(bf(CNS.1 ~ mi(ML.1) + mi(diet.breadth) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                 bf(diet.breadth | mi() + resp_trunc(lb=-2.02) ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
               family=gaussian,
               data=st.cephdat,
               data2=list(cormat=ceph100_cormat[[1]]),
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               backend="cmdstanr", cores=4)

#loop model over list of matrices
m4.diet_loops <- vector("list", 100) 
for (i in seq_along(m4.diet_loops)) {
  m4.diet_loops[[i]] <- update(m4.diet,
                               data2 = list(cormat = ceph100_cormat[[i]],
                                            backend = "cmdstanr",
                                            cores = 4)
  )
}

#combine models
m4.diet_comb <- combine_models(m4.diet_loops[[i]])
save(m4.diet_comb, file="m4diet100.rda")
summary(m4.diet_comb)

## sensitivity check: diet habitat interaction----
m4.db <- brm(bf(CNS.1 ~ mi(ML.1) + mi(diet.breadth)*benthic + depth.mean + pos.latmean + articles.read + (1|gr(phy.species,cov=cormat))) +
               bf(diet.breadth | mi() + resp_trunc(lb=-2.02) ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
               bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
             family=gaussian,
             data=st.cephdat,
             data2=list(cormat=cormat),
             prior = c(prior(normal(0,1), class = Intercept),
                       prior(normal(0,0.5), class = b)),
             backend="cmdstanr", cores=4)
summary(m4.db)

## main model: predator breadth ----
adjustmentSets(cephdag, outcome="CNS", exposure="predators")
# { ML, WoS, benthic, depth, latitude }
min(na.omit(st.cephdat$predator.breadth))
m5.preds <- brm(bf(CNS.1 ~ mi(ML.1) + mi(predator.breadth) + depth.mean + pos.latmean + benthic + articles.read + (1 | gr(phy.species, cov = cormat))) +
                  bf(predator.breadth | mi() + resp_trunc(lb=-1.55) ~ 1 + mi(ML.1) + depth.mean + pos.latmean + benthic + articles.read + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))),
                family = gaussian,
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = st.cephdat,
                data2 = list(cormat = ceph100_cormat[[1]]),
                backend = "cmdstanr", cores = 4)

#loop model over list of matrices
m5.preds_loops <- vector("list", 100) 
for (i in seq_along(m5.preds_loops)) {
  m5.preds_loops[[i]] <- update(m5.preds,
                                data2 = list(cormat = ceph100_cormat[[i]],
                                             backend = "cmdstanr",
                                             cores = 4)
  )
}

#combine models
m5.preds_comb <- combine_models(m5.preds_loops[[i]])
summary(m5.preds_comb)
plot(conditional_effects(m5.preds, effects = "predator.breadth", resp = "CNS1"), points = TRUE)

## sensitivity check: predator breadth habitat interaction----
m5.phab <- brm(bf(CNS.1 ~ mi(ML.1) + mi(predator.breadth)*benthic + articles.read + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))) +
                 bf(predator.breadth | mi() + resp_trunc(lb=-1.55) ~ 1 + mi(ML.1) + articles.read + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))),
               family = gaussian,
               prior = c(prior(normal(0,0.5), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               iter = 2000, warmup = 1000, chains = 4,
               data = st.cephdat,
               data2 = list(cormat=cormat),
               backend = "cmdstanr", cores = 4)
summary(m5.phab) 

## main model: habitat----
m5.hab3 <- brm(bf(CNS.1 ~ mi(ML.1) + habitat3 + WoS + (1 | gr(phy.species, cov = cormat))) +
                 bf(ML.1 | mi() ~ 1 + WoS + (1 | gr(phy.species, cov = cormat))),
               prior = c(prior(normal(0, 1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               data = st.cephdat,
               data2 = list(cormat=ceph100_cormat[[1]]),
               backend = "cmdstanr",
               cores = 4)
summary(m5.hab3)

#loop model over list of matrices
m5.hab3_loops <- vector("list", 100) 
for (i in seq_along(m5.hab3_loops)) {
  m5.hab3_loops[[i]] <- update(m5.hab3,
                               data2 = list(cormat = ceph100_cormat[[i]],
                                            backend = "cmdstanr",
                                            cores = 4)
  )
}

#combine models
m5.hab3_comb <- combine_models(m5.hab3_loops[[i]])
summary(m5.hab3_comb)
conditional_effects(m5.hab3_comb, effects = "habitat3", resp = "CNS1")
save(m5.hab3_comb, file="m5.hab100.rda")

## main model: latitude range ----
adjustmentSets(cephdag, exposure="latitude", outcome="CNS") 
# { WoS, benthic, phylogeny }
m5.latrange <- brm(bf(CNS.1 ~ mi(ML.1) + lat.range + WoS + benthic + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML.1 | mi() ~ 1 + lat.range + benthic + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = st.cephdat,
                   data2 = list(cormat=ceph100_cormat[[1]]),
                   backend = "cmdstanr",
                   cores = 4)

#loop model over list of matrices
m5.latrange_loops <- vector("list", 100) 
for (i in seq_along(m5.latrange_loops)) {
  m5.latrange_loops[[i]] <- update(m5.latrange,
                                   data2 = list(cormat = ceph100_cormat[[i]],
                                                backend = "cmdstanr",
                                                cores = 4)
  )
}

#combine models
m5.latrange_comb <- combine_models(m5.latrange_loops[[i]])
save(m5.latrange_comb, file="m5latrange100.rda")
summary(m5.latrange_comb)

#sensitvity check: mean latitude
m5.latmean <- brm(bf(CNS.1 ~ mi(ML.1) + pos.latmean + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                    bf(ML.1 | mi() ~ 1 + pos.latmean + benthic + (1 | gr(phy.species, cov = cormat))),
                  prior = c(prior(normal(0,1), class = Intercept),
                            prior(normal(0,0.5), class = b)),
                  data = st.cephdat,
                  data2 = list(cormat=cormat),
                  backend = "cmdstanr",
                  cores = 4)
summary(m5.latmean)

## main model: minimum depth----
m5.mindepth <- brm(bf(CNS.1 ~ mi(ML.1) + depth.min + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML.1 | mi() ~ 1 + depth.min  + benthic + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = st.cephdat,
                   data2 = list(cormat=ceph100_cormat[[1]]),
                   backend = "cmdstanr",
                   cores = 4)

#loop model over list of matrices
m5.min_loops <- vector("list", 100) 
for (i in seq_along(m5.min_loops)) {
  m5.min_loops[[i]] <- update(m5.mindepth,
                              data2 = list(cormat = ceph100_cormat[[i]],
                                           backend = "cmdstanr",
                                           cores = 4)
  )
}

#combine models
m5mindepth_comb <- combine_models(m5.min_loops[[i]])
summary(m5mindepth_comb)
save(m5mindepth_comb, "m5mindepth100.rda")

## main model: minimum depth interaction w/ habitat----
m5.mindb <- brm(bf(CNS.1 ~ mi(ML.1) + depth.min*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                  bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data = st.cephdat,
                data2 = list(cormat=ceph100_cormat[[1]]),
                backend = "cmdstanr",
                cores = 4)

#loop model over list of matrices
m5.mindb_loops <- vector("list", 100) 
for (i in seq_along(m5.mindb_loops)) {
  m5.mindb_loops[[i]] <- update(m5.mindb,
                                data2 = list(cormat = ceph100_cormat[[i]],
                                             backend = "cmdstanr",
                                             cores = 4)
  )
}

#combine models
m5mindb_comb <- combine_models(m5.mindb_loops[[i]])
save(m5mindb_comb, file="m5mindb100.rda")
summary(m5mindb_comb)

## main model: maximum depth interaction w/ habitat----
m5.maxdb <- brm(bf(CNS.1 ~ mi(ML.1) + depth.max*benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                     bf(ML.1 | mi() ~ 1 + (1 | gr(phy.species, cov = cormat))),
                   prior = c(prior(normal(0,1), class = Intercept),
                             prior(normal(0,0.5), class = b)),
                   data = st.cephdat,
                   data2 = list(cormat = ceph100_cormat[[i]]),
                   backend = "cmdstanr",
                   cores = 4)

#loop model over list of matrices
m5.maxdb_loops <- vector("list", 100) 
for (i in seq_along(m5.maxdb_loops)) {
  m5.maxdb_loops[[i]] <- update(m5.maxdb,
                                data2 = list(cormat = ceph100_cormat[[i]],
                                             backend = "cmdstanr",
                                             cores = 4)
  )
}

#combine models
m5maxdb_comb <- combine_models(m5.maxdb_loops[[i]])
save(m5maxdb_comb, file="m5maxdb100.rda")
summary(m5maxdb_comb)


