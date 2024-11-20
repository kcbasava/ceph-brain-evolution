#load packages----
library(brms)
library(ape)
library(tidyverse)
library(cmdstanr)
library(ggdag)
library(dagitty)

options(scipen=999) #turn off scientific notation

#read in data and phylogenies
cephtreeMCC <- read.nexus("cephtreeMCC.tree") #consensus (maximum clade credibility) tree
cormat <- vcv(cephtreeMCC, corr=TRUE) #correlation matrix for consensus phylogeny

#the following code can be run to reproduce the standardized version of the dataset with additional CNS and ML measures----
#or can be read in as a csv below

#cephdat2 <- read.csv("cephdat2.csv")

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
# View(st.cephdat2[c(9:39)])
# st.cephdat2[c(9:39)] <- lapply(st.cephdat2[c(9:39)], function(x) round(standardize(log(x)), digits=3))
# st.cephdat2 <- st.cephdat2[,-c(13,14)] #remove original latitude columns with negative values
# 
# df <- st.cephdat2
# 
# #if CNS.2 and CNS.3 are both present, take the mean; else fill in CNS.2
# df$CNS.23 <- apply(df[c(34,35)], 1, function(row) {
#   if (is.na(row[1]) & !is.na(row[2])) {
#     return(row[2])
#   } else if (!is.na(row[1]) & is.na(row[2])) {
#     return(row[1])
#   } else if (is.na(row[1]) & is.na(row[2])) {
#     return(NA)
#   } else {
#     return(mean(c(row[1], row[2])))
#   }
# })
# 
# df$ML.23 <- apply(df[c(36,37)], 1, function(row) {
#   if (is.na(row[1]) & !is.na(row[2])) {
#     return(row[2])
#   } else if (!is.na(row[1]) & is.na(row[2])) {
#     return(row[1])
#   } else if (is.na(row[1]) & is.na(row[2])) {
#     return(NA)
#   } else {
#     return(mean(c(row[1], row[2])))
#   }
# })
# write.csv(df, file="st.cephdat2.csv", row.names=FALSE)

st.cephdat2 <- read.csv("st.cephdat2.csv")
st.cephdat2 <- st.cephdat2[!is.na(st.cephdat2$CNS.23),] #remove species without additional brain measures

st.cephdat2$benthic <- factor(st.cephdat2$benthic)
st.cephdat2$sociality.bin <- factor(st.cephdat2$sociality.bin)
st.cephdat2$sociality.3 <- factor(st.cephdat2$sociality.3)
st.cephdat2$habitat3 <- factor(st.cephdat2$habitat3)
st.cephdat2$depth_cat <- factor(st.cephdat2$depth_cat)


#with alternate CNS measures----

## maximum age of sexual maturity ----
m1.smmax2 <- brm(bf(CNS.23 ~ ML.23 + mi(lifespan.max) + mi(matage.max) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                   bf(matage.max | mi() ~ 1 + ML.23 + mi(lifespan.max) + WoS + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) + 
                   bf(lifespan.max | mi() ~ 1 + ML.23 + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))),
                 data=st.cephdat2, 
                 data2=list(cormat=cormat), 
                 prior = c(prior(normal(0,1), class = Intercept),
                           prior(normal(0,0.5), class = b)),
                 backend="cmdstanr", cores=4)
summary(m1.smmax2)
# somewhat positive 

## sociality----
m2.soc2 <- brm(bf(CNS.23 ~ ML.23 + sociality.bin + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))),              data=st.cephdat2, 
               data2=list(cormat=cormat), 
               prior = c(prior(normal(0,1), class = Intercept),
                         prior(normal(0,0.5), class = b)),
               backend="cmdstanr", cores=4)
summary(m2.soc2)
# still mostly negative

## habitat----
m5.benth2 <- brm(bf(CNS.23 ~ ML.23 + benthic + WoS + (1 | gr(phy.species, cov = cormat))),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0,0.5), class = b)),
                 data = st.cephdat2,
                 iter=4000, chains=4,
                 data2 = list(cormat = cormat),
                 backend = "cmdstanr",
                 cores = 4)
summary(m5.benth2)
#less but still mostly positive

m5.preds2 <- brm(bf(CNS.23 ~ ML.23 + mi(predator.breadth)*benthic + depth.mean + pos.latmean + benthic + WoS + (1 | gr(phy.species, cov = cormat))) +
                   bf(predator.breadth | mi() ~ 1 + ML.23 + depth.mean + pos.latmean + benthic + WoS + (1 | gr(phy.species, cov = cormat))),
                 family = gaussian,
                 prior = c(prior(normal(0,1), class = Intercept),
                           prior(normal(0,0.5), class = b)),
                 data = st.cephdat2,
                 data2 = list(cormat = cormat),
                 backend = "cmdstanr", cores = 4)
summary(m5.preds2)
# negative instead of positive

#maximum depth----
m5.maxdepth2 <- brm(bf(CNS.23 ~ ML.23 + depth.max + benthic + WoS + (1 | gr(phy.species, cov = cormat))),                   
                    prior = c(prior(normal(0,1), class = Intercept),                                                                                                                  prior(normal(0,0.5), class = b)),
                    data = st.cephdat2,
                    data2 = list(cormat = cormat),
                    backend = "cmdstanr",
                    iter=4000, chains=4,
                    cores = 4)
summary(m5.maxdepth2)
#still pretty negative

#only with complete cases for focal variable----
## foraging/hunting repertoire----
mcc.hunt <- brm(bf(CNS.1 ~ mi(ML.1) + foraging.repertoire + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
                family=gaussian,
                data=st.cephdat,
                data2=list(cormat=cormat),
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                backend="cmdstanr", cores=4)
summary(mcc.hunt)
#still negative 

## defense repertoire----
mcc.antipr <- brm(bf(CNS.1 ~ mi(ML.1) +  defense.repertoire + articles.read + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                    bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
                  family=gaussian,
                  data=st.cephdat,
                  data2=list(cormat=cormat),
                  prior = c(prior(normal(0,1), class = Intercept),
                            prior(normal(0,0.5), class = b)),
                  backend="cmdstanr", cores=4)
# still a little positive

## diet----
mcc.diet <- brm(bf(CNS.1 ~ mi(ML.1) + diet.breadth*benthic + articles.read + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))) +
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1|gr(phy.species,cov=cormat))),
                family=gaussian,
                data=st.cephdat,
                data2=list(cormat=cormat),
                iter=4000, chains=4,
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                backend="cmdstanr", cores=4)
summary(mcc.diet)
#slightly negative for pelagic + slightly positive for benthic but basically 0 for both

mcc.preds <- brm(bf(CNS.1 ~ mi(ML.1) + predator.breadth*benthic + depth.mean + pos.latmean + benthic + articles.read + (1 | gr(phy.species, cov = cormat))) +
                   bf(ML.1 | mi() ~ 1 + benthic + depth.mean + pos.latmean + (1 | gr(phy.species, cov = cormat))),
                 family = gaussian,
                 prior = c(prior(normal(0,1), class = Intercept),
                           prior(normal(0,0.5), class = b)),
                 data = st.cephdat,
                 data2 = list(cormat = cormat),
                 backend = "cmdstanr",
                 cores = 4)
summary(mcc.preds)
#slightly positive for both, more so with benthic 

#sociality within decapodiformes----
alltaxa <- read.csv("/Users/kiranbasava/nonhumans/di_cephproject/phylos/alltaxa.csv")

decadat <- merge(st.cephdat, alltaxa, by="phy.species")
decadat <- decadat[decadat$Superorder=="Decapodiformes",]

m.socdec <- brm(bf(CNS.1 ~ mi(ML.1) + sociality.bin + WoS + benthic + depth.mean + (1|gr(phy.species,cov=cormat))) + 
                  bf(ML.1 | mi() ~ 1 + benthic + depth.mean + (1|gr(phy.species, cov=cormat))), 
                prior = c(prior(normal(0,1), class = Intercept),
                          prior(normal(0,0.5), class = b)),
                data=decadat, 
                data2=list(cormat=cormat), 
                backend="cmdstanr", cores=4)
summary(m.socdec)
#still negative
