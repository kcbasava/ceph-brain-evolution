library(ape)
library(dagitty)

#updating with new depth data november 19, 2024
cephdat <- read.csv("/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/cephdat.csv")
View(cephdat)

#standardize function from rethinking (https://rdrr.io/github/rmcelreath/rethinking/src/R/utilities.r)
standardize <- function(x) {
  x <- scale(x)
  z <- as.numeric(x)
  attr(z,"scaled:center") <- attr(x,"scaled:center")
  attr(z,"scaled:scale") <- attr(x,"scaled:scale")
  return(z)
}

#log and standardize continuous vars----
#I manually changed Argonaut, Chiroteuthis, and Megalocranchia eq.dist to 0.000001 instead of 0 in cephdat.csv so can log
#and to updated depth measures of 0 
st.cephdat <- cephdat
View(st.cephdat[c(7:30)])
st.cephdat[c(7:30)] <- lapply(st.cephdat[c(7:30)], function(x) round(standardize(log(x)), digits=3))
st.cephdat <- st.cephdat[,-c(12,13)] #remove original latitude columns with negative values

write.csv(st.cephdat, "/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution/st.cephdat.csv", row.names=FALSE)

cephtrees100 <- read.nexus("cephtrees100.trees")

#list of correlation matrices for 100 trees----
ceph100cormat <-list()
for (i in 1:length(cephtrees100) ) {
  ceph100cormat[[i]] <- vcv.phylo(cephtrees100[[i]], corr=TRUE)
}

save(ceph100cormat, file="ceph100_cormat.rda")

#DAG code----
cephdag <- dagitty("dag {
ML -> CNS
behavior -> CNS
benthic -> CNS
depth -> CNS
diet -> CNS
latitude -> CNS
lifespan -> CNS
maturity -> CNS
phylogeny -> CNS
predators -> CNS
sociality -> CNS
latitude -> ML
phylogeny -> ML
depth -> ML
latitude -> behavior
WoS -> behavior
phylogeny -> behavior
depth -> behavior
phylogeny -> lifespan
predators -> lifespan
WoS -> lifespan
ML -> lifespan
benthic -> lifespan
depth -> maturity
benthic -> maturity
WoS -> maturity
latitude -> maturity
lifespan -> maturity
phylogeny -> maturity
ML -> maturity
WoS -> diet
ML -> diet
benthic -> diet
depth -> diet
latitude -> diet
phylogeny -> diet
ML -> predators
WoS -> predators
depth -> predators
latitude -> predators
benthic -> predators
benthic -> sociality
depth -> sociality
WoS -> sociality
phylogeny -> sociality
phylogeny -> latitude
WoS -> latitude
benthic -> latitude
WoS -> depth
benthic -> depth
phylogeny -> depth
phylogeny -> benthic
}")
save(cephdag, file="cephdag.rda")


cephdag
adjustmentSets(cephdag, outcome="CNS", exposure="depth")
#{ WoS, benthic, phylogeny }

adjustmentSets(cephdag, outcome="CNS", exposure="benthic")
#{phylogeny}

adjustmentSets(cephdag, outcome="CNS", exposure="latitude")
#{ WoS, benthic, phylogeny }

adjustmentSets(cephdag, outcome="CNS", exposure="diet")
#{ ML, WoS, benthic, depth, latitude, phylogeny }

adjustmentSets(cephdag, outcome="CNS", exposure="sociality")
#{ WoS, benthic, depth, phylogeny }

adjustmentSets(cephdag, outcome="CNS", exposure="maturity")
#{ ML, WoS, benthic, depth, latitude, lifespan, phylogeny }

adjustmentSets(cephdag, outcome="CNS", exposure="behavior") #encompasses cognition, defense, and foraging
#{ WoS, depth, latitude, phylogeny }

