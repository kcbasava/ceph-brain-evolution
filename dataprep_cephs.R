library(ape)
library(ggdag)
library(dagitty)

setwd("/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses")
cephdat <- read.csv("cephdat.csv")
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
st.cephdat <- cephdat
View(st.cephdat[c(9:35)])
st.cephdat[c(9:35)] <- lapply(st.cephdat2[c(9:35)], function(x) round(standardize(log(x)), digits=3))
st.cephdat <- st.cephdat[,-c(13,14)] #remove original latitude columns with negative values

write.csv(st.cephdat, "st.cephdat.csv", row.names=FALSE)

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