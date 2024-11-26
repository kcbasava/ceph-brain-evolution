#load packages
library(dplyr)
library(tibble)
library(stringr)
library(modelsummary)
library(viridis)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(ggplot2)
library(brms)
library(ape)
library(mice)

getwd()
setwd("/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/ceph-brain-evolution")
st.cephdat <- read.csv("st.cephdat.csv") #logged and standardized data
st.cephdat$benthic <- factor(st.cephdat$benthic)
st.cephdat$sociality.bin <- factor(st.cephdat$sociality.bin)
st.cephdat$sociality.3 <- factor(st.cephdat$sociality.3)
st.cephdat$habitat3 <- factor(st.cephdat$habitat3)
st.cephdat$depth_cat <- factor(st.cephdat$depth_cat)

cephdat <- read.csv("cephdat.csv")
cephdat$benthic <- as.factor(cephdat$benthic)
cephdat$sociality.bin <- as.factor(cephdat$sociality.bin)
cephdat$sociality.3 = as.factor(cephdat$sociality.3)
cephdat$habitat3 = as.factor(cephdat$habitat3)

options(scipen=999) #turn off scientific notation

load(file="nov2024_fits/m.hab.rda") #3-category habitat
load(file="nov2024_fits/m.benth.rda") #binary habitat

load(file="nov2024_fits/m.dhab.rda") #diet-habitat
load(file="nov2024_fits/m.phab.rda") #predator breadth-habitat
load(file="nov2024_fits/m.mindb.rda") 
load(file="nov2024_fits/m.maxdb.rda")

load(file="nov2024_fits/m.smmax.rda") #max age maturity
load(file="nov2024_fits/m.smmin.rda") #min age maturity
load(file="nov2024_fits/m.soc.rda") #sociality binary 
load(file="nov2024_fits/m.soc3.rda") #sociality 3-category

#tb_imp plot for interactions----
#imputation function
impute <- function(data, vars, model) {
  
  # Match brms var names
  ymi_vars <- gsub("[\\._]", "", vars)
  
  for (i in seq_along(vars)) {
    # Extract imputed values from brms model
    imputed_values <- brms::as_draws_df(model, variable = paste0("^Ymi_", ymi_vars[i]), regex = TRUE)
    
    # Remove unnecessary columns
    imputed_values$.chain <- NULL
    imputed_values$.iteration <- NULL
    imputed_values$.draw <- NULL
    
    # Identify missing indices in the dataset
    missing_indices <- which(is.na(data[[vars[i]]]))
    
    # Compute median of imputed values
    imputed_median <- apply(imputed_values[missing_indices, ], 2, median)
    
    # Create column indicating missing values, useful for plotting
    data[[paste0(vars[i], "_imputed")]] <- ifelse(is.na(data[[vars[i]]]), 1, 0)
    
    # Impute missing values in the dataset
    data[missing_indices, ][[vars[i]]] <- imputed_median
    
  }
  return(data)
}

#Figure 4 (main paper) flag for TB unable to remake (points shifted to the right)----
## 3 category habitat----
h_imp <- st.cephdat

# imputed mantle length: extract, compute and plug in posterior median
h_impute <- brms::as_draws_df(m.hab, variable = "^Ymi_ML1", regex = TRUE)
h_impute$.chain <- NULL; h_impute$.iteration <- NULL; h_impute$.draw <- NULL
h_miss_idx <- which(is.na(h_imp$ML.1))
h_imp[h_miss_idx,]$ML.1 <- apply(h_impute[h_miss_idx,], 2, median)

h_grid <- rbind(transform(h_imp, habitat3 = 0),
              transform(h_imp, habitat3 = 1),
              transform(h_imp, habitat3 = 2))

# predict full posterior; set covariates to 0...
mhab_imp.epred <- h_grid %>%
  mutate(ML.1 = 0, WoS = 0, depth.mean=0) %>%
  add_epred_draws(m.hab, resp = "CNS1")

hab_postpts <- 
  mhab_imp.epred |>
  ggplot(aes(x = habitat3,
             y = .epred)) +
  stat_pointinterval(point_interval=median_hdci) + #by default showing 66% and 95% intervals
  geom_point(aes(x=as.numeric(habitat3)-1, y = CNS.1), data = h_imp, position = position_jitter(width = 0.15), colour='blue') +
  labs(x="Habitat", y = "CNS") +
  scale_x_continuous(breaks=c(1,2,3), labels=(c("pelagic", "variable", "benthic"))) + theme(axis.text.x=element_text(size=11)) + theme_classic()
hab_postpts <- hab_postpts + theme_classic()
hab_postpts
save(hab_postpts, file="nov24_hab_postpts.rda")

##binary benthic----
bp_imp <- st.cephdat
# imputed mantle length: extract, compute and plug in posterior median
bp_impute <- brms::as_draws_df(m.benth, variable = "^Ymi_ML1", regex = TRUE)
bp_impute$.chain <- NULL; bp_impute$.iteration <- NULL; bp_impute$.draw <- NULL
bp_miss_idx <- which(is.na(bp_imp$ML.1))
bp_imp[bp_miss_idx,]$ML.1 <- apply(bp_impute[bp_miss_idx,], 2, median)

View(bp_imp)

# prepare prediction grid; could be more fine-grained, e.g. 0.5 SDs
bppredgrid <- rbind(transform(bp_imp, benthic = 0),
                    transform(bp_imp, benthic = 1))

# predict full posterior with covariates at 0...
bp_imp.epred <- bppredgrid %>%
  mutate(ML.1 = 0, WoS = 0) %>%
  add_epred_draws(m.benth, resp = "CNS1")

# ... and plot posterior with with point intervals
bp_postpts <- bp_imp.epred |>
  ggplot(aes(x = benthic,
             y = .epred)) +
  stat_pointinterval(point_interval=median_hdci) +
  geom_point(aes(x = as.numeric(benthic)-1, y = CNS.1), data = bp_imp, position=position_jitter(width = 0.15), colour='blue') +
  labs(x="Habitat", y = "CNS") +
  scale_x_continuous(breaks=c(0,1), labels=(c("pelagic","benthic"))) + theme(axis.text.x=element_text(size=11)) + theme_classic()
bp_postpts
save(bp_postpts, file="nov24_bp_postpts.rda")

## stack habitat plots----
(hab_postpts / bp_postpts) + plot_layout(axes="collect_y")

#Figure S2----
## a) depth*benthic/pelagic----
### minimum depth: benthic interaction----
mindb.imp <- impute(data = st.cephdat, var = c("ML.1"), model = m.mindb)

# Generate combinations of values for depth.min and benthic
combinations <- expand.grid(depth.min = -1:2, benthic = 0:1)

# Create prediction grid for all combinations of depth.min and benthic
minnd <- do.call(rbind, lapply(1:nrow(combinations), function(i) {
  transform(mindb.imp, 
            depth.min = combinations$depth.min[i],
            benthic = combinations$benthic[i])
}))

# Alternatively: predict full posterior...
mindb.imp.epred <- minnd |>
  # unblock mutate(...) to get CIs similar to 
  # brms conditional_effects() defaults 
  # i.e. a conditional estimate with covariates at 0:
  mutate(ML.1 = 0, WoS = 0) |> 
  add_epred_draws(m.mindb, resp = "CNS1", re_formula = NA)

range(mindb.imp.epred$.epred)
hist(mindb.imp.epred$depth.min)
hist(st.cephdat$depth.min)

# ... and plot posterior with varying intervals!
mindbplot <- mindb.imp.epred |> 
  ggplot(aes(x = depth.min,
             y = .epred)) + 
  stat_lineribbon(aes(color = as.factor(benthic), fill = as.factor(benthic)), 
                  .width = c(.95),
                  point_interval = "median_hdci", # set point estimate and interval type method
                  alpha=0.5) + 
  geom_point(aes(x = depth.min, 
                 y = CNS.1, 
                 color = as.factor(benthic)),
             alpha = 0.75,
             data = mindb.imp) + 
  labs(x="Minimum depth", y = "CNS") +
  scale_fill_manual(values=c("darkblue","grey"), guide="none") +
  scale_color_manual(name="Habitat",values=c("0"= "blue","1"="black")) +
  theme(axis.text.x=element_text(size=13)) + theme_classic()
mindbplot <- mindbplot + theme(legend.position="none")
mindbplot

### maximum depth: benthic interaction----
maxdb.imp <- impute(data = st.cephdat, var = c("ML.1"), model = m.maxdb)

# Generate combinations of values for depth.max and benthic
combinationsmax <- expand.grid(depth.max = -3:2, benthic = 0:1)

# Create prediction grid for all combinations of depth.max and benthic
ndmax <- do.call(rbind, lapply(1:nrow(combinationsmax), function(i) {
  transform(maxdb.imp, 
            depth.max = combinationsmax$depth.max[i],
            benthic = combinationsmax$benthic[i])
}))

# Alternatively: predict full posterior...
maxdb.imp.epred <- ndmax |>
  mutate(ML.1 = 0, WoS = 0) |> 
  add_epred_draws(m.maxdb, resp = "CNS1", re_formula = NA)

# ... and plot posterior with varying intervals!
maxdbplot <- maxdb.imp.epred |> 
  ggplot(aes(x = depth.max,
             y = .epred)) + 
  stat_lineribbon(aes(color = as.factor(benthic), fill = as.factor(benthic)), 
                  .width = c(.95),
                  point_interval = "median_hdci", # set point estimate and interval type method
                  alpha=0.5) + 
  geom_point(aes(x = depth.max, 
                 y = CNS.1, 
                 color = as.factor(benthic)),
             alpha = 0.75,
             data = maxdb.imp) + 
  labs(x="Maximum depth", y = "CNS") +
  scale_fill_manual(values=c("darkblue","grey"), guide="none") +
  scale_color_manual(name="Habitat",values=c("0"= "blue","1"="black"), labels=c("Pelagic","Benthic")) +
  theme(axis.text.x=element_text(size=13)) + theme_classic()
maxdbplot 

## b)----
### diet.breadth:benthic----
dhab.imp <- impute(data = st.cephdat, var = c("ML.1","diet.breadth"), model = m.dhab.constrained)

dhabgrid <- expand.grid(diet.breadth = -2:2, benthic = 0:1)

nddhab <- do.call(rbind, lapply(1:nrow(dhabgrid), function(i) {
  transform(dhab.imp, 
            diet.breadth = dhabgrid$diet.breadth[i],
            benthic = dhabgrid$benthic[i])
}))

m.dhab.imp.epred <- nddhab |>
  mutate(ML.1 = 0, pos.latmean=0, depth.mean=0, articles.read=0) |> 
  add_epred_draws(m.dhab.constrained, resp = "CNS1", re_formula = NA)

dhabplot <- m.dhab.imp.epred |> 
  ggplot(aes(x = diet.breadth,
             y = .epred)) + 
  stat_lineribbon(aes(color = as.factor(benthic), fill = as.factor(benthic)), 
                  .width = c(.95),
                  point_interval = "median_hdci", # set point estimate and interval type method
                  alpha=0.5) + 
  geom_point(aes(x = diet.breadth, 
                 y = CNS.1, 
                 color = as.factor(benthic)),
             alpha = 0.75,
             data = dhab.imp) + 
  labs(x="Diet breadth", y = "CNS") +
  scale_fill_manual(values=c("darkblue","grey"), guide="none") +
  scale_color_manual(name="Habitat",values=c("0"= "blue","1"="black"), guide="none") +
  theme(axis.text.x=element_text(size=13))
dhabplot <- dhabplot + theme_classic()
dhabplot

### predator.breadth:benthic----
phab.imp <- impute(data = st.cephdat, var = c("ML.1","predator.breadth"), model = m.phab.constrained)

phabgrid <- expand.grid(predator.breadth = -2:2, benthic = 0:1)

nd.phab <- do.call(rbind, lapply(1:nrow(phabgrid), function(i) {
  transform(phab.imp, 
            predator.breadth = phabgrid$predator.breadth[i],
            benthic = phabgrid$benthic[i])
}))

phab.imp.epred <- nd.phab |>
  mutate(ML.1 = 0, articles.read = 0, depth.mean=0, pos.latmean=0) |> 
  add_epred_draws(m.phab.constrained, resp = "CNS1", re_formula = NA)

phabplot <- phab.imp.epred |> 
  ggplot(aes(x = predator.breadth,
             y = .epred)) + 
  stat_lineribbon(aes(color = as.factor(benthic), fill = as.factor(benthic)), 
                  .width = c(.95),
                  point_interval = "median_hdci", # set point estimate and interval type method
                  alpha=0.5) + 
  geom_point(aes(x = predator.breadth, 
                 y = CNS.1, 
                 color = as.factor(benthic)),
             alpha = 0.75,
             data = phab.imp) + 
  labs(x="Predator breadth", y = "CNS") +
  scale_fill_manual(values=c("darkblue","grey"), guide="none") +
  scale_color_manual(name="Habitat",values=c("0"= "blue","1"="black"), guide="none") +
  theme(axis.text.x=element_text(size=13)) + theme_classic()
phabplot <- phabplot + theme_classic()
phabplot

## c) depth ranges, habitat, and EQ----
depthplot_dat <- cephdat[,c("phy.species","habitat3","depth.max","depth.min","CNS.1","ML.1")]
depthplot_dat <- complete(mice(depthplot_dat, maxit=10))
depthplot_dat$CNS.1 <- log(depthplot_dat$CNS.1)
depthplot_dat$ML.1 <- log(depthplot_dat$ML.1)
depthplot_dat$EQ <- depthplot_dat$CNS.1/depthplot_dat$ML.1
cephCNSord <- depthplot_dat[order(depthplot_dat$EQ),]
cephCNSord$depth.max <- -1*cephCNSord$depth.max
cephCNSord$depth.min <- -1*cephCNSord$depth.min
cephCNSord$depth.mean <- (cephCNSord$depth.max+cephCNSord$depth.min)/2
EQ.depthplot_new <- ggplot(cephCNSord, aes(x=factor(EQ, levels=unique(EQ)), ymin=depth.max, ymax=depth.min, color=factor(habitat3))) + 
  geom_errorbar(width=0.3, size=1, aes(color=factor(habitat3))) +
  scale_color_manual(name="Habitat type", values=c("0"="#2a0593", "1"="#d6556d", "2" = '#feba2c'), labels=c("Pelagic","Varying","Benthic")) +
  scale_x_discrete(breaks=factor(cephCNSord$EQ, levels=unique(cephCNSord$EQ)), labels=(gsub("_", " ", cephCNSord$phy.species)), guide=guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(breaks = seq(-7000, 0, by = 500)) +
  xlab("Species left to right by ascending relative brain size (CNS/ML)") +
  ylab("Depth ranges") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
EQ.depthplot_new

## combine minimum, maximum, predator, and diet plots----
maxdbplot <- maxdbplot + theme(legend.position="none")

pdhab_plots <- (phabplot / dhabplot) + plot_layout(guides="collect", axes="collect_y")

minmaxdb_plots <- (mindbplot / maxdbplot) + plot_layout(axes='collect_y', guides='collect') 

combineddepthnov24 <- ( pdhab_plots | minmaxdb_plots) / EQ.depthplot_new
save(combineddepthnov24, file="combineddepthnov24.rda")

#Figure S3 flag for TB unable to remake (points shifted to the right)-----
## a) sociality----
### binary sociality plot----
s2_imp <- st.cephdat
# imputed mantle length: extract, compute and plug in posterior median
s2_impute <- brms::as_draws_df(m.soc, variable = "^Ymi_ML1", regex = TRUE)
s2_impute$.chain <- NULL; s2_impute$.iteration <- NULL; s2_impute$.draw <- NULL
s2_miss_idx <- which(is.na(s2_imp$ML.1))
soc_imp <- apply(s2_impute[s2_miss_idx,], 2, median)
s2_imp[s2_miss_idx,]$ML.1 <- soc_imp

View(s2_imp)

# prepare prediction grid; could be more fine-grained, e.g. 0.5 SDs
nds2 <- rbind(transform(s2_imp, sociality.bin = 0),
              transform(s2_imp, sociality.bin = 1))

# predict full posterior with covariates at 0...
msoc_imp.epred <- nds2 %>%
  mutate(ML.1 = 0, WoS = 0, depth.mean=0) %>%
  add_epred_draws(m.soc, resp = "CNS1")

# ... and plot posterior with with point intervals
soc2_postpts <- msoc_imp.epred |>
  ggplot(aes(x = sociality.bin,
             y = .epred)) +
  stat_pointinterval(point_interval=median_hdci) +
  geom_point(aes(x = as.numeric(sociality.bin)-1, y = CNS.1), data = s2_imp, position=position_jitter(width = 0.15), colour='blue') +
  labs(x="Sociality type", y = "CNS") +
  scale_x_continuous(breaks=c(0,1), labels=(c("solitary","gregarious"))) + theme(axis.text.x=element_text(size=11)) + theme_classic()

soc2_postpts <- soc2_postpts + theme_classic()
soc2_postpts
save(soc2_postpts, file="nov24_soc2_postpts.rda")

### sociality 3 categories----
s3_imp <- st.cephdat

# imputed mantle length: extract, compute and plug in posterior median
s3_impute <- brms::as_draws_df(m.soc3, variable = "^Ymi_ML1", regex = TRUE)
s3_impute$.chain <- NULL; s3_impute$.iteration <- NULL; s3_impute$.draw <- NULL
s3_miss_idx <- which(is.na(s3_imp$ML.1))
soc3_imp <- apply(s3_impute[s3_miss_idx,], 2, median)
s3_imp[s3_miss_idx,]$ML.1 <- soc3_imp

nds3 <- rbind(transform(s3_imp, sociality.3 = 1),
              transform(s3_imp, sociality.3 = 2),
              transform(s3_imp, sociality.3 = 3))

# predict full posterior; set covariates to 0...
msoc3_imp.epred <- nds3 %>%
  mutate(ML.1 = 0, WoS = 0, depth.mean=0) %>%
  add_epred_draws(m.soc3, resp = "CNS1")

soc3_postpts <- 
  msoc3_imp.epred |>
  ggplot(aes(x = sociality.3,
             y = .epred)) +
  stat_pointinterval(point_interval=median_hdci) + #by default showing 66% and 95% intervals
  geom_point(aes(x=as.numeric(sociality.3), y = CNS.1), data = s3_imp, position = position_jitter(width = 0.15), colour='blue') +
  labs(x="Sociality type", y = "CNS") +
  scale_x_continuous(breaks=c(1,2,3), labels=(c("solitary", "tolerant", "gregarious"))) + theme(axis.text.x=element_text(size=11)) + theme_classic()
soc3_postpts <- soc3_postpts + theme_classic()
soc3_postpts
save(soc3_postpts, file="soc3_postpts.rda")

## b) maturity----
### maximum age----
var_imp1 <- brms::as_draws_df(m.smmax, variable = "^Ymi_matagemax", regex = TRUE)
var_imp1$.chain <- NULL; var_imp1$.iteration <- NULL; var_imp1$.draw <- NULL
View(var_imp1)
# where imputed variable is missing
var_miss_idx1 <- which(is.na(st.cephdat$matage.max))

# compute posterior median and HPDI for imputed species
# NB: the posterior distributions are multimodal (see below); 
# therefore, specify highest central posterior density, hdci()
var_imp_mu1 <- apply(var_imp1[var_miss_idx1,], 2, median)
var_imp_ci1 <- apply(var_imp1[var_miss_idx1,], 2, tidybayes::hdci)

View(var_imp_ci1)

# get CNS size for species with imputed predictor
Ymi1 <- st.cephdat$CNS.1[var_miss_idx1]
View(Ymi1)

# collect in df for plotting
Xmiss1 <- data.frame(Ymi1=Ymi1, mu = var_imp_mu1, cl = var_imp_ci1[1,], cu = var_imp_ci1[2,])

all_matage <- c(st.cephdat$matage.max, Xmiss1$mu)
x_range <- range(all_matage, na.rm = TRUE)
x_vals <- seq(from = x_range[1], to = x_range[2], length.out = 100)
ce.matur1 <- conditional_effects(m.smmax, effects = "matage.max", resp = "CNS1", int_conditions = list(matage.max=x_vals))

plot_matage.max <- ggplot() +
  geom_point(data = st.cephdat, aes(x = matage.max, y = CNS.1), color = "blue") +
  geom_line(data = ce.matur1$CNS1.CNS1_matage.max, aes(x = matage.max, y = estimate__)) +
  geom_ribbon(data = ce.matur1$CNS1.CNS1_matage.max, aes(x = matage.max, ymin =lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data=Xmiss1, aes(x=mu,y=Ymi1), color="black", alpha=0.4) +
  labs(x="Maximum age of sexual maturity", y = "CNS") +
  theme(axis.text.x=element_text(size=11)) + theme_classic()
save(plot_matage.max, file="plot_matage.max.rda")

### minimum age-----
tidybayes::get_variables(m.smmin)
var_imp5 <- brms::as_draws_df(m.smmin, variable = "^Ymi_matagemin", regex = TRUE)
var_imp5$.chain <- NULL; var_imp5$.iteration <- NULL; var_imp5$.draw <- NULL

# where imputed variable is missing
var_miss_idx5 <- which(is.na(st.cephdat$matage.min))

# compute posterior median and HPDI for imputed species
# NB: the posterior distributions are multimodal (see below); 
# therefore, specify highest central posterior density, hdci()
var_imp_mu5 <- apply(var_imp5[var_miss_idx5,], 2, median)
var_imp_ci5 <- apply(var_imp5[var_miss_idx5,], 2, tidybayes::hdci)

# get CNS size for species with imputed predictor
Ymi5 <- st.cephdat$CNS.1[var_miss_idx5]

# collect in df for plotting
Xmiss5 <- data.frame(Ymi=Ymi5, mu = var_imp_mu5, cl = var_imp_ci5[1,], cu = var_imp_ci5[2,])
all_minage <- c(st.cephdat$matage.min, Xmiss5$mu)
range(all_minage, na.rm = TRUE)
xmin_vals <- seq(from = -1.984, to = 2.544825, length.out = 100)
ce.matur2 <- conditional_effects(m.smmin, effects = "matage.min", resp = "CNS1", int_conditions = list(matage.min=xmin_vals))

plot_matage.min <- ggplot() +
  geom_point(data = st.cephdat, aes(x = matage.min, y = CNS.1), color = "blue") +
  geom_line(data = ce.matur2$CNS1.CNS1_matage.min, aes(x = matage.min, y = estimate__)) +
  geom_ribbon(data = ce.matur2$CNS1.CNS1_matage.min, aes(x = matage.min, ymin =lower__, ymax = upper__), alpha = 0.2) +
  geom_point(data=Xmiss5, aes(x=mu,y=Ymi),color="black", alpha=0.4)  +
  labs(x="Minimum age of sexual maturity", y = "CNS") +
  theme(axis.text.x=element_text(size=11)) + theme_classic()
plot_matage.min
save(plot_matage.min, file="/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/plot_matage.min.rda")

## combine sociality and maturity plots----
load("/Users/kiranbasava/nonhumans/di_cephproject/analyses/cephalopod_analyses/plot_matage.max.rda")
socmatplot <- ((soc2_postpts / soc3_postpts) + plot_layout(axes='collect_y', guides='collect') | (plot_matage.max / plot_matage.min)) 
#manually edited to remove redundant CNS labels 
