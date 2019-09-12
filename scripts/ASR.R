### Likelihood Ancestral State Inference And Model Selection
# 
#
#
# Author Jesus Martinez-Gomez, Cheslea Specht
#
library(phytools)
library(corHMM)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(xtable)


### Set Working Directory and Read in Data
setwd("/Volumes/GoogleDrive/My\ Drive/Projects/Chapter1")
tree <- read.tree("data/TraitDataCollection/output/trimed_big_constraint_fixed_dated.tre")
data <-  read.nexus.data("data/TraitDataCollection/output/cleanTraitData.nex") %>%
  unlist() %>%
  as.data.frame() %>%
  rownames_to_column()
states <- c("Cyme","Panicle","Raceme","Umbel")

# Settings #####################
# Eventually I will try to make this script work with everything
# num_models =  # Then number of models 
# num_sim = # number of simulations
# num_simASR = # number of ASR simulations
# states <- c("Cyme","Panicle","Raceme","Umbel")

# Speficy Model Here
#
#



# Model Fitting #####################

# Goal: Determine best fit model from 15 derivates

# I fit the standard models (ER,SYM,ARD) as well as the hypothesized model raceme-derived and cyme-derived.
# As the latter models only make explict assumption about the evolution of raceme -> umbel and cyme -> umbel. Therefore
# to accurate explore all possibilities I fit revirsible ( raceme <-> umbel, cyme <-> umbel ) and irreversible models 
# (raceme -> umbel, cyme -> umbel ) versions. Reversible and irreversible models were tested in ER, SYM and ARD backgrounds.
# A total of 15 models were tested. 
# Models
# 
 
### Fit Standard Models ####

#ER
fit.ER <- rayDISC(phy=tree, data = data,  model="ER",
                node.states="marginal", diagn=FALSE)
#SYM
fit.SYM <- rayDISC(phy=tree, data = data,  model="SYM",
              node.states="marginal", diagn=FALSE)
#ARD
fit.ARD <- rayDISC(phy=tree, data = data,  model="ARD",
               node.states="marginal", diagn=FALSE)

### Fit Raceme -> Umbel Models ###

# Fit Raceme_Reversible_ARD
rate.mat.Raceme_Rev_ARD<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="ARD") %>%
  rate.par.drop(drop.par=c(10,11)) 

fit.Raceme_Rev_ARD <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Raceme_Rev_ARD,
                         node.states="marginal", diagn=FALSE)

# Fit Raceme_Irreversible_ARD
rate.mat.Raceme_Irr_ARD<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="ARD") %>%
  rate.par.drop(drop.par=c(10,11,3,6,9))

fit.Raceme_Irr_ARD <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Raceme_Irr_ARD,
                            node.states="marginal", diagn=FALSE)

# Fit Raceme_Reversible_SYM
rate.mat.Raceme_Rev_SYM<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="SYM") 
   rate.mat.Raceme_Rev_SYM[1,4] <- NA 
   rate.mat.Raceme_Rev_SYM[2,4] <- NA 
  
fit.Raceme_Rev_SYM <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Raceme_Rev_SYM,
                                node.states="marginal", diagn=FALSE)
  
# Fit Raceme_Irreversible_SYM
rate.mat.Raceme_Irr_SYM<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="SYM") 
rate.mat.Raceme_Irr_SYM[1,4] <- NA 
rate.mat.Raceme_Irr_SYM[2,4] <- NA 
rate.mat.Raceme_Irr_SYM[4,1] <- NA 
rate.mat.Raceme_Irr_SYM[4,2] <- NA 
rate.mat.Raceme_Irr_SYM[4,3] <- NA 

fit.Raceme_Irr_SYM <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Raceme_Irr_SYM,
                              node.states="marginal", diagn=FALSE)                              

# Fit Raceme_Reversible_ER
rate.mat.Raceme_Rev_ER<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="ER") 
  rate.mat.Raceme_Rev_ER[1,4] <- NA 
  rate.mat.Raceme_Rev_ER[2,4] <- NA 

fit.Raceme_Rev_ER <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Raceme_Rev_ER,
                                node.states="marginal", diagn=FALSE)
  
# Fit Raceme_Irreversible_ER
rate.mat.Raceme_Irr_ER<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="ER") 
rate.mat.Raceme_Irr_ER[1,4] <- NA 
rate.mat.Raceme_Irr_ER[2,4] <- NA 
rate.mat.Raceme_Irr_ER[4,1] <- NA 
rate.mat.Raceme_Irr_ER[4,2] <- NA 
rate.mat.Raceme_Irr_ER[4,3] <- NA 

fit.Raceme_Irr_ER <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Raceme_Irr_ER,
                             node.states="marginal", diagn=FALSE)

### Fit cyme -> Umbel Models ###
# Fit Cyme_Reversible_ARD
rate.mat.Cyme_Rev_ARD<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="ARD") %>%
  rate.par.drop(drop.par=c(11,12)) 

fit.Cyme_Rev_ARD <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Cyme_Rev_ARD,
                              node.states="marginal", diagn=FALSE)

# Fit Cyme_Irreversible_ARD
rate.mat.Cyme_Irr_ARD<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="ARD") %>%
  rate.par.drop(drop.par=c(11,12,3,6,9))

fit.Cyme_Irr_ARD <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Cyme_Irr_ARD,
                              node.states="marginal", diagn=FALSE)

# Fit Cyme_Reversible_SYM
rate.mat.Cyme_Rev_SYM<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="SYM") 
rate.mat.Cyme_Rev_SYM[3,4] <- NA 
rate.mat.Cyme_Rev_SYM[2,4] <- NA 

fit.Cyme_Rev_SYM <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Cyme_Rev_SYM,
                              node.states="marginal", diagn=FALSE)

# Fit Cyme_Irreversible_SYM
rate.mat.Cyme_Irr_SYM<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="SYM") 
rate.mat.Cyme_Irr_SYM[3,4] <- NA 
rate.mat.Cyme_Irr_SYM[2,4] <- NA 
rate.mat.Cyme_Irr_SYM[4,1] <- NA 
rate.mat.Cyme_Irr_SYM[4,2] <- NA 
rate.mat.Cyme_Irr_SYM[4,3] <- NA 

fit.Cyme_Irr_SYM <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Cyme_Irr_SYM,
                              node.states="marginal", diagn=FALSE)   

# Fit Cyme_Reversible_ER
rate.mat.Cyme_Rev_ER<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="ER") 
rate.mat.Cyme_Rev_ER[3,4] <- NA 
rate.mat.Cyme_Rev_ER[2,4] <- NA 

fit.Cyme_Rev_ER <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Cyme_Rev_ER,
                             node.states="marginal", diagn=FALSE)

# Fit Cyme_Irreversible_ER
rate.mat.Cyme_Irr_ER<-rate.mat.maker(hrm=FALSE,ntraits=1, nstates=4, model="ER") 
rate.mat.Cyme_Irr_ER[3,4] <- NA 
rate.mat.Cyme_Irr_ER[2,4] <- NA 
rate.mat.Cyme_Irr_ER[4,1] <- NA 
rate.mat.Cyme_Irr_ER[4,2] <- NA 
rate.mat.Cyme_Irr_ER[4,3] <- NA 

fit.Cyme_Irr_ER <- rayDISC(phy=tree, data = data, rate.mat = rate.mat.Cyme_Irr_ER,
                             node.states="marginal", diagn=FALSE)

### Compare AIC scores ####

model_selecton <- data.frame("ER"  <- c(fit.ER$loglik,   fit.ER$AIC,  fit.ER$AICc),
                      "Raceme_Rev_ER" <-c(fit.Raceme_Rev_ER$loglik, fit.Raceme_Rev_ER$AIC, fit.Raceme_Rev_ER$AICc),
                      "Raceme_Irr_ER" <-c(fit.Raceme_Irr_ER$loglik, fit.Raceme_Irr_ER$AIC, fit.Raceme_Irr_ER$AICc),
                      "Cyme_Rev_ER" <- c(fit.Cyme_Rev_ER$loglik,fit.Cyme_Rev_ER$AIC, fit.Cyme_Rev_ER$AICc),
                      "Cyme_Irr_ER" <- c(fit.Cyme_Irr_ER$loglik,fit.Cyme_Irr_ER$AIC, fit.Cyme_Irr_ER$AICc),
                      "SYM" <- c(fit.SYM$loglik, fit.SYM$AIC, fit.SYM$AICc),
                      "Raceme_Rev_SYM"<-c(fit.Raceme_Rev_SYM$loglik, fit.Raceme_Rev_SYM$AIC, fit.Raceme_Rev_SYM$AICc),
                      "Raceme_Irr_SYM"<-c(fit.Raceme_Irr_SYM$loglik, fit.Raceme_Irr_SYM$AIC, fit.Raceme_Irr_SYM$AICc),
                      "Cyme_Rev_SYM"<- c(fit.Cyme_Rev_SYM$loglik,fit.Cyme_Rev_SYM$AIC,fit.Cyme_Rev_SYM$AICc),
                      "Cyme_Irr_SYM"<- c(fit.Cyme_Irr_SYM$loglik,fit.Cyme_Irr_SYM$AIC,fit.Cyme_Irr_SYM$AICc),
                      "ARD" <- c(fit.ARD$loglik, fit.ARD$AIC, fit.ARD$AICc),
                      "Raceme_Rev_ARD" <-c(fit.Raceme_Rev_ARD$loglik, fit.Raceme_Rev_ARD$AIC, fit.Raceme_Rev_ARD$AICc),
                      "Raceme_Irr_ARD" <-c(fit.Raceme_Irr_ARD$loglik, fit.Raceme_Irr_ARD$AIC, fit.Raceme_Irr_ARD$AICc),
                      "Cyme_Rev_ARD"<- c(fit.Cyme_Rev_ARD$loglik, fit.Cyme_Rev_ARD$AIC, fit.Cyme_Rev_ARD$AICc),
                      "Cyme_Irr_ARD"<- c(fit.Cyme_Irr_ARD$loglik, fit.Cyme_Irr_ARD$AIC, fit.Cyme_Irr_ARD$AICc),
                      row.names = c("-lnL","AIC","AICc")                      ) 
model_names <- c("fit.ER",
                 "fit.Raceme_Rev_ER",
                 "fit.Raceme_Irr_ER",
                 "fit.Cyme_Rev_ER",
                 "fit.Cyme_Irr_ER",
                 "fit.SYM",
                 "fit.Raceme_Rev_SYM",
                 "fit.Raceme_Irr_SYM",
                 "fit.Cyme_Rev_SYM",
                 "fit.Cyme_Irr_SYM",
                 "fit.ARD",
                 "fit.Raceme_Rev_ARD",
                 "fit.Raceme_Irr_ARD",
                 "fit.Cyme_Rev_ARD",
                 "fit.Cyme_Irr_ARD") 
colnames(model_selecton) <- model_names


### Generate Latex Table
latex_table <- model_selecton %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  arrange(AICc) 
  
latex_table$deltaAICc<- sapply(latex_table$AICc,  function(x) x-min(latex_table$AICc))

latex_table <- xtable(latex_table)

print.xtable(latex_table, type="latex", 
             file="output/corHMM/ModelSelection.tex")



### Make Plot - with Delta AICc < 10
cols = c("#FEFA78","#46A299","#DFBC75","#CF65B5")
pdf("output/corHMM/EvolutionPrezASR.pdf",width=8, height=11, useDingbats=FALSE)

plot(tree, type = "fan", cex = .225,label.offset = 2, main = "Raceme_Rev_ARD",
     show.tip.label = FALSE)
  nodelabels(pie=fit.Raceme_Rev_ARD$states, cex=0.3,piecol = cols )
  tiplabels( pie=fit.Raceme_Rev_ARD$tip.states, cex=0.2,piecol = cols )
  
plot(tree, type = "fan", cex = .225,label.offse = 2, main = "Raceme_Irr_ARD",
     show.tip.label = FALSE)
  nodelabels(pie=fit.Raceme_Irr_ARD$states, cex=0.3,piecol = cols )
  tiplabels( pie=fit.Raceme_Irr_ARD$tip.states, cex=0.2,piecol = cols )

plot(tree, type = "fan", cex = .225,label.offset = 2, main = "Raceme_Irr_SYM",
     show.tip.label = FALSE)
  nodelabels(pie=fit.Raceme_Irr_SYM$states, cex=0.3,piecol = cols ) 
  tiplabels( pie=fit.Raceme_Irr_SYM$tip.states, cex=0.2,piecol = cols )
  
plot(tree, type = "fan", cex = .225,label.offset = 2, main = "ARD",
     show.tip.label = FALSE)
  nodelabels(pie=fit.ARD$states, cex=0.3,piecol = cols ) 
  tiplabels( pie=fit.ARD$tip.states, cex=0.2,piecol = cols )
  
dev.off()
 
#### Plot Phylogeny with Tips mapped on Using ggtree this was for Evolution Meeting
# DFBC75 = BROWN
# CF65B5 = umbel
# FEFA78 = yell cyme
# 46A299 = blue panicle

library(ggtree)
data1 <- data
rownames(data1) <- data$rowname
data1$rowname  <- NULL
 
pdf("output/corHMM/EvolutionPrezPhylo.pdf",width=8, height=11, useDingbats=FALSE)
  
  p1 <- ggtree(tree, layout="rectangular") 
  p10 <- gheatmap(p1, data1,
                  offset=0.0, 
                  width=0.1,hjust=0) + 
    scale_fill_manual(values=c("#FEFA78","#DFBC75","#46A299","#CF65B5"))
plot(p10)
    
  
dev.off()


### Save Results

results <- list(fit.ER,
                fit.Raceme_Rev_ER, 
                fit.Raceme_Irr_ER,
                fit.Cyme_Rev_ER,
                fit.Cyme_Irr_ER,
                fit.SYM,
                fit.Raceme_Rev_SYM,
                fit.Raceme_Irr_SYM,
                fit.Cyme_Rev_SYM,
                fit.Cyme_Irr_SYM,
                fit.ARD,
                fit.Raceme_Rev_ARD,
                fit.Raceme_Irr_ARD,
                fit.Cyme_Rev_ARD,
                fit.Cyme_Irr_ARD,
                model_names)

save(results, file = "output/corHMM/corHMM_results.RData")


