###########################################################################################################
# Meta-analysis of plant responses to mycorrhizal fungi for Hoeksema et al. (2018, Communications Biology)#
###########################################################################################################

# This R code was compiled by Jason Hoeksema, with substantial contributions from others, especially Wolfgang Viechtbauer, James Meadow, Sounak Chakraborty, and Marc Lajeunesse
# If you have questions, please contact Jason at hoeksema@olemiss.edu, although
# if your questions are about the mechanics of using the metafor package, please first consult the metafor website, which
# has extensive documentation here: http://www.metafor-project.org/doku.php

rm(list=ls())
setwd('C:/Users/Hoeksema/Documents/All Me/Current_projects/NESCent/Big_Analysis')
# setwd('/home/hoeksema/Documents/Hoeksema/mycodb_analysis') # edit to set desired working directory, which should contain
# the following files: MycoDB_version4.csv, PlantTree_version2.tre, and FungalTree_version2.tre (available at https://datadryad.org//resource/doi:10.5061/dryad.723m1.4)

# read in and subset MycoDB to create data frame with which to replicate analyses from Hoeksema et al. (2018)
mycodb.v4 <- data.frame(read.csv("MycoDB_version4.csv",header=T)) 
dim(mycodb.v4) # a data frame with 50 columns and 4581 rows
mycodb.v4 <- subset(mycodb.v4, Hoeksema2018 == "YES")
dim(mycodb.v4) # should be a dataframe with 50 columns and 4008 rows

# create new SD columns in which to impute missing SD values using 2 alternative methods, including 10 replicates for HotDeck_NN method
mycodb.v4$ctrl_sd_Bracken <- mycodb.v4$ctrl_sd
mycodb.v4$ctrl_sd_HotDeck.1 <- mycodb.v4$ctrl_sd
mycodb.v4$ctrl_sd_HotDeck.2 <- mycodb.v4$ctrl_sd
mycodb.v4$ctrl_sd_HotDeck.3 <- mycodb.v4$ctrl_sd
mycodb.v4$ctrl_sd_HotDeck.4 <- mycodb.v4$ctrl_sd
mycodb.v4$ctrl_sd_HotDeck.5 <- mycodb.v4$ctrl_sd
mycodb.v4$ctrl_sd_HotDeck.6 <- mycodb.v4$ctrl_sd
mycodb.v4$ctrl_sd_HotDeck.7 <- mycodb.v4$ctrl_sd
mycodb.v4$ctrl_sd_HotDeck.8 <- mycodb.v4$ctrl_sd
mycodb.v4$ctrl_sd_HotDeck.9 <- mycodb.v4$ctrl_sd
mycodb.v4$ctrl_sd_HotDeck.10 <- mycodb.v4$ctrl_sd

mycodb.v4$trt_sd_Bracken <- mycodb.v4$trt_sd
mycodb.v4$trt_sd_HotDeck.1 <- mycodb.v4$trt_sd
mycodb.v4$trt_sd_HotDeck.2 <- mycodb.v4$trt_sd
mycodb.v4$trt_sd_HotDeck.3 <- mycodb.v4$trt_sd
mycodb.v4$trt_sd_HotDeck.4 <- mycodb.v4$trt_sd
mycodb.v4$trt_sd_HotDeck.5 <- mycodb.v4$trt_sd
mycodb.v4$trt_sd_HotDeck.6 <- mycodb.v4$trt_sd
mycodb.v4$trt_sd_HotDeck.7 <- mycodb.v4$trt_sd
mycodb.v4$trt_sd_HotDeck.8 <- mycodb.v4$trt_sd
mycodb.v4$trt_sd_HotDeck.9 <- mycodb.v4$trt_sd
mycodb.v4$trt_sd_HotDeck.10 <- mycodb.v4$trt_sd

# re-order some levels to make reference level more intuitive
mycodb.v4$NONMYCOCONTROL2 <- relevel(as.factor(mycodb.v4$NONMYCOCONTROL2), ref="mics_not_added")
mycodb.v4$INOC.COMPLEXITY <- relevel(as.factor(mycodb.v4$INOC.COMPLEXITY), ref="Single")

# create some simpler variables and variable names for subsequent use in models
mycodb.v4$Id <- mycodb.v4$NONCTLTRTSETID
mycodb.v4$PlantSpecies <- mycodb.v4$PlantSpecies2018
mycodb.v4$FungalGenus2018 <- as.factor(mycodb.v4$FungalGenus2018)
mycodb.v4$FungalGenus <- as.factor(mycodb.v4$FungalGenus2018)
mycodb.v4$PAPERTITLE <- as.factor(mycodb.v4$PAPERTITLE)
mycodb.v4$PaperID <- as.integer(mycodb.v4$PAPERTITLE)
mycodb.v4$PaperID <- as.factor(mycodb.v4$PaperID)

# create separate 'full' frames for AM and EM symbiosis
# these 'full' datasets contain some rows with multiple fungi and thus without a specific fungal genus name
AMfull.dat <- subset(mycodb.v4, MYCORRHIZAETYPE=='AM')
dim(AMfull.dat)
EMfull.dat <- subset(mycodb.v4, MYCORRHIZAETYPE=='EM')
dim(EMfull.dat) # EMfull.dat was not used in further analyses, as it was barely larger than EMsub.dat

# create 'sub' datasets that contain only rows with single fungi and thus all rows have fungal genus names
withfungalnames.dat <- subset(mycodb.v4, !is.na(FungalGenus))
dim(withfungalnames.dat) # a data frame with 74 columns and 3399 rows
withfungalnames.dat$PlantSpecies <- factor(withfungalnames.dat$PlantSpecies)
withfungalnames.dat$FungalGenus <- factor(withfungalnames.dat$FungalGenus)
AMsub.dat <- subset(withfungalnames.dat, MYCORRHIZAETYPE=='AM')
dim(AMsub.dat)
EMsub.dat <- subset(withfungalnames.dat, MYCORRHIZAETYPE=='EM')
dim(EMsub.dat)

# create 2 new variables for evolutionary origins of ectomycorrhizal (EM) lifestyle in plants,
# one with 2 origins (Pinaceae and Rosids), another with 7 origins (Pinaceae plus 6 clades within the Rosids)
EMsub.dat$EMPLANT_2_ORIGINS <- ifelse(EMsub.dat$PlantFamily=='pinaceae','pinaceae', 'rosids')
EMsub.dat$EMPLANT_2_ORIGINS <- factor(EMsub.dat$EMPLANT_2_ORIGINS)
levels(EMsub.dat$EMPLANT_2_ORIGINS)
table(EMsub.dat$EMPLANT_2_ORIGINS) # ultimately, this variable was not used in analyses presented in the paper
EMsub.dat$PlantFamily <- factor(EMsub.dat$PlantFamily)
levels(EMsub.dat$PlantFamily)
table(EMsub.dat$PlantFamily)
EMsub.dat$EMPLANT_7_ORIGINS <- EMsub.dat$PlantFamily
# collapse the 9 EM plant families in EMsub.dat into 7 clades:
levels(EMsub.dat$EMPLANT_7_ORIGINS) <- c("fagales","dipterocarpaceae","fabaceae","fagales","myrtaceae","fagales","phyllanthaceae","pinaceae","salicaceae")
EMsub.dat$EMPLANT_7_ORIGINS <- factor(EMsub.dat$EMPLANT_7_ORIGINS)
levels(EMsub.dat$EMPLANT_7_ORIGINS)
table(EMsub.dat$EMPLANT_7_ORIGINS) 
# EMPLANT_7_ORIGINS corresponds to the "Plant Origin" variable in Hoeksema et al. (2018), and was used in all final EM analyses

# factor and print levels of the 2 variables coding for evolutionary origins of EM lifestyle in fungi
EMsub.dat$EMF_ORIGIN_TED <- as.factor(EMsub.dat$EMF_ORIGIN_TED)
levels(EMsub.dat$EMF_ORIGIN_TED)
EMsub.dat$EMF_ORIGIN_ALT <- as.factor(EMsub.dat$EMF_ORIGIN_ALT)
levels(EMsub.dat$EMF_ORIGIN_ALT)
# Note: EMF_ORIGIN_TED and ALT are identical except ALT only uses one EM origin for Hydnotrya and Tuber
# EMF_ORIGIN_TED correspond to the "Fungal Origin" variable in Hoeksema et al. (2018), and was used for all final EM analyses

# create factors coding for evolutionary origin specificity in which the levels are unique combinations of plant & fungal evolutionary origins
EMsub.dat$origin_combo_TED2 <- paste(EMsub.dat$EMF_ORIGIN_TED, EMsub.dat$EMPLANT_2_ORIGINS, sep="_x_")
EMsub.dat$origin_combo_TED7 <- paste(EMsub.dat$EMF_ORIGIN_TED, EMsub.dat$EMPLANT_7_ORIGINS, sep="_x_")
EMsub.dat$origin_combo_ALT2 <- paste(EMsub.dat$EMF_ORIGIN_ALT, EMsub.dat$EMPLANT_2_ORIGINS, sep="_x_")
EMsub.dat$origin_combo_ALT7 <- paste(EMsub.dat$EMF_ORIGIN_ALT, EMsub.dat$EMPLANT_7_ORIGINS, sep="_x_")
# The factor "origin_combo_TED7" corresponds to the Plant x Fungal Origin variable in Hoeksema et al. (2018), and was used for all final EM analyses
EMsub.dat$origin_combo_TED2 <- factor(EMsub.dat$origin_combo_TED2)
EMsub.dat$origin_combo_TED7 <- factor(EMsub.dat$origin_combo_TED7)
EMsub.dat$origin_combo_ALT2 <- factor(EMsub.dat$origin_combo_ALT2)
EMsub.dat$origin_combo_ALT7 <- factor(EMsub.dat$origin_combo_ALT7)

############################################
# Calculate ESTVAR4 = the primary variance estimator used in the primary analyses of the Hoeksema et al. (2018) paper
# It is based on Eq. 1 from Hedges et al. 1999, but replacing (imputing) CV (SD/mean) in observations
# in which SD is unknown (due to missing SD values) with its median from those in which it is known. 
# We perform this imputation separately for AMsub, AMfull, and EMsub.

library(metafor) # needed here for replmiss() function

AMfull.median_ctrl_cv <- with(AMfull.dat, median(ctrl_sd / ctrl_mass, na.rm=TRUE))
AMfull.median_trt_cv <- with(AMfull.dat, median(trt_sd / trt_mass, na.rm=TRUE))
AMfull.dat$ESTVAR4 <- with(AMfull.dat, trt_sd^2 / (trt_mass^2 * trt_reps) + ctrl_sd^2 / (ctrl_mass^2 * ctrl_reps))
AMfull.dat$ESTVAR4 <- with(AMfull.dat, replmiss(ESTVAR4, median(trt_sd^2  / trt_mass^2,  na.rm=TRUE) / trt_reps +
                                                  median(ctrl_sd^2 / ctrl_mass^2, na.rm=TRUE) / ctrl_reps))
AMsub.median_ctrl_cv <- with(AMsub.dat, median(ctrl_sd / ctrl_mass, na.rm=TRUE))
AMsub.median_trt_cv <- with(AMsub.dat, median(trt_sd / trt_mass, na.rm=TRUE))
AMsub.dat$ESTVAR4 <- with(AMsub.dat, trt_sd^2 / (trt_mass^2 * trt_reps) + ctrl_sd^2 / (ctrl_mass^2 * ctrl_reps))
AMsub.dat$ESTVAR4 <- with(AMsub.dat, replmiss(ESTVAR4, median(trt_sd^2  / trt_mass^2,  na.rm=TRUE) / trt_reps +
                                                  median(ctrl_sd^2 / ctrl_mass^2, na.rm=TRUE) / ctrl_reps))
EMsub.median_ctrl_cv <- with(EMsub.dat, median(ctrl_sd / ctrl_mass, na.rm=TRUE))
EMsub.median_trt_cv <- with(EMsub.dat, median(trt_sd / trt_mass, na.rm=TRUE))
EMsub.dat$ESTVAR4 <- with(EMsub.dat, trt_sd^2 / (trt_mass^2 * trt_reps) + ctrl_sd^2 / (ctrl_mass^2 * ctrl_reps))
EMsub.dat$ESTVAR4 <- with(EMsub.dat, replmiss(ESTVAR4, median(trt_sd^2  / trt_mass^2,  na.rm=TRUE) / trt_reps +
                                                median(ctrl_sd^2 / ctrl_mass^2, na.rm=TRUE) / ctrl_reps))
###########################################
# Impute missing values of SD using two alternative methods: 
# Bracken1992 and HotDeck_NN from impute_SD() function in metagear package, 
# for comparing results with those obtained using ESTVAR4 as the variance estimator
# Do this separately for AMsub and EMsub subsets
# For HotDeck_NN, do this separately for trt and ctrl SD values, with 10 replicate data sets

library(metagear)
AMsub.dat <- impute_SD(AMsub.dat, "ctrl_sd_Bracken", "ctrl_mass", method="Bracken1992")
AMsub.dat <- impute_SD(AMsub.dat, "trt_sd_Bracken", "trt_mass", method="Bracken1992")
EMsub.dat <- impute_SD(EMsub.dat, "ctrl_sd_Bracken", "ctrl_mass", method="Bracken1992")
EMsub.dat <- impute_SD(EMsub.dat, "trt_sd_Bracken", "trt_mass", method="Bracken1992")

AMsub.dat.1a <- impute_SD(AMsub.dat, "ctrl_sd_HotDeck.1", "ctrl_mass", method="HotDeck_NN", range=3)
AMsub.dat.1a <- data.frame(AMsub.dat.1a)
AMsub.dat.1b <- impute_SD(AMsub.dat, "trt_sd_HotDeck.1", "trt_mass", method="HotDeck_NN", range=3)
AMsub.dat.1b <- data.frame(AMsub.dat.1b)

AMsub.dat.2a <- impute_SD(AMsub.dat, "ctrl_sd_HotDeck.2", "ctrl_mass", method="HotDeck_NN", range=3)
AMsub.dat.2a <- data.frame(AMsub.dat.2a)
AMsub.dat.2b <- impute_SD(AMsub.dat, "trt_sd_HotDeck.2", "trt_mass", method="HotDeck_NN", range=3)
AMsub.dat.2b <- data.frame(AMsub.dat.2b)

AMsub.dat.3a <- impute_SD(AMsub.dat, "ctrl_sd_HotDeck.3", "ctrl_mass", method="HotDeck_NN", range=3)
AMsub.dat.3a <- data.frame(AMsub.dat.3a)
AMsub.dat.3b <- impute_SD(AMsub.dat, "trt_sd_HotDeck.3", "trt_mass", method="HotDeck_NN", range=3)
AMsub.dat.3b <- data.frame(AMsub.dat.3b)

AMsub.dat.4a <- impute_SD(AMsub.dat, "ctrl_sd_HotDeck.4", "ctrl_mass", method="HotDeck_NN", range=3)
AMsub.dat.4a <- data.frame(AMsub.dat.4a)
AMsub.dat.4b <- impute_SD(AMsub.dat, "trt_sd_HotDeck.4", "trt_mass", method="HotDeck_NN", range=3)
AMsub.dat.4b <- data.frame(AMsub.dat.4b)

AMsub.dat.5a <- impute_SD(AMsub.dat, "ctrl_sd_HotDeck.5", "ctrl_mass", method="HotDeck_NN", range=3)
AMsub.dat.5a <- data.frame(AMsub.dat.5a)
AMsub.dat.5b <- impute_SD(AMsub.dat, "trt_sd_HotDeck.5", "trt_mass", method="HotDeck_NN", range=3)
AMsub.dat.5b <- data.frame(AMsub.dat.5b)

AMsub.dat.6a <- impute_SD(AMsub.dat, "ctrl_sd_HotDeck.6", "ctrl_mass", method="HotDeck_NN", range=3)
AMsub.dat.6a <- data.frame(AMsub.dat.6a)
AMsub.dat.6b <- impute_SD(AMsub.dat, "trt_sd_HotDeck.6", "trt_mass", method="HotDeck_NN", range=3)
AMsub.dat.6b <- data.frame(AMsub.dat.6b)

AMsub.dat.7a <- impute_SD(AMsub.dat, "ctrl_sd_HotDeck.7", "ctrl_mass", method="HotDeck_NN", range=3)
AMsub.dat.7a <- data.frame(AMsub.dat.7a)
AMsub.dat.7b <- impute_SD(AMsub.dat, "trt_sd_HotDeck.7", "trt_mass", method="HotDeck_NN", range=3)
AMsub.dat.7b <- data.frame(AMsub.dat.7b)

AMsub.dat.8a <- impute_SD(AMsub.dat, "ctrl_sd_HotDeck.8", "ctrl_mass", method="HotDeck_NN", range=3)
AMsub.dat.8a <- data.frame(AMsub.dat.8a)
AMsub.dat.8b <- impute_SD(AMsub.dat, "trt_sd_HotDeck.8", "trt_mass", method="HotDeck_NN", range=3)
AMsub.dat.8b <- data.frame(AMsub.dat.8b)

AMsub.dat.9a <- impute_SD(AMsub.dat, "ctrl_sd_HotDeck.9", "ctrl_mass", method="HotDeck_NN", range=3)
AMsub.dat.9a <- data.frame(AMsub.dat.9a)
AMsub.dat.9b <- impute_SD(AMsub.dat, "trt_sd_HotDeck.9", "trt_mass", method="HotDeck_NN", range=3)
AMsub.dat.9b <- data.frame(AMsub.dat.9b)

AMsub.dat.10a <- impute_SD(AMsub.dat, "ctrl_sd_HotDeck.10", "ctrl_mass", method="HotDeck_NN", range=3)
AMsub.dat.10a <- data.frame(AMsub.dat.10a)
AMsub.dat.10b <- impute_SD(AMsub.dat, "trt_sd_HotDeck.10", "trt_mass", method="HotDeck_NN", range=3)
AMsub.dat.10b <- data.frame(AMsub.dat.10b)

EMsub.dat.1a <- impute_SD(EMsub.dat, "ctrl_sd_HotDeck.1", "ctrl_mass", method="HotDeck_NN", range=3)
EMsub.dat.1a <- data.frame(EMsub.dat.1a)
EMsub.dat.1b <- impute_SD(EMsub.dat, "trt_sd_HotDeck.1", "trt_mass", method="HotDeck_NN", range=3)
EMsub.dat.1b <- data.frame(EMsub.dat.1b)

EMsub.dat.2a <- impute_SD(EMsub.dat, "ctrl_sd_HotDeck.2", "ctrl_mass", method="HotDeck_NN", range=3)
EMsub.dat.2a <- data.frame(EMsub.dat.2a)
EMsub.dat.2b <- impute_SD(EMsub.dat, "trt_sd_HotDeck.2", "trt_mass", method="HotDeck_NN", range=3)
EMsub.dat.2b <- data.frame(EMsub.dat.2b)

EMsub.dat.3a <- impute_SD(EMsub.dat, "ctrl_sd_HotDeck.3", "ctrl_mass", method="HotDeck_NN", range=3)
EMsub.dat.3a <- data.frame(EMsub.dat.3a)
EMsub.dat.3b <- impute_SD(EMsub.dat, "trt_sd_HotDeck.3", "trt_mass", method="HotDeck_NN", range=3)
EMsub.dat.3b <- data.frame(EMsub.dat.3b)

EMsub.dat.4a <- impute_SD(EMsub.dat, "ctrl_sd_HotDeck.4", "ctrl_mass", method="HotDeck_NN", range=3)
EMsub.dat.4a <- data.frame(EMsub.dat.4a)
EMsub.dat.4b <- impute_SD(EMsub.dat, "trt_sd_HotDeck.4", "trt_mass", method="HotDeck_NN", range=3)
EMsub.dat.4b <- data.frame(EMsub.dat.4b)

EMsub.dat.5a <- impute_SD(EMsub.dat, "ctrl_sd_HotDeck.5", "ctrl_mass", method="HotDeck_NN", range=3)
EMsub.dat.5a <- data.frame(EMsub.dat.5a)
EMsub.dat.5b <- impute_SD(EMsub.dat, "trt_sd_HotDeck.5", "trt_mass", method="HotDeck_NN", range=3)
EMsub.dat.5b <- data.frame(EMsub.dat.5b)

EMsub.dat.6a <- impute_SD(EMsub.dat, "ctrl_sd_HotDeck.6", "ctrl_mass", method="HotDeck_NN", range=3)
EMsub.dat.6a <- data.frame(EMsub.dat.6a)
EMsub.dat.6b <- impute_SD(EMsub.dat, "trt_sd_HotDeck.6", "trt_mass", method="HotDeck_NN", range=3)
EMsub.dat.6b <- data.frame(EMsub.dat.6b)

EMsub.dat.7a <- impute_SD(EMsub.dat, "ctrl_sd_HotDeck.7", "ctrl_mass", method="HotDeck_NN", range=3)
EMsub.dat.7a <- data.frame(EMsub.dat.7a)
EMsub.dat.7b <- impute_SD(EMsub.dat, "trt_sd_HotDeck.7", "trt_mass", method="HotDeck_NN", range=3)
EMsub.dat.7b <- data.frame(EMsub.dat.7b)

EMsub.dat.8a <- impute_SD(EMsub.dat, "ctrl_sd_HotDeck.8", "ctrl_mass", method="HotDeck_NN", range=3)
EMsub.dat.8a <- data.frame(EMsub.dat.8a)
EMsub.dat.8b <- impute_SD(EMsub.dat, "trt_sd_HotDeck.8", "trt_mass", method="HotDeck_NN", range=3)
EMsub.dat.8b <- data.frame(EMsub.dat.8b)

EMsub.dat.9a <- impute_SD(EMsub.dat, "ctrl_sd_HotDeck.9", "ctrl_mass", method="HotDeck_NN", range=3)
EMsub.dat.9a <- data.frame(EMsub.dat.9a)
EMsub.dat.9b <- impute_SD(EMsub.dat, "trt_sd_HotDeck.9", "trt_mass", method="HotDeck_NN", range=3)
EMsub.dat.9b <- data.frame(EMsub.dat.9b)

EMsub.dat.10a <- impute_SD(EMsub.dat, "ctrl_sd_HotDeck.10", "ctrl_mass", method="HotDeck_NN", range=3)
EMsub.dat.10a <- data.frame(EMsub.dat.10a)
EMsub.dat.10b <- impute_SD(EMsub.dat, "trt_sd_HotDeck.10", "trt_mass", method="HotDeck_NN", range=3)
EMsub.dat.10b <- data.frame(EMsub.dat.10b)

# Calculate estimated variance of log response ratio, using Eq. 1 from Hedges et al. 1999, 
# with imputed values of SD from both the Bracken1992 and HotDeck_NN methods (Bracken1992 = ESTVAR5, HotDeck_NN = ESTVAR6.1-ESTVAR6.10)
AMsub.dat$ESTVAR5<-((AMsub.dat$trt_sd_Bracken)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat$ctrl_sd_Bracken)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR5<-as.numeric(AMsub.dat$ESTVAR5)

AMsub.dat$ESTVAR6.1<-((AMsub.dat.1b$trt_sd_HotDeck.1)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat.1a$ctrl_sd_HotDeck.1)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR6.1<-as.numeric(AMsub.dat$ESTVAR6.1)
AMsub.dat$ESTVAR6.2<-((AMsub.dat.2b$trt_sd_HotDeck.2)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat.2a$ctrl_sd_HotDeck.2)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR6.2<-as.numeric(AMsub.dat$ESTVAR6.2)
AMsub.dat$ESTVAR6.3<-((AMsub.dat.3b$trt_sd_HotDeck.3)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat.3a$ctrl_sd_HotDeck.3)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR6.3<-as.numeric(AMsub.dat$ESTVAR6.3)
AMsub.dat$ESTVAR6.4<-((AMsub.dat.4b$trt_sd_HotDeck.4)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat.4a$ctrl_sd_HotDeck.4)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR6.4<-as.numeric(AMsub.dat$ESTVAR6.4)
AMsub.dat$ESTVAR6.5<-((AMsub.dat.5b$trt_sd_HotDeck.5)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat.5a$ctrl_sd_HotDeck.5)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR6.5<-as.numeric(AMsub.dat$ESTVAR6.5)
AMsub.dat$ESTVAR6.6<-((AMsub.dat.6b$trt_sd_HotDeck.6)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat.6a$ctrl_sd_HotDeck.6)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR6.6<-as.numeric(AMsub.dat$ESTVAR6.6)
AMsub.dat$ESTVAR6.7<-((AMsub.dat.7b$trt_sd_HotDeck.7)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat.7a$ctrl_sd_HotDeck.7)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR6.7<-as.numeric(AMsub.dat$ESTVAR6.7)
AMsub.dat$ESTVAR6.8<-((AMsub.dat.8b$trt_sd_HotDeck.8)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat.8a$ctrl_sd_HotDeck.8)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR6.8<-as.numeric(AMsub.dat$ESTVAR6.8)
AMsub.dat$ESTVAR6.9<-((AMsub.dat.9b$trt_sd_HotDeck.9)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat.9a$ctrl_sd_HotDeck.9)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR6.9<-as.numeric(AMsub.dat$ESTVAR6.9)
AMsub.dat$ESTVAR6.10<-((AMsub.dat.10b$trt_sd_HotDeck.10)^2/(AMsub.dat$trt_reps*(AMsub.dat$trt_mass)^2))+ ((AMsub.dat.10a$ctrl_sd_HotDeck.10)^2/(AMsub.dat$ctrl_reps*(AMsub.dat$ctrl_mass)^2))
AMsub.dat$ESTVAR6.10<-as.numeric(AMsub.dat$ESTVAR6.10)

EMsub.dat$ESTVAR5<-((EMsub.dat$trt_sd_Bracken)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat$ctrl_sd_Bracken)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR5<-as.numeric(EMsub.dat$ESTVAR5)

EMsub.dat$ESTVAR6.1<-((EMsub.dat.1b$trt_sd_HotDeck.1)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat.1a$ctrl_sd_HotDeck.1)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR6.1<-as.numeric(EMsub.dat$ESTVAR6.1)
EMsub.dat$ESTVAR6.2<-((EMsub.dat.2b$trt_sd_HotDeck.2)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat.2a$ctrl_sd_HotDeck.2)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR6.2<-as.numeric(EMsub.dat$ESTVAR6.2)
EMsub.dat$ESTVAR6.3<-((EMsub.dat.3b$trt_sd_HotDeck.3)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat.3a$ctrl_sd_HotDeck.3)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR6.3<-as.numeric(EMsub.dat$ESTVAR6.3)
EMsub.dat$ESTVAR6.4<-((EMsub.dat.4b$trt_sd_HotDeck.4)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat.4a$ctrl_sd_HotDeck.4)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR6.4<-as.numeric(EMsub.dat$ESTVAR6.4)
EMsub.dat$ESTVAR6.5<-((EMsub.dat.5b$trt_sd_HotDeck.5)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat.5a$ctrl_sd_HotDeck.5)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR6.5<-as.numeric(EMsub.dat$ESTVAR6.5)
EMsub.dat$ESTVAR6.6<-((EMsub.dat.6b$trt_sd_HotDeck.6)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat.6a$ctrl_sd_HotDeck.6)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR6.6<-as.numeric(EMsub.dat$ESTVAR6.6)
EMsub.dat$ESTVAR6.7<-((EMsub.dat.7b$trt_sd_HotDeck.7)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat.7a$ctrl_sd_HotDeck.7)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR6.7<-as.numeric(EMsub.dat$ESTVAR6.7)
EMsub.dat$ESTVAR6.8<-((EMsub.dat.8b$trt_sd_HotDeck.8)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat.8a$ctrl_sd_HotDeck.8)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR6.8<-as.numeric(EMsub.dat$ESTVAR6.8)
EMsub.dat$ESTVAR6.9<-((EMsub.dat.9b$trt_sd_HotDeck.9)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat.9a$ctrl_sd_HotDeck.9)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR6.9<-as.numeric(EMsub.dat$ESTVAR6.9)
EMsub.dat$ESTVAR6.10<-((EMsub.dat.10b$trt_sd_HotDeck.10)^2/(EMsub.dat$trt_reps*(EMsub.dat$trt_mass)^2))+ ((EMsub.dat.10a$ctrl_sd_HotDeck.10)^2/(EMsub.dat$ctrl_reps*(EMsub.dat$ctrl_mass)^2))
EMsub.dat$ESTVAR6.10<-as.numeric(EMsub.dat$ESTVAR6.10)

# filter AMsub and EMsub dataframes to contain only columns needed for subsequent analyses
AMsub.vars <- c("EXPERIMENTID","CTLTRTSETID","Id","LASTNAME1","LASTNAME2","PAPERYEAR","JOURNALNAME",
"PaperID","EFFECTSIZE1", "PlantFamily","PlantSpecies","FungalGenus","PLANTLIFEHISTORY","DOMESTICATED",
"FUNGROUP","NONMYCOCONTROL2","FERTP","FERTN","STERILIZED","LOCATION","ESTVAR4","ESTVAR5","ESTVAR6.1","ESTVAR6.2",
"ESTVAR6.3","ESTVAR6.4","ESTVAR6.5","ESTVAR6.6","ESTVAR6.7","ESTVAR6.8","ESTVAR6.9","ESTVAR6.10")
AMsub.dat <- AMsub.dat[AMsub.vars]

EMsub.vars <- c("EXPERIMENTID","CTLTRTSETID","Id","LASTNAME1","LASTNAME2","PAPERYEAR","JOURNALNAME",
"PaperID","EFFECTSIZE1","PlantFamily","PlantSpecies","FungalGenus","EMF_ORIGIN_TED","EMPLANT_7_ORIGINS","origin_combo_TED7",
"FUNGROUP","NONMYCOCONTROL2","FERTP","FERTN","STERILIZED","LOCATION","ESTVAR4","ESTVAR5","ESTVAR6.1","ESTVAR6.2","ESTVAR6.3","ESTVAR6.4","ESTVAR6.5",
"ESTVAR6.6","ESTVAR6.7","ESTVAR6.8","ESTVAR6.9","ESTVAR6.10")
EMsub.dat <- EMsub.dat[EMsub.vars]

# count how many unique publications are in each data set
levels(EMfull.dat$PaperID)
levels(AMfull.dat$PaperID)
levels(EMsub.dat$PaperID)
levels(AMsub.dat$PaperID)

#### Factor some categorical variables so that they do not contain un-used levels, as
# the presence of un-used levels will give an error in rma.mv() about "Model matrix not of full rank."
AMfull.dat$PlantFamily<-factor(AMfull.dat$PlantFamily)
AMfull.dat$PlantSpecies<-factor(AMfull.dat$PlantSpecies)
AMfull.dat$FungalGenus<-factor(AMfull.dat$FungalGenus)
AMfull.dat$FUNGROUP<-factor(AMfull.dat$FUNGROUP)
AMfull.dat$NONMYCOCONTROL2<-factor(AMfull.dat$NONMYCOCONTROL2)
AMfull.dat$FERTP<-factor(AMfull.dat$FERTP)
AMfull.dat$FERTN<-factor(AMfull.dat$FERTN)
AMfull.dat$INOC.COMPLEXITY<-factor(AMfull.dat$INOC.COMPLEXITY)
AMfull.dat$STERILIZED<-factor(AMfull.dat$STERILIZED)
AMfull.dat$LOCATION<-factor(AMfull.dat$LOCATION)
AMfull.dat$DOMESTICATED<-factor(AMfull.dat$DOMESTICATED)
AMfull.dat$PLANTLIFEHISTORY<-factor(AMfull.dat$PLANTLIFEHISTORY)
AMfull.dat$PaperID<-factor(AMfull.dat$PaperID)

AMsub.dat$PlantFamily<-factor(AMsub.dat$PlantFamily)
AMsub.dat$PlantSpecies<-factor(AMsub.dat$PlantSpecies)
AMsub.dat$FungalGenus<-factor(AMsub.dat$FungalGenus)
AMsub.dat$FUNGROUP<-factor(AMsub.dat$FUNGROUP)
AMsub.dat$NONMYCOCONTROL2<-factor(AMsub.dat$NONMYCOCONTROL2)
AMsub.dat$FERTP<-factor(AMsub.dat$FERTP)
AMsub.dat$FERTN<-factor(AMsub.dat$FERTN)
AMsub.dat$STERILIZED<-factor(AMsub.dat$STERILIZED)
AMsub.dat$LOCATION<-factor(AMsub.dat$LOCATION)
AMsub.dat$DOMESTICATED<-factor(AMsub.dat$DOMESTICATED)
AMsub.dat$PLANTLIFEHISTORY<-factor(AMsub.dat$PLANTLIFEHISTORY)
AMsub.dat$PaperID<-factor(AMsub.dat$PaperID)

EMsub.dat$PlantFamily<-factor(EMsub.dat$PlantFamily)
EMsub.dat$PlantSpecies<-factor(EMsub.dat$PlantSpecies)
EMsub.dat$FungalGenus<-factor(EMsub.dat$FungalGenus)
EMsub.dat$FUNGROUP<-factor(EMsub.dat$FUNGROUP)
EMsub.dat$NONMYCOCONTROL2<-factor(EMsub.dat$NONMYCOCONTROL2)
EMsub.dat$FERTP<-factor(EMsub.dat$FERTP)
EMsub.dat$FERTN<-factor(EMsub.dat$FERTN)
EMsub.dat$STERILIZED<-factor(EMsub.dat$STERILIZED)
EMsub.dat$LOCATION<-factor(EMsub.dat$LOCATION)
EMsub.dat$PaperID<-factor(EMsub.dat$PaperID)

# make tables showing levels of predictors (and their counts) for data sets
table(AMfull.dat$PlantFamily)
table(AMfull.dat$FungalGenus)
table(AMfull.dat$FUNGROUP)
table(AMfull.dat$NONMYCOCONTROL2)
table(AMfull.dat$FERTP)
table(AMfull.dat$FERTN)
table(AMfull.dat$INOC.COMPLEXITY)
table(AMfull.dat$STERILIZED)
table(AMfull.dat$LOCATION)
table(AMfull.dat$DOMESTICATED)
table(AMfull.dat$PLANTLIFEHISTORY)

table(AMsub.dat$PlantFamily)
table(AMsub.dat$FungalGenus)
table(AMsub.dat$FUNGROUP)
table(AMsub.dat$NONMYCOCONTROL2)
table(AMsub.dat$FERTP)
table(AMsub.dat$FERTN)
table(AMsub.dat$STERILIZED)
table(AMsub.dat$LOCATION)
table(AMsub.dat$DOMESTICATED)
table(AMsub.dat$PLANTLIFEHISTORY)

table(EMsub.dat$PlantFamily)
table(EMsub.dat$FungalGenus)
table(EMsub.dat$FUNGROUP)
table(EMsub.dat$NONMYCOCONTROL2)
table(EMsub.dat$FERTP)
table(EMsub.dat$FERTN)
table(EMsub.dat$STERILIZED)
table(EMsub.dat$LOCATION)

### create additional factors needed below, for plant x fungus interaction in datasets with fungal names
AMsub.dat$plantxfungus <- paste(AMsub.dat$PlantSpecies, AMsub.dat$FungalGenus, sep="XX")
EMsub.dat$plantxfungus <- paste(EMsub.dat$PlantSpecies, EMsub.dat$FungalGenus, sep="XX")
AMsub.dat$plantxfungus <- factor(AMsub.dat$plantxfungus)
EMsub.dat$plantxfungus <- factor(EMsub.dat$plantxfungus)

### count the numbers of plant species, fungal genera, and plant x fungus combos in AMsub, AMfull, and EMsub datasets
nlevels(AMfull.dat$PlantSpecies)
nlevels(AMfull.dat$FungalGenus)
nlevels(AMsub.dat$PlantSpecies)
nlevels(AMsub.dat$FungalGenus)
nlevels(AMsub.dat$plantxfungus)
nlevels(EMsub.dat$PlantSpecies)
nlevels(EMsub.dat$FungalGenus)
nlevels(EMsub.dat$plantxfungus)
nlevels(EMsub.dat$origin_combo_TED7)

### create identical copies of plant and fungal taxonomic factors, for associating with phylogenies
AMfull.dat$PlantSpecies.phyl <- AMfull.dat$PlantSpecies

AMsub.dat$PlantSpecies.phyl <- AMsub.dat$PlantSpecies
AMsub.dat$FungalGenus.phyl  <- AMsub.dat$FungalGenus
AMsub.dat$plant.phylxfungus.phyl <- AMsub.dat$plantxfungus
AMsub.dat$plant_phylxfungus <- AMsub.dat$plantxfungus
AMsub.dat$plantxfungus_phyl <- AMsub.dat$plantxfungus

EMsub.dat$PlantSpecies.phyl <- EMsub.dat$PlantSpecies
EMsub.dat$FungalGenus.phyl  <- EMsub.dat$FungalGenus
EMsub.dat$plant.phylxfungus.phyl <- EMsub.dat$plantxfungus
EMsub.dat$plant_phylxfungus <- EMsub.dat$plantxfungus
EMsub.dat$plantxfungus_phyl<- EMsub.dat$plantxfungus

##########################################################
# use AMfull.dat for full analysis of largest dataset for AM fungi with data for all fixed factors (include INOC.COMPLEXITY, but leave out fungal taxonomy/phylogeny)
# use AMsub.dat for full analysis of largest dataset for AM fungi that has complete fungal names, and data for all fixed factors (except INOC.COMPLEXITY)
# use EMsub.dat for full analysis of largest dataset for EM fungi that has complete fungal names, and data for all fixed factors (except INOC.COMPLEXITY)
# EMfull was determined to be not worthy of analysis, since it is barely larger than EMsub and does not allow inclusion of fungal taxonomy/phylogeny
##########################################################

##########################################################
# Both AMsub and EMsub have observations (1 in AMsub, 3 in EMsub) that result in major outliers of ESTVAR6 when using HotDeck_NN imputation
# (4-5 orders of magnitude larger than the median value of ESTVAR6)
# Calculations to show ratio of outlier to median:
max(AMsub.dat$ESTVAR6.5)/median(AMsub.dat$ESTVAR6.5) # smallest ratio of outlier to median value in AMsub
max(AMsub.dat$ESTVAR6.8)/median(AMsub.dat$ESTVAR6.8) # largest ratio of outlier to median value in AMsub
max(EMsub.dat$ESTVAR6.7)/median(EMsub.dat$ESTVAR6.7) # smallest ratio of outlier to median value in EMsub
max(EMsub.dat$ESTVAR6.9)/median(EMsub.dat$ESTVAR6.9) # largest ratio of outlier to median value in EMsub

# These outliers prevent model fitting using likelihood, and result in biased error estimates
# It is always the same observations resulting in these outliers, across multiple replicate imputed datasets,
# so drop these observations from data before model fitting (1 from AMsub, 3 from EMsub) when testing alternative imputation methods:
AMsub2.dat <- AMsub.dat[-c(which.max(AMsub.dat$ESTVAR6.1)), ]
EMsub2.dat <- EMsub.dat[-c(which.max(EMsub.dat$ESTVAR6.1)), ]
EMsub2.dat <- EMsub2.dat[-c(which.max(EMsub2.dat$ESTVAR6.1)), ]
EMsub2.dat <- EMsub2.dat[-c(which.max(EMsub2.dat$ESTVAR6.1)), ]

##########################################################
# Read in and manipulate plant and fungal phylogenies
library(ape)

# read in plant and fungal phylogenetic trees in Newick format
tree.F <- read.tree(file="FungalTree_version2.tre") # available on Dryad at http://
R.F <- vcv(phy=tree.F, corr=TRUE) # corr=TRUE gives a correlation matrix
round(R.F[1:10,1:10],3) #view part of the R.F file

tree.P <- read.tree(file="PlantTree_version2.tre") # available on Dryad at http://
R.P <- vcv(phy=tree.P, corr=TRUE)
round(R.P[1:10,1:10],3) # view part of the R.P file

# make non-phylogenetic versions of plant and fungal correlation matrices for use below
R.F.nophyl <- R.F
R.F.nophyl[upper.tri(R.F.nophyl) | lower.tri(R.F.nophyl)] = 0
R.P.nophyl <- R.P
R.P.nophyl[upper.tri(R.P.nophyl) | lower.tri(R.P.nophyl)] = 0

### check that all of the plants/fungi in the AMsub data are actually in the trees
is.element(unique(AMsub.dat$PlantSpecies), tree.P$tip)
unique(AMsub.dat$PlantSpecies)[!is.element(unique(AMsub.dat$PlantSpecies), tree.P$tip)]
is.element(unique(AMsub.dat$FungalGenus),  tree.F$tip)
unique(AMsub.dat$FungalGenus)[!is.element(unique(AMsub.dat$FungalGenus),  tree.F$tip)]

### check that all of the plants/fungi in the EMsub data are actually in the trees
is.element(unique(EMsub.dat$PlantSpecies), tree.P$tip)
unique(EMsub.dat$PlantSpecies)[!is.element(unique(EMsub.dat$PlantSpecies), tree.P$tip)]
is.element(unique(EMsub.dat$FungalGenus),  tree.F$tip)
unique(EMsub.dat$FungalGenus)[!is.element(unique(EMsub.dat$FungalGenus),  tree.F$tip)]

### check that all of the plants in the AMfull data are actually in the plant tree
is.element(unique(AMfull.dat$PlantSpecies), tree.P$tip)
unique(AMfull.dat$PlantSpecies)[!is.element(unique(AMfull.dat$PlantSpecies), tree.P$tip)]

# Calculate tensor product of phylogenetic vcv matrices, for observed species in AMsub dataset
# for use with plant x fungus interaction effects
observedAM <- data.frame(levels(AMsub.dat$plantxfungus))
colnames(observedAM) <- "PFcombo"
library(reshape2)
observedAM_split<-colsplit(observedAM$PFcombo, "XX", c("plant","fungus"))
RO.AM <- matrix(NA, nrow=length(observedAM$PFcombo), ncol=length(observedAM$PFcombo))
rownames(RO.AM) <- levels(AMsub.dat$plantxfungus)
colnames(RO.AM) <- levels(AMsub.dat$plantxfungus)

for (i in 1:length(observedAM_split$plant)) {
   for (j in 1:length(observedAM_split$fungus)) {
      if (i <= j)
         next
            RO.AM[i,j] <- R.P[observedAM_split$plant[i], observedAM_split$plant[j]] * R.F[observedAM_split$fungus[i], observedAM_split$fungus[j]]
   }
}
diag(RO.AM) <- 1
RO.AM[upper.tri(RO.AM)] <- t(RO.AM)[upper.tri(RO.AM)]

# calculate tensor product of plant phylogeny x fungal taxonomy for AM fungi
R_PP_FT.AM <- matrix(NA, nrow=length(observedAM$PFcombo), ncol=length(observedAM$PFcombo))
rownames(R_PP_FT.AM) <- levels(AMsub.dat$plantxfungus)
colnames(R_PP_FT.AM) <- levels(AMsub.dat$plantxfungus)
for (i in 1:length(observedAM_split$plant)) {
   for (j in 1:length(observedAM_split$fungus)) {
      if (i <= j)
         next
         R_PP_FT.AM[i,j] <- R.P[observedAM_split$plant[i], observedAM_split$plant[j]] * R.F.nophyl[observedAM_split$fungus[i], observedAM_split$fungus[j]]
   }
}
diag(R_PP_FT.AM) <- 1
R_PP_FT.AM[upper.tri(R_PP_FT.AM)] <- t(R_PP_FT.AM)[upper.tri(R_PP_FT.AM)]

# calculate tensor product of plant taxonomy x fungal phylogeny for AM fungi
R_PT_FP.AM <- matrix(NA, nrow=length(observedAM$PFcombo), ncol=length(observedAM$PFcombo))
rownames(R_PT_FP.AM) <- levels(AMsub.dat$plantxfungus)
colnames(R_PT_FP.AM) <- levels(AMsub.dat$plantxfungus)
for (i in 1:length(observedAM_split$plant)) {
   for (j in 1:length(observedAM_split$fungus)) {
      if (i <= j)
         next
       R_PT_FP.AM[i,j] <- R.P.nophyl[observedAM_split$plant[i], observedAM_split$plant[j]] * R.F[observedAM_split$fungus[i], observedAM_split$fungus[j]]
   }
}
diag(R_PT_FP.AM) <- 1
R_PT_FP.AM[upper.tri(R_PT_FP.AM)] <- t(R_PT_FP.AM)[upper.tri(R_PT_FP.AM)]

## do the same for EM fungi
observedEM <- data.frame(levels(EMsub.dat$plantxfungus))
colnames(observedEM) <- "PFcombo"
library(reshape2)
observedEM_split<-colsplit(observedEM$PFcombo, "XX", c("plant","fungus"))
RO.EM <- matrix(NA, nrow=length(observedEM$PFcombo), ncol=length(observedEM$PFcombo))
rownames(RO.EM) <- levels(EMsub.dat$plantxfungus)
colnames(RO.EM) <- levels(EMsub.dat$plantxfungus)
for (i in 1:length(observedEM_split$plant)) {
   for (j in 1:length(observedEM_split$fungus)) {
      if (i <= j)
         next
       RO.EM[i,j] <- R.P[observedEM_split$plant[i], observedEM_split$plant[j]] * R.F[observedEM_split$fungus[i], observedEM_split$fungus[j]]
   }
}
diag(RO.EM) <- 1
RO.EM[upper.tri(RO.EM)] <- t(RO.EM)[upper.tri(RO.EM)]

# calculate tensor product of plant phylogeny x fungal taxonomy for EM fungi
R_PP_FT.EM <- matrix(NA, nrow=length(observedEM$PFcombo), ncol=length(observedEM$PFcombo))
rownames(R_PP_FT.EM) <- levels(EMsub.dat$plantxfungus)
colnames(R_PP_FT.EM) <- levels(EMsub.dat$plantxfungus)
for (i in 1:length(observedEM_split$plant)) {
   for (j in 1:length(observedEM_split$fungus)) {
      if (i <= j)
         next
            R_PP_FT.EM[i,j] <- R.P[observedEM_split$plant[i], observedEM_split$plant[j]] * R.F.nophyl[observedEM_split$fungus[i], observedEM_split$fungus[j]]
   }
}
diag(R_PP_FT.EM) <- 1
R_PP_FT.EM[upper.tri(R_PP_FT.EM)] <- t(R_PP_FT.EM)[upper.tri(R_PP_FT.EM)]

# calculate tensor product of plant taxonomy x fungal phylogeny for EM fungi
R_PT_FP.EM <- matrix(NA, nrow=length(observedEM$PFcombo), ncol=length(observedEM$PFcombo))
rownames(R_PT_FP.EM) <- levels(EMsub.dat$plantxfungus)
colnames(R_PT_FP.EM) <- levels(EMsub.dat$plantxfungus)
for (i in 1:length(observedEM_split$plant)) {
   for (j in 1:length(observedEM_split$fungus)) {
      if (i <= j)
         next
          R_PT_FP.EM[i,j] <- R.P.nophyl[observedEM_split$plant[i], observedEM_split$plant[j]] * R.F[observedEM_split$fungus[i], observedEM_split$fungus[j]]
   }
}
diag(R_PT_FP.EM) <- 1
R_PT_FP.EM[upper.tri(R_PT_FP.EM)] <- t(R_PT_FP.EM)[upper.tri(R_PT_FP.EM)]

# optionally, save workspace up to this point:
save.image(file="workspace_phylometa_data.RData")
# load a saved workspace image:
# load(file.choose())

#
####
#######
###############
######################### Model fitting
###############
#######
####
#

# load packages needed for model fitting
library(metafor)
library(multcomp)
library(glmulti)

# optional code for setting likelihood model-fitting optimizer for metafor functions
# optimizer <- "nlminb"
# control=list(optimizer="nlminb", rel.tol=1e-8, eval.max=10000, iter.max=10000)
# nlminb is often fastest, but occasionally fails to converge
optimizer <- "optim" # this is the default
optmethod <- "BFGS"  # (only relevant when optimizer="optim")
# BFGS almost always converges, but is usually a bit slower than nlminb

# a custom function we will use below to extract relative variable importance (RVI) values from glmulti objects:
extractRVI <- function(x) {
  ww = exp(-(x@crits - x@crits[1])/2)
  ww = ww/sum(ww)
  clartou = function(x) {
    pieces <- sort(strsplit(x, ":")[[1]])
    if (length(pieces) > 1)
      paste(pieces[1], ":", pieces[2], sep = "")
    else x
  }
  tet = lapply(x@formulas, function(x) sapply(attr(delete.response(terms(x)), "term.labels"), clartou))
  allt <- unique(unlist(tet))
  imp <- sapply(allt, function(x) sum(ww[sapply(tet, function(t) x %in% t)]))
  return(sort(imp))
}

##################
####################################### ANALYSIS OF EMsub data
##################
################## MODEL SELECTION, using both REMLr and ML to compare results
###################

# Warning: The model selection procedures below may take from hours to days or weeks, depending on processor speed and the number of processors being used
# It is recommended to fit individual models first and estimate per-model run time, then calculate number of total models using method=d in glmulti, and estimate total run time.

sink('Output_01 -- GLMULTI model selection, EMsub exhaustive, ML.txt') # combined with sink() below, this will put output from enclosed code in a .txt file in the working directory
optmethod <- "BFGS"
rma.glmulti.EM <- function(formula, data, ...) {
  rma.mv(as.formula(paste(deparse(formula))), V=ESTVAR4, method="ML",
         random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | FungalGenus.phyl, ~ 1 | plant.phylxfungus.phyl,
                       ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus, ~ 1 | plantxfungus,
                       ~ 1 | EMF_ORIGIN_TED, ~ 1 | EMPLANT_7_ORIGINS, ~ 1 | origin_combo_TED7), 
         R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl=R_PT_FP.EM), 
         data=EMsub.dat, control=list(optimizer=optimizer, optmethod=optmethod, 
                                      sigma2.init=c(.04,.16,.6,.08,.0001,.0001,.0001,.13,.0001,.0001,.0001,.0001,.18,0.0001)),outlist="minimal") # initial values to speed model fitting
}
# outlist = "minimal" reduces the size of rma.mv model fit objects so that the final glmulti object is of manageable size
time.start <- proc.time()
options(warn=-1)
modelSelEM1 <- glmulti(EFFECTSIZE1 ~ FERTN + FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, method="h", 
                       data=EMsub.dat, level=1, fitfunction=rma.glmulti.EM, crit="aicc", confsetsize=1000, plotty=FALSE)
options(warn=1)
time.end <- proc.time()
cat("Minutes:", ((time.end - time.start)/60)[3], "\n")

summary(modelSelEM1)
print(modelSelEM1)
extractRVI(modelSelEM1)
weightable(modelSelEM1)

sink()

sink('Output_02 -- GLMULTI model selection, EMsub exhaustive, REMLr.txt')
optmethod <- "BFGS"
rma.glmulti.EM <- function(formula, data, ...) {
  rma.mv(as.formula(paste(deparse(formula))), V=ESTVAR4, method="REML",
         random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | FungalGenus.phyl, ~ 1 | plant.phylxfungus.phyl,
                       ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus, ~ 1 | plantxfungus,
                       ~ 1 | EMF_ORIGIN_TED, ~ 1 | EMPLANT_7_ORIGINS, ~ 1 | origin_combo_TED7), 
         R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl=R_PT_FP.EM), 
         data=EMsub.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, 
                                      sigma2.init=c(.04,.16,.6,.08,.0001,.0001,.0001,.13,.0001,.0001,.0001,.0001,.18,0.0001)),outlist="minimal")
}
time.start <- proc.time()
options(warn=-1)
modelSelEM2 <- glmulti(EFFECTSIZE1 ~ FERTN + FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, method="h", 
                       data=EMsub.dat, level=1, fitfunction=rma.glmulti.EM, crit="aicc", confsetsize=1000, plotty=FALSE)
options(warn=1)
time.end <- proc.time()
cat("Minutes:", ((time.end - time.start)/60)[3], "\n")

summary(modelSelEM2)
print(modelSelEM2)
extractRVI(modelSelEM2)
weightable(modelSelEM2)

sink()

####################
#################### Fitting (using REML) of best EM model from model selection, for parameter estimation
####################

sink('Output_03 --EMsub BEST model fit with REML.txt')
# Note: Here, we use the model with 4 fixed effect predictors: STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP
# Those four fixed effects were the ones identified as important (with RVI >0.5) by both ML and REMLr model selection
# The model with all four of those was the AIC-best model according to ML model selection, and according to REML model selection
# was within 0.35 AIC units of the best model (which had only FERTN and FERTP).
cat("==EMbest REML ==\n")
optimizer <- "optim" 
optmethod <- "BFGS"
print(system.time(result_EMbest_REML <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR4, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, 
                                                             ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                             ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                             ~ 1 | EMF_ORIGIN_TED, ~ 1 | EMPLANT_7_ORIGINS, ~ 1 | origin_combo_TED7), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.04,0.16, 0.65, 0.03, 0.0001, 0.0001,0.0001, 0.13, 0.0001, 0.0001, 0.0001,0.0001,0.18,0.0001)))))
summary(result_EMbest_REML, digits=3)

cat("Marginal mean for STERILIZED = STERno  \n")
predict(result_EMbest_REML, newmods=rbind(c(0,0.5,0.5,0.5)))
cat("Marginal mean for STERILIZED = STERyes  \n")
predict(result_EMbest_REML, newmods=rbind(c(1,0.5,0.5,0.5)))
cat("Marginal mean for NONMYCOCONTROL2 = mics_not_added  \n")
predict(result_EMbest_REML, newmods=rbind(c(0.5,0,0.5,0.5)))
cat("Marginal mean for NONMYCOCONTROL2 = microbes_added  \n")
predict(result_EMbest_REML, newmods=rbind(c(0.5,1,0.5,0.5)))
cat("Marginal mean for FERTN = Nno  \n")
predict(result_EMbest_REML, newmods=rbind(c(0.5,0.5,0,0.5)))
cat("Marginal mean for FERTN = Nyes  \n")
predict(result_EMbest_REML, newmods=rbind(c(0.5,0.5,1,0.5)))
cat("Marginal mean for FERTP = Pno  \n")
predict(result_EMbest_REML, newmods=rbind(c(0.5,0.5,0.5,0)))
cat("Marginal mean for FERTP = Pyes  \n")
predict(result_EMbest_REML, newmods=rbind(c(0.5,0.5,0.5,1)))

### compute 'typical sampling variance' according to Higgins and Thompson (2002)
### see: http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
k <- result_EMbest_REML$k
wi <- 1/result_EMbest_REML$vi
s2 <- (k-1) * sum(wi) / (sum(wi)^2 - sum(wi^2))

cat("### marginal R^2 (fixed effects alone)\n")
R2GLMM.mar.EM <- var(fitted(result_EMbest_REML)) / (var(fitted(result_EMbest_REML)) + sum(result_EMbest_REML$sigma2) + s2)
R2GLMM.mar.EM

cat("### conditional R^2 (fixed plus random effects)\n")
R2GLMM.con.EM <- (var(fitted(result_EMbest_REML)) + sum(result_EMbest_REML$sigma2)) / (var(fitted(result_EMbest_REML)) + sum(result_EMbest_REML$sigma2) + s2)
R2GLMM.con.EM

cat("### R^2 attributable to random effects\n")
R2GLMM.con.EM-R2GLMM.mar.EM

cat("### partial R^2 for Plant Phylogeny\n")
R2GLMM.con.EM.plantphylo <- (result_EMbest_REML$sigma2[4]) / (var(fitted(result_EMbest_REML)) + sum(result_EMbest_REML$sigma2) + s2)
R2GLMM.con.EM.plantphylo

cat("### partial R^2 for Plant Phylogeny x Fungal Phylogeny interaction\n")
R2GLMM.con.EM.phybyphy <- (result_EMbest_REML$sigma2[8]) / (var(fitted(result_EMbest_REML)) + sum(result_EMbest_REML$sigma2) + s2)
R2GLMM.con.EM.phybyphy

cat("### partial R^2 for Plant Origin\n")
R2GLMM.con.EM.plantorigin <- (result_EMbest_REML$sigma2[13]) / (var(fitted(result_EMbest_REML)) + sum(result_EMbest_REML$sigma2) + s2)
R2GLMM.con.EM.plantorigin

cat("### partial R^2 for Plant Origin x Fungal Origin combinations\n")
R2GLMM.con.EM.origincombos <- (result_EMbest_REML$sigma2[14]) / (var(fitted(result_EMbest_REML)) + sum(result_EMbest_REML$sigma2) + s2)
R2GLMM.con.EM.origincombos

cat("### partial R^2 for Study ID \n")
R2GLMM.con.EM.ID <- (result_EMbest_REML$sigma2[1]) / (var(fitted(result_EMbest_REML)) + sum(result_EMbest_REML$sigma2) + s2)
R2GLMM.con.EM.ID

cat("### partial R^2 for Control Set \n")
R2GLMM.con.EM.CTRLSET <- (result_EMbest_REML$sigma2[2]) / (var(fitted(result_EMbest_REML)) + sum(result_EMbest_REML$sigma2) + s2)
R2GLMM.con.EM.CTRLSET

cat("### partial R^2 for Paper \n")
R2GLMM.con.EM.PAPER <- (result_EMbest_REML$sigma2[3]) / (var(fitted(result_EMbest_REML)) + sum(result_EMbest_REML$sigma2) + s2)
R2GLMM.con.EM.PAPER

sink()

####################
#################### Fitting (using REML) of random-only EM model, for estimation of overall weighted average effect size
####################
sink('Output_04_EMsub random-only model.txt')
cat("==EMsub_random_only ==\n")
optimizer <- "optim" 
optmethod <- "BFGS"
print(system.time(randonly_EM.REML <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR4, Rscale="cor0", method="REML", 
                                          random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, 
                                                        ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                        ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                        ~ 1 | EMF_ORIGIN_TED, ~ 1 | EMPLANT_7_ORIGINS, ~ 1 | origin_combo_TED7), 
                                          R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                          data=EMsub.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.0001,0.081, 0.511, 0.069, 0.0001, 0.0001,0.009, 0.0001, 0.0001, 0.0001, 0.0001,.0001,.18,.0001)))))
summary(randonly_EM.REML, digits=3)
cat("==Sum of variance components in random-only model for EM: ==\n")
varcomp.sum.EM.REML <- sum(randonly_EM.REML$sigma2)
varcomp.sum.EM.REML

overall.EM <- predict(randonly_EM.REML) # calculate overall weighted mean effect size from random-only model
overall.EM.pcnt <- 100*(exp(overall.EM$pred)-1)# transform to percent improvement in plant growth
overall.EM.pcnt.SE <- 100*((exp(overall.EM$se))-1)
overall.EM.pcnt
overall.EM.pcnt.SE 

sink()

###############################
# Obtain likelihood profiles and profile likelihood confidence intervals for variance components from BEST EM model, fit with REMLr
###############################
sink('Output_05 -- confint output for EMsub BEST REML model.txt')
par(mfrow=c(4,4))
for (i in 1:14) {
confint(result_EMbest_REML, sigma2=i, control=list(vc.min=0, vc.max=1.75))
profile(result_EMbest_REML, sigma2=i, xlim=c(0,1))
}
sink()

##################################################################################
############################ Calculate BLUPS for random effects in AIC-best models
##################################################################################
# library(devtools) # needed for install_github
# options(unzip = "unzip") # may be needed on Linux / RStudio before installng from github, but may cause problems on other platforms
# devtools::install_github("wviechtb/metafor") # installs development version of metafor
# get BLUPS for best EM model using REML
blups_EM <- ranef(result_EMbest_REML, verbose=TRUE)  # verbose=TRUE gives progress info
names(blups_EM) # returns a list of data frames with element names equal to the random effects
# For example, the 'plant phylogeny x fungal phylogeny' BLUPs are here:
# blups_EM$plant.phylxfungus.phyl
# First column for the BLUPs, then the SEs of the BLUPs, and then prediction interval bounds
# A histogram of the 'plant phylogeny x fungal phylogeny' BLUPs:
hist(blups_EM$plant.phylxfungus.phyl$intrcpt, breaks=50, col="lightgray")
# Or the standardized BLUPs:
hist(blups_EM$plant.phylxfungus.phyl$intrcpt / blups_EM$plant.phylxfungus.phyl$se, breaks=50, col="lightgray")

####################
#################### Create Figure 1, illustrating plant phylogeny x fungal phylogeny effect in the EM data
####################
# Make a plant tree containing only the EM plants
EMplants <- levels(EMsub.dat$PlantSpecies)
library(ape)
EMplant.tree<-drop.tip(tree.P,tree.P$tip.label[-match(EMplants, tree.P$tip.label)])
# edit EM plant tree to change node label 'eucalyptus' to 'myrtaceae':
EMplant.tree$node.label[19] <- "myrtaceae"

# Make a fungal tree containing only the EM fungi
EMfungi <- levels(EMsub.dat$FungalGenus)
EMfungi.tree<-drop.tip(tree.F,tree.F$tip.label[-match(EMfungi, tree.F$tip.label)])

# put the plant phylogeny x fungus phylogeny EM blups into a separate data frame:
plantxfungus_blups_EM <- data.frame(blups_EM$plant.phylxfungus.phyl)
plantxfungus_blups_EM$combo <- rownames(plantxfungus_blups_EM)
plantxfungus_blups_EM$pcnt <- 100*(exp(plantxfungus_blups_EM$intrcpt)-1) # transform blups to have value of percent benefit to plant growth
plantxfungus_blups_EM$value <- plantxfungus_blups_EM$pcnt+overall.EM.pcnt # center blups around overall weighted mean effect size
library(reshape2)
plantxfungus_blups_EM <- cbind(plantxfungus_blups_EM, colsplit(plantxfungus_blups_EM$combo, "XX", c("plant","fungus")))

# put EM Plant Origin blups into a separate data frame
plantorigin_blups_EM <- data.frame(blups_EM$EMPLANT_7_ORIGINS)
plantorigin_blups_EM$origin <- rownames(plantorigin_blups_EM)
plantorigin_blups_EM$pcnt <- 100*(exp(plantorigin_blups_EM$intrcpt)-1) # transform blups to have value of percent benefit to plant growth
plantorigin_blups_EM$pcnt_centered <- plantorigin_blups_EM$pcnt+overall.EM.pcnt
plantorigin_blups_EM

# install package for creating figure illustrating plant phylogeny x fungal phylogeny interaction for EM symbiosis:
# library(devtools) # needed for installing dualingTrees from github
# options(unzip = "unzip") # may be needed on Linux / RStudio before installng from github, but may cause problems on other platforms
# install_github("jfmeadow/dualingTrees-pkg", subdir="dualingTrees")
# optionally, save objects for making EM figure on another machine:
# save(list=c("EMplant.tree", "EMfungi.tree", "plantxfungus_blups_EM","plantorigin_blups_EM"),file = "EM_dualplot_objects.RData")

library(dualingTrees)

blup_trees_list1 <-
  input_trees(
    x_tree = EMfungi.tree,
    y_tree = EMplant.tree,
    x_key = plantxfungus_blups_EM$fungus,
    y_key = plantxfungus_blups_EM$plant,
    response = plantxfungus_blups_EM$pcnt,
    bubble_scale=1.5,
   #bubble_transform='sqrt',
    y_lab_cutoff = -2,
    pn_color_cutoff =  0,
    y_node_labs = c('dipterocarpaceae','phyllanthaceae','fabaceae','fagales','myrtaceae','pinaceae','salicaceae'),
    y_node_cex = plantorigin_blups_EM) # accepts a data.frame with row.names exactly matching a subset of tree internal nodes, and first column = values for point plotting (positive and negative colors as > or < 0)

# pdf('~/EM_dual_plot.pdf',
# width = 8, height = 8)
# par(new=FALSE)

plot_trees(
  trees_list = blup_trees_list1,
  pn_cols = c('#409ab1', '#dd1d18'),
  x_tree_col = 'gray40',
  x_bar_axis_offset = .7,
  y_bar_axis_offset = .8,
  y_tip_connect_length = 0.085,
  x_space = .3,
  y_space = .5,
  png_filename = 'EM_dual_plot.png',  # creates .png file in the working directory
  w_inches = 8,
  h_inches = 8,
  x_edge_width = 1,      # width of branches on left tree
  y_edge_width = 1,      # width of branches on top tree
  x_tip_width = 1,       # width of bars at left tree tips
  y_tip_width = 1,       # width of bars at top tree tips
  leg_text_pos = .35,
  leg_title = 'Plant response\nto EM fungi',
  pn_leg_labels = c('above 80.3%','below 80.3%'))

#optional text for labeling the Phyllanthaceae node:
#par(mar=c(0,0,0,0), mfrow=c(1,1), new=TRUE)
#plot(0, xlim=c(0,1), ylim=c(0,1), ann=FALSE, axes=FALSE, type='n')
#text(x=0.06, y=0.46, 'Phyllanthaceae', col=1, cex=.68)
#dev.off()

##############
##############
############## Fit EM models using variance estimators calculated using alternative imputation methods
##############
##############
# ESTVAR5 & ESTVAR6 were calculated using imputed SDs, which were 
# imputed using 2 different imputation methods (Bracken1999 and HotDeck_NN) in the R package metagear

sink('Output_06 --EMsub BEST model fit using ESTVAR5 and ESTVAR6.txt')
optimizer <- "optim" 
optmethod <- "BFGS"
cat("==EM reduced model, REML, Bracken imputation ==\n")
EMbest_REML_ESTVAR5 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR5, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR5, digits=3)

cat("==\n")
cat("==EM reduced model, REML, HotDeck_NN imputation ==\n")
# note that model fitting with ESTVAR6 uses EMsub2.dat, which has 998 observations (3 fewer than EMsub.dat, after dropping 3 ESTVAR6 outliers)

EMbest_REML_ESTVAR6.1 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.1, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR6.1, digits=3)

EMbest_REML_ESTVAR6.2 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.2, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7,
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR6.2, digits=3)

EMbest_REML_ESTVAR6.3 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.3, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR6.3, digits=3)

EMbest_REML_ESTVAR6.4 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.4, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7,
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR6.4, digits=3)

EMbest_REML_ESTVAR6.5 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.5, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7,
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR6.5, digits=3)

EMbest_REML_ESTVAR6.6 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.6, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7,
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR6.6, digits=3)

EMbest_REML_ESTVAR6.7 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.7, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR6.7, digits=3)

EMbest_REML_ESTVAR6.8 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.8, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7,
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR6.8, digits=3)

EMbest_REML_ESTVAR6.9 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.9, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR6.9, digits=3)

EMbest_REML_ESTVAR6.10 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.10, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP ,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7,
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMbest_REML_ESTVAR6.10, digits=3)

sink()


sink('Output_07 --EMsub SATURATED model fit using ESTVAR5 and ESTVAR6.txt')
cat("==EM saturated REML, Bracken imputation ==\n")
optimizer <- "optim" 
optmethod <- "BFGS"
EMsaturated_REML_ESTVAR5 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR5, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR5, digits=3)

cat("==\n")
cat("==EM saturated REML, HotDeck_NN imputation, 10 iterations ==\n")

EMsaturated_REML_ESTVAR6.1 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.1, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7,
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR6.1, digits=3)

EMsaturated_REML_ESTVAR6.2 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.2, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR6.2, digits=3)

EMsaturated_REML_ESTVAR6.3 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.3, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR6.3, digits=3)

EMsaturated_REML_ESTVAR6.4 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.4, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7,
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR6.4, digits=3)

EMsaturated_REML_ESTVAR6.5 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.5, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7,
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR6.5, digits=3)

EMsaturated_REML_ESTVAR6.6 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.6, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7,
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR6.6, digits=3)

EMsaturated_REML_ESTVAR6.7 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.7, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR6.7, digits=3)

EMsaturated_REML_ESTVAR6.8 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.8, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR6.8, digits=3)

EMsaturated_REML_ESTVAR6.9 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.9, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR6.9, digits=3)

EMsaturated_REML_ESTVAR6.10 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.10, Rscale="cor0", method="REML", 
                                               mods = ~ STERILIZED + NONMYCOCONTROL2 + FERTN + FERTP + LOCATION + FUNGROUP,
                                               random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, 
											~ 1 | EMPLANT_7_ORIGINS, ~ 1 | EMF_ORIGIN_TED, ~ 1 | origin_combo_TED7, 
											~ 1 | plant.phylxfungus.phyl, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus,
                                                            	 ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),
                                               R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.EM, plant_phylxfungus=R_PP_FT.EM, plantxfungus_phyl= R_PT_FP.EM), 
                                               data=EMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F,sigma2.init=c(0.0001,0.0001,0.0001,0.2,0.0001,0.0001,0.0001,0.1,0.0001,0.0001,0.0001,0.04,0.16,0.65)))
summary(EMsaturated_REML_ESTVAR6.10, digits=3)

sink()


#################
################# Bayesian (MCMCglmm) analyses of EM
#################
# Quoting Horvathova et al. 2011: " For the random effects of all the models, we used V = 1 and nu = 0.002, 
# which are equivalent to an inverse gamma prior, widely used in the statistical literature (Gelman & Hill 2007). "
# Gelman, A. & Hill, J. 2007 Data analysis using regression and multilevel/hierarchical models. Cambridge: Cambridge University Press.
# This same approach was used by Cornwallis et al. (Nature 466, no. 7309 (2010): 969), and is also followed here.

library(MCMCglmm)
library(MASS)

# construct inverse phylogeny matrices
R.P.inv <- as(ginv(R.P), "dgCMatrix")
dimnames(R.P.inv) <- dimnames(R.P)
R.F.inv <- as(ginv(R.F), "dgCMatrix")
dimnames(R.F.inv) <- dimnames(R.F)
RO.EM.inv <- as(ginv(RO.EM), "dgCMatrix")        
dimnames(RO.EM.inv) <- dimnames(RO.EM)
R_PP_FT.EM.inv <- as(ginv(R_PP_FT.EM), "dgCMatrix")  
dimnames(R_PP_FT.EM.inv) <- dimnames(R_PP_FT.EM)
R_PT_FP.EM.inv <- as(ginv(R_PT_FP.EM), "dgCMatrix")  
dimnames(R_PT_FP.EM.inv) <- dimnames(R_PT_FP.EM)

sink('Output_08_MCMCglmm EM models.txt')

cat("### EM saturated model, used for estimating fixed effects\n")
####################################################### EM Model with all fixed effects
prior.EM <- list(R=list(V=diag(1), nu=0.002), G=list(G1=list(V=diag(1),nu=0.002), G2=list(V=diag(1),nu=0.002), G3=list(V=diag(1),nu=0.002), 
                                                     G4=list(V=diag(1),nu=0.002), G5=list(V=diag(1),nu=0.002), G6=list(V=diag(1),nu=0.002), 
                                                     G7=list(V=diag(1),nu=0.002), G8=list(V=diag(1),nu=0.002), G9=list(V=diag(1),nu=0.002), 
                                                     G10=list(V=diag(1),nu=0.002), G11=list(V=diag(1),nu=0.002), G12=list(V=diag(1),nu=0.002),
                                                     G13=list(V=diag(1),nu=0.002), G14=list(V=diag(1),nu=0.002))) 


mcmcModel.EM <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~ Id + CTLTRTSETID + PaperID + PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus  + plantxfungus +
                       plant_phylxfungus + plantxfungus_phyl + plant.phylxfungus.phyl + EMF_ORIGIN_TED + EMPLANT_7_ORIGINS + origin_combo_TED7, 
                     mev = EMsub.dat$ESTVAR4, 
                     data=EMsub.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 
summary(mcmcModel.EM)

cat(" \n")
cat("### EM best model, used for estimating random effects comparable to those from maximum likelihood\n")
####################################################### EM Model with four fixed effects found important in likelihood model selection
mcmcModel_best.EM <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                         random = ~ Id + CTLTRTSETID + PaperID + PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus  + plantxfungus +
                           plant_phylxfungus + plantxfungus_phyl + plant.phylxfungus.phyl + EMF_ORIGIN_TED + EMPLANT_7_ORIGINS + origin_combo_TED7, 
                         mev = EMsub.dat$ESTVAR4, 
                         data=EMsub.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 
summary(mcmcModel_best.EM)

cat(" \n")
cat("### EM random-only model\n")
####################################################### EM Model with no fixed effects
mcmcModel.EM.nofixed <- MCMCglmm(fixed= EFFECTSIZE1 ~  1, 
                                 random = ~ Id + CTLTRTSETID + PaperID + PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus  + plantxfungus + 
                                   plant_phylxfungus + plantxfungus_phyl + plant.phylxfungus.phyl + EMF_ORIGIN_TED + EMPLANT_7_ORIGINS + origin_combo_TED7, 
                                 mev = EMsub.dat$ESTVAR4, 
                                 data=EMsub.dat, 
                                 nitt=100000, thin=5, burnin=25000, 
                                 prior = prior.EM,
                                 verbose = TRUE, singular.ok=TRUE,
                                 ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                                 plantxfungus_phyl=R_PT_FP.EM.inv , plant.phylxfungus.phyl = RO.EM.inv)) 
summary(mcmcModel.EM.nofixed)

sink()

######################################################## EM models using alternative imputations of SD
sink('Output_09_MCMCglmm EM models with alternative imputed SDs.txt')

prior.EM <- list(R=list(V=diag(1), nu=0.002), G=list(G1=list(V=diag(1),nu=0.002), G2=list(V=diag(1),nu=0.002), G3=list(V=diag(1),nu=0.002), 
                                                     G4=list(V=diag(1),nu=0.002), G5=list(V=diag(1),nu=0.002), G6=list(V=diag(1),nu=0.002), 
                                                     G7=list(V=diag(1),nu=0.002), G8=list(V=diag(1),nu=0.002), G9=list(V=diag(1),nu=0.002), 
                                                     G10=list(V=diag(1),nu=0.002), G11=list(V=diag(1),nu=0.002), G12=list(V=diag(1),nu=0.002),
                                                     G13=list(V=diag(1),nu=0.002), G14=list(V=diag(1),nu=0.002))) 

cat("### EM saturated model, ESTVAR5 (Bracken)\n")

EMsaturated_mcmc_Bracken <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                     mev = EMsub.dat$ESTVAR5, 
                     data=EMsub.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_Bracken)

cat(" \n")
cat("### EM reduced model, ESTVAR5 (Bracken)\n")
####################################################### EM Model with four fixed effects found important in likelihood model selection
EMbest_mcmc_Bracken <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                         mev = EMsub.dat$ESTVAR5, 
                         data=EMsub.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_Bracken)

cat(" \n")
cat("### EM saturated model, ESTVAR6 (HotDeck_NN), with 10 replicate imputed data sets \n") # note that these analyses are conducted with EMsub2.dat, which has 998 observations

EMsaturated_mcmc_ESTVAR6.1 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                     mev = EMsub2.dat$ESTVAR6.1, 
                     data=EMsub2.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_ESTVAR6.1)

EMsaturated_mcmc_ESTVAR6.2 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                     mev = EMsub2.dat$ESTVAR6.2, 
                     data=EMsub2.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_ESTVAR6.2)

EMsaturated_mcmc_ESTVAR6.3 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                     mev = EMsub2.dat$ESTVAR6.3, 
                     data=EMsub2.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_ESTVAR6.3)

EMsaturated_mcmc_ESTVAR6.4 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                     mev = EMsub2.dat$ESTVAR6.4, 
                     data=EMsub2.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_ESTVAR6.4)

EMsaturated_mcmc_ESTVAR6.5 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                     mev = EMsub2.dat$ESTVAR6.5, 
                     data=EMsub2.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_ESTVAR6.5)

EMsaturated_mcmc_ESTVAR6.6 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                     mev = EMsub2.dat$ESTVAR6.6, 
                     data=EMsub2.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_ESTVAR6.6)

EMsaturated_mcmc_ESTVAR6.7 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                     mev = EMsub2.dat$ESTVAR6.7, 
                     data=EMsub2.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_ESTVAR6.7)

EMsaturated_mcmc_ESTVAR6.8 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                     mev = EMsub2.dat$ESTVAR6.8, 
                     data=EMsub2.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_ESTVAR6.8)

EMsaturated_mcmc_ESTVAR6.9 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                     mev = EMsub2.dat$ESTVAR6.9, 
                     data=EMsub2.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_ESTVAR6.9)

EMsaturated_mcmc_ESTVAR6.10 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2 + LOCATION + FUNGROUP, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                     mev = EMsub2.dat$ESTVAR6.10, 
                     data=EMsub2.dat, 
                     nitt=100000, thin=5, burnin=25000, 
                     prior = prior.EM,
                     verbose = TRUE, singular.ok=TRUE,
                     ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                     plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMsaturated_mcmc_ESTVAR6.10)

cat(" \n")
cat("### EM reduced model, ESTVAR6 (HotDeck_NN), with 10 replicate imputed data sets\n")
####################################################### EM Model with four fixed effects found important in likelihood model selection
EMbest_mcmc_ESTVAR6.1 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                         mev = EMsub2.dat$ESTVAR6.1, 
                         data=EMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_ESTVAR6.1)

EMbest_mcmc_ESTVAR6.2 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                         mev = EMsub2.dat$ESTVAR6.2, 
                         data=EMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_ESTVAR6.2)

EMbest_mcmc_ESTVAR6.3 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                         mev = EMsub2.dat$ESTVAR6.3, 
                         data=EMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_ESTVAR6.3)

EMbest_mcmc_ESTVAR6.4 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                         mev = EMsub2.dat$ESTVAR6.4, 
                         data=EMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_ESTVAR6.4)

EMbest_mcmc_ESTVAR6.5 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                         mev = EMsub2.dat$ESTVAR6.5, 
                         data=EMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_ESTVAR6.5)

EMbest_mcmc_ESTVAR6.6 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                         mev = EMsub2.dat$ESTVAR6.6, 
                         data=EMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_ESTVAR6.6)

EMbest_mcmc_ESTVAR6.7 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                         mev = EMsub2.dat$ESTVAR6.7, 
                         data=EMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_ESTVAR6.7)

EMbest_mcmc_ESTVAR6.8 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                         mev = EMsub2.dat$ESTVAR6.8, 
                         data=EMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_ESTVAR6.8)

EMbest_mcmc_ESTVAR6.9 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                         mev = EMsub2.dat$ESTVAR6.9, 
                         data=EMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_ESTVAR6.9)

EMbest_mcmc_ESTVAR6.10 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN +  FERTP + STERILIZED + NONMYCOCONTROL2, 
                     random = ~  PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + EMPLANT_7_ORIGINS + EMF_ORIGIN_TED + origin_combo_TED7 + 
						plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,
                         mev = EMsub2.dat$ESTVAR6.10, 
                         data=EMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.EM,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.EM.inv, 
                                         plantxfungus_phyl=R_PT_FP.EM.inv, plant.phylxfungus.phyl = RO.EM.inv)) 

summary(EMbest_mcmc_ESTVAR6.10)

# calculate posterior SE (and print it with corresponding posterior mean) for variance components in models fit to data imputed by Bracken imputation:

for ( i in 1:14) {
  print(paste0("EMsaturated_Bracken_mean.",i))
  print(mean(EMsaturated_mcmc_Bracken$VCV[,i]))
  print(paste0("SE.",i))
  print(sqrt(var(EMsaturated_mcmc_Bracken$VCV[,i])))
}
for ( i in 1:14) {
  print(paste0("EMbest_Bracken_mean.",i))
  print(mean(EMbest_mcmc_Bracken$VCV[,i]))
  print(paste0("SE.",i))
  print(sqrt(var(EMbest_mcmc_Bracken$VCV[,i])))
}

# from Bayesian model fits of 10 replicate data sets imputed with HotDeck_NN, calculate aggregate
# mean and SE of all 14 variance component estimates across the 10 replicates using Rubin's rules, e.g., equations from pg. 96
# of Nakagawa, 2015, Chapter 4 from Ecological Statistics: Contemporary Theory and Application

for ( i in 1:14) {
	print(paste0("EMsaturated_HotDeck_NN_mean.",i))
	print(mean(c(mean(EMsaturated_mcmc_ESTVAR6.1$VCV[,i]), mean(EMsaturated_mcmc_ESTVAR6.2$VCV[,i]),mean(EMsaturated_mcmc_ESTVAR6.3$VCV[,i]),
			mean(EMsaturated_mcmc_ESTVAR6.4$VCV[,i]),mean(EMsaturated_mcmc_ESTVAR6.5$VCV[,i]),mean(EMsaturated_mcmc_ESTVAR6.6$VCV[,i]),
			mean(EMsaturated_mcmc_ESTVAR6.7$VCV[,i]),mean(EMsaturated_mcmc_ESTVAR6.8$VCV[,i]),mean(EMsaturated_mcmc_ESTVAR6.9$VCV[,i]),
			mean(EMsaturated_mcmc_ESTVAR6.10$VCV[,i]))))
	vW <- mean(c(var(EMsaturated_mcmc_ESTVAR6.1$VCV[,i]), var(EMsaturated_mcmc_ESTVAR6.2$VCV[,i]),var(EMsaturated_mcmc_ESTVAR6.3$VCV[,i]),
			var(EMsaturated_mcmc_ESTVAR6.4$VCV[,i]),var(EMsaturated_mcmc_ESTVAR6.5$VCV[,i]),var(EMsaturated_mcmc_ESTVAR6.6$VCV[,i]),
			var(EMsaturated_mcmc_ESTVAR6.7$VCV[,i]),var(EMsaturated_mcmc_ESTVAR6.8$VCV[,i]),var(EMsaturated_mcmc_ESTVAR6.9$VCV[,i]),
			var(EMsaturated_mcmc_ESTVAR6.10$VCV[,i])))
	vB <- var(c(mean(EMsaturated_mcmc_ESTVAR6.1$VCV[,i]), mean(EMsaturated_mcmc_ESTVAR6.2$VCV[,i]),mean(EMsaturated_mcmc_ESTVAR6.3$VCV[,i]),
			mean(EMsaturated_mcmc_ESTVAR6.4$VCV[,i]),mean(EMsaturated_mcmc_ESTVAR6.5$VCV[,i]),mean(EMsaturated_mcmc_ESTVAR6.6$VCV[,i]),
			mean(EMsaturated_mcmc_ESTVAR6.7$VCV[,i]),mean(EMsaturated_mcmc_ESTVAR6.8$VCV[,i]),mean(EMsaturated_mcmc_ESTVAR6.9$VCV[,i]),
			mean(EMsaturated_mcmc_ESTVAR6.10$VCV[,i])))
	vT <- vW + vB + vB/10
	SE <- sqrt(vT)
	print(paste0("SE.",i))
	print(SE)
}

for ( i in 1:14) {
	print(paste0("EMbest_HotDeck_NN_mean.",i))
	print(mean(c(mean(EMbest_mcmc_ESTVAR6.1$VCV[,i]), mean(EMbest_mcmc_ESTVAR6.2$VCV[,i]),mean(EMbest_mcmc_ESTVAR6.3$VCV[,i]),
			mean(EMbest_mcmc_ESTVAR6.4$VCV[,i]),mean(EMbest_mcmc_ESTVAR6.5$VCV[,i]),mean(EMbest_mcmc_ESTVAR6.6$VCV[,i]),
			mean(EMbest_mcmc_ESTVAR6.7$VCV[,i]),mean(EMbest_mcmc_ESTVAR6.8$VCV[,i]),mean(EMbest_mcmc_ESTVAR6.9$VCV[,i]),
			mean(EMbest_mcmc_ESTVAR6.10$VCV[,i]))))
	vW <- mean(c(var(EMbest_mcmc_ESTVAR6.1$VCV[,i]), var(EMbest_mcmc_ESTVAR6.2$VCV[,i]),var(EMbest_mcmc_ESTVAR6.3$VCV[,i]),
			var(EMbest_mcmc_ESTVAR6.4$VCV[,i]),var(EMbest_mcmc_ESTVAR6.5$VCV[,i]),var(EMbest_mcmc_ESTVAR6.6$VCV[,i]),
			var(EMbest_mcmc_ESTVAR6.7$VCV[,i]),var(EMbest_mcmc_ESTVAR6.8$VCV[,i]),var(EMbest_mcmc_ESTVAR6.9$VCV[,i]),
			var(EMbest_mcmc_ESTVAR6.10$VCV[,i])))
	vB <- var(c(mean(EMbest_mcmc_ESTVAR6.1$VCV[,i]), mean(EMbest_mcmc_ESTVAR6.2$VCV[,i]),mean(EMbest_mcmc_ESTVAR6.3$VCV[,i]),
			mean(EMbest_mcmc_ESTVAR6.4$VCV[,i]),mean(EMbest_mcmc_ESTVAR6.5$VCV[,i]),mean(EMbest_mcmc_ESTVAR6.6$VCV[,i]),
			mean(EMbest_mcmc_ESTVAR6.7$VCV[,i]),mean(EMbest_mcmc_ESTVAR6.8$VCV[,i]),mean(EMbest_mcmc_ESTVAR6.9$VCV[,i]),
			mean(EMbest_mcmc_ESTVAR6.10$VCV[,i])))
	vT <- vW + vB + vB/10
	SE <- sqrt(vT)
	print(paste0("SE.",i))
	print(SE)
}

sink()

###################
################### Funnel plots, diagnosing potential publication bias, using best EM model
###################

### set up 2x2 array for plotting
par(mfrow=c(2,2))

### draw funnel plots
funnel(result_EMbest_REML, pch=16, cex=.7, col=rgb(0,0,0,.1), main="Standard Error for EMsub")
funnel(result_EMbest_REML, pch=16, cex=.7, col=rgb(0,0,0,.1), yaxis="vi", main="Sampling Variance for EMsub")
funnel(result_EMbest_REML, pch=16, cex=.7, col=rgb(0,0,0,.1), yaxis="seinv", main="Inverse Standard Error for EMsub")
funnel(result_EMbest_REML, pch=16, cex=.7, col=rgb(0,0,0,.1),yaxis="vinv", main="Inverse Sampling Variance for EMsub")


##################
####################################### ANALYSIS OF AM data
##################

##################
################## MODEL SELECTION, using both REMLr and ML to compare results
##################
# Warning: The model selection procedures below may take from days to weeks, depending on processor speed and the number of processors being used
# It is recommended to fit individual models first and estimate per-model run time, then calculate number of total models using method=d in glmulti, and estimate total run time.

sink("Output_10 -- GLMULTI model selection, AMfull, exhaustive screening, ML.txt")
rma.glmulti.AMfull.1 <- function(formula, data, ...) {
  
  ### track progress via out file
  options(deparse.cutoff = 10000)
  cat(paste(capture.output(print(formula))[1], "\n"), file="tmp_exhaustive_AMfull_ML.out", append=TRUE)
  options(deparse.cutoff = 60)
  
  control <- list(optimizer="optim", optmethod="BFGS", sigma2.init=c(.12,.14,.18,.10,.10))
  res <- try(rma.mv(formula, V=ESTVAR4, method=ML,
                    random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies),
                    R = list(PlantSpecies.phyl=R.P), Rscale="cor0",
                    data = AMfull.dat, control = control, outlist = "minimal"), silent = TRUE)
  
  if (inherits(res, "try-error")) {
    control$optimizer <- "nlminb"
    res <- try(rma.mv(formula, V=ESTVAR4, method=ML,
                      random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies),
                      R = list(PlantSpecies.phyl=R.P), Rscale="cor0",
                      data = AMfull.dat, control = control, outlist = "minimal"), silent = TRUE)
  }
    return(res)
}
time.start <- proc.time()
options(warn=-1)
modelSel.AMfull.ml <- glmulti(y="EFFECTSIZE1",xr=c("FERTN","FERTP","STERILIZED","NONMYCOCONTROL2","FUNGROUP","LOCATION","DOMESTICATED","PLANTLIFEHISTORY","INOC.COMPLEXITY"),
                 exclude=c("FERTN:STERILIZED","FERTN:NONMYCOCONTROL2","FERTN:LOCATION","FERTN:DOMESTICATED","FERTN:PLANTLIFEHISTORY","FERTN:INOC.COMPLEXITY",
                           "FERTP:STERILIZED","FERTP:NONMYCOCONTROL2","FERTP:LOCATION","FERTP:DOMESTICATED","FERTP:PLANTLIFEHISTORY","FERTP:INOC.COMPLEXITY",
                           "STERILIZED:FUNGROUP","STERILIZED:LOCATION","STERILIZED:DOMESTICATED","STERILIZED:PLANTLIFEHISTORY","STERILIZED:INOC.COMPLEXITY",
                           "NONMYCOCONTROL2:FUNGROUP","NONMYCOCONTROL2:LOCATION","NONMYCOCONTROL2:DOMESTICATED","NONMYCOCONTROL2:PLANTLIFEHISTORY","NONMYCOCONTROL2:INOC.COMPLEXITY",
                           "FUNGROUP:LOCATION","FUNGROUP:DOMESTICATED","FUNGROUP:PLANTLIFEHISTORY","FUNGROUP:INOC.COMPLEXITY","LOCATION:DOMESTICATED",
                           "LOCATION:PLANTLIFEHISTORY","LOCATION:INOC.COMPLEXITY","DOMESTICATED:PLANTLIFEHISTORY","DOMESTICATED:INOC.COMPLEXITY","PLANTLIFEHISTORY:INOC.COMPLEXITY"),
                 maxsize=-1, method="h", data=AMfull.dat, level=2, fitfunction=rma.glmulti.AMfull.1, crit="aicc", confsetsize=1440, plotty=FALSE, estmethod=ML)
options(warn=1)
time.end <- proc.time()
cat("Minutes:", ((time.end - time.start)/60)[3], "\n")

print(modelSel.AMfull.ml)
extractRVI(modelSel.AMfull.ml)
weights10 <- weightable(modelSel.AMfull.ml)
weights10[1:20,] # list the best 20 models

sink()

save(modelSel.AMfull.ml, file="exhaustive_AMfull_ml.rdata"))


sink("Output_11 -- GLMULTI model selection, AMfull, exhaustive screening, REML.txt")
rma.glmulti.AMfull.2 <- function(formula, data, ...) {
  
  ### track progress via out file
  options(deparse.cutoff = 10000)
  cat(paste(capture.output(print(formula))[1], "\n"), file="tmp_exhaustive_AMfull_REML.out", append=TRUE)
  options(deparse.cutoff = 60)
  
  control <- list(optimizer="optim", optmethod="BFGS", REMLf=FALSE, sigma2.init=c(.12,.14,.18,.10,.10))
  res <- try(rma.mv(formula, V=ESTVAR4, method=REML,
                    random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies),
                    R = list(PlantSpecies.phyl=R.P), Rscale="cor0",
                    data = AMfull.dat, control = control, outlist = "minimal"), silent = TRUE)
  
  if (inherits(res, "try-error")) {
    control$optimizer <- "nlminb"
    res <- try(rma.mv(formula, V=ESTVAR4, method=REML,
                      random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies),
                      R = list(PlantSpecies.phyl=R.P), Rscale="cor0",
                      data = AMfull.dat, control = control, outlist = "minimal"), silent = TRUE)
  }
    return(res)
}
time.start <- proc.time()
options(warn=-1)
modelSel.AMfull.reml <- glmulti(y="EFFECTSIZE1",xr=c("FERTN","FERTP","STERILIZED","NONMYCOCONTROL2","FUNGROUP","LOCATION","DOMESTICATED","PLANTLIFEHISTORY","INOC.COMPLEXITY"),
                 exclude=c("FERTN:STERILIZED","FERTN:NONMYCOCONTROL2","FERTN:LOCATION","FERTN:DOMESTICATED","FERTN:PLANTLIFEHISTORY","FERTN:INOC.COMPLEXITY",
                           "FERTP:STERILIZED","FERTP:NONMYCOCONTROL2","FERTP:LOCATION","FERTP:DOMESTICATED","FERTP:PLANTLIFEHISTORY","FERTP:INOC.COMPLEXITY",
                           "STERILIZED:FUNGROUP","STERILIZED:LOCATION","STERILIZED:DOMESTICATED","STERILIZED:PLANTLIFEHISTORY","STERILIZED:INOC.COMPLEXITY",
                           "NONMYCOCONTROL2:FUNGROUP","NONMYCOCONTROL2:LOCATION","NONMYCOCONTROL2:DOMESTICATED","NONMYCOCONTROL2:PLANTLIFEHISTORY","NONMYCOCONTROL2:INOC.COMPLEXITY",
                           "FUNGROUP:LOCATION","FUNGROUP:DOMESTICATED","FUNGROUP:PLANTLIFEHISTORY","FUNGROUP:INOC.COMPLEXITY","LOCATION:DOMESTICATED",
                           "LOCATION:PLANTLIFEHISTORY","LOCATION:INOC.COMPLEXITY","DOMESTICATED:PLANTLIFEHISTORY","DOMESTICATED:INOC.COMPLEXITY","PLANTLIFEHISTORY:INOC.COMPLEXITY"),
                 maxsize=-1, method="h", data=AMfull.dat, level=2, fitfunction=rma.glmulti.AMfull.2, crit="aicc", confsetsize=1440, plotty=FALSE, estmethod=REML)
options(warn=1)
time.end <- proc.time()
cat("Minutes:", ((time.end - time.start)/60)[3], "\n")

print(modelSel.AMfull.reml)
extractRVI(modelSel.AMfull.reml)
weights11 <- weightable(modelSel.AMfull.reml)
weights11[1:20,] # list best top 20 models

sink()

save(modelSel.AMfull.reml, file="exhaustive_AMfull_reml.rdata"))


############################################################################

sink("Output_12 -- GLMULTI model selection, AMsub, exhaustive screening, ML.txt")

rma.glmulti.AMsub.1 <- function(formula, data, ...) {
  
  ### track progress via out file
  cat(paste(capture.output(print(formula))[1], "\n"), file="tmp_exhaustive_AMsub_ML.out", append=TRUE)
  control <- list(optimizer="optim", optmethod="BFGS", sigma2.init=c(.09,.16,.15,.02,.15,.035,.045))
  
  res <- try(rma.mv(formula, V=ESTVAR4, method=ML,
                    random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus),
                    R = list(PlantSpecies.phyl=R.P, plant_phylxfungus=R_PP_FT.AM), Rscale="cor0",
                    data = AMsub.dat, control = control, outlist = "minimal"), silent = TRUE)
  
  if (inherits(res, "try-error")) {
    control$optimizer <- "nlminb"
    res <- try(rma.mv(formula, V=ESTVAR4, method=ML,
                      random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus),
                      R = list(PlantSpecies.phyl=R.P, plant_phylxfungus=R_PP_FT.AM), Rscale="cor0",
                      data = data, control = control, outlist = "minimal"), silent = TRUE)
  }
  return(res)
}

time.start <- proc.time()
options(warn=-1)
modelSel.AMsub.ml <- glmulti(y="EFFECTSIZE1",xr=c("FERTN","FERTP","STERILIZED","NONMYCOCONTROL2","FUNGROUP","LOCATION","DOMESTICATED","PLANTLIFEHISTORY"),
                 exclude=c("FERTN:STERILIZED","FERTN:NONMYCOCONTROL2","FERTN:LOCATION","FERTN:DOMESTICATED","FERTN:PLANTLIFEHISTORY","FERTP:STERILIZED",
                           "FERTP:NONMYCOCONTROL2","FERTP:LOCATION","FERTP:DOMESTICATED","FERTP:PLANTLIFEHISTORY","STERILIZED:FUNGROUP",
                           "STERILIZED:LOCATION","STERILIZED:DOMESTICATED","STERILIZED:PLANTLIFEHISTORY","NONMYCOCONTROL2:FUNGROUP","NONMYCOCONTROL2:LOCATION",
                           "NONMYCOCONTROL2:DOMESTICATED","NONMYCOCONTROL2:PLANTLIFEHISTORY","FUNGROUP:LOCATION","FUNGROUP:DOMESTICATED",
                           "FUNGROUP:PLANTLIFEHISTORY","LOCATION:DOMESTICATED","LOCATION:PLANTLIFEHISTORY","DOMESTICATED:PLANTLIFEHISTORY"),
                 maxsize=-1, method="h", data=AMsub.dat, level=2, fitfunction=rma.glmulti.AMsub.1, crit="aicc", confsetsize=720, plotty=FALSE, estmethod=ML)
options(warn=1)
time.end <- proc.time()
cat("Minutes:", ((time.end - time.start)/60)[3], "\n")

print(modelSel.AMsub.ml)
extractRVI(modelSel.AMsub.ml)
weights12 <- weightable(modelSel.AMsub.ml)
weights12[1:21,] # list the best 21 models

sink()

save(modelSel.AMsub.ml, file="exhaustive_AMsub_ml.rdata"))


sink("Output_13 -- GLMULTI model selection, AMsub, exhaustive screening, REML.txt")

rma.glmulti.AMsub.2 <- function(formula, data, ...) {
  
  ### track progress via out file
  cat(paste(capture.output(print(formula))[1], "\n"), file="tmp_exhaustive_AMsub_REML.out", append=TRUE)
  control <- list(optimizer="optim", optmethod="BFGS", REMLf=FALSE, sigma2.init=c(.09,.16,.15,.02,.15,.035,.045))
  
  res <- try(rma.mv(formula, V=ESTVAR4, method=REML,
                    random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus),
                    R = list(PlantSpecies.phyl=R.P, plant_phylxfungus=R_PP_FT.AM), Rscale="cor0",
                    data = AMsub.dat, control = control, outlist = "minimal"), silent = TRUE)
  
  if (inherits(res, "try-error")) {
    control$optimizer <- "nlminb"
    res <- try(rma.mv(formula, V=ESTVAR4, method=REML,
                      random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus),
                      R = list(PlantSpecies.phyl=R.P, plant_phylxfungus=R_PP_FT.AM), Rscale="cor0",
                      data = data, control = control, outlist = "minimal"), silent = TRUE)
  }
  return(res)
}

time.start <- proc.time()
options(warn=-1)
modelSel.AMsub.reml <- glmulti(y="EFFECTSIZE1",xr=c("FERTN","FERTP","STERILIZED","NONMYCOCONTROL2","FUNGROUP","LOCATION","DOMESTICATED","PLANTLIFEHISTORY"),
                 exclude=c("FERTN:STERILIZED","FERTN:NONMYCOCONTROL2","FERTN:LOCATION","FERTN:DOMESTICATED","FERTN:PLANTLIFEHISTORY","FERTP:STERILIZED",
                           "FERTP:NONMYCOCONTROL2","FERTP:LOCATION","FERTP:DOMESTICATED","FERTP:PLANTLIFEHISTORY","STERILIZED:FUNGROUP",
                           "STERILIZED:LOCATION","STERILIZED:DOMESTICATED","STERILIZED:PLANTLIFEHISTORY","NONMYCOCONTROL2:FUNGROUP","NONMYCOCONTROL2:LOCATION",
                           "NONMYCOCONTROL2:DOMESTICATED","NONMYCOCONTROL2:PLANTLIFEHISTORY","FUNGROUP:LOCATION","FUNGROUP:DOMESTICATED",
                           "FUNGROUP:PLANTLIFEHISTORY","LOCATION:DOMESTICATED","LOCATION:PLANTLIFEHISTORY","DOMESTICATED:PLANTLIFEHISTORY"),
                 maxsize=-1, method="h", data=AMsub.dat, level=2, fitfunction=rma.glmulti.AMsub.2, crit="aicc", confsetsize=720, plotty=FALSE, estmethod=REML)
options(warn=1)
time.end <- proc.time()
cat("Minutes:", ((time.end - time.start)/60)[3], "\n")

print(modelSel.AMsub.reml)
extractRVI(modelSel.AMsub.reml)
weights13 <- weightable(modelSel.AMsub.reml)
weights13[1:21,] # list the best 21 models

sink()

save(modelSel.AMsub.reml, file="exhaustive_AMsub_reml.rdata"))

############################################################################

####################
#################### Fitting (using REML) of best AMfull models from model selection, for parameter estimation
####################

sink('Output_14--AMfull models, REML.txt')

cat("==Model AM_full_random-only ==\n")
cat("==Note: This was the AIC-best model according to REMLr model selection ==\n")
optimizer <- "optim" 
optmethod <- "BFGS"
print(system.time(randonly_AMfull.REML <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR4, Rscale="cor0", method="REML", 
                                        random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies),   
                                        R=list(PlantSpecies.phyl=R.P), 
                                        data=AMfull.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(.12,0.14, 0.18, 0.10, 0.16)))))
summary(randonly_AMfull.REML, digits=3)

overall.AMfull <- predict(randonly_AMfull.REML) # calculate overall weighted mean effect size, from random-only AMfull model
overall.AMfull.pcnt <- 100*((exp(overall.AMfull$pred))-1) # transform to percent improvement in plant growth
overall.AMfull.pcnt.SE <- 100*((exp(overall.AMfull$se))-1)
overall.AMfull
overall.AMfull.pcnt
overall.AMfull.pcnt.SE

### compute 'typical sampling variance' according to Higgins and Thompson (2002)
### see: http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate

k3 <- randonly_AMfull.REML$k
wi3 <- 1/randonly_AMfull.REML$vi
s23 <- (k3-1) * sum(wi3) / (sum(wi3)^2 - sum(wi3^2))

cat("### marginal R^2 (fixed effects alone)\n")
R2GLMM.mar.AMfull <- var(fitted(randonly_AMfull.REML)) / (var(fitted(randonly_AMfull.REML)) + sum(randonly_AMfull.REML$sigma2) + s23)
R2GLMM.mar.AMfull

cat("### conditional R^2 (fixed plus random effects)\n")
R2GLMM.con.AMfull <- (var(fitted(randonly_AMfull.REML)) + sum(randonly_AMfull.REML$sigma2)) / (var(fitted(randonly_AMfull.REML)) + sum(randonly_AMfull.REML$sigma2) + s23)
R2GLMM.con.AMfull

cat("### R^2 attributable to random effects\n")
R2GLMM.con.AMfull-R2GLMM.mar.AMfull

cat("### partial conditional R^2 for Plant Phylogeny\n")
R2GLMM.con.AMfull.plantphylo <- (randonly_AMfull.REML$sigma2[4]) / (var(fitted(randonly_AMfull.REML)) + sum(randonly_AMfull.REML$sigma2) + s23)
R2GLMM.con.AMfull.plantphylo

cat("### partial conditional R^2 for Plant Species\n")
R2GLMM.con.AMfull.plantsp <- (randonly_AMfull.REML$sigma2[5]) / (var(fitted(randonly_AMfull.REML)) + sum(randonly_AMfull.REML$sigma2) + s23)
R2GLMM.con.AMfull.plantsp

cat("### partial conditional R^2 for Id \n")
R2GLMM.con.AMfull.ID <- (randonly_AMfull.REML$sigma2[1]) / (var(fitted(randonly_AMfull.REML)) + sum(randonly_AMfull.REML$sigma2) + s23)
R2GLMM.con.AMfull.ID

cat("### partial conditional R^2 for Control Set \n")
R2GLMM.con.AMfull.CTRLSET <- (randonly_AMfull.REML$sigma2[2]) / (var(fitted(randonly_AMfull.REML)) + sum(randonly_AMfull.REML$sigma2) + s23)
R2GLMM.con.AMfull.CTRLSET

cat("### partial conditional R^2 for Paper \n")
R2GLMM.con.AMfull.PAPER <- (randonly_AMfull.REML$sigma2[3]) / (var(fitted(randonly_AMfull.REML)) + sum(randonly_AMfull.REML$sigma2) + s23)
R2GLMM.con.AMfull.PAPER

### shorthand code for splitting up conditional R^2 into its components:
c(var(fitted(randonly_AMfull.REML)), randonly_AMfull.REML$sigma2) / (var(fitted(randonly_AMfull.REML)) + sum(randonly_AMfull.REML$sigma2) + s23)

cat("\n")
cat("==AMfull BEST model from ML model selection fit here with REML ==\n")
optimizer <- "optim" 
optmethod <- "BFGS"
print(system.time(result_AMfull_best_ML <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR4, Rscale="cor0", method="REML", 
                                                    mods = ~ FERTP + STERILIZED + INOC.COMPLEXITY,
                                                    random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies), 
                                                    R=list(PlantSpecies.phyl=R.P), 
                                                    data=AMfull.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(.12,0.14, 0.18, 0.10, 0.16)))))
summary(result_AMfull_best_ML, digits=3)

k5 <- result_AMfull_best_ML$k
wi5 <- 1/result_AMfull_best_ML$vi
s25 <- (k5-1) * sum(wi5) / (sum(wi5)^2 - sum(wi5^2))

cat("### marginal R^2 (fixed effects alone)\n")
R2GLMM.mar.AMfull.LM <- var(fitted(result_AMfull_best_ML)) / (var(fitted(result_AMfull_best_ML)) + sum(result_AMfull_best_ML$sigma2) + s25)
R2GLMM.mar.AMfull.LM

cat("Marginal mean for FERTP = NO  \n")
predict(result_AMfull_best_ML, newmods=rbind(c(0,0.5,0.5,0.5)))
cat("Marginal mean for FERTP = YES  \n")
predict(result_AMfull_best_ML, newmods=rbind(c(1,0.5,0.5,0.5)))

cat("Marginal mean for STERILIZED = NO  \n")
predict(result_AMfull_best_ML, newmods=rbind(c(0.5,0,0.5,0.5)))
cat("Marginal mean for STERILIZED = YES  \n")
predict(result_AMfull_best_ML, newmods=rbind(c(0.5,1,0.5,0.5)))

cat("Marginal mean for INOC.COMPLEXITY = Single  \n")
predict(result_AMfull_best_ML, newmods=rbind(c(0.5,0.5,0,0)))
cat("Marginal mean for INOC.COMPLEXITY = Multi  \n")
predict(result_AMfull_best_ML, newmods=rbind(c(0.5,0.5,1,0)))
cat("Marginal mean for INOC.COMPLEXITY = Whole  \n")
predict(result_AMfull_best_ML, newmods=rbind(c(0.5,0.5,0,1)))

cat("==AMfull BEST model, plus FUNGROUP, fit with REML==\n")
# FUNGROUP had an RVI slightly greater than 0.5 in ML model selection
optimizer <- "optim" 
optmethod <- "BFGS"
print(system.time(result_AMfull_best_with_FUNGROUP_ML <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR4, Rscale="cor0", method="REML", 
                                                                  mods = ~ FUNGROUP + FERTP + STERILIZED + INOC.COMPLEXITY,
                                                                  random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies), 
                                                                  R=list(PlantSpecies.phyl=R.P), 
                                                                  data=AMfull.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(.12,0.14, 0.18, 0.10, 0.16)))))
summary(result_AMfull_best_with_FUNGROUP_REML, digits=3)

cat("Marginal mean for FUNGROUP = C3grass  \n")
predict(result_AMfull_best_with_FUNGROUP_ML, newmods=rbind(c(0,0,0,0,0,0.5,0.5,0.5,0.5)))
cat("Marginal mean for FUNGROUP = C4grass  \n")
predict(result_AMfull_best_with_FUNGROUP_ML, newmods=rbind(c(1,0,0,0,0,0.5,0.5,0.5,0.5)))
cat("Marginal mean for FUNGROUP = Nfixforb  \n")
predict(result_AMfull_best_with_FUNGROUP_ML, newmods=rbind(c(0,1,0,0,0,0.5,0.5,0.5,0.5)))
cat("Marginal mean for FUNGROUP = Nfixwood  \n")
predict(result_AMfull_best_with_FUNGROUP_ML, newmods=rbind(c(0,0,1,0,0,0.5,0.5,0.5,0.5)))
cat("Marginal mean for FUNGROUP = nonNforb  \n")
predict(result_AMfull_best_with_FUNGROUP_ML, newmods=rbind(c(0,0,0,1,0,0.5,0.5,0.5,0.5)))
cat("Marginal mean for FUNGROUP = nonNwood  \n")
predict(result_AMfull_best_with_FUNGROUP_ML, newmods=rbind(c(0,0,0,0,1,0.5,0.5,0.5,0.5)))

sink()


########################################

sink('Output_15--AMsub models, REML.txt')

cat("==Model AMsub_random_only ==\n")
# this was also the best AMsub model according to model selection using REMLr
optimizer <- "optim" 
optmethod <- "BFGS"
print(system.time(randonly_AMsub.REML <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR4, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, 
                                                     ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus), 
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.1,0.16, 0.16, 0.003, 0.16, 0.0001,0.0001, 0.0001, 0.06, 0.0001, 0.0001)))))
summary(randonly_AMsub.REML, digits=3)

overall.AMsub <- predict(randonly_AMsub.REML) # calculate overall weighted mean effect size, from random-only AMsub model
overall.AMsub.pcnt <- 100*((exp(overall.AMsub$pred))-1) # transform to percent increase in plant growth
overall.AMsub.pcnt.SE <- 100*((exp(overall.AMsub$se))-1)
overall.AMsub.pcnt
overall.AMsub.pcnt.SE

### compute 'typical sampling variance' according to Higgins and Thompson (2002)
### see: http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
k6 <- randonly_AMsub.REML$k
wi6 <- 1/randonly_AMsub.REML$vi
s26 <- (k6-1) * sum(wi6) / (sum(wi6)^2 - sum(wi6^2))

cat("### marginal R^2 (fixed effects alone)\n")
R2GLMM.mar.AMsub <- var(fitted(randonly_AMsub.REML)) / (var(fitted(randonly_AMsub.REML)) + sum(randonly_AMsub.REML$sigma2) + s26)
R2GLMM.mar.AMsub

cat("### conditional R^2 (fixed plus random effects)\n")
R2GLMM.con.AMsub <- (var(fitted(randonly_AMsub.REML)) + sum(randonly_AMsub.REML$sigma2)) / (var(fitted(randonly_AMsub.REML)) + sum(randonly_AMsub.REML$sigma2) + s26)
R2GLMM.con.AMsub

cat("### R^2 attributable to random effects\n")
R2GLMM.con.AMsub-R2GLMM.mar.AMsub

cat("### partial conditional R^2 for Plant Phylogeny\n")
R2GLMM.con.AMsub.plantphylo <- (randonly_AMsub.REML$sigma2[4]) / (var(fitted(randonly_AMsub.REML)) + sum(randonly_AMsub.REML$sigma2) + s26)
R2GLMM.con.AMsub.plantphylo

cat("### partial conditional R^2 for Plant Species\n")
R2GLMM.con.AMsub.plantsp <- (randonly_AMsub.REML$sigma2[5]) / (var(fitted(randonly_AMsub.REML)) + sum(randonly_AMsub.REML$sigma2) + s26)
R2GLMM.con.AMsub.plantsp

cat("### partial conditional R^2 for Plant Phylogeny x Fungal Genus Interaction\n")
R2GLMM.con.AMsub.plantsp <- (randonly_AMsub.REML$sigma2[9]) / (var(fitted(randonly_AMsub.REML)) + sum(randonly_AMsub.REML$sigma2) + s26)
R2GLMM.con.AMsub.plantsp

cat("### partial conditional R^2 for Id \n")
R2GLMM.con.AMsub.ID <- (randonly_AMsub.REML$sigma2[1]) / (var(fitted(randonly_AMsub.REML)) + sum(randonly_AMsub.REML$sigma2) + s26)
R2GLMM.con.AMsub.ID

cat("### partial conditional R^2 for Control Set \n")
R2GLMM.con.AMsub.CTRLSET <- (randonly_AMsub.REML$sigma2[2]) / (var(fitted(randonly_AMsub.REML)) + sum(randonly_AMsub.REML$sigma2) + s26)
R2GLMM.con.AMsub.CTRLSET

cat("### partial conditional R^2 for Paper \n")
R2GLMM.con.AMsub.PAPER <- (randonly_AMsub.REML$sigma2[3]) / (var(fitted(randonly_AMsub.REML)) + sum(randonly_AMsub.REML$sigma2) + s26)
R2GLMM.con.AMsub.PAPER


cat("==Model AMsub saturated, i.e., with all fixed effects ==\n")
optimizer <- "optim" 
optmethod <- "BFGS"
print(system.time(saturated_AMsub.REML <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR4, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                        random = list(~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID, ~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, 
                                                      ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                      ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus),   
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.1,0.16, 0.16, 0.003, 0.16, 0.0001,0.0001, 0.0001, 0.06, 0.0001, 0.0001)))))
summary(saturated_AMsub.REML, digits=3)
cat("==Sum of residual variance components in saturated mixed model for AMsub: ==\n")
varcomp.resid.AMsub.REML <- sum(saturated_AMsub.REML$sigma2)
varcomp.resid.AMsub.REML

sink()

##############
##############
############## Fit AMsub models using variance estimators calculated with SDs from alternative imputation methods
##############
##############

sink('Output_16 --AMsub BEST model fit using Bracken and HotDeck_NN imputation.txt')
optimizer <- "optim" 
optmethod <- "BFGS"
cat("==Model AMsub random-only, Bracken imputation ==\n")
AMbest_REML_ESTVAR5 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR5, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR5, digits=3)

cat("==\n")
cat("==Model AMsub random-only, HotDeck_NN imputation, 10 replicates ==\n") 
# note, these models are fit to the AMsub2.dat data, which had an outlier value of ESTVAR6 removed
# When model fitting was attempted with the full AMsub.dat data set, model fitting typically failed with 
# the error: "Ratio of largest to smallest sampling variance extremely large. Cannot obtain stable results."

AMbest_REML_ESTVAR6.1 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.1, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR6.1, digits=3)

AMbest_REML_ESTVAR6.2 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.2, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR6.2, digits=3)

AMbest_REML_ESTVAR6.3 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.3, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR6.3, digits=3)

AMbest_REML_ESTVAR6.4 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.4, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR6.4, digits=3)

AMbest_REML_ESTVAR6.5 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.5, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR6.5, digits=3)

AMbest_REML_ESTVAR6.6 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.6, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR6.6, digits=3)

AMbest_REML_ESTVAR6.7 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.7, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),  
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR6.7, digits=3)

AMbest_REML_ESTVAR6.8 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.8, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR6.8, digits=3)

AMbest_REML_ESTVAR6.9 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.9, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),  
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR6.9, digits=3)

AMbest_REML_ESTVAR6.10 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.10, Rscale="cor0", method="REML", 
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID), 
                                       R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                       data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMbest_REML_ESTVAR6.10, digits=3)

sink()


sink('Output_17 --AMsub SATURATED model fit using Bracken and HotDeck_NN imputation.txt')
optimizer <- "optim" 
optmethod <- "BFGS"
cat("==Model AMsub saturated, Bracken imputation ==\n")
AMsaturated_REML_ESTVAR5 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR5, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),  
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR5, digits=3)

cat("==\n")
cat("==Model AMsub saturated, HotDeck_NN imputation with 10 replicates ==\n")
# note that these models use AMsub2.dat data, which lacks the ESTVAR6 outlier
AMsaturated_REML_ESTVAR6.1 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.1, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),  
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR6.1, digits=3)

AMsaturated_REML_ESTVAR6.2 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.2, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),   
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR6.2, digits=3)

AMsaturated_REML_ESTVAR6.3 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.3, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),  
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR6.3, digits=3)

AMsaturated_REML_ESTVAR6.4 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.4, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),    
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR6.4, digits=3)

AMsaturated_REML_ESTVAR6.5 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.5, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),   
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR6.5, digits=3)

AMsaturated_REML_ESTVAR6.6 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.6, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),  
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR6.6, digits=3)

AMsaturated_REML_ESTVAR6.7 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.7, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),    
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR6.7, digits=3)

AMsaturated_REML_ESTVAR6.8 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.8, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),   
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR6.8, digits=3)

AMsaturated_REML_ESTVAR6.9 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.9, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),   
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR6.9, digits=3)

AMsaturated_REML_ESTVAR6.10 <- rma.mv(yi=EFFECTSIZE1, V=ESTVAR6.10, Rscale="cor0", method="REML", 
                                        mods = ~ FERTP*FERTN + FUNGROUP*FERTN + FUNGROUP*FERTP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,
                                       random = list(~ 1 | PlantSpecies.phyl, ~ 1 | PlantSpecies, ~ 1 | FungalGenus.phyl, ~ 1 | FungalGenus, ~ 1 | plant.phylxfungus.phyl, 
                                                     ~ 1 | plant_phylxfungus, ~ 1 | plantxfungus_phyl, ~ 1 | plantxfungus, ~ 1 | Id, ~ 1 | CTLTRTSETID, ~ 1 | PaperID),   
                                        R=list(PlantSpecies.phyl=R.P, FungalGenus.phyl=R.F, plant.phylxfungus.phyl=RO.AM, plant_phylxfungus=R_PP_FT.AM, plantxfungus_phyl= R_PT_FP.AM), 
                                        data=AMsub2.dat, control=list(optimizer=optimizer, optmethod=optmethod, REMLf=F, sigma2.init=c(0.01,0.1,0.0001,0.0001,0.0001,0.06,0.0001,0.0001,0.1,0.15,0.15)))
summary(AMsaturated_REML_ESTVAR6.10, digits=3)

sink()


###############################
# Obtain likelihood profiles and profile likelihood confidence intervals for variance components from BEST AM models, 
# according to REMLr model selection, fit with REMLr
###############################

sink('Output_18 -- confint output for AMsub BEST REML model.txt')
cat("AMsub ==\n")
par(mfrow=c(3,4))
for (i in 1:11) {
confint(randonly_AMsub.REML, sigma2=i, control=list(vc.min=0, vc.max=2))
profile(randonly_AMsub.REML, sigma2=i, xlim=c(0,1))}
sink()

sink('Output_19 -- confint output for AMfull BEST REML model.txt')
cat("AMfull ==\n")
par(mfrow=c(2,3))
for (i in 1:5) {
confint(randonly_AMfull.REML, sigma2=i, control=list(vc.min=0, vc.max=2))
profile(randonly_AMfull.REML, sigma2=i, xlim=c(0,1))}
sink()

cat("Alternative code for profiling variance components  \n")
print(system.time(profile.out1 <- profile(randonly_AMfull.REML, startmethod="prev", plot=FALSE)))
# this approach profiles all variance components, and the resulting object is a list where each element corresponds to one component
# Then, if we want the profile plot for the first component (for example), use:
plot(profile.out1[[1]]$sigma2, profile.out1[[1]]$ll, type="o", pch=19)
# Or if we want all profile plots at once, we can use:
par(mfrow=c(2,6))
invisible(lapply(profile.out1, function(x) plot(x$sigma2, x$ll, type="o", pch=19)))
sink()

############################
############################ Calculate BLUPS for random effects in AIC-best models
###########################
# get BLUPS for best (according to REMLr model selection) AMfull REML model
blups_AMfull.REML <- ranef(randonly_AMfull.REML, verbose=TRUE)  

# get BLUPS for best (according to REMLr model selection) AMsub REML model
blups_AMsub.REML <- ranef(randonly_AMsub.REML, verbose=TRUE)

####################
#################### Create Figure 2, illustrating plant phylogeny x fungal genus effect in the AMsub data
####################
# Make a plant tree containing only the AMsub plants
AMsub.plants <- levels(AMsub.dat$PlantSpecies)
library(ape)
AMsub.plant.tree<-drop.tip(tree.P,tree.P$tip.label[-match(AMsub.plants, tree.P$tip.label)])

# Make a fungal tree containing only the AMsub fungi
AMfungi <- levels(AMsub.dat$FungalGenus)
AMfungi.tree<-drop.tip(tree.F,tree.F$tip.label[-match(AMfungi, tree.F$tip.label)])

# put the plant phylogeny x fungal genus AMsub blups into a separate data frame:
plantxfungus_blups_AM <- data.frame(blups_AMsub.REML$plant_phylxfungus)
plantxfungus_blups_AM$combo <- rownames(plantxfungus_blups_AM)
plantxfungus_blups_AM$pcnt <- 100*(exp(plantxfungus_blups_AM$intrcpt)-1) # transform blups to have value of percent benefit to plant growth
plantxfungus_blups_AM$pcnt_centered <- plantxfungus_blups_AM$pcnt+overall.AMsub.pcnt # create column with blups centered around overall weighted mean effect size
library(reshape2)
plantxfungus_blups_AM <- cbind(plantxfungus_blups_AM, colsplit(plantxfungus_blups_AM$combo, "XX", c("plant","fungus")))

# save objects for making this figure on another machine:
# save(list=c("AMsub.plant.tree", "AMfungi.tree", "plantxfungus_blups_AM"),file = "AMsub_dualplot_objects.RData")

library(dualingTrees)

blup_trees_list2 <-
  input_trees(
    x_tree = AMfungi.tree,
    y_tree = AMsub.plant.tree,
    x_key = plantxfungus_blups_AM$fungus,
    y_key = plantxfungus_blups_AM$plant,
    response = plantxfungus_blups_AM$pcnt,
    bubble_scale=1.5,
   #bubble_transform='sqrt',
    y_lab_cutoff = 0,
    pn_color_cutoff = 0,
    y_node_labs = c('anacardiaceae','asteraceae','fabaceae','myrtaceae','poaceae','rosaceae'))

plot_trees(
  trees_list = blup_trees_list2,
  pn_cols = c('#409ab1', '#dd1d18'),
  x_tree_col = 'gray40',
  x_bar_axis_offset = .8,
  y_bar_axis_offset = .7,
  y_tip_connect_length = 0.075,
  x_space = .3,
  y_space = .5,
  y_lab_cex = 0.4,
  png_filename = 'AMsub_dual_plot.png',  # creates .png file in the working directory
  w_inches = 8,
  h_inches = 12,
  x_edge_width = .75, # width of branches on left tree
  y_edge_width = .75, # width of branches on top tree
  x_tip_width = .75,  # width of bars at left tree tips
  y_tip_width = .75,  # width of bars at top tree tips
  tall_layout = TRUE,
  leg_text_pos = .35,
  leg_title = 'Plant response\nto AM fungi',
  pn_leg_labels = c('above 62.0%','below 62.0%'))

#################
#################
################# Bayesian (MCMCglmm) analyses of AM symbiosis
#################
#################

library(MCMCglmm)
library(MASS)

# construct inverse phylogeny matrices
R.P.inv <- as(ginv(R.P), "dgCMatrix")
dimnames(R.P.inv) <- dimnames(R.P)
R.F.inv <- as(ginv(R.F), "dgCMatrix")
dimnames(R.F.inv) <- dimnames(R.F)
RO.AM.inv <- as(ginv(RO.AM), "dgCMatrix")        
dimnames(RO.AM.inv) <- dimnames(RO.AM)
R_PP_FT.AM.inv <- as(ginv(R_PP_FT.AM), "dgCMatrix")  
dimnames(R_PP_FT.AM.inv) <- dimnames(R_PP_FT.AM)
R_PT_FP.AM.inv <- as(ginv(R_PT_FP.AM), "dgCMatrix")  
dimnames(R_PT_FP.AM.inv) <- dimnames(R_PT_FP.AM)

sink('Output_20_MCMCglmm AM models.txt')

####################################################### AM-full Model saturated
cat("### \n")
cat("### AM-full saturated model\n")

prior.AMfull<-list(R=list(V=diag(1), nu=0.002), G=list(G1=list(V=diag(1),nu=0.002), G2=list(V=diag(1),nu=0.002), G3=list(V=diag(1),nu=0.002), 
                                           G4=list(V=diag(1),nu=0.002), G5=list(V=diag(1),nu=0.002))) 
Model.AMfull <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY + INOC.COMPLEXITY,  
                          random = ~ Id + CTLTRTSETID + PaperID + PlantSpecies.phyl + PlantSpecies, 
                          mev = AMfull.dat$ESTVAR4, 
                          data=AMfull.dat, 
                          nitt=100000, thin=5, burnin=25000, 
                          prior = prior.AMfull,
                          verbose = TRUE, singular.ok=TRUE,
                          ginverse = list(PlantSpecies.phyl  = R.P.inv)) 
summary(Model.AMfull)

cat("### \n")
cat("### AM-full random-only model\n")
####################################################### AM-full Model with no fixed effects
Model.AM.nofixed <- MCMCglmm(EFFECTSIZE1 ~  1,  
                             random = ~ Id + CTLTRTSETID + PaperID + PlantSpecies.phyl + PlantSpecies, 
                             mev = AMfull.dat$ESTVAR4, 
                             data=AMfull.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMfull,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv)) 
summary(Model.AM.nofixed)

cat("### \n")
cat("### AM-sub saturated model\n")
####################################################### AM-sub Model saturated
prior.AMsub <-list(R=list(V=diag(1), nu=0.002), G=list(G1=list(V=diag(1),nu=0.002), G2=list(V=diag(1),nu=0.002), G3=list(V=diag(1),nu=0.002), 
                                                            G4=list(V=diag(1),nu=0.002), G5=list(V=diag(1),nu=0.002), G6=list(V=diag(1),nu=0.002), 
                                                            G7=list(V=diag(1),nu=0.002), G8=list(V=diag(1),nu=0.002), G9=list(V=diag(1),nu=0.002), 
                                                            G10=list(V=diag(1),nu=0.002), G11=list(V=diag(1),nu=0.002))) 
Model.AMsub <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ Id + CTLTRTSETID + PaperID + PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plantxfungus + plant_phylxfungus + plantxfungus_phyl, 
                         mev = AMsub.dat$ESTVAR4, 
                         data=AMsub.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(Model.AMsub)

cat("### \n")
cat("### AM-sub random-only model\n")
####################################################### AM-sub Model with no fixed effects
Model.AMsub.nofixed <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                             random = ~ Id + CTLTRTSETID + PaperID + PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plantxfungus + plant_phylxfungus + plantxfungus_phyl  + plant.phylxfungus.phyl, 
                             mev = AMsub.dat$ESTVAR4, 
                             data=AMsub.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(Model.AMsub.nofixed)

sink()


#################
############################ AM-sub models using alternative imputations of SD
#################

sink('Output_21_MCMCglmm AM-sub models with alternative imputed SDs.txt')

prior.AMsub <-list(R=list(V=diag(1), nu=0.002), G=list(G1=list(V=diag(1),nu=0.002), G2=list(V=diag(1),nu=0.002), G3=list(V=diag(1),nu=0.002), 
                                                            G4=list(V=diag(1),nu=0.002), G5=list(V=diag(1),nu=0.002), G6=list(V=diag(1),nu=0.002), 
                                                            G7=list(V=diag(1),nu=0.002), G8=list(V=diag(1),nu=0.002), G9=list(V=diag(1),nu=0.002), 
                                                            G10=list(V=diag(1),nu=0.002), G11=list(V=diag(1),nu=0.002))) 
cat("### \n")
cat("### AM-sub saturated model, Bracken\n")
AMsub.saturated.mcmc.Bracken <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                         mev = AMsub.dat$ESTVAR5, 
                         data=AMsub.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.Bracken)

cat("### \n")
cat("### AM-sub random-only model, Bracken\n")
AMsub.nofixed.mcmc.Bracken <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                             mev = AMsub.dat$ESTVAR5, 
                             data=AMsub.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.Bracken)

cat("### \n")
cat("### AM-sub saturated model, HotDeck_NN imputation, 10 replicates\n")
AMsub.saturated.mcmc.ESTVAR6.1 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                         mev = AMsub2.dat$ESTVAR6.1, 
                         data=AMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.ESTVAR6.1)

AMsub.saturated.mcmc.ESTVAR6.2 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                         mev = AMsub2.dat$ESTVAR6.2, 
                         data=AMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.ESTVAR6.2)

AMsub.saturated.mcmc.ESTVAR6.3 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                         mev = AMsub2.dat$ESTVAR6.3, 
                         data=AMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.ESTVAR6.3)

AMsub.saturated.mcmc.ESTVAR6.4 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,  
                         mev = AMsub2.dat$ESTVAR6.4, 
                         data=AMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.ESTVAR6.4)

AMsub.saturated.mcmc.ESTVAR6.5 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                         mev = AMsub2.dat$ESTVAR6.5, 
                         data=AMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.ESTVAR6.5)

AMsub.saturated.mcmc.ESTVAR6.6 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,  
                         mev = AMsub2.dat$ESTVAR6.6, 
                         data=AMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.ESTVAR6.6)

AMsub.saturated.mcmc.ESTVAR6.7 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,  
                         mev = AMsub2.dat$ESTVAR6.7, 
                         data=AMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.ESTVAR6.7)

AMsub.saturated.mcmc.ESTVAR6.8 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                         mev = AMsub2.dat$ESTVAR6.8, 
                         data=AMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.ESTVAR6.8)

AMsub.saturated.mcmc.ESTVAR6.9 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                         mev = AMsub2.dat$ESTVAR6.9, 
                         data=AMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.ESTVAR6.9)

AMsub.saturated.mcmc.ESTVAR6.10 <- MCMCglmm(fixed= EFFECTSIZE1 ~  FERTN*FERTP + FERTN*FUNGROUP +  FERTP*FUNGROUP + STERILIZED*NONMYCOCONTROL2 + LOCATION + DOMESTICATED + PLANTLIFEHISTORY,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,  
                         mev = AMsub2.dat$ESTVAR6.10, 
                         data=AMsub2.dat, 
                         nitt=100000, thin=5, burnin=25000, 
                         prior = prior.AMsub,
                         verbose = TRUE, singular.ok=TRUE,
                         ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv, plant.phylxfungus.phyl = RO.AM.inv, plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv)) 
summary(AMsub.saturated.mcmc.ESTVAR6.10)

cat("### \n")
cat("### AM-sub random-only model, HotDeck_NN imputation, 10 replicates\n")
AMsub.nofixed.mcmc.ESTVAR6.1 <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,  
                             mev = AMsub2.dat$ESTVAR6.1, 
                             data=AMsub2.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.ESTVAR6.1)

AMsub.nofixed.mcmc.ESTVAR6.2 <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                             mev = AMsub2.dat$ESTVAR6.2, 
                             data=AMsub2.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.ESTVAR6.2)

AMsub.nofixed.mcmc.ESTVAR6.3 <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                             mev = AMsub2.dat$ESTVAR6.3, 
                             data=AMsub2.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.ESTVAR6.3)

AMsub.nofixed.mcmc.ESTVAR6.4 <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,  
                             mev = AMsub2.dat$ESTVAR6.4, 
                             data=AMsub2.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.ESTVAR6.4)

AMsub.nofixed.mcmc.ESTVAR6.5 <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,  
                             mev = AMsub2.dat$ESTVAR6.5, 
                             data=AMsub2.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.ESTVAR6.5)

AMsub.nofixed.mcmc.ESTVAR6.6 <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                             mev = AMsub2.dat$ESTVAR6.6, 
                             data=AMsub2.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.ESTVAR6.6)

AMsub.nofixed.mcmc.ESTVAR6.7 <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,  
                             mev = AMsub2.dat$ESTVAR6.7, 
                             data=AMsub2.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.ESTVAR6.7)

AMsub.nofixed.mcmc.ESTVAR6.8 <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,  
                             mev = AMsub2.dat$ESTVAR6.8, 
                             data=AMsub2.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.ESTVAR6.8)

AMsub.nofixed.mcmc.ESTVAR6.9 <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID, 
                             mev = AMsub2.dat$ESTVAR6.9, 
                             data=AMsub2.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.ESTVAR6.9)

AMsub.nofixed.mcmc.ESTVAR6.10 <- MCMCglmm(fixed= EFFECTSIZE1 ~  1,  
                         random = ~ PlantSpecies.phyl + PlantSpecies + FungalGenus.phyl + FungalGenus + plant.phylxfungus.phyl + plant_phylxfungus + plantxfungus_phyl + plantxfungus + Id + CTLTRTSETID + PaperID,  
                             mev = AMsub2.dat$ESTVAR6.10, 
                             data=AMsub2.dat, 
                             nitt=100000, thin=5, burnin=25000, 
                             prior = prior.AMsub,
                             verbose = TRUE, singular.ok=TRUE,
                             ginverse = list(PlantSpecies.phyl  = R.P.inv, FungalGenus.phyl  = R.F.inv,  plant_phylxfungus = R_PP_FT.AM.inv, plantxfungus_phyl=R_PT_FP.AM.inv, plant.phylxfungus.phyl = RO.AM.inv)) 
summary(AMsub.nofixed.mcmc.ESTVAR6.10)

################### calculate posterior SE (and print it with corresponding posterior mean) for variance components 
# in models fit to data imputed by Bracken imputation:

for ( i in 1:11) {
  print(paste0("AMsub.saturated.Bracken_mean.",i))
  print(mean(AMsub.saturated.mcmc.Bracken$VCV[,i]))
  print(paste0("SE.",i))
  print(sqrt(var(AMsub.saturated.mcmc.Bracken$VCV[,i])))
}
for ( i in 1:11) {
  print(paste0("AMsub.nofixed.Bracken_mean.",i))
  print(mean(AMsub.nofixed.mcmc.Bracken$VCV[,i]))
  print(paste0("SE.",i))
  print(sqrt(var(AMsub.nofixed.mcmc.Bracken$VCV[,i])))
}

# from Bayesian model fits of 10 replicate data sets imputed with HotDeck_NN, calculate aggregate
# mean and SE of all 11 variance component estimates across the 10 replicates using Rubin's rules, e.g., equations from pg. 96
# of Nakagawa, 2015, Chapter 4 from Ecological Statistics: Contemporary Theory and Application
for ( i in 1:11) {
	print(paste0("AMsub.saturated.HotDeck_NN.mean.",i))
	print(mean(c(mean(AMsub.saturated.mcmc.ESTVAR6.1$VCV[,i]), mean(AMsub.saturated.mcmc.ESTVAR6.2$VCV[,i]),mean(AMsub.saturated.mcmc.ESTVAR6.3$VCV[,i]),
			mean(AMsub.saturated.mcmc.ESTVAR6.4$VCV[,i]),mean(AMsub.saturated.mcmc.ESTVAR6.5$VCV[,i]),mean(AMsub.saturated.mcmc.ESTVAR6.6$VCV[,i]),
			mean(AMsub.saturated.mcmc.ESTVAR6.7$VCV[,i]),mean(AMsub.saturated.mcmc.ESTVAR6.8$VCV[,i]),mean(AMsub.saturated.mcmc.ESTVAR6.9$VCV[,i]),
			mean(AMsub.saturated.mcmc.ESTVAR6.10$VCV[,i]))))
	vW <- mean(c(var(AMsub.saturated.mcmc.ESTVAR6.1$VCV[,i]), var(AMsub.saturated.mcmc.ESTVAR6.2$VCV[,i]),var(AMsub.saturated.mcmc.ESTVAR6.3$VCV[,i]),
			var(AMsub.saturated.mcmc.ESTVAR6.4$VCV[,i]),var(AMsub.saturated.mcmc.ESTVAR6.5$VCV[,i]),var(AMsub.saturated.mcmc.ESTVAR6.6$VCV[,i]),
			var(AMsub.saturated.mcmc.ESTVAR6.7$VCV[,i]),var(AMsub.saturated.mcmc.ESTVAR6.8$VCV[,i]),var(AMsub.saturated.mcmc.ESTVAR6.9$VCV[,i]),
			var(AMsub.saturated.mcmc.ESTVAR6.10$VCV[,i])))
	vB <- var(c(mean(AMsub.saturated.mcmc.ESTVAR6.1$VCV[,i]), mean(AMsub.saturated.mcmc.ESTVAR6.2$VCV[,i]),mean(AMsub.saturated.mcmc.ESTVAR6.3$VCV[,i]),
			mean(AMsub.saturated.mcmc.ESTVAR6.4$VCV[,i]),mean(AMsub.saturated.mcmc.ESTVAR6.5$VCV[,i]),mean(AMsub.saturated.mcmc.ESTVAR6.6$VCV[,i]),
			mean(AMsub.saturated.mcmc.ESTVAR6.7$VCV[,i]),mean(AMsub.saturated.mcmc.ESTVAR6.8$VCV[,i]),mean(AMsub.saturated.mcmc.ESTVAR6.9$VCV[,i]),
			mean(AMsub.saturated.mcmc.ESTVAR6.10$VCV[,i])))
	vT <- vW + vB + vB/10
	SE <- sqrt(vT)
	print(paste0("SE.",i))
	print(SE)
}

for ( i in 1:11) {
	print(paste0("AMsub.nofixed.HotDeck_NN.mean.",i))
	print(mean(c(mean(AMsub.nofixed.mcmc.ESTVAR6.1$VCV[,i]), mean(AMsub.nofixed.mcmc.ESTVAR6.2$VCV[,i]),mean(AMsub.nofixed.mcmc.ESTVAR6.3$VCV[,i]),
			mean(AMsub.nofixed.mcmc.ESTVAR6.4$VCV[,i]),mean(AMsub.nofixed.mcmc.ESTVAR6.5$VCV[,i]),mean(AMsub.nofixed.mcmc.ESTVAR6.6$VCV[,i]),
			mean(AMsub.nofixed.mcmc.ESTVAR6.7$VCV[,i]),mean(AMsub.nofixed.mcmc.ESTVAR6.8$VCV[,i]),mean(AMsub.nofixed.mcmc.ESTVAR6.9$VCV[,i]),
			mean(AMsub.nofixed.mcmc.ESTVAR6.10$VCV[,i]))))
	vW <- mean(c(var(AMsub.nofixed.mcmc.ESTVAR6.1$VCV[,i]), var(AMsub.nofixed.mcmc.ESTVAR6.2$VCV[,i]),var(AMsub.nofixed.mcmc.ESTVAR6.3$VCV[,i]),
			var(AMsub.nofixed.mcmc.ESTVAR6.4$VCV[,i]),var(AMsub.nofixed.mcmc.ESTVAR6.5$VCV[,i]),var(AMsub.nofixed.mcmc.ESTVAR6.6$VCV[,i]),
			var(AMsub.nofixed.mcmc.ESTVAR6.7$VCV[,i]),var(AMsub.nofixed.mcmc.ESTVAR6.8$VCV[,i]),var(AMsub.nofixed.mcmc.ESTVAR6.9$VCV[,i]),
			var(AMsub.nofixed.mcmc.ESTVAR6.10$VCV[,i])))
	vB <- var(c(mean(AMsub.nofixed.mcmc.ESTVAR6.1$VCV[,i]), mean(AMsub.nofixed.mcmc.ESTVAR6.2$VCV[,i]),mean(AMsub.nofixed.mcmc.ESTVAR6.3$VCV[,i]),
			mean(AMsub.nofixed.mcmc.ESTVAR6.4$VCV[,i]),mean(AMsub.nofixed.mcmc.ESTVAR6.5$VCV[,i]),mean(AMsub.nofixed.mcmc.ESTVAR6.6$VCV[,i]),
			mean(AMsub.nofixed.mcmc.ESTVAR6.7$VCV[,i]),mean(AMsub.nofixed.mcmc.ESTVAR6.8$VCV[,i]),mean(AMsub.nofixed.mcmc.ESTVAR6.9$VCV[,i]),
			mean(AMsub.nofixed.mcmc.ESTVAR6.10$VCV[,i])))
	vT <- vW + vB + vB/10
	SE <- sqrt(vT)
	print(paste0("SE.",i))
	print(SE)
}

sink()

###########################################################################################
########### FUNNEL PLOTS USING BEST REML MODEL for AMfull
### set up 2x2 array for plotting
par(mfrow=c(2,2))

### draw funnel plots
funnel(randonly_AMfull.REML, pch=16, cex=.7, col=rgb(0,0,0,.1), main="Standard Error for AMfull")
funnel(randonly_AMfull.REML, pch=16, cex=.7, col=rgb(0,0,0,.1), yaxis="vi", main="Sampling Variance for AMfull")
funnel(randonly_AMfull.REML, pch=16, cex=.7, col=rgb(0,0,0,.1), yaxis="seinv", main="Inverse Standard Error for AMfull")
funnel(randonly_AMfull.REML, pch=16, cex=.7, col=rgb(0,0,0,.1), yaxis="vinv", main="Inverse Sampling Variance for AMfull")


