################################################################################
############################# SETTING THE SCENE ################################
################################################################################

# ------------------------------ EXPLANATION  ----------------------------------
# This is an R script explaining part of the analysis done in the manuscript:

# Rehren, J., Coll, M., Jiddawi,N., Kluger, L.C., Omar, O., Christensen, V.,  
# Pennino, M.G., Wolff, M. Evaluating ecosystem impacts of gear regulations in 
# a data-limited fishery using a temporal food web model. 
# ICES Journal of Marine Science. (submitted December 2021).

# The R script contains the analysis of the questionnaire data obtained 
# in February 2021 through interviews with Chwaka Bay fishers about their 
# fisheries catches in the past and the present. The result of this analysis was
# used for the parameterization of the vulnerability parameter in the Chwaka Bay 
# Ecosim model.

# ------------------------------------------------------------------------------

# --------------------------- LOAD PACKAGES AND FUNCTIONS ----------------------
# 1.) CLEAR SCREEN AND SET WORKING DIRECTORY ----
ls()
rm(list=ls())

# 2.) LOAD PACKAGES ----

library(dplyr)             # Data wrangling
library(spatstat)          # Calculating a weighted median

# ------------------------------------------------------------------------------

################################################################################
############################## LOAD AND PREP DATA ##############################
################################################################################

# 1.) LOAD THE DATA FRAME DATA ----

# This is the information about the demographics and fishing behaviour of each 
# participant
pDat <- read.table("Data/Information_participants.csv", sep=";", 
                   header=T)

# This data frame contains the catch values in 2021 and when the participant 
# started fishing
spcDat <- read.table("Data/Current and historical catch per species.csv", 
                  sep=";", header=T)

# This data frame contains information on total catch and effort in 2021 and
# when the participant started fishing

cDat <- read.table("Data/Total current and historical catch.csv", 
                   sep=";", header=T)

# These two data frames contain the proportions of taxonomic groups in the 
# functional groups 'Other herbivorous fish' and 'Other carnivorous fish'

carniProp <- read.table("Data/Other_carnivorous_fish.csv", sep=";", header=T)

herbiProp <- read.table("Data/Other_herbivorous_fish.csv", sep=";", header=T)

# 2.) SUBSET FOR EXPERT FISHERS ----

# We subset both data frames to keep only participants that have been fishing 
# for a minimum of 15 years.

pDatExp <- subset(pDat, pDat$Years_fishing_Chwaka>=15)
idExp <- pDatExp$Id
spCExp <- subset(spcDat, spcDat$Id %in% idExp)

# 3.) CORRECT FOR CHANGES IN FISHING EFFORT ----

# We correct the catch changes in fishing effort betwen current and historical 
# fishing periods

cDatL <- split(cDat, cDat$Time)                            # split for handling

cAdjust <- as.numeric(cDatL$Current$Duration_fishing_trip)/
  as.numeric(cDatL$Historical$Duration_fishing_trip)       # calculate the ratio
                                                           # between current and
                                                           # historical effort

cDatL$Current$Effort_change <- cAdjust                     # will be adjusted

# create Effort change variable in spCExp data frame
spCExp$Effort_change <- 1

# match the information of effort
spCExp$Effort_change <- cDatL$Current$Effort_change[match(spCExp$Id,
                                                          cDatL$Current$Id)]

# Set the Effort change for historical catches to 1 again as these will not be
# adjusted

spCExp$Effort_change[spCExp$Time=="Historical"] <- 1

# Correct the current catch by the effort change
spCExp$CatchNum <- as.numeric(spCExp$Catch)
spCExp$Catch_corrected <- (spCExp$CatchNum/spCExp$Effort_change)


################################################################################
################################## DATA ANALYSIS ###############################
################################################################################

# 1.) CALCULATE THE RATIO BETWEEN HISTORICAL AND CURRENT CATCH ----

spCExpL <- split(spCExp, spCExp$Time)

spCExpFin <- spCExpL$Current

spCExpFin$Ratio <- spCExpL$Historical$Catch_corrected/spCExpL$Current$Catch_corrected

# we add back the Historical catch and the fishing years

spCExpFin$Catch_Historical <- spCExpL$Historical$Catch_corrected
spCExpFin$Years_fishing <- 2021-spCExpL$Historical$Time_year


# We take out the two interviews in which the fisher has reported a change in 
# fishing effort (Effort_change!=1) because by correcting for it, we
# have altered the direction of catch change 

ind5 <- which(spCExpFin$Ratio<1 & spCExpFin$Id!="NJ0214C10")

removInt <- unique(spCExpFin[ind5,]$Id)

spCExpFin2 <- spCExpFin[!(spCExpFin$Id %in% removInt),]

# 2.) CALCULATE THE MEDIAN AND QUANTILES ----

# Because of the few entries for crabs and lobsters, we group them directly
# into the group "Crabs_Lobsters"

spCExpFin2$SpeciesCorr <- spCExpFin2$Species

spCExpFin2$SpeciesCorr[spCExpFin2$SpeciesCorr=="Lobsters"] <- "Crabs_Lobsters"
spCExpFin2$SpeciesCorr[spCExpFin2$SpeciesCorr=="Crabs"] <- "Crabs_Lobsters"

# and we calculate the median and quantiles for each taxonomic group.
spCExpFin3 <- spCExpFin2 %>% 
  group_by(SpeciesCorr) %>% 
  mutate(Ratio_median=median(Ratio, na.rm = T)) %>% 
  mutate(Q25Ratio=as.numeric(quantile(Ratio, na.rm=T, type = 1)[[2]])) %>% 
  mutate(Q75Ratio=as.numeric(quantile(Ratio, na.rm=T, type = 1)[[4]]))

spCExpFinAll <- spCExpFin3

# 3.) GROUP BY FUNCTIONAL GROUPS ----
 ## 3.1.) Other carnivorous fish ----

# taxonomic groups fromthe functional group other carnivorous fish for which
# information was collected:
carniNam <- c("Lethrinus_harak","Lethrinus_variegatus", "Lethrinus_mahsena",
              "Serranidae", "Mullidae", "Anguilliformes", "Cheilinus_trilobatus")

# we subset the data frame to contain only those
carniPropSub <- subset(carniProp, carniProp$Latin %in% carniNam)

# recalculate the proportion
carniPropSub$Group_contribution_new <- 1
carniPropSub$Group_contribution_new <- carniPropSub$Total/sum(carniPropSub$Total)

# and we subset the interview data frame
carni <- subset(spCExpFinAll, spCExpFinAll$SpeciesCorr %in% carniNam)

 ## 3.2.) Other herbivorous fish ----
# taxonomic groups fromthe functional group other carnivorous fish for which
# information was collected:
herbNam <- c("Acanthurus_spp", "Hipposcarus_harid")

# we subset the data frame to contain only those
herbiPropSub <- subset(herbiProp, herbiProp$Latin %in% herbNam)

# recalculate the proportion
herbiPropSub$Group_contribution_new <- 1
herbiPropSub$Group_contribution_new <- herbiPropSub$Total/sum(herbiPropSub$Total)

# and we subset the interview data frame
herb <- subset(spCExpFinAll, spCExpFinAll$SpeciesCorr %in% herbNam)

# 4.) CALCULATE MEDIAN AND QUANTILES PER FUNCTIONAL GROUP ----
 ## 4.1.) Other carnivorous fish ----

carni$contribution <- carniPropSub$Group_contribution[match(carni$SpeciesCorr, carniPropSub$Latin)]

colInt <- c("SpeciesCorr", "Ratio_median", "contribution")
uni <- unique(carni[,colInt])

weightedMedian <- weighted.median(uni$Ratio_median, uni$contribution)

carnivorous_fish <- data.frame(SpeciesCorr="Carnivorous_fish", 
                               Ratio_median=weightedMedian, 
                               Q25Ratio=quantile(uni$Ratio_median, na.rm=T, 
                                                 type = 1)[2], 
                               Q75Ratio=quantile(uni$Ratio_median, na.rm=T, 
                                                 type = 1)[4])

 ## 4.2.) Other herbivorous fish ----

herb$contribution <- herbiPropSub$Group_contribution[match(herb$SpeciesCorr, 
                                                           herbiPropSub$Latin)]

uniherbi <- unique(herb[,colInt])

weightedMedianHerbi <- weighted.median(uniherbi$Ratio_median, 
                                       uniherbi$contribution)

herbivorous_fish <- data.frame(SpeciesCorr="Herbivorous_fish", 
                               Ratio_median=weightedMedianHerbi, 
                               Q25Ratio=quantile(uniherbi$Ratio_median, 
                                                 na.rm=T, type = 1)[2], 
                               Q75Ratio=quantile(uniherbi$Ratio_median, 
                                                 na.rm=T, type = 1)[4])
 
# 5.) BIND ALL FUNCTIONAL GROUPS ----

# get the data for the other functional groups
spNam <- c("Siganus_sutor", "Leptoscarus_vaigiensis", "Lethrinus_lentjan", 
              "Lethrinus_borbonicus", "Lutjanus_fulviflamma", "Scarus_ghobban",
           "Crabs_Lobsters", "Octopus", "Squids")
sp <- subset(spCExpFinAll, spCExpFinAll$SpeciesCorr %in% spNam)

# we subset them
colInt2 <- c("SpeciesCorr", "Ratio_median", "Q25Ratio", "Q75Ratio")
sp <- unique(sp[,colInt2])

# and we bind all functional groups in one df
dfFinal <- rbind(sp, carnivorous_fish, herbivorous_fish)
