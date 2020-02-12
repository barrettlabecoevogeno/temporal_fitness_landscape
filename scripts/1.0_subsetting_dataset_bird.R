##########################################################################################
# Created by Marc-Olivier Beausoleil
# McGill University 
# Created Tuesday, February 12, 2020 
# Why: 
# Requires 
# NOTES: 
##########################################################################################

# Define variables --------------------------------------------------------
pdf = FALSE 
misidentification.fixing = TRUE 

# Load libraries ----------------------------------------------------------
library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)

# Source scripts ----------------------------------------------------------
source('scripts/custom_PCA.R')
source('scripts/0.1_misc.R')

# New dataset, or sorting informations ------------------------------------
whole.data = read.csv(file = "data/raw_data/DF_1999-2019.csv", header = TRUE, sep = ",")
col.numeric = c("Beak.LengthA","Beak.DepthA","Beak.WidthA",
                "Tarsus","Wing.Chord","Mass",
                "MedianBeakLength","MedianBeakWidth","MedianBeakDepth")
# Make the columns as numeric (it's going to generate NAs)
whole.data[, col.numeric] = apply(whole.data[,col.numeric], 2, 
                                  function(x) as.numeric(as.character(x)))

# Define AGE in finches as the sex (which was used by the team as an indication of maturity so age)
whole.data$age = as.character(whole.data$Sex0)
whole.data$age <- plyr::revalue(x = whole.data$age, 
                                replace = c(#""=NA, 
                                  "f"="f", # Adults 
                                  "m"="m", # Adults 
                                  "j"="j", "m;j"="j","f;j"="j","f?;j"="j", # Juvenile 
                                  "m?"=NA,"f?"=NA, # NA 
                                  "fledgling"="fldg", 
                                  "nestling"="nestl"))
# There are some weird birds that should be considered in the analysis, here is a list of some of them: 
# Misidentification -------------------------------------------------------
if (misidentification.fixing) {
  # Correct species name 
  whole.data[whole.data$Species1 %in% "fortis/magnirostris","Species1"] <- "fortis;magnirostris"
  ## Seems to be an outlier of the fortis group
  whole.data = whole.data[whole.data$BANDFINAL!="KGSK2106",]
  # I decided to change the species for these guys because they don't cluster properly otherwise 
  whole.data[whole.data$BANDFINAL %in% c("JP1087","SH157_S","UM174","SH123_S"),"Species1"] <- "fortis"
  whole.data[whole.data$BANDFINAL %in% c("JP3846","JP4935","JP737","SH1137"),  "Species1"] <- "fuliginosa"
  whole.data[whole.data$BANDFINAL == "KGSK2147" & whole.data$Species1 == "fuliginosa","Species1"] <- "fortis"
  whole.data[whole.data$BANDFINAL == "JP819","Species1"]<-"fuliginosa"
  whole.data[whole.data$BANDFINAL == "JP3765","Species1"]<-"scandens"
  list.misidentified.magnirostris.fortis = c("JH1A459", "JP1589", "JP2528", "JP2541", "JP3400", "JP3423", "JP3497", "JP3532", 
                                             "JP3598", "JP3607", "JP3626", "JP3656", "JP3665", "JP3818", "JP3860", "JP3865", 
                                             "JP3890", "JP4431", "JP4439", "JP4444", "JP4453", "JP4517", "JP4548", "JP4744", 
                                             "JP4958", "P739", "PP1A257", "SH1069_L","SH1102_L", "SH1110", "SH1143", "SH1145_L", 
                                             "SH1154_L", "SH192_L", "SH602_L", "SH978_L")
  whole.data[whole.data$BANDFINAL %in% list.misidentified.magnirostris.fortis,"Species1"] <- "fortis"
  # This individual is peculiar! It is a migrant from AB (2004) to EG (2006)
  # whole.data[whole.data$BANDFINAL == "JP1014",]
  whole.data = whole.data[whole.data$BANDFINAL != "JP1014",]
  
  # Notes on weird birds ----------------------------------------------------
  # what do you want to see
  # diagnose = c("BANDFINAL","Species1","Sex0","Year","Notes")
  # The list:
  # list.weird = c("JP9001", "P339", "JP714", "JP3114", "P585", "JP5168", "JP1781",
  #                "JP3042", "SH118", "JP3053", "JP5029", "JP2462", "JP1093",
  #                "P560", "JP5167", "203", "SH675", "JP1724", "SH410", "JP3051",
  #                "JP1484", "JP3574", "SH543", "SH405", "JP000", "JP023",
  #                "JP090", "JP008", "JP101", "SH424", "JP3472",
  #                "JP1386bis", "JP3410", "JP3401", "JP2363", "JP3392",
  #                "JP3351", "SH099", "JP3431", "JP3432", "SH896",
  #                "SH602", "JP3112")
}

# Subset Database for analysis --------------------------------------------
datata = whole.data %>%
  dplyr::filter(Sex0 %in% c("", "f", "f;j", "f?", 
                            "f?;j", "fledgling", "j", "m", "m;j",
                            "m?", "nestling"), # For the survival, I should also consider the individuals that were juveniles but became adults, but only taking the adult traits
                Year %in% c(1999:2019),
                !(BANDFINAL %in% ""), # Remove empty bands 
                Site %in% c("El Garrapatero"), # Keep only 1 site
                Species1 %in% c("fortis", 
                                # "fortis;fuliginosa","fortis;magnirostris","fortis;scandens", 
                                "fuliginosa", "magnirostris","scandens")) %>% # Keep ground finches and species limits 
  droplevels(.)
datata$BANDFINAL = as.character(datata$BANDFINAL)
datata$BANDFINAL_Year = paste(datata$BANDFINAL,datata$Year, sep = "_")

# Remove NAs of median beak traits --------------------------------------------------------------
# This line of code will remove ALL individuals that have AT LEAST 1 NA OR MORE FOR THE BEAK TRAITS 
datata = datata[-which(apply(is.na(datata[,c("MedianBeakLength", "MedianBeakWidth", "MedianBeakDepth")]),1,sum) >= 1),]
# The commands below will find the NAs in the selected columns and I'll remove them... 
# I won't run this at the moment, but might be needed in the future 
# data.find.na = apply(datata[,c("Tarsus", "Wing.Chord", "Mass")],2,function(x) which(is.na(x)))
# datata.no.na = datata[-unique(unlist(data.find.na)),cols]


# Recapture matrix --------------------------------------------------------
recap = datata %>%  # Make a DB 
  group_by(Year) %>% # Group by year to count the captures in various years
  dplyr::select(BANDFINAL) %>% # Select BANDINFAL (UNIQUE ID of finches)
  table() %>% # Make a table of occurance 
  t() %>% # Transpose the matrix
  decostand(method = "pa") %>% # Make 0s and 1s
  as.data.frame.matrix()

ch = paste('y', unique(sort(datata$Year)), sep='.')

# rename the columns of the recapture history 
colnames(recap) <- ch

# Create a Maxseen column 
recap$maxseen <- apply(recap[,ch],1,sum)

# Create a column with the band numbers (used to merge the data)
recap$BANDFINAL = row.names(recap)

# Create a capture history (ch) column
recap$ch = tidyr::unite(data = recap[,ch], 
                        col = ch, sep = "")
# Creates column with the first instance of capture
recap$first = apply(recap[,ch], 1, which.max) 
## update "first" so that it correctly reflects selected years (other way of coding it )
# recap$first <- apply(recap[,ch], 1, function(x) min(which(x==1)))

cbind(test1,test2)
# Creates column with the last instance of capture
recap$last = apply(recap[,ch], 1, function(x) max(which(x==1))) 


# Merge the recapture history to the main data  ---------------------------
datata2 = merge(x = datata,
                y = recap, 
                by.x = "BANDFINAL",
                by.y = "BANDFINAL")

# Apparent survival -------------------------------------------------------
# Correct the capture history (APPARENT SURVIVAL). May take several seconds 
ch.corr = known.state.cjs(datata2[,ch])

# Make recapture history with one column name 
datata2$X = as.matrix(datata2[,ch])
datata2$X.corr = as.matrix(ch.corr)

# Create variable that shows the APPARENT SURVIVAL 
datata2$maxseen.corr = apply(ch.corr,1,sum)
max(datata2$maxseen.corr)

# Removing duplicated entries in the original data ------------------------
# Select columns 
cols = c("BANDFINAL", "BANDFINAL_Year", "Date", "Year","Site", "Sex0", "Species1", 
         "Tarsus", "Wing.Chord", "Mass", "age","Notes","By.",
         "MedianBeakLength", "MedianBeakWidth","MedianBeakDepth")
# Before computing the PCA, you should remove the duplicates 
# But for the analysis of the fitness landscapes in various years, the duplicates will be reused 
datata3 = (datata2[!duplicated(datata2$BANDFINAL), ]) 

# PCA Beak  ---------------------------------------------------------------
# Calculate the PCA on Median Beak traits 
pca.all = vegan::rda(datata3[,c("MedianBeakLength", "MedianBeakWidth", "MedianBeakDepth")])

# Record the PCA scores 
pca.all.sco = scores(pca.all)

# Add PCA scores to the data 
datata3$PC1 = pca.all.sco$sites[,1]
datata3$PC2 = pca.all.sco$sites[,2]

# Calculate standard deviation of the pca scores 
sd.pc1 = sd(datata3$PC1)
sd.pc2 = sd(datata3$PC2)

# Standardized the pca scores 
datata3$PC1.std = datata3$PC1/sd.pc1
datata3$PC2.std = datata3$PC2/sd.pc2

# The datata and datata2 contains DUPLICATED BAND NUMBERS 
dim(datata) # Contains duplicated BANDFINAL
dim(datata2)# Contains duplicated BANDFINAL
dim(datata3)# NO duplicated BANDFINAL


# Generate output data ----------------------------------------------------
bird.data = datata3 %>% droplevels()
yr.list.subset = unique(sort(bird.data$Year))
save(bird.data,yr.list.subset,file = "output/bird.data.RData")
write.csv(bird.data,"output/bird.data.csv")

table(bird.data$Species1,
      bird.data$Year)

