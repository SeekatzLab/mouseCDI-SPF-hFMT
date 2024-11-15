### HUMOFMT/allhumos6B: Combining all SCFA and BA data (in-lab HPLC data and UM Metabolomics core data)
### 5.6.24
### Anna M. Seekatz

### based on following analyses done previously:
	- HPLC data:
		- PreValidation (includes normalization of 2 most recent plates, HF9 (2 plates), HF10/GF3 (4 plates)): /Users/annaseekatz/Box Sync/Projects/HUMO.w.RM/allhumos6/allhumos6B/SCFA_PreValidation/SCFA_all_combining.R
		- PostVal: (same, but after validating some of the peaks in A2): /Users/annaseekatz/Box Sync/Projects/HUMO.w.RM/allhumos6/allhumos6B/SCFA_PostVal
		- notes on plates: /Users/annaseekatz/Box Sync/Projects/HUMO.w.RM/HPLC/HPLC_Dec2019Edited.R
		- PREVIOUS combinining of plates HF34, HF5, HF7/GF2: /Users/annaseekatz/Box Sync/Projects/HUMO.w.RM/allhumos6/allhumos6B/flux/allhumos6B_FileCheck.R (lines 846+)
	- UM core (cecal) data:
		- allhumos4/metabcore_analysis/humofmt.metabcore_analysis.R
	
### inputs and original file information, per plate:
	- /Data/allhf_runs/allhf_seq.meta.merge.txt: used for metadata
	- HPLC data: /Data/metabolomics/HPLC/HPLC_RawData
		- /scfa1: HF34 data
			- lcm_samplekey.txt
			- lcm.run1_data_05.27.16.txt
			- lcm.run2_data_05.27.16.txt
		- /scfa2: HF5 data (part of another project)
			- HF5plA_all.compounds_10.05.17.txt
			- scfa2_samplekey.txt
		- /scfa3: GF2.HF7.HF8 (contains data for current and other projects)
			- scfa3_samplekey.txt
			- HF7.GF2_all.compounds_08.01.18.txt
		- /scf4: H9 (2 plates) and HF10/GF3 (4 plates) (also contains data for other project)
			- HF9.HF10_all.compounds.raw.txt
			- scfa4.5_samplekey.txt
	- UM core data (cecal): /Data/metabolomics/UM-MetaboCore
		- humofmt_metabcore_scfa.data.txt: SCFA measurements
		- BAs_weightnorm.modified.txt: BA measurements
		- this data was normalized by the core

### output information:
	- /Users/annaseekatz/Box Sync/Projects/HUMO.w.RM/Manuscript/Code/SCFA_ReAnalysis:
	- /Data/metabolomics/HPLC/HPLC_RawData: normalized data output
		- scfa1/scfa1_all.concentrations.txt
		- scfa2/scfa2_all.concentrations.txt (note: these were part of another project)
		- scfa3/scfa3_all.concentrations.txt
		- scfa4.5/scfa4.5_all.concentrations.txt
	- Data/metabolomics/HPLC/:
		- allhf_scfa.seq.meta.merge.txt: final product used for fecal scfa analysis (merged with meta)
	- Data/metabolomics/UM-MetabCore:
		- HFm_scfa.ba.cecal.seq.meta.merge.txt: all UM cecal measurements
	
		
# wd: /Data/metabolomics/HPLC
		
``{r}

####-------------------#### HF3.4 plate:

library(dplyr)
library(tidyverse)
library(tidyr)
library(data.table)
#library(purrr)
library(readr)

### read in and recopy rawdata:
# samplekey:
key<-read.table(file="HPLC_RawData/scfa1/lcm_samplekey.txt", header=TRUE)
# data from run1:
lcm1<-read.table(file="HPLC_RawData/scfa1/lcm.run1_data_05.27.16.txt", header=TRUE)
lcm2<-read.table(file="HPLC_RawData/scfa1/lcm.run2_data_05.27.16.txt", header=TRUE)

### re-normalize by run:
## run1:
key <- read.table(file="HPLC_RawData/scfa1/scfa1_samplekey.txt", header=TRUE)
lcm1 <- read.table(file="HPLC_RawData/scfa1/lcm.run1_data_05.27.16.txt", header=TRUE)

# define standards and merge with data to get standard curve:
standards<-key[key$sample_type==c("standard"), c("lcm_sampleID", "sample_name", "concentration_mM")]
standards<-standards[!is.na(standards$concentration_mM), ]
st.data<-merge(standards, lcm1, by.x="lcm_sampleID", by.y="Sample_ID")

# test out one compound:
butyrate<-st.data[st.data$compound==c("butyrate"), ]
lm1.butyrate<-lm(concentration_mM~Area, data=butyrate)
summary(lm1.butyrate)
plot(concentration_mM~Area, data=butyrate)
intercept1<-coef(summary(lm1.butyrate))["(Intercept)","Estimate"]
slope1<-coef(summary(lm1.butyrate))["Area","Estimate"]		
	# this calls out the slope
summary(lm1.butyrate)$adj.r.squared
	# this calls out the R-squared
intercept1+slope1*butyrate$Area
	# this gives you the calculated concentration of the standards...
# looks good

# let's do this for every compound, together:
dat <- data.table(x=st.data$Area, y=st.data$concentration_mM, grp=st.data$compound)
compounds<-dat[,list(r2=summary(lm(y~x))$r.squared , slope=summary(lm(y~x))$coef[2], intercept=summary(lm(y~x))$coef[1] ), by=grp]
# this gives you the R squared for all compounds, by standard
#write.table(compounds, file="HPLC_RawData/scfa1/lcm1_standards.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# let's apply this to our data, and calculate the concentrations for the unknowns:
adj.data<-merge(lcm1, compounds, by.x="compound", by.y="grp")
#adj.data$Area[is.na(adj.data$Area)]<-0  #since including the intercept, this does not work
adj.data$calc_conc<-adj.data$intercept+adj.data$slope*adj.data$Area
	# now, we also want to adjust for weight
	# there are two ways: 
		# fecal amount/PBS_added x calc_conc: technically, this is what the Schmidt lab uses
			# however, I would also need to account for the dilution done at the end (create another variable)
		# calc_conc/dilution: since we did our 1/10 dilutions based on weight/vol, ours are technically adjusted
	# let's do both, then you can choose at the end
# first, let's get the important variables and merge with weight data:
results<-adj.data[, c("compound", "run", "Data_Filename", "Sample_ID", "Ret_Time", "Area", "Height", "Peak_Start", "Peak_End", "Relative_Retention_Time", "r2", "slope", "intercept", "calc_conc")]
sample.key<-key[key$sample_type==c("unknown"), ]
adj.results<-merge(sample.key, results, by.x="lcm_sampleID", by.y="Sample_ID")
adj.results$calc_conc[is.na(adj.results$calc_conc)]<-0
	# now, take into consideration the dilution
	# first, the Schmidt way:
adj.results$dilfactor<-10 * adj.results$dilution
adj.results$calc_conc_adj.weight<-abs((adj.results$calc_conc * adj.results$dilfactor * (adj.results$PBS_added/1000)) / adj.results$fecal_weight)
	# using the dilution instead:
adj.results$calc_conc_adj.df<-adj.results$calc_conc / adj.results$dilution
adj.results1<-adj.results
# done!!!
#write.table(adj.results1, file="HPLC_RawData/scfa1/scfa1_run1_concentrations.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

## run2:
lcm2<-read.table(file="HPLC_RawData/scfa1/lcm.run2_data_05.27.16.txt", header=TRUE)
key<-read.table(file="HPLC_RawData/scfa1/scfa1_samplekey.txt", header=TRUE)
standards<-key[key$sample_type==c("standard"), c("lcm_sampleID", "sample_name", "concentration_mM")]
standards<-standards[!is.na(standards$concentration_mM), ]

# merge with the data to get standards curve
st.data<-merge(standards, lcm2, by.x="lcm_sampleID", by.y="Sample_ID")

# calculate standard curves:
dat <- data.table(x=st.data$Area, y=st.data$concentration_mM, grp=st.data$compound)
compounds<-dat[,list(r2=summary(lm(y~x))$r.squared , slope=summary(lm(y~x))$coef[2], intercept=summary(lm(y~x))$coef[1] ), by=grp]
#write.table(compounds, file="HPLC_RawData/scfa1/lcm2_standards.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# merge with data:
adj.data<-merge(lcm2, compounds, by.x="compound", by.y="grp")
adj.data$calc_conc<-adj.data$intercept+adj.data$slope*adj.data$Area

results<-adj.data[, c("compound", "run", "Data_Filename", "Sample_ID", "Ret_Time", "Area", "Height", "Peak_Start", "Peak_End", "Relative_Retention_Time", "r2", "slope", "intercept", "calc_conc")]
sample.key<-key[key$sample_type==c("unknown"), ]
adj.results<-merge(sample.key, results, by.x="lcm_sampleID", by.y="Sample_ID")
adj.results$calc_conc[is.na(adj.results$calc_conc)]<-0
	# by weight:
adj.results$dilfactor<-10 * adj.results$dilution
adj.results$calc_conc_adj.weight<-abs((adj.results$calc_conc * adj.results$dilfactor * (adj.results$PBS_added/1000)) / adj.results$fecal_weight)
	# by dilution:
adj.results$calc_conc_adj.df<-adj.results$calc_conc / adj.results$dilution
adj.results2<-adj.results
#write.table(adj.results2, file="HPLC_RawData/scfa1/scfa_run2_concentrations.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# combined data together and merge with sample info:
all.results<-rbind(adj.results1, adj.results2)
#write.table(all.results, file="HPLC_RawData/scfa1/scfa1_all.concentrations.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#meta<-read.table(file="../allhumos2/allhumos_metadata_temp.txt", header=TRUE)
#combined<-merge(meta, all.results, by.x="sampleID", by.y="sample_name")

####-------------------#### HF5 plate:

### note: for this plate, raw data files were separated by compound already
# code differs a bit from above ()

### read in and combine data with the sample key info
# look at standard curve:
key<-read.table(file="HPLC_RawData/scfa2/scfa2_samplekey.txt", header=TRUE)
data<-read.csv(file="HPLC_RawData/scfa2/HF5plA_all.compounds_10.05.17.txt", header=TRUE, na.strings = c("NA", "<NA>", "missed", "", "*"), sep="\t")

# add column that matches SampleID on HPLC sheet:
key$SampleID[key$sampletype=="unknown"]<-paste(key$DNA_Plate[key$sampletype=="unknown"], key$Well[key$sampletype=="unknown"], sep="-")
key$SampleID[key$sampletype=="standard"]<-as.character(key$Well[key$sampletype=="standard"])
standards<-key[key$sampletype==c("standard"), c("SampleID", "sampleID", "concentration_mM")]
standards<-standards[!is.na(standards$concentration_mM), ]

# merge with the data to get standards curve
st.data<-merge(standards, data, by.x="SampleID", by.y="SampleID")

### test out one compound:
butyrate<-st.data[st.data$compound==c("butyrate"), ]
lm1.butyrate<-lm(concentration_mM~Area, data=butyrate)
summary(lm1.butyrate)
plot(concentration_mM~Area, data=butyrate)

# can calculate slope etc:
intercept1<-coef(summary(lm1.butyrate))["(Intercept)","Estimate"]
slope1<-coef(summary(lm1.butyrate))["Area","Estimate"]		
summary(lm1.butyrate)$adj.r.squared
intercept1+slope1*butyrate$Area
	# this gives you the calculated concentration of the standards...
	
### let's do it for all the compounds:
dat <- data.table(x=st.data$Area, y=st.data$concentration_mM, grp=st.data$compound)
compounds<-dat[,list(r2=summary(lm(y~x))$r.squared , slope=summary(lm(y~x))$coef[2], intercept=summary(lm(y~x))$coef[1] ), by=grp]
#write.table(compounds, file="HPLC_RawData/scfa2/HF5plA_standards.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# let's apply this to our data, and calculate the concentrations for the unknowns:
adj.data<-merge(data, compounds, by.x="compound", by.y="grp")
#adj.data$Area[is.na(adj.data$Area)]<-0  #since including the intercept, this does not work
adj.data$calc_conc<-adj.data$intercept+adj.data$slope*adj.data$Area
	# now, we also want to adjust for weight
	# there are two ways: 
		# fecal amount/PBS_added x calc_conc: technically, this is what the Schmidt lab uses
			# however, I would also need to account for the dilution done at the end (create another variable)
		# calc_conc/dilution: since we did our 1/10 dilutions based on weight/vol, ours are technically adjusted
	# let's do both, then you can choose at the end
# first, let's get the important variables and merge with weight data:
results<-adj.data[, c("compound", "DataFilename", "SampleID", "Ret.Time", "Area", "Height", "PeakStart", "PeakEnd", "RelativeRetentionTime", "r2", "slope", "intercept", "calc_conc")]
sample.key<-key[key$sampletype==c("unknown"), ]
adj.results<-merge(sample.key, results, by.x="SampleID", by.y="SampleID")
adj.results$calc_conc[is.na(adj.results$calc_conc)]<-0
	# now, take into consideration the dilution
	# if you want the Schmidt way (note: would need to add weights, PBS from previous sheets):
#adj.results$dilfactor<-10 * adj.results$dilution
#adj.results$calc_conc_adj.weight<-abs((adj.results$calc_conc * adj.results$dilfactor * (adj.results$PBS_added/1000)) / adj.results$fecal_weight)
	# using the dilution instead:
adj.results$calc_conc_adj.df<-adj.results$calc_conc / adj.results$dilution
adj.results1<-adj.results
# done!!!
#write.table(adj.results1, file="HPLC_RawData/scfa2/scfa2_all.concentrations.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

####-------------------#### HF7.GF2 plate:

### note: for this data set, there are MANY MORE samples
# code is maximized to read in multiples of data files
# have to also account for differences across plates, so normalization is also more complicated


### read in files and normalize data:
data<-read.csv(file="HPLC_RawData/scfa3/HF7.GF2_all.compounds_08.01.18.txt", header=TRUE, na.strings = c("NA", "<NA>", "missed", "", "*"), sep="\t")
keydf <- read.table(file="HPLC_RawData/scfa3/scfa3_samplekey.txt", header=TRUE)

### normalize (by plate!):
# let's also pick what values we need
df<-as.data.frame(data[, c("plateID", "compound", "SampleName", "SampleID", "SampleType", "Ret.Time", "Area", "PeakStart", "PeakEnd")])

### we also need to get a range of acceptable retention times: by both plate and compound, for each standard group (within a plate)
# let's parse our table of ALL compounds into a vector of tables, by compound
# then, when we apply the function to get a minimum/maximum per plate, and parse out info that way
# note: that WOULD be the easiest, but I got sick of trying to figure it out
# so instead, I just subset the standards, calculated the min/max per plate and compound, then added these as a column to the full dataframe
# these can then be used to define whether a value for an unknown should be accepted or not

# first, combine plateID and standards to categorize each standard value correctly in case we need to parse out some data
# ID which ones are standards, unknowns, controls
df$SampleType <- as.character(df$SampleType)
df$SampleType[grep("STD", df$SampleID)]<-"STD"
controls<-c("BLANK", "Water", "PBS")
df$SampleType[df$SampleID %in% controls]<-"control"		
df$STDID[df$SampleType=="STD"]<-paste(df$plateID[df$SampleType=="STD"], df$SampleID[df$SampleType=="STD"], sep="_")

# split into a list of dataframes
#dflist<-split(df, df$compound)

# technically, this works to get a minimum by a certain value
#test %>% 
#  group_by(plateID, compound) %>% 
#  summarise(min = min(Ret.Time, na.rm=T), max = max(Ret.Time, na.rm=T))
  
# get a standard table
stdf<-droplevels(df[df$SampleType=="STD", ])

# calculate standard min/max per plate and compound:
stable<-stdf %>% 
  group_by(plateID, compound) %>% 
  summarise(minRT = min(Ret.Time, na.rm=T), maxRT = max(Ret.Time, na.rm=T))
stable$PlCID<-paste(stable$plateID, stable$compound, sep="_")  

# add these to your regular df:
df$PlCID<-paste(df$plateID, df$compound, sep="_")
dfe<-merge(df, stable[, 3:5], by="PlCID")   

# now, reject any peaks for each compound that are outside of the range:
dfe$PeakQual[dfe$Ret.Time >= dfe$minRT & dfe$Ret.Time <= dfe$maxRT]<-"OK"
dfe$PeakQual[dfe$Ret.Time < dfe$minRT | dfe$Ret.Time > dfe$maxRT]<-"Reject"
dfe$PeakQual[is.na(dfe$Ret.Time)]<-NA
summary(as.factor(dfe$PeakQual))
# this is probably a bit too conservative, and rejecting actual peaks.
# based on manual review, the retention time can probably be expanded a bit
# doing this will NOT replace the "Undetected" values--those remain "0", as no peak was identified within that range


###--------------###
### fixing values:

# let's give a little bit more of a cushion for accepting / rejecting peaks, and assume that an NA peak = 0, now that a manual review has been done:
dfe$PeakQual[dfe$Ret.Time >= dfe$minRT - 0.04 & dfe$Ret.Time <= dfe$maxRT + 0.04]<-"OK"
dfe$PeakQual[dfe$Ret.Time < dfe$minRT - 0.04| dfe$Ret.Time > dfe$maxRT + 0.04]<-"Reject"
dfe$PeakQual[is.na(dfe$Ret.Time)]<-"Undetected"
summary(as.factor(dfe$PeakQual))
	# this produced only 73 'reject' values--seems appropriate

# double check that the STD standard names have the mix concentrations in each of them (in the past, I have forgotten to change these)
dim(dfe[grep("STD", dfe$SampleID),])
dim(dfe[grep("Mix", dfe$SampleName),])		#these should match

###----### IF you need to add 'x mM Mix' to the standards' SampleID, do the following (not necessary for this example):
# let's do this with a merge instead of defining each one individually
temp <- dfe[grep("Mix", dfe$SampleName), c("SampleName", "SampleID")] %>%
		distinct() %>%
		left_join(dfe, ., by = "SampleID")
temp$SampleName.x[temp$SampleType %in% c("STD")] <- temp$SampleName.y[temp$SampleType %in% c("STD")]
temp$SampleName.y <- NULL
colnames(temp)[colnames(temp)=="SampleName.x"] <- "SampleName"
dfe <- temp
###----###

# Finally, designate 'Undetected' peaks from the raw data files (NA's previously) as 0:
dfe$Area[dfe$PeakQual == "Undetected"] <- 0

### now, let's combine this data with the sample key info
# with this information, we can model the concentration for all unknowns

# merge data with key information
d<-merge(dfe, keydf, by.x="SampleID", by.y="WellID")

### Now, we can check standard curves!
### and then calculate concentration based on the area for all compounds

### test out standard curve one compound:
standard<-d[d$compound==c("butyrate") & d$SampleType=="STD", ]
lm1.standard<-lm(concentration_mM~Area, data=standard)
summary(lm1.standard)
plot(concentration_mM~Area, data=standard)
	## all of them look really good!

# can calculate slope etc:
but<-d[d$compound==c("butyrate") & d$SampleType=="STD", ]
lm1.but<-lm(concentration_mM~Area, data=but)
intercept1<-coef(summary(lm1.but))["(Intercept)","Estimate"]
slope1<-coef(summary(lm1.but))["Area","Estimate"]		
summary(lm1.but)$sadj.r.squared
intercept1+slope1*but$Area
	# this gives you the calculated concentration of the standards
	# BUT, we need to do it per compound, possibly per plate?
	
### let's calculate the actual concentrations for our data based on these values	
# identify standards
std<-d[d$SampleType=="STD", ]
st.table <- data.table(x=std$Area, y=std$concentration_mM, grp=std$compound)

# calculate slope, r2, and intercept for each compound:
# note: in this case, we are just combining the area data for ALL standards (we eliminated peaks in the data that we do not want already)
compounds<-st.table[,list(r2=summary(lm(y~x))$r.squared , slope=summary(lm(y~x))$coef[2], intercept=summary(lm(y~x))$coef[1] ), by=grp]
#write.table(compounds, file="HPLC_RawData/scfa3/HF7.GF2_standards.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# let's apply this to our data, and calculate the concentrations for the unknowns:
adj.data<-merge(d, compounds, by.x="compound", by.y="grp")
adj.data$calc_conc<-adj.data$intercept+adj.data$slope*adj.data$Area

# still need to adjust for weight/dilution:
# fecal amount/PBS_added x calc_conc: technically, this is what the Schmidt lab uses
	# however, I would also need to account for the dilution done at the end (create another variable)
	# Schmidt way would be:  adj.results$calc_conc_adj.weight<-abs((adj.results$calc_conc * adj.results$dilfactor * (adj.results$PBS_added/1000)) / adj.results$fecal_weight)
	# if you want the Schmidt way (note: would need to add weights, PBS from previous sheets):
# calc_conc/dilution: since we did our 1/10 dilutions based on weight/vol, ours are technically adjusted
# let's do this instead
results<-adj.data[adj.data$SampleType=="Unknown",]
results$calc_conc[is.na(results$calc_conc)]<-0		#anything NA means it was not there
results$calc_conc[results$PeakQual=="Reject"]<-NA		#anything rejected means undetermined (NA)
	# now, take into consideration the dilution	
#adj.results$dilfactor<-10 * adj.results$dilution
results$dilution<-as.numeric(results$dilution)
results$calc_conc_adj<-results$calc_conc / results$dilution

# remove some extraction/control samples from final product:
results$SampleType[grep("extraction", results$sampleID)]<-"control"

adj.results<-results[results$SampleType=="Unknown",]

#write.table(adj.results, file="HPLC_RawData/scfa3/scfa3_all.concentrations.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

``


####-------------------#### 

### SCFA4 plates
	
``{r}``

### note: this is combining 2 different runs of 4 plates each...
# code has been optimized to do this for both of the runs at the same time!


### Before we normalize, we need to make sure that peaks outside of the RT range are rejected
	# - a 'rejected' peak is a peak that cannot be resolved, generally because another compound is obscuring its detection, making it impossible to detect its presence or measure its presence
	# - the machine software automatically assigns a nearby peak to this
	# - during manual curation, I try to resolve the 'real' peak if I can
		# if I cannot, I either leave it there with incorrect RT values to be resolved later (here), or I assign it NA if it really looks like it is not there
	# - during manual curation, I make sure 'NA' values are actually 'NA', or re-assign it as 'NA' in the above scenario
		# - an 'NA' in the raw data for me indicates that it is 'undetected', which gets a 0 peak area assignment (it is not present)

### we need to get a range of acceptable retention times: by both plate and compound, for each standard group (within a plate)
# let's parse our table of ALL compounds into a vector of tables, by compound
# then, when we apply the function to get a minimum/maximum per plate, and parse out info that way
# note: that WOULD be the easiest, but I got sick of trying to figure it out
# so instead, I just subset the standards, calculated the min/max per plate and compound, then added these as a column to the full dataframe
# these can then be used to define whether a value for an unknown should be accepted or not for all of your data files at the same time

df <- read.table("HPLC_RawData/scfa4.5/HF9.HF10_all.compounds.raw.txt", header=T, sep="\t")

### first, combine plateID and standards to categorize each standard value correctly in case we need to parse out some data
# ID which ones are standards, unknowns, controls
df$SampleType[grep("STD", df$SampleID)]<-"STD"
controls<-c("BLANK", "Water", "PBS")
df$SampleType[df$SampleID %in% controls]<-"control"

### then, we need to add each plate up into their respective experiment/plate/runs from the HPLC:
df$Exp <- as.factor(gsub("_.*", "", df$DataFileOrigin) %>%
			gsub(".*-", "", .)) 
summary(df$Exp)
# you should have as many categories here as you have different plates/experiments

### all the standards need to be compared within their own compound, as well as within their own run/plate
### FIRST, let's make sure our peaks are within the appropriate RT values
stdf <- filter(df, SampleType == "STD")

# redefine some variables:
stdf$compound <- as.factor(stdf$compound)

# calculate standard min/max per plate and compound:
stable<-stdf %>% 
  group_by(Exp, compound) %>% 
  summarise(minRT = min(Ret.Time, na.rm=T), maxRT = max(Ret.Time, na.rm=T))
stable$PlCID<-paste(stable$Exp, stable$compound, sep="_")  
as.data.frame(stable)
# you should have a minRT/maxRT for each compound, per plate/run

# now, apply these values to the FULL dataframe:
df$PlCID<-paste(df$Exp, df$compound, sep="_")
dfe<-merge(df, stable[, 3:5], by="PlCID")

### now, reject any peaks for each compound that are outside of the range:
dfe$PeakQual[dfe$Ret.Time >= dfe$minRT & dfe$Ret.Time <= dfe$maxRT]<-"OK"
dfe$PeakQual[dfe$Ret.Time < dfe$minRT | dfe$Ret.Time > dfe$maxRT]<-"Reject"
dfe$PeakQual[is.na(dfe$Ret.Time)]<-NA
summary(dfe$PeakQual)
# this sets which peaks you are setting as acceptable ('OK'), vs rejecting ('Reject')

###--------------###
### When I manually went through each compound within each plate--there are some anomalies. I work with the following conventions:
	# NAs in my dataset are assumed to be a concentration of 0 in this case (will be classified as 'Undetected' in the next step)
		# this could be bc it is not there but also that it is swamped by another compound, and we cannot determine it's area (rare, but does happen)
	# in looking at the data, some peaks in the data set look to be acceptable despite not COMPLETELY filling in these min/max RT values
		# we will increase the filtering for reject/OK to allow for a little bit more room (+/-0.05 from the min and max)
		# this should still reject the peaks that are REALLY off
	# If I am rejecting a peak, this is bc there is a peak of another compound swamping the compound I am looking at, or just another compound
		# these peaks cannot be determined to be there or not, so I just classify them as 'Reject', and basically is treated as missing data (a true 'NA')

###--------------###
### fixing values:

# let's give a little bit more of a cushion for accepting / rejecting peaks, and assume that an NA peak = 0, now that a manual review has been done:
dfe$PeakQual[dfe$Ret.Time >= dfe$minRT - 0.08 & dfe$Ret.Time <= dfe$maxRT + 0.08]<-"OK"
dfe$PeakQual[dfe$Ret.Time < dfe$minRT - 0.08| dfe$Ret.Time > dfe$maxRT + 0.08]<-"Reject"
dfe$PeakQual[is.na(dfe$Ret.Time)]<-"Undetected"
summary(as.factor(dfe$PeakQual))
	# this produced only 51 'reject' values--pretty good

# double check that the STD standard names have the mix concentrations in each of them (in the past, I have forgotten to change these)
dim(dfe[grep("STD", dfe$SampleID),])
dim(dfe[grep("Mix", dfe$SampleName),])		#these should match

###----### IF you need to add 'x mM Mix' to the standards' SampleID, do the following (not necessary for this example):
# let's do this with a merge instead of defining each one individually
temp <- dfe[grep("Mix", dfe$SampleName), c("SampleName", "SampleID")] %>%
		distinct() %>%
		left_join(dfe, ., by = "SampleID")
temp$SampleName.x[temp$SampleType %in% c("STD")] <- temp$SampleName.y[temp$SampleType %in% c("STD")]
temp$SampleName.y <- NULL
colnames(temp)[colnames(temp)=="SampleName.x"] <- "SampleName"
dfe <- temp
###----###

# Finally, designate 'Undetected' peaks from the raw data files (NA's previously) as 0:
dfe$Area[dfe$PeakQual == "Undetected"] <- "0"

# overwrite file, if you want:
#write.table(dfe, file="HPLC_RawData/scfa4.5/HF9.HF10_all.compounds.raw.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)	


###--------------###
### combine with sample info:
dfe <- read.table("HPLC_RawData/scfa4.5/HF9.HF10_all.compounds.raw.txt", header=T, sep="\t")
key <- read.table("HPLC_RawData/scfa4.5/scfa4.5_samplekey.txt", header=T, sep="\t")

# this file has all the info for all samples, including the standards
# BUT, key to merging these correctly with the data is that you are merging by wellID, PER EXPERIMENT
	# you need to make sure the experiment, plate, and well are matched up correctly

# you will need to add an Exp column here that matches what you have defined in the rawdata
# if you already ensured that everything is matching, you are good to go!
# if you used 'p1' and 'plateA', for example, interchangeably, no worries, you can change it here
key$Exp <- gsub("_samplekey.txt", "", key$DataFileOrigin) %>%
			gsub("HPLC_", "", .) %>%
			gsub("_", ".", .) %>%
			gsub("plate", "Plate", .) #%>%
			#gsub("plateA", "p1", .) %>%	# in case you need to do this
			#gsub("plateB", "p2", .) %>%
			#gsub("plateC", "p3", .) %>%
			#gsub("plateD", "p4", .)
			
# Create a wellkey to combine each well-exp-plate appropriately
key$wellkey <- paste(key$Exp, key$Well, sep="-") %>%
				gsub("[.]", "-", .) %>%
				gsub("p", "P", .)
				
# Create wellkey for merge for standards:
key$WellID <- as.character(key$WellID)
key$sampletype <- as.character(key$sampletype)
key$wellkey[key$sampletype =="standard"] <- key$WellID[key$sampletype =="standard"]

# Create a wellkey in the data file to match this (with experiment, plate, and well, but ONLY in samples, not STD)--the wellID in this file is SampleID
# note: for this one, there are some duplicate samples that were ran in PlateB, named as PlateB2: for this, we will rename it PlateB again so that it merges properly
dfe$SampleType <- as.character(dfe$SampleType)
dfe$Exp <- as.character(dfe$Exp)
dfe$wellkey[is.na(dfe$SampleType)] <- dfe$Exp[is.na(dfe$SampleType)] %>%
			gsub("PlateB2", "PlateB", .) %>%
			gsub("\\..*", "", .) %>%
			paste(., dfe$SampleID[is.na(dfe$SampleType)], sep = "-") %>%
			gsub("P1", "PlateA", .) %>%
			gsub("P2", "PlateB", .) %>%
			gsub("P3", "PlateC", .) %>%
			gsub("P4", "PlateD", .) 
# note: sometimes this does not work (even in tidyr) because wellkey gets replaced as an integer? unclear why
# also add control / standard wellkey to data file:
dfe$SampleID <- as.character(dfe$SampleID)
dfe$wellkey[is.na(dfe$wellkey)] <- dfe$SampleID[is.na(dfe$wellkey)]

d<-merge(dfe, key[, c("wellkey", "sampleID", "DNA_Plate", "concentration_mM", "dilution")], by.x="wellkey", by.y="wellkey", all.x=T)
d$SampleType[is.na(d$SampleType)] <- "unknown"
# all EXCEPT standards should join nicely

# looks like everything joined, except the standards and control! let's fix their sampleID, and add the ACTUAL concentration according to our sample key data
d$SampleName <- as.character(d$SampleName)
d$SampleID <- as.character(d$SampleID) 
d$sampleID[is.na(d$sampleID)] <- d$SampleName[is.na(d$sampleID)]
d$concentration_mM[d$SampleType %in% c("STD")] <- d$SampleName[d$SampleType %in% c("STD")] %>%
	gsub("mMMix", "", .)

# write out results if you want, now with metadata
#write.table(d, file="HPLC_RawData/scfa4.5/HF9.HF10_all.compounds.raw_w.samplekey.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	## check that ALL samples have a 0.1 in their dilution column--this means that they joined correctly!

###--------------###
### Calculate concentration and standard curve!
d <- read.table("HPLC_RawData/scfa4.5/HF9.HF10_all.compounds.raw_w.samplekey.txt", header=T, sep="\t")

# do for one compound
standard <- d %>%
			filter(compound == "butyrate", SampleType == "STD")

lm1.standard<-lm(concentration_mM~Area, data=standard)
summary(lm1.standard)
plot(concentration_mM~Area, data=standard)
	## you can do this for all the compounds, if you want
	## this is ALL runs combined--all of them look really good!
	## of course, will have to apply normalization by PLATE (next)

#----
# now, let's calculate per run:

standard <- d %>%
			filter(compound == "butyrate", SampleType == "STD", Exp == "HF10.PlateB")

lm1.standard<-lm(concentration_mM~Area, data=standard)
summary(lm1.standard)
plot(concentration_mM~Area, data=standard)
	## LOOKS EVEN BETTER! Is there anything more lovely than a standard curve with an R-squared value of 0.9999? I think not
	
# these are the pieces f information we need to extract per compound, per run:
intercept1<-coef(summary(lm1.standard))["(Intercept)","Estimate"]
slope1<-coef(summary(lm1.standard))["Area","Estimate"]		
# ... to use this to calculate the concentration:
intercept1+slope1*Area  # the Area will be taken care of later, but we need to figure out how to add columns of intercept1 and slope1 to the file, per run, per compound


### let's do this for each compound and run to make a table of standard values to asses to our table
# note: once again, there is likely to be a more efficient, fancier way of doing this, but here we are

### LOOP IT UP for the standards

# let's rename a variable to account for the compound x run
d$rc <- paste(d$Exp, d$compound, sep="_") 

# split your df into a list of dfs, by your variable, separating ONLY the standards:
standard <- d %>%
		filter(SampleType == "STD") %>%
		split(., f = .$rc)

# run the loop over the data		
std.table <- data.frame(NULL) 
	for(n in 1:length(standard)) {
  		lm.std<-lm(concentration_mM~Area, data=standard[[n]])
  		print(summary(lm.std))
  		standard[[n]]$intercept <- coef(summary(lm.std))["(Intercept)","Estimate"]   
  		standard[[n]]$slope <- coef(summary(lm.std))["Area","Estimate"]
  		standard[[n]]$Rsq <- summary(lm.std)$adj.r.squared
  		}
head(std.table)  # note: this did not print a table, just added to the 'standard' table--whatever, it works!

# now, just flatten the standard table and extract what you need
std.table <- standard %>%
	reduce(rbind) %>%
	select(rc, intercept, slope, Rsq) %>%
	distinct() %>%
	as.data.frame()
#write.table(std.table, file="HPLC_RawData/scfa4.5/HF9.HF10_standards.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# finally, combine THIS data to your data table (d) above:
df <- merge(d, std.table, by="rc")


### now, you can calculate the concentration_mM for each compound

df$calc_conc <- df$intercept + df$slope*df$Area
# note: I am not sure if I should just keep the y-intercept as 0--technically, 0 concentration should = 0 in the Area! 
# will let you make this call--for now, I will keep the calculated y-intercept

# if keeping the calculated y-intercept, need to re-assign 0 for the undetected peaks
# also note that some of the y-intercepts are <0; will reassign these to 0 as well
df$calc_conc <- as.numeric(df$calc_conc)
df$calc_conc[df$PeakQual == "Undetected"] <- 0
df$calc_conc[df$calc_conc < 0] <- 0

# let's also add NA for the peaks that were rejected bc they couldn't be resolved:
df$calc_conc[df$PeakQual == "Reject"] <- NA


#----
### The FINAL normalization step is to take into account the weight/dilution! 
	# - if you prepped the samples with a 1 in 10 dilution, all samples should already be standardized to weight
	# - now, all you need to do is re-assign the concentration to get mmol/kg stool
	# - note: for some samples, I had to dilute further since there was not enough--I always update this in my sample meta file
	
# let's finally reassign the SampleType for our samples (unknowns):
df$SampleType <- as.character(df$SampleType)
df$SampleType[is.na(df$SampleType)] <- "unknown"
df$SampleType <- as.factor(df$SampleType)

# now, for the unknowns, adjust the concentration_mM to reflect the dilution
# you also need to 
df$calc_conc_adj[df$SampleType =="unknown"] <- df$calc_conc[df$SampleType =="unknown"] / df$dilution[df$SampleType =="unknown"]


# select only your samples:
df <- df %>%
	filter(SampleType == "unknown")
	
# let's add "NA" to the calculated area
	
#write.table(df, file="HPLC_RawData/scfa4.5/scfa4.5_all.concentrations.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)	

```

####-------------------#### 

### Step 2: combine all data with metadata of current project
	
``{r}

### need to look at all the past SCFA data
library(reshape2)
library(tidyverse)
library(tidyr)
library(dplyr)

# let's merge all the current data in a readable way
# order of columns: sampleID, seqID, compound, DataFilename, Ret_Time, Area, Peak_Start, Peak_End, RelativeRetentionTime, r2, slope, intercept, calc_conc, calc_conc_adj.df
scfa1 <- read.table(file="HPLC_RawData/scfa1/scfa1_all.concentrations.txt", header=TRUE) %>%
	select(seqID, compound, DataFileName = Data_Filename, Ret_Time, Area, Peak_Start, Peak_End, r2, slope, intercept, calc_conc, calc_conc_adj.df) %>%
	mutate(PeakQual = "OK")
scfa2 <- read.table("HPLC_RawData/scfa2/scfa2_all.concentrations.txt", header=T) %>%
	select(seqID, compound, DataFileName = DataFilename, Ret_Time = Ret.Time, Area, Peak_Start = PeakStart, Peak_End = PeakEnd, r2, slope, intercept, calc_conc, calc_conc_adj.df) %>%
	mutate(PeakQual = "OK")
scfa3 <- read.table("HPLC_RawData/scfa3/scfa3_all.concentrations.txt", header=T) %>%
	select(seqID = sampleID, compound, DataFileName = PlCID, Ret_Time = Ret.Time, Area, Peak_Start = PeakStart, Peak_End = PeakEnd, r2, slope, intercept, calc_conc, calc_conc_adj.df = calc_conc_adj, PeakQual)
scfa4 <- read.table("HPLC_RawData/scfa4.5/scfa4.5_all.concentrations.txt", header=T) %>%
	select(seqID = sampleID, compound, DataFileName = DataFileOrigin, Ret_Time = Ret.Time, Area, Peak_Start = PeakStart, Peak_End = PeakEnd, r2 = Rsq, slope, intercept, calc_conc, calc_conc_adj.df = calc_conc_adj, PeakQual)

# bind together and reshape into wide data structure
scfa <- rbind(scfa1, scfa2, scfa3, scfa4) %>%
	#write_tsv("all.scfa.concentrations_combined.txt")
	mutate(calc_conc_adj.df = abs(calc_conc_adj.df)) %>%
	filter(PeakQual == "OK") %>%
	drop_na(seqID) %>%
	filter(compound %in% c("acetate", "propionate", "butyrate", "succinate")) %>%
	select(seqID, compound, calc_conc_adj.df) %>%
	pivot_wider(names_from = compound, values_from = calc_conc_adj.df) #%>%
	# write_tsv("all.scfa.concentrations_combined_filtered.wide.txt")

# note: some seqIDs were missing; found them using this code:
#test <- scfa %>% dplyr::summarise(n = dplyr::n(), .by = c(seqID, compound)) %>% dplyr::filter(n > 1L)
# in looking at the samplekey maps from scfa4.5, these are Kwi/Steve samples; just gt rid of the seqID with NA

# other reshape option to get each variable in their own column
#scfa_filtered <- dcast(scfa, sampleID ~ compound) %>%
#				.[, c("sampleID", "acetate", "butyrate", "propionate", "lactate", "succinate")]

## combine with meta, ensuring names match:	
meta <- read.table("../Data/16S/allhf_runs/allhf_seq.meta.merge.txt", header=T)

# test merge to see how many match:
scfa.meta <- meta %>%
	right_join(., scfa, by = "seqID") %>%
	#write_tsv("TEST_scfa.meta.txt")
	# when compared to inner_join, missing 5 seqIDs
	
scfa.meta %>% filter(is.na(seqID) | is.na(sampleID)) # lists names that don't match
seqID
#         GF2_hFMT1B
# HF8_myFMT1_regular
# HF8_myFMT2_regular
#          16_1903_2:  found it; 
#          d9_1835_q
meta %>% filter(seqID == "d16_1903_2")			# exists; need to change this name in scfa key (HF10, plate C)
meta %>% filter(str_detect(seqID, "d9_1835"))	# exists; options include, _1, _2, or _3
scfa %>% filter(str_detect(seqID, "d9_1835"))	# _2 and _3 are covered, so must be _1!!
scfa %>% filter(str_detect(seqID, "GF2_hFMT1"))	# exists; looks like there is an extra 'd' on the seqID side of META...will change here...
# for other 2 samples, these seem to be FMT samples that were not extracted; will not worry about them for now 

# correct sample names:
fixed <- scfa %>%
	mutate(seqID = str_replace_all(seqID, pattern = "16_1903_2", replacement = "d16_1903_2")) %>%
	mutate(seqID = str_replace_all(seqID, pattern = "d9_1835_q", replacement = "d9_1835_1"))
# re-read in file, and merge again:
meta <- read.table("../Data/16S/allhf_runs/allhf_seq.meta.merge.txt", header=T) %>%
	inner_join(., fixed, by ="seqID")
	
write_tsv(meta, "/Data/metabolomics/HPLC/allhf_scfa.seq.meta.merge.txt")

### all data is merged!!


### some last looks at the data
# in looking at how many cecal samples there are, there are some anomalies...
all.scfa %>% filter(type %in% c("cecal")) %>% select(seqID, Experiment, FMT_type, FMT_input, day, acetate, butyrate, type) %>% as.data.frame()
all.scfa %>% filter(Experiment %in% c("HF3", "HF4")) %>% select(seqID, Experiment, FMT_type, FMT_input, day, acetate, butyrate, type) %>% as.data.frame()
all.scfa %>% filter(Experiment %in% c("HF8", "HF9", "HF10")) %>% select(seqID, Experiment, FMT_type, FMT_input, day, acetate, butyrate, type) %>% as.data.frame()



```

####--------------------------------

### Step 3: for UM data, only need to SCFA and BA files with metadata
	- /Data/metabolomics/UM-MetaboCore/humofmt_metabcore_scfa.data.txt
		
```{r}

library(dplyr)
library(tidyverse)
library(tidyr)
library(data.table)
#library(purrr)
library(readr)

## wd: /Data/metabolomics

## read in core data and metadata:
corescfa <- read_tsv("UM-MetaboCore/humofmt_metabcore_scfa.data.txt") #%>%
			#write_tsv("humofmt_metabcore_scfa.data.txt")
coreba <- read_tsv("UM-MetaboCore/BAs_weightnorm.modified.txt") #%>%
			#write_tsv("BAs_weightnorm.modified.txt")
meta <- read.table("../16S/allhf_runs/allhf_seq.meta.merge.txt", header=T)

## merge core data together, then meta:
merged <- corescfa %>%
			select(SampleID, Acetate_nmol, Propionate_nmol, Butyrate_nmol, Isovalerate_nmol, Valerate_nmol, Hexanoate_nmol, Heptanoate_nmol, Octanoate_nmol) %>%
			inner_join(coreba, ., by="SampleID") %>%
			select(-CoreID, -CoreName, -wet_weight) %>%
			right_join(meta, ., by = c("seqID" = "SampleID")) %>%
			write_tsv("TEST_merged.core.txt")

# only 19 match...let's look at the full metadata sheet and see what names are not correctly aligned
allmeta <- read.table("../16S/allhf_runs/allhf_seq.meta.merge.txt", header=T)

merged %>% filter(is.na(seqID) | is.na(sampleID)) # 31 names don't match

# let's see what they might be:
allmeta %>% filter(str_detect(seqID, "1510")) 
allmeta %>% filter(str_detect(seqID, "1511"))
allmeta %>% filter(str_detect(seqID, "1512"))
allmeta %>% filter(str_detect(seqID, "975"))
allmeta %>% filter(str_detect(seqID, "907"))
allmeta %>% filter(str_detect(seqID, "908"))
	# these mice were not sequenced, at least not in this batch--all pre-FMT
allmeta %>% filter(str_detect(seqID, "1211")) %>% select(seqID, sampleID, Experiment, day, FMT, type, Cdiff_infection)
allmeta %>% filter(str_detect(seqID, "1212")) %>% select(seqID, sampleID, Experiment, day, FMT, type, Cdiff_infection)
allmeta %>% filter(str_detect(seqID, "1213")) %>% select(seqID, sampleID, Experiment, day, FMT, type, Cdiff_infection)
allmeta %>% filter(str_detect(seqID, "1206")) %>% select(seqID, sampleID, Experiment, day, FMT, type, Cdiff_infection)
allmeta %>% filter(str_detect(seqID, "1208")) %>% select(seqID, sampleID, Experiment, day, FMT, type, Cdiff_infection)
	# these mice (1211_1, _2, _3; 1210_1, _2; 1212_1, _2, _3; 1206_1, _2; 1208_1, _2) will be under seqID of d19_1211_* rather than cecal_d19_121*_*
	# they are also noCDI mice in this cage
allmeta %>% filter(str_detect(seqID, "1191")) %>% select(seqID, sampleID, Experiment, day, FMT, type, Cdiff_infection)
	# d42_1191_1 instead of cecal_d42_1191_1
allmeta %>% filter(str_detect(seqID, "1270")) %>% select(seqID, sampleID, Experiment, day, FMT, type, Cdiff_infection)
	# cecal_1270_4 instead of cecal_d42_1270_4
	# however, this is also a yWT_spore, and there is only one of them...
	# so, will exclude this one from the mFMT samples
allmeta %>% filter(str_detect(seqID, "1275")) %>% select(seqID, sampleID, Experiment, day, FMT, type, Cdiff_infection)
	# cecal_1275_1 instead of cecal_d42_1275_1
	# cecal_1275_2 instead of cecal_d42_1275_2
allmeta %>% filter(str_detect(seqID, "1294")) %>% select(seqID, sampleID, Experiment, day, FMT, type, Cdiff_infection)
	# cecal_1294_1 instead of cecal_d42_1945_1
	
## we want to keep all samples except the noCDI (any sample with a d19) day for now
## but, we do need to add some metadata to cecal sampels that were not sequenced (and thus not in our normal meta data)
## these samples do, however, have the metadata information in their seqID
#cages <- c("1191", "1270", "1275", "1294") # couldn't get this to work with case_when
#healthy <- c("cecal_dn7_1510_1", "cecal_dn7_1510_2", "cecal_dn7_1510_3")
corescfa <- read_tsv("UM-MetaboCore/humofmt_metabcore_scfa.data.txt") #%>%
			#write_tsv("humofmt_metabcore_scfa.data.txt")
coreba <- read_tsv("UM-MetaboCore/BAs_weightnorm.modified.txt") #%>%
			#write_tsv("BAs_weightnorm.modified.txt")
meta <- read.table("../16S/allhf_runs/allhf_seq.meta.merge.txt", header=T)

## add meta data to NA columns based on sample name:
meta_correct <- corescfa %>%
			select(SampleID, Acetate_nmol, Propionate_nmol, Butyrate_nmol, Isovalerate_nmol, Valerate_nmol, Hexanoate_nmol, Heptanoate_nmol, Octanoate_nmol) %>%
			inner_join(coreba, ., by="SampleID") %>%
			select(-CoreID, -CoreName, -wet_weight) %>%
			mutate(SampleID = str_replace_all(SampleID, pattern = "cecal_d42_1191", replacement = "cecal_1191")) %>%
			mutate(SampleID = str_replace_all(SampleID, pattern = "cecal_d42_127", replacement = "cecal_127")) %>%
			mutate(SampleID = str_replace_all(SampleID, pattern = "cecal_d42_1294", replacement = "cecal_1294")) %>%
			right_join(meta, ., by = c("seqID" = "SampleID")) %>%
			filter(!seqID == "cecal_1270_4") %>%	#this one was treated with spore prep
			filter(!str_detect(seqID, "d19")) %>%	#get rid of d19 samples, as these were noCDI
				# let's split the NAs away and add some metadata based on the seqID:
			filter(is.na(sampleID)) %>%
			separate(seqID, into = c("type", "day", "CageID", "mouse"), sep = "_", remove = FALSE) %>%
			mutate(day = str_remove_all(day, paste(c("d", "n"), collapse = "|")), 
					day = str_replace(day, "7", "-7"),
					day = as.integer(day))

# then, join these back together with the samples that DO match metadata
merged <- corescfa %>%
			select(SampleID, Acetate_nmol, Propionate_nmol, Butyrate_nmol, Isovalerate_nmol, Valerate_nmol, Hexanoate_nmol, Heptanoate_nmol, Octanoate_nmol) %>%
			inner_join(coreba, ., by="SampleID") %>%
			select(-CoreID, -CoreName, -wet_weight) %>%
			filter(str_detect(SampleID, 'd21|d42')) %>% #only keep post-FMT samples
			mutate(SampleID = str_replace_all(SampleID, pattern = "cecal_d42_1191", replacement = "cecal_1191")) %>%
			mutate(SampleID = str_replace_all(SampleID, pattern = "cecal_d42_127", replacement = "cecal_127")) %>%
			mutate(SampleID = str_replace_all(SampleID, pattern = "cecal_d42_1294", replacement = "cecal_1294")) %>%
			right_join(meta, ., by = c("seqID" = "SampleID")) %>%
			filter(!seqID == "cecal_1270_4") %>%
			bind_rows(., meta_correct) %>%
			mutate(Mouse_status = case_when(is.na(FMT_input) ~ "SPF-Young", TRUE ~ Mouse_status)) %>%
			mutate(Mouse_genotype = case_when(is.na(FMT_input) ~ "BL6", TRUE ~ Mouse_genotype)) %>%
			mutate(group = case_when(
				day == -7 ~ "healthy",
				day %in% c(0, 1, 4, 9) ~ "preFMT",
				FMT_input == "yWT" ~ "mFMT", 
				FMT_input == "none" ~ "noFMT",
				FMT_input == "D2_A" ~ "hFMT",
				TRUE ~ NA_character_
				))
write_tsv(merged, "UM-MetaboCore/HFm_scfa.ba.cecal.seq.meta.merge.txt")	
			
```
