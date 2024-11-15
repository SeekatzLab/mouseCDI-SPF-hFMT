### HUMOFMT/allhumos6B: Additional File prep
### 9.30.24
### Anna M. Seekatz

### Inputs: 
	- /Data/16S/allhf_runs (contains meta data for all samples processed from mothur outside of manuscript project)
		- allhf_seq.meta.merge.txt (meta for all samples)
		- seqID_key.txt (key for interpreting sampleID and seqID)
		- meta_controls.txt (meta data for controls per sequencing run)
	- Data/16S/allhf_runs/mothurfiles (raw data files processed from mothur; see mothur logfile for steps used)
		- allhumos6B.final.groups.summary (mothur-generated summary parameters for all samples)
		- allhumos6B.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary (mothur-generated 'phylotype' or reference-based taxonomy)
		- allhumos6B.final.shared: mothur-generated OTU abundances
		- allhumos6B.final.0.03.cons.taxonomy: taxonomic classification for mothur-generated OTUs

### Outputs: into /Data
	- /Data/16S/allhf_runs: combined data for all samples in runs
		- allhf_seq.meta.merge.txt: contains summary parameters and meta
		- allhf_all.genera.txt: genus-level relative abundance of all detected genera for all samples
		- allhf_genfrac2p.txt: genus-level relative abundance of top 98% of genera for all samples
		- allhf_genfrac2p.meta.txt: same as above, only includes meta data (transformed data)
	- /Data/16S: input files used for ALL CURRENT ANALYSES
		- HFm_summary.txt: summary parameters from mothur
		- HFm_summary.nmds.new.groups.txt: same as above, but includes additional meta data groups for analysis (including FMT inputs)
		- HFm.otus.shared: OTU abundances
		- taxonomy_meta.txt: additional color scheme information for OTU taxonomy


###-------------###

### Step 1:  combine mothur-generated sequencing data with meta data, weights, and cfu for multiple runs (includes samples not part of current study)
	- mothur was originally ran on a more comprehensive data set that included samples beyond those used in this study --> allhf pre-name
	- later file prep will subset samples for the current project --> HFm prefix

```{r}
library(tidyr)
library(reshape2)
library(tidyverse)
library(readr)
library(ggplot2)

### Step 1.1: combine meta data with summary metrics calculated from OTUs in mothur (see mothur file)

# cd /Data/16S/allhf_runs

meta <-read.table("allhf_runs/allhf_meta.wts.cfu.txt", header=T, na.strings = c("NA", "<NA>", "missed"))	# no duplications here!!
sums <- read.table("allhf_runs/mothurfiles/allhumos6B.final.groups.summary", header=T)

# let's merge the data with the seqID key made before:
key <- read.table("allhf_runs/seqID_key.txt", header=T)			# no duplications here!! 
metakey <- merge(meta, key, by.x="sampleID", by.y="sampleID", all.x=T, all.y=T)
which(duplicated(metakey$sampleID))						# check for dupes
#write_tsv(metakey, "testkey.txt")
	# check this file to make sure there are no missing seqID's
	# if there are, add to seqID_key.txt and repeat above
#which(is.na(metakey$seqID))		# FIX; add to seqID_key

merged <- merge(metakey, sums[, 2:ncol(sums)], by.x="seqID", by.y="group", all.x=T, all.y=T)
dim(merged)		# lots of mismatched data points...sigh
#write_tsv(merged,"seq.meta.merge.txt")	
which(duplicated(merged$sampleID))			# only ones duplicated are the "controls"...let's add those to the key!

# looks good! (see last notes below)

# last task: let's add some of the random experimental data to the human, mouse_survey, and control_ samples:
control_meta <- read.table("allhf_runs/meta_controls.txt", header=T)
dim(control_meta)	#217
head(control_meta)

dim(merged)
rownames(merged) <- merged$seqID
merged.matrix <- merged[, 31:ncol(merged)]
mouse.meta <- merged[merged$sampleType %in% c("mouse"), 1:30]

meta.combined <- rbind(control_meta, mouse.meta)
alldata <- merge(meta.combined, merged.matrix, by.x="seqID", by.y="row.names", all.x=T, all.y=T)
dim(alldata)		#8816--let's get rid of the straggler sample we do not want:
alldata <- merge(meta.combined, merged.matrix, by.x="seqID", by.y="row.names")
dim(alldata)		#8815

#write_tsv(alldata,"allhf_runs/allhf_seq.meta.merge.txt")	
	
# last step: adding "fecal" to the sequencing samples
# note: we are assuming that these are all fecals, since all cecals are generally plated, EXCEPT the following:
	# "healthy_cecal_1912_1", "healthy_cecal_1912_2", "healthy_cecal_1913_1", "healthy_cecal_1913_2"
# all others with sequencing data will be named fecal (if this is incorrect, we can fix later)

alldata$type[alldata$sampleID %in% c("healthy_cecal_1912_1", "healthy_cecal_1912_2", "healthy_cecal_1913_1", "healthy_cecal_1913_2")] <- "cecal"
alldata$type[grep("cecum_", alldata$sampleID)] <- "cecal"
alldata$type[grep("colon_", alldata$sampleID)] <- "colonic"
alldata$type[!is.na(alldata$nseqs) & is.na(alldata$type)] <- "fecal"

#write_tsv(alldata,"allhf_runs/allhf_seq.meta.merge.txt")

# 7/18/23: if you wanted to add more FMTinput representatives:
alldata[alldata$sampleType %in% c('FMTinput'), c("seqID")]
alldata[alldata$sampleType %in% c('mouse_survey'), c("seqID")]
df.inputs3 <- alldata[alldata$sampleType %in% c('mouse_survey'), ]	#need to incude some of these as baseline FMTinputs
## the regular sporeprep that I have used in the past does not have enough seqs to pass QC...let's add the "unfiltered_sporeprep1" "unfiltered_sporeprep2"
## let's also add a few more from the mouse_survey:
alldata$sampleType[alldata$seqID %in% c("unfiltered_sporeprep1","unfiltered_sporeprep2", "R_WT_1", "Rwt_2", "RAG205", "RAG207")] <- "FMTinput"
alldata$FMT_type[alldata$seqID %in% c("unfiltered_sporeprep1","unfiltered_sporeprep2", "R_WT_1", "Rwt_2", "RAG205", "RAG207")] <- "mFMT"
alldata$FMT[alldata$seqID %in% c("unfiltered_sporeprep1","unfiltered_sporeprep2")] <- "mFMT2"
alldata$FMT[alldata$seqID %in% c("R_WT_1", "Rwt_2")] <- "mFMT3"
alldata$FMT[alldata$seqID %in% c("RAG205", "RAG207")] <- "mFMT4"
alldata$FMT_input[alldata$seqID %in% c("unfiltered_sporeprep1","unfiltered_sporeprep2")] <- "yWT_spore"
alldata$FMT_input[alldata$seqID %in% c("R_WT_1", "Rwt_2")] <- "rWT"
alldata$FMT_input[alldata$seqID %in% c("RAG205", "RAG207")] <- "yRag"
#write_tsv(alldata,"allhf_runs/allhf_seq.meta.merge.txt")

# 7/15/24: found another typo in original meta file:
test <- read.table(file="allhf_runs/allhf_seq.meta.merge.txt", header=T)
test %>% filter(str_detect(seqID, "GF2_hFMT1")) #GF2_hFMT1B has a 'd' in front of it
fixed <- test %>%
	mutate(seqID = str_replace_all(seqID, pattern = "dGF2_hFMT1", replacement = "GF2_hFMT1")) %>%
	mutate(sampleID = str_replace_all(sampleID, pattern = "dGF2_hFMT1", replacement = "GF2_hFMT1")) %>%
	write_tsv(.,"allhf_runs/allhf_seq.meta.merge.txt")
	
### step 1.2: combine metadata with phylotype (i.e., reference-based taxonomy):

### Converting taxonomic classifications to phylotype:

taxo<-read.table(file="allhf_runs/mothurfiles/allhumos6B.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary", header=TRUE)
taxo.samples<-taxo[, 6:ncol(taxo)]

sums <- read.table("allhf_runs/allhf_seq.meta.merge.txt", header=T)
meta <- sums[, 1:30]

# if you were to filter based on some type of name:
#cx<-read.table(file="mesal_otufrac.w.meta.txt", header=TRUE, sep="\t")
#names<-as.character(cx$seqID)
#taxo.names<-taxo.samples[, colnames(taxo.samples) %in% names]	#get only mesalamine samples
#taxo.filtered<-cbind(taxo[1:5], taxo.names)
#taxo.present<-taxo.filtered[rowSums(taxo.filtered[6:ncol(taxo.filtered)]) > 0, ]
#write.table(taxo.present, file="mothurfiles/mesal.only_wang.tax.summary.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)	

# filter out genus-level assignments only, and assign a phylum:
tax<-taxo
tax2<-tax[which(tax$taxlevel==2), ]
tax2[, c("rankID", "taxon")]
tax6<-tax[which(tax$taxlevel==6), ]
tax6$rankID<-gsub("^0.1.1.*", "20_Archaea_unclassified", tax6$rankID)
tax6$rankID<-gsub("^0.1.2.*", "20_Euryarchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.1.3.*", "20_Pacearchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.1.4.*", "20_Thaumarchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.1.5.*", "20_Woesearchaeota", tax6$rankID)
tax6$rankID<-gsub("^0.2.1\\..*", "10_Acidobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.2\\..*", "04_Actinobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.3\\..*", "20_Atribacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.4\\..*", "20_BRC1", tax6$rankID)
tax6$rankID<-gsub("^0.2.5\\..*", "11_Bacteria_unclassified", tax6$rankID)
tax6$rankID<-gsub("^0.2.6\\..*", "01_Bacteroidetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.7\\..*", "20_Candidatus_Saccharibacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.8\\..*", "20_Chlamydiae", tax6$rankID)
tax6$rankID<-gsub("^0.2.9\\..*", "20_Chloroflexi", tax6$rankID)
tax6$rankID<-gsub("^0.2.10..*", "20_Cloacimonetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.11..*", "20_Cyanobacteria/Chloroplast", tax6$rankID)
tax6$rankID<-gsub("^0.2.12..*", "20_Deferribacteres", tax6$rankID)
tax6$rankID<-gsub("^0.2.13..*", "20_Deinococcus-Thermus", tax6$rankID)
tax6$rankID<-gsub("^0.2.14..*", "20_Fibrobacteres", tax6$rankID)
tax6$rankID<-gsub("^0.2.15..*", "02_Firmicutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.16..*", "06_Fusobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.17..*", "20_Gemmatimonadetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.18..*", "20_Ignavibacteriae", tax6$rankID)
tax6$rankID<-gsub("^0.2.19..*", "20_Latescibacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.20..*", "20_Lentisphaerae", tax6$rankID)
tax6$rankID<-gsub("^0.2.21..*", "20_Microgenomates", tax6$rankID)
tax6$rankID<-gsub("^0.2.22..*", "20_Nitrospinae", tax6$rankID)
tax6$rankID<-gsub("^0.2.23..*", "20_Nitrospirae", tax6$rankID)
tax6$rankID<-gsub("^0.2.24..*", "20_Parcubacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.25..*", "20_Planctomycetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.26..*", "03_Proteobacteria", tax6$rankID)
tax6$rankID<-gsub("^0.2.27..*", "09_Spirochaetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.28..*", "08_Synergistetes", tax6$rankID)
tax6$rankID<-gsub("^0.2.29..*", "07_Tenericutes", tax6$rankID)
tax6$rankID<-gsub("^0.2.30..*", "20_Thermotogae", tax6$rankID)
tax6$rankID<-gsub("^0.2.31..*", "05_Verrucomicrobia", tax6$rankID)
tax6$rankID<-gsub("^0.2.33..*", "20_candidate_division_WPS-1", tax6$rankID)
tax6$rankID<-gsub("^0.3.1..*", "11_unknown_unclassified", tax6$rankID)
colnames(tax6)[2]<-"phylum"

# filter some columns, turn into a matrix and check for duplication:
subtax6<-subset(tax6, select=-c(taxlevel, daughterlevels))
subtax6<-subtax6[order(subtax6$phylum, -subtax6$total), ]
taxmatrix<-subtax6[, c(4:ncol(subtax6))]
which(duplicated(subtax6$taxon))
	# none are duplicated, but follow trhough with remainder if they were:
	# fix the duplicated row:
	#subtax6$taxon<-as.character(subtax6$taxon)
	#subtax6$taxon[487]<-"Actinobacteria_unclassified2"
	#subtax6$taxon<-as.factor(subtax6$taxon)
rownames(taxmatrix)<-subtax6$taxon
genera<- taxmatrix[, colSums(taxmatrix)>2500,]
	# get rel. abund fraction:
genmatrix<-as.data.frame(t(genera))
genera.fr<-genmatrix/rowSums(genmatrix)*100
genus.fr<-t(genera.fr)
all.genera<-cbind(subtax6[1:3], genus.fr)
#write.table(all.genera, file="allhf_runs/allhf_all.genera.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
	# get top 1%:
phyla<-subtax6[1:3]
genus1<- genus.fr[rowSums(genus.fr>=1)>=1,]
namelist<-as.character(rownames(genus1))
phyla1p<-phyla[phyla$taxon %in% namelist, ]
genera1<-cbind(phyla1p, genus1)
	# get top 2%
genus2<- genus.fr[rowSums(genus.fr>=2)>=2,]
namelist<-as.character(rownames(genus2))
phyla2p<-phyla[phyla$taxon %in% namelist, ]
genera2<-cbind(phyla2p, genus2)
#write.table(genera2, file="allhf_runs/allhf_genfrac2p.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# let's also add metadata and colors for later:
summary(as.factor(genera2$phylum))
	bac<-colorRampPalette(c("darkgreen","green3","lightgreen","seagreen", "lightgreen"))(n=10)
	firm<-colorRampPalette(c("midnightblue","mediumblue","blue3","blue","dodgerblue4","dodgerblue1","deepskyblue4","deepskyblue1","skyblue3","skyblue","steelblue4","steelblue1","royalblue4","royalblue1","slateblue4","purple3","orchid3","plum4","plum1","pink3","pink","lightpink1","lightpink3","palevioletred4","palevioletred1","magenta4","deeppink4","mediumvioletred","magenta3","magenta1","thistle"))(n=52)
	pro<-colorRampPalette(c("yellow2","darkgoldenrod3","goldenrod2","orange2","yellow4"))(n=17)
	actino<-colorRampPalette(c("tan", "brown", "darkred"))(n=5)
	verruco<-c("hotpink")
	#fuso<-colorRampPalette(c("tan", "red", "darkred"))(n=3)
	fuso<-colorRampPalette(c("red"))(n=1)
	synten<-colorRampPalette(c("cyan", "darkcyan"))(n=2)
	uncl<-c("black")
	other<-colorRampPalette(c("grey50"))(n=4)
	#other<-colorRampPalette(c("grey50"))(n=1)
color<-c(bac, firm, pro, actino, verruco, fuso, uncl, synten, other)
genera2<-cbind(genera2[1:3], color, genera2[4:ncol(genera2)])
#write.table(genera2, file="allhf_runs/allhf_genfrac2p.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# read in file and combine with meta:
#genera2<-read.table(file="mesal_genfrac2p.txt", header=TRUE, sep="\t", comment.char = "")
genbar<-genera2
rownames(genbar)<-genbar$taxon
	rm_g<-subset(genbar, select =-c(phylum, taxon, color, total) )
	barg<-as.data.frame(t(rm_g))
	barg$other<-100-rowSums(barg)
	others<-100-colSums(barg)
	barg$sampleID<-rownames(barg)
	col.gen<-c(as.character(genbar$color), "grey80")
	#barg$sampleID<-gsub("X19_", "19_", barg$sampleID)
bar<-merge(meta, barg, by.x=c("seqID"), by.y=c("sampleID"))
#write.table(bar, 'allhf_runs/allhf_genfrac2p.meta.txt',quote=FALSE,sep="\t", col.names=TRUE, row.names=FALSE)


```

###-------------###

### Step 2: subset metadata / mothur-generated raw data for current project
	- contains info for all samples used in current study (prefix = HFm)
	- these were input files used for analyses

```{r}
library(tidyr)
library(reshape2)
library(tidyverse)

# cd /Data/16S

### Step 2.1: subset mothur-generated summary statistics; add some additional analysis-specific groups

df<-read.table(file="allhf_runs/allhf_seq.meta.merge.txt", header=T)

## For this experiment we are using the following experimental groups:
	# only SPF experiments; exclude HF5 and HF6, since those were contaminated experiments
	# all mFMT groups (no filtrate; no tributyrin groups, though)
	# the following hFMT groups: AMS001, AMS005, D2_A, D4_A, D5_A, R4_F, rWT, yWT, yRag, yWT_spore

# first, let's subset the data for all mice:
# we do not want the contaminated OR the tributyrin treatments here
df_nocontam <- df[!df$Experiment %in% c("HF5", "HF5-contam", "HF6", "HF6-contam") , ] %>%
					 .[!.$FMT %in% c("hFMT6_tri"), ]
d <- df_nocontam[df_nocontam$Mouse_status %in% c("SPF-Young") & 
				df_nocontam$FMT_input %in% c("AMS001", "AMS005", "D2_A", "D4_A", "D5_A", "R4_F", "rWT", "yWT", "yRag", "yWT_spore", "none") &
				df_nocontam$day > -8 & 
				df_nocontam$Cdiff_infection %in% c("CDI"), ] %>%
				droplevels(.)

### note: let's also add the inputs that are relevant to these, and spit that out for future use:
summary(df$sampleType)		# all inputs should be under "FMTinput"
df[df$seqID %in% c('AMS01'), ]			# need to somehow also include these...
df.inputs <-df[df$sampleType %in% c('FMTinput'), ]		# these are all the types possible...in the WHOLE data set

# pull out inputs relevant to this experiment
df.inputs[df.inputs$Experiment %in% c('HF1', 'HF2', 'HF3', 'HF4', 'HF7', 'HF8', 'HF9', 'HF10'), c("seqID", "FMT_input", "Experiment")]
	# looks like this doesn't exactly match--let's get a couple extra replicates for some of them, and remove some of them

# filtering:
summary(d$FMT_input)		# these are the ones you want, plus some extra from the human samples
input.list <- as.character(unique(d$FMT_input))
df.inputs <- df.inputs[df.inputs$FMT_input %in% input.list & df.inputs$Experiment %in% c('HF1', 'HF2', 'HF3', 'HF4', 'HF7', 'HF8', 'HF9', 'HF10'), ]	
df.inputs[, c("seqID", "FMT_input", "Experiment")]
summary(as.factor(df.inputs$FMT_input))

# let's also add the actual human samples into this:
df[df$sampleType %in% c('human'),  c("seqID", "FMT_input", "Experiment")]		# select the human samples that match our inputs
df.inputs2 <- df[df$sampleType %in% c('human') & df$seqID %in% c('AMS01', 'AMS05', 'D2_A', 'D4_A', 'D5_A', 'R4_F'), ]

# finally, if you wanted to add mouse_survey samples in here too (includes untreated mouse cecal samples, as well)
#df.inputs3 <- df[df$sampleType %in% c('mouse_survey'), ]	# let's make the relevant mouse samples here 'FMT_input', and add the appropriate FMT types
	# note: not adding these right now, EXCEPT for the healthy cecal samples:
healthy.cecal <- df[df$seqID %in% c("healthy_cecal_1913_1",  "healthy_cecal_1913_2", "healthy_cecal_1912_1", "healthy_cecal_1912_2"), ]
	
# add these together, and to the filtered dataset that also includes the mice:
d.w.inputs <- rbind(d, df.inputs, df.inputs2, healthy.cecal)

d <- d.w.inputs

### set a group that you want to compare:
# in this case, let's keep all the hFMT's in the same group, but split the mouse into the following:
	# yWT = regular mFMT
	# rWT, yRag, yWT_spore = other mice (in case someone asks about other mice)
d$group[d$FMT_type %in% c("hFMT")] <- "hFMT"
d$group[d$FMT_input %in% c("yWT")] <- "mFMT"
d$group[d$FMT_input %in% c("none")] <- "noFMT"
d$group[d$FMT_input %in% c("rWT", "yRag", "yWT_spore")] <- paste(d$FMT_type[d$FMT_input %in% c("rWT", "yRag", "yWT_spore")], "other", sep="_")
summary(as.factor(d$group))
d <- d[order(d$group, d$day),]
levels(as.factor(d$group))

# let's modify the FMTinput and group for the human samples, as well:
d$group[d$sampleType %in% c('human')] <- "hFMT"
d$FMT_input[d$seqID =='AMS01'] <- 'AMS001'
d$FMT_input[d$seqID =='AMS05'] <- 'AMS005'
d$FMT_input[d$seqID =='D2_A'] <- 'D2_A'
d$FMT_input[d$seqID =='D4_A'] <- 'D4_A'
d$FMT_input[d$seqID =='R4_F'] <- 'R4_F'
d$FMT_input[d$seqID =='D5_A'] <- 'D5_A'
d$group[d$sampleType %in% c('mouse_survey')] <- "mFMT"
# dim(d) 4039

#write.table(d, file="HFm_summary.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)	

### also added some additional analysis groups:

# ok, let's actually do the NMDS now
data <- read.table("HFm_summary.txt", header=T, sep="\t")
levels(as.factor(data$group))		#we want to match the same levels as above, EXCEPT, also match the poop donors

# we need to add another group-type that separates out pre-post-FMT
data$GROUP <- data$FMT_input
data$GROUP[data$sampleType %in% c("mouse_survey") & data$type=='cecal'] <- "pre-abx"
data$GROUP[data$day %in% c(-7) & data$type=='fecal'] <- "pre-abx"
data$GROUP[data$day %in% c(0) & data$type=='fecal'] <- "post-cef"
data$GROUP[data$day %in% c(1,4) & data$type=='fecal'] <- "post-cdi"
data$GROUP[data$day %in% c(7,9,11) & data$type=='fecal'] <- "post-vanco"
data$GROUP <- as.factor(factor(data$GROUP, levels=c('pre-abx', 'post-cef', 'post-cdi', 'post-vanco', 
								'D2_A', 'D4_A', 'D5_A', 'R4_F', 'AMS001', 'AMS005', 
								'none',
								'rWT', 'yRag', 'yWT_spore', 'yWT')))
# set colors as close to decided key as possible:
# these are what you have set
group.levels = c("healthy", "hFMT1", "hFMT2", "hFMT3", "mFMT", "mFMT_other", "noFMT" )
colors.group = c("#3B99B1", "#EACB2B", "#E87700", "#E8A419", "#56B29E", "#9FC095", "#F5191C")
# designing from this:
#hcl.colors(15, palette="Zissou 1")		# this will list the colors
#par(mar=c(0.5, 0.5, 1, 0.2), mfrow=c(2,1), oma = c(4, 4, 0.5, 0.5))
#barplot(1:length(hcl.colors(15, palette="Zissou 1")), col=hcl.colors(15, palette="Zissou 1"), las=2, names=hcl.colors(15, palette="Zissou 1"))
#barplot(1:length(group.levels), col=colors.group, las=2, names=group.levels)
# for abx (pre-FMT), want a grey color scheme:
#hcl.colors(5, palette="Light Grays")
#barplot(1:length(hcl.colors(5, palette="Light Grays")), col=hcl.colors(5, palette="Light Grays"), las=2, names=hcl.colors(5, palette="Light Grays"))
nmds.col <- c('#3B99B1', "#474747", "#A8A8A8", "#E2E2E2", 
						'#EACB2B', '#E79812', '#E8A91B', '#ED5300', '#E87700', '#E8A419', 
						'#F5191C',
						'#8BBD94', '#A6C293', '#C1C88C', '#56B29E')
# define which ones are inputs for pch
data$TYPE <- data$type
data$TYPE[data$sampleType %in% c("human", "FMTinput")] <- "FMTinput"
data$TYPE <- as.factor(data$TYPE)

dpch <- c(17, 19, 8) #cecal, fecal, FMTinput; colored circles
dpch2 <- c(24, 21, 8) #cecal, fecal, FMTinput
dcex <- c(1,1,2)

# to save the groups we just made:
#write_tsv(data, "Data/HFm_summary.nmds.new.groups.txt")



###----------------###

### Step 2.2: subset OTU abundances for project samples:

#otus <- read.table("~/Box\ Sync/Projects/HUMO.w.RM/allhumos6/allhumos6B/flux/mothurfiles/allhumos6B.final.shared", header=T, sep="\t", row.names=2) %>%
#		select(., -c(label, numOtus))		# use this f you do NOT need to rerun in mothur!
otus <- read.table("allhf_runs/mothurfiles/allhumos6B.final.shared", header=T, sep="\t")
		# use this if you need to rerun as-is in mothur
# read in what samples we are using, plus get rid of samples < 2500 sequences
d <- read.table("HFm_summary.txt", header=T, sep="\t") %>%
		.[!is.na(.$simpson), ] %>%
		.[.$nseqs > 2500, ]

# filter out the samples we do not need from OTU file; remove 0 lines from this; and remove samples with < 2500 sequences total
RemoveZeros <- function(x) {
    if(is.numeric(x)) {
        sum(x) > 0
    } else {
        TRUE
    }
}

otus.matrix <- otus %>% 
	select(., -c(label, numOtus, Group))
rownames(otus.matrix) <- otus$Group
otus.filtered <- otus.matrix[rownames(otus.matrix) %in% as.character(d$seqID), ] %>%
	.[, sapply(.,  RemoveZeros)] %>%
	cbind(otus[otus$Group %in% as.character(d$seqID), c(1:3)], .)
otus.filtered$numOtus <- length(otus.filtered[4:length(otus.filtered)])	#if you wanted to keep the numOtus...

#write_tsv(otus.filtered, "HFm.otus.shared")



```

###-------------###

### Step 3: organize OTU taxonomy file to match phylotype color scheme


```
###-------------####

# organizing OTU taxonomy file:
OTUids <- read.table(file="Data/16S/allhf_runs/mothurfiles/allhumos6B.final.0.03.cons.taxonomy", header=T)

# split taxonomy into groups:
otus <- OTUids %>% 
	mutate(Taxonomy = gsub("([0-9])", "", Taxonomy)) %>%
	separate_wider_delim(Taxonomy, delim = "();", names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
	select(-Species)

# organize and add colors by phylum:
genus_counts <- otus %>%
	group_by(Phylum, Genus) %>% 
	summarize(count=n()) %>% 
	arrange(Phylum, desc(count))
		# too many!
	
# there are a lot, so let's match this up with the previous taxonomy colors and fill in the blanks after:
taxcolors <- read.table("Data/16S/allhf_runs/allhf_genfrac2p.txt", header=T, comment.char = "") %>%
	#write_tsv(., "Data/HFm_genfrac2p_colors.txt") 
	select(taxon, color)

# define phyum colors; replace with that if not found from color palette above
otus_colors <- otus %>%
	left_join(., taxcolors, join_by("Genus" == "taxon")) %>%
	mutate(phylum_color = case_when(
		Phylum == "Firmicutes" ~ "lightcyan", 
		Phylum == "Bacteroidetes" ~ "darkseagreen1", 
		Phylum == "Verrucomicrobia" ~ "lightpink", 
		Phylum == "Actinobacteria" ~ "seashell", 
		Phylum == "Proteobacteria" ~ "lightgoldenyellow", 
		Phylum == "Fusobacteria" ~ "salmon", 
		Phylum == "Bacteria_unclassified" ~ "black", 
		TRUE ~ "grey50")
			) %>%
	mutate(color = coalesce(color,phylum_color)) #%>%
	#write_tsv(., "Data/16S/taxonomy_meta.txt")

```

```

###-------------###

### Step 4: combine all sample meta data together for Table S1 / S2
	- all 16S metadata (with CFUs, weights)
	- combined with all HPLC and UM-Metabolomics core data
	- Metabolon data: Table S2

```{r}

library(dplyr)


# wd: /Data

# Step 4.1: 16S rRNA seq metadata, Table S1:

d <- read.table("16S/HFm_summary.txt", header=T, sep="\t") %>%
		.[!is.na(.$simpson), ] %>%
		.[.$nseqs > 2500, ] %>%
		select(seqID, sampleID, mouse, percent_weight, sex, FMT_input, FMT_type, Experiment, day, type, CFU, sampleType, group, shannon, nseqs) #%>%
		#write_tsv(., "TableXX_16S.forSRA.txt")
# this information was also used for SRA / BioSample info
# added project information from BioProject / Biosample:
biosample <- read.table("../../SRA/BioSampleObjects_16S.txt", header=T, sep="\t") %>%
			select(Accession, SPUID) %>%
			mutate(across('SPUID', str_replace, 'unclassified_', '')) %>%
			mutate(across('SPUID', str_replace, '_metagenome', ''))
biosample2 <- read.table("../../SRA/BioSampleObjects_UMFMT.txt", header=T, sep="\t") %>%
			select(Accession, Sample.Name) %>%
			mutate(across('Sample.Name', str_replace, 'unclassified_', '')) %>%
			mutate(across('Sample.Name', str_replace, '_metagenome', ''))
d <- read.table("TableXX_16S.forSRA.txt", header=T, sep="\t") %>%
		mutate(BioProjectID = case_when(
							seqID %in% c("D2_A", "D4_A", "D5_A", "R4_F") ~ "PRJNA384621", 			
							TRUE ~ "PRJNA1168499")
						) %>%
		left_join(., biosample, join_by("seqID" == "SPUID")) %>%
		left_join(., biosample2, join_by("seqID" == "Sample.Name")) %>%
		mutate(BioSampleID = coalesce(Accession.x, Accession.y)) %>%
		select(seqID, BioProjectID, BioSampleID, sampleID, mouse, percent_weight, sex, FMT_input, FMT_type, Experiment, day, type, CFU, sampleType, group, shannon, nseqs) #%>%
		#write_tsv(., "TableXX_16S.forSRA.txt")

		
# Step 4.2, read in metabolon data, Table S2:
d.metab <- read.csv('Data/metabolomics/Metabolon/metabData_metabolon.csv', header = TRUE, sep = ",", row.names=1) %>%
					write_tsv(., "TableS2_metabolon.txt")

# Step 4.3, HPLC (fecal) SCFA data used in project, Table S3:
d.scfa <- read.table("metabolomics/HPLC/HFm_scfa.seq.meta.merge.txt", header=T, sep="\t") %>%
					filter(type == "fecal" & sampleType =="mouse") %>%
					filter(!(group =="mFMT_other" & day > 11)) %>%
					filter(day > -8) %>%
					mutate(group = case_when(
							day %in% c(-7, 0, 1, 4, 9, 11) ~ "preFMT", 			
							TRUE ~ group),
						comp = factor(group, levels=c("mFMT", "mFMT_other", "hFMT", "noFMT", "preFMT"))
						)	%>%
					#select(seqID, propionate:butyrate) %>%
					#coalesce_join(d, ., by = "seqID") %>%
					select(seqID, sampleID, mouse, percent_weight, sex, FMT_input, FMT_type, Experiment, day, type, CFU, sampleType, group, propionate:butyrate) #%>%
					#write_tsv(., "TableS3_HPLC.scfa.txt")
					
					
# Step 4.4, output ba/scfa from cecal content done by UM core, Table S4:
d.ba <- read_tsv("metabolomics/UM-MetaboCore/HFm_scfa.ba.cecal.seq.meta.merge.txt") %>%
			filter(!group=="preFMT") %>%
			unite("group_c", day,group, sep ="_", remove=FALSE) %>%
			mutate(group_c = factor(group_c, levels=c('-7_healthy', '21_mFMT', '21_hFMT', '21_noFMT', 
										'42_mFMT', '42_hFMT', '42_noFMT'))) %>%
			select(seqID, sampleID, mouse, percent_weight, sex, FMT_input, FMT_type, Experiment, day, type, CFU, sampleType, group, DCA:Octanoate_nmol) #%>%
			#write_tsv(., "TableS4_UM.scfa.ba.txt")
			
										
# let's combine tables 1, 3, and 4 together, as they contain mostly the same info
# note: could not get coalesce_join to work
all.tables <- d %>% 
			full_join(d.scfa, by = "seqID") %>%
			mutate(across(
    			.cols = grep("\\.x$", names(.), value = TRUE),   # Select columns ending with .x
   				 .fns = ~ coalesce(.x, get(sub("\\.x$", ".y", cur_column())))  # Coalesce with corresponding .y column
  				)) %>%
  			rename_with(~ gsub("\\.x$", "", .), grep("\\.x$", names(.), value = TRUE)) %>%  # Rename .x columns
  			select(-grep("\\.y$", names(.), value = TRUE)) %>%
  				# now, join and coalesce with d.ba table
  			full_join(d.ba, by = "seqID") %>%	
  			mutate(across(
    			.cols = grep("\\.x$", names(.), value = TRUE),   # Select columns ending with .x
   				 .fns = ~ coalesce(.x, get(sub("\\.x$", ".y", cur_column())))  # Coalesce with corresponding .y column
  				)) %>%
  			rename_with(~ gsub("\\.x$", "", .), grep("\\.x$", names(.), value = TRUE)) %>%  # Rename .x columns
  			select(-grep("\\.y$", names(.), value = TRUE)) %>%
  			write_tsv(., "TableS1_16S.scfa.ba.txt")
  