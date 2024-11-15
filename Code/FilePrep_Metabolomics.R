### HUMOFMT/allhumos6B: Additional File prep for metabolomics data (Figure 5)
### 9.30.24
### Anna M. Seekatz

### Inputs: 
	- /Data/16S/HFm_summary.txt (contains most of the necessary metadata; OTU summary metrics)
	- /Data/metabolomics/Metabolon:
		- metadata_metabolon.csv
		- metabData_metabolon.csv

### Outputs: into /Data
	- /Data/metabolomics:
		- 
	- /Data/metabolomics/Metabolon:
		- sample_metadata_metabolon.csv: metadata for samples sent to Metabolon
		- raw_scaled_metabolon.csv: raw Metabolon data (select information)
		- compound_meta_metabolon.csv: used for compound information


###-------------###

### Step 1: clean up some Metabolon generated data for analysis and comine with metadata


```{r}
library(tidyr)
library(reshape2)
library(tidyverse)
library(readr)
library(ggplot2)


# filter data for manuscript:
meta <- read.csv('/Data/metabolomics/Metabolon/metadata_metabolon.csv', header = TRUE, sep = ",") %>%
	filter(!Experiment %in% c("GF2", "GF3") & !status %in% c("noCDI","human")) %>%
	mutate(status = factor(status, levels = c("cleared", "colonized")))
#write.csv(meta, file="ForPub/sample_metadata_metabolon.csv", row.names=F)

data <- read.csv('/Data/metabolomics/Metabolon/metabData_metabolon.csv', header = TRUE, sep = ",", row.names=1) %>%
	select(-COMP_ID, -PATHWAY.SORTORDER, -SUPERPATHWAY, -SUBPATHWAY, -KEGG, -PLATFORM, -RI, -MASS, -PUBCHEM, -CAS, -HMDB, -CHEMICAL_ID) %>%
	select(BIOCHEMICAL, meta$metaboID) %>%
	drop_na()
#write.csv(data, file="ForPub/raw_scaled_metabolon.csv", row.names=F)

# setup compound metadata, adding direction of compound amounts based on colonized/cleared:
meta <- read.csv('/Data/metabolomics/Metabolon/sample_metadata_metabolon.csv', header = TRUE, sep = ",") %>%
	#mutate_at(., vars(type:status), as.factor) # if only a few conversions to factors
	mutate_if(., is.character, as.factor)
# select compounds up/down with cleared v colonized
data_dir <- data %>%
  gather(key="metaboID", value="ScaledImpData", -BIOCHEMICAL) %>%
  mutate(metaboID=factor(metaboID, levels=levels(meta$metaboID))) %>%
  left_join(., meta, by="metaboID") %>%
  select(BIOCHEMICAL, ScaledImpData, status) %>%
  group_by(status, BIOCHEMICAL) %>%
  mutate(ScaledImpData=log10(ScaledImpData)) %>%
  summarize(mean=mean(ScaledImpData)) %>%
  spread(key=status, value=mean) %>%
  mutate(dir=factor(ifelse(colonized<cleared,"Up in cleared","Down in cleared"),levels=c("Up in cleared","Down in cleared"))) %>%
  select(BIOCHEMICAL, dir)
  
# select compounds up/down with mFMT v noFMT
data_dir_groups <- data %>%
  gather(key="metaboID", value="ScaledImpData", -BIOCHEMICAL) %>%
  mutate(metaboID=factor(metaboID, levels=levels(meta$metaboID))) %>%
  left_join(., meta, by="metaboID") %>%
  select(BIOCHEMICAL, ScaledImpData, group) %>%
  group_by(group, BIOCHEMICAL) %>%
  mutate(ScaledImpData=log10(ScaledImpData)) %>%
  summarize(mean=mean(ScaledImpData)) %>%
  spread(key=group, value=mean) %>%
  mutate(dir_noFMT=factor(ifelse(noFMT<mFMT,"Up in mFMT","Down in mFMT"),levels=c("Up in mFMT","Down in mFMT"))) %>%
  mutate(dir_hFMT=factor(ifelse(hFMT<mFMT,"Up in mFMT","Down in mFMT"),levels=c("Up in mFMT","Down in mFMT"))) %>%
  mutate(dir_healthy=factor(ifelse(healthy<mFMT,"Up in mFMT","Down in mFMT"),levels=c("Up in mFMT","Down in mFMT"))) %>%
  select(BIOCHEMICAL, dir_noFMT, dir_hFMT, dir_healthy)
  
# then, combine with metadata that relates to compounds
compound_meta <- read.csv('metabData_metabolon.csv', header = TRUE, sep = ",", row.names=1) %>%
	drop_na(c21_1835_3) %>%		#just getting rid of human samples here
	select(COMP_ID:HMDB) %>%
	left_join(., data_dir, by="BIOCHEMICAL") %>%
	left_join(., data_dir_groups, by="BIOCHEMICAL")
#write.csv(compound_meta, file="/Data/metabolomics/Metabolon/compound_meta_metabolon.csv", row.names=F)

```

###-------------###

### Step 2: combining data HPLC (SCFA, fecal) sample info


```{r}

```