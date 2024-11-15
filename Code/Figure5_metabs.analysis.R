#### Figure 5: metabolomics analysis
### 7.28.24
### Anna M. Seekatz

### Figures:
	- Fig. 5A: NMDS of all metabolon data
	- Fig. 5B: RF of metabolon
	
### Input files: located in Data/metabolomics/Metabolon
	- metadata_metabolon.csv
	- metabData_metabolon.csv (based off of 'ScaledImpData' tab directly from '/CLEM-01-20VW+/CLEM-01-20VW+ CDT CECAL CONTENT.XLSX')
	- sample_metadata_metabolon.csv (same as a above; just with sample metadata in a different format)
	- raw_scaled_metabolon.csv  (same as above; just with values as a dataframe)
	- compound_meta_metabolon.csv (same as above; just with compound info as a separate dataframe)
	


###-------------------------
### Figure 5A: NMDS based off of all metabolon values

### inputs:
	- metadata_metabolon.csv
	- metabData_metabolon.csv
	- ForPub/sample_metadata_metabolon.csv
	- ForPub/raw_scaled_metabolon.csv
	- HFm_scfa.fecal.seq.meta.merge.txt
	- HFm_scfa.ba.cecal.seq.meta.merge.txt

```{r}
###-------------###
### Figure 5A, NMDS:

library(vegan)
library(tidyr)
library(tidyverse)
library(dplyr)
library(janitor)

# read in data, selecting only SPF experiments and relevant data for metab info:
meta <- read.csv('Data/metabolomics/Metabolon/metadata_metabolon.csv', header = TRUE, sep = ",") %>%
	filter(!Experiment %in% c("GF2", "GF3") & !status %in% c("noCDI","human"))
data <- read.csv('Data/metabolomics/Metabolon/metabData_metabolon.csv', header = TRUE, sep = ",", row.names=1) %>%
	select(-COMP_ID, -PATHWAY.SORTORDER, -PLATFORM, -RI, -MASS, -PUBCHEM, -CAS, -HMDB, -CHEMICAL_ID)

# set dataframe for calculating distance
# for now, eliminated all NAs, since those were compounds different between mice and humans included in datasets
d <- data %>% select(-SUPERPATHWAY, -SUBPATHWAY, -KEGG) %>%
		column_to_rownames("BIOCHEMICAL") %>%
		drop_na() %>%
		t() %>%
  		as.data.frame() %>%
  		rownames_to_column("metaboID") %>%
  		left_join(meta, ., by="metaboID") %>%
  		select(!(seqID:status)) %>%
  		column_to_rownames("metaboID") %>%
  		clean_names() #779 compounds left

# look at dispersion of homogeneity:
set.seed(13)
d.dist <- vegdist(d, method="bray")
d.div <- adonis2(d.dist ~ group, data = meta, permutations = 999, method="bray")
dispersion <- betadisper(d.dist, group=meta$group)
anova(dispersion)
TukeyHSD(dispersion)
# also tried method = horn; similar pattern, just more emphasized
# means that groups within themselves, at least compared to other groups, are very homogeneous? not variable within their own group?

# look at adonis for group differences:
d.div <- adonis2(d.dist ~ group, data = meta, permutations = 999, method="bray")
#d.div <- adonis2(d.dist ~ group, data = meta, permutations = 999, method="horn")
	# cant do posthoc, but whatever

# calculate NMDS
set.seed(13)
metab_nmds <- metaMDS(d.dist) %>%
	scores() %>%
	as_tibble(rownames="metaboID")
# low stress, so good

# rejoin with metadata for graphing:
meta_nmds <- inner_join(meta, metab_nmds, by="metaboID")
#write_tsv(meta_nmds, "Data/HFm_bray.nmds.meta.txt")

# adonis test on samples (overall group):
test <- adonis2(as.dist(d.dist)~meta_nmds$group)

# visualize
data <- meta_nmds					
data <- data[order(data$subgroup), ]
data$group <- as.factor(data$group)
data$subgroup <- as.factor(data$subgroup)
data$status <- as.factor(data$status)
data$subgroup <- as.factor(factor(data$subgroup, levels=c(
								'hFMT1', 'hFMT2', 'hFMT3',
								'noFMT',
								'mFMT', 'healthy')))
# colors from before:
#group.levels = c("healthy", "hFMT1", "hFMT2", "hFMT3", "mFMT", "mFMT_other", "noFMT" )
#colors.group = c("#3B99B1", "#EACB2B", "#E87700", "#E8A419", "#56B29E", "#9FC095", "#F5191C")
nmds.col <- c("#EACB2B", "#E87700", "#E8A419", 
						"#F5191C",
						"#56B29E", "#3B99B1")
#barplot(1:length(nmds.col), col=nmds.col, las=2, names=paste(levels(data$subgroup), nmds.col, sep="\n"), cex.names=0.8)

# define some shapes around other factors:
dpch <- c(17, 19) #cleared v colonized
dpch2 <- c(24, 21) #cleared v colonized
#dcex <- c(1,1,2)

# test out plotting by itself:
par(mar = c(5, 4, 4, 8), xpd = TRUE)
plot(data$NMDS1, data$NMDS2, xaxt = "n", yaxt = "n",
			#xlim=c(min(data$NMDS1, na.rm=T), max(data$NMDS1, na.rm=T) ), ylim=c(min(data$NMDS2, na.rm=TRUE), max(data$NMDS2, na.rm=TRUE)), 
			xlim=c(-0.9,.9), ylim=c(-.9,.9), 
			ylab="NMDS2", xlab="NMDS1", 
			col=adjustcolor(nmds.col[data$subgroup], alpha=0.8), 
			pch=dpch[data$status])
# add tick marks:
axis(1, at = c(-.9, -.5, 0, .5, .9))
axis(2, at = c(-.9, -.5, 0, .5, .9))
# add black points around cleared only
points(data$NMDS1[data$status %in% c("cleared")], data$NMDS2[data$status %in% c("cleared")], col='grey20', 
	pch=dpch2[data$status[data$status %in% c("cleared")]])
# legend:
legend("topright", inset = c(- 0.5, 0),                   # Create legend outside of plot
       legend = as.character(levels(data$subgroup)),
       col=adjustcolor(nmds.col, alpha=0.8),
       cex=1, pt.cex=c(1), pch=c(19,19,19,19,17,17))
#or
legend("bottomleft", #inset = c(- 0.5, 0),                   # Create legend outside of plot
       legend = as.character(levels(data$subgroup)),
       col=adjustcolor(nmds.col, alpha=0.8),
       cex=1, pt.cex=c(1), pch=c(19,19,19,19,17,17))
       
### Figure 5B: pairwise distances based on Bray-Curtis

require(readr)
require(tidyr) 
library(tidyverse)
library(ggplot2)
library(FSA)
library(janitor)

# read in data:
sample_meta <- read.csv('Data/metabolomics/Metabolon/sample_metadata_metabolon.csv', header = TRUE, sep = ",") %>%
	mutate_if(., is.character, as.factor)
metabs <- read.csv('Data/metabolomics/Metabolon/raw_scaled_metabolon.csv', header=TRUE, row.names=1, sep = ",") %>%
			t(.)

# calculate bray-curtis distance:
set.seed(13)
metabs.dist <- metabs %>%
	vegdist(method="bray")	%>% #default is bray
	as.matrix() %>%
	as_tibble(rownames = "samples") %>%
	pivot_longer(-samples) %>%
	filter(samples < name) 

### step 1: modify distance matrix; merge with meta data; filter out relevant comparisons
meta.dist <- metabs.dist %>% 
			rename(sampleID1 = samples) %>%
			rename(sampleID2 = name) %>%
			inner_join(., sample_meta, by = c("sampleID1" = "metaboID")) %>%
			inner_join(., sample_meta, by = c("sampleID2" = "metaboID")) %>%
			select(sampleID1, sampleID2, group.x, group.y, status.x, status.y, value) %>%
				# merge both sampleIDs with metadata
			rename_with(~ c("group1", "status1", "group2", "status2"), 
				all_of(c("group.x", "status.x", "group.y", "status.y"))
				) %>%
			filter(!sampleID1 == sampleID2) %>%	
				# get rid of self-comparisons (which are 0)
				# define comparison
			mutate(comparison = case_when(
			group1 == "hFMT" & group2 == "mFMT" | group2 == "hFMT" & group1 == "mFMT" ~ "mFMT-hFMT",
			group1 == "noFMT" & group2 == "mFMT" | group2 == "noFMT" & group1 == "mFMT" ~ "mFMT-noFMT",
			group1 == "hFMT" & group2 == "noFMT" | group2 == "hFMT" & group1 == "noFMT" ~ "hFMT-noFMT",
			group1 == "hFMT" & group2 == "healthy" | group2 == "hFMT" & group1 == "healthy" ~ "healthy-hFMT",
			group1 == "mFMT" & group2 == "healthy" | group2 == "mFMT" & group1 == "healthy" ~ "healthy-mFMT",
			group1 == "noFMT" & group2 == "healthy" | group2 == "noFMT" & group1 == "healthy" ~ "healthy-noFMT",
			group1 == "mFMT" & group1 == group2 ~ "within-mFMT",
			group1 == "hFMT" & group1 == group2 ~ "within-hFMT",
			group1 == "noFMT" & group1 == group2 ~ "within-noFMT",
			group1 == "healthy" & group1 == group2 ~ "within-healthy",
			TRUE ~ NA_character_),
			comparison = factor(comparison, levels=c("healthy-mFMT", "healthy-hFMT", "healthy-noFMT", 
												"mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT", 
												"within-healthy", "within-mFMT", "within-hFMT", "within-noFMT"))
				) %>% 	
			#filter(!is.na(comparison)) 
			
# graph
level_orders <- c("healthy-mFMT", "healthy-hFMT", "healthy-noFMT", 
					"mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT", 
					"within-healthy", "within-mFMT", "within-hFMT", "within-noFMT")

meta.dist %>% ggplot(aes(x=comparison, y=value)) +
	geom_jitter(width=0.25, color="grey", alpha=0.6, size = 0.5) +
	stat_summary(fun.data=median_hilow, color="red", size=0.5,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="Bray-Curtis distances") +
	scale_x_discrete(breaks=level_orders, 
			labels=str_replace(level_orders, "-", ":")) +
	scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.2)) +
	geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey20") +
	geom_vline(xintercept = 6.5, linetype = "dashed", color = "grey20") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

## stats:
library(FSA)
# inputs:
meta.dist %>% subset(., comparison %in% c("healthy-mFMT", "healthy-hFMT", "healthy-noFMT")) %>%
	dunnTest(value ~ comparison, data = .)
#                    Comparison         Z      P.unadj        P.adj
#1  healthy-hFMT - healthy-mFMT  5.748791 8.988368e-09 1.797674e-08
#2 healthy-hFMT - healthy-noFMT -4.449065 8.624509e-06 8.624509e-06
#3 healthy-mFMT - healthy-noFMT -9.073609 1.151390e-19 3.454169e-19
meta.dist %>% subset(., comparison %in% c("mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT")) %>%
	dunnTest(value ~ comparison, data = .)
#               Comparison         Z      P.unadj        P.adj
#1  hFMT-noFMT - mFMT-hFMT  3.359513 7.807986e-04 7.807986e-04
#2 hFMT-noFMT - mFMT-noFMT -4.056111 4.989655e-05 9.979310e-05
#3  mFMT-hFMT - mFMT-noFMT -6.939901 3.923754e-12 1.177126e-11
meta.dist %>% subset(., comparison %in% c("within-healthy", "within-mFMT", "within-hFMT", "within-noFMT")) %>%
	dunnTest(value ~ comparison, data = .)
#                     Comparison          Z      P.unadj        P.adj
#1  within-healthy - within-hFMT -1.8076161 7.066628e-02 2.119988e-01
#2  within-healthy - within-mFMT  2.9490599 3.187422e-03 1.274969e-02
3     within-hFMT - within-mFMT  6.0429840 1.512897e-09 9.077381e-09
#4 within-healthy - within-noFMT -1.0256912 3.050371e-01 6.100743e-01
#5    within-hFMT - within-noFMT  0.6797357 4.966718e-01 4.966718e-01
#6    within-mFMT - within-noFMT -4.3541217 1.336015e-05 6.680075e-05	
       
       
```

###-------------------------
### Figure 5C: Random Forest analysis

### inputs:
	- Data/metabolomics/Metabolon/sample_metadata_metabolon.csv
	- Data/metabolomics/Metabolon/raw_scaled_metabolon.csv
	- Data/metabolomics/Metabolon/compound_meta_metabolon.csv

```{r}

# random forest:
# adapted from: https://github.com/asmcmill/FMT-Manuscript/blob/main/All%20Metabolites.Rmd

library(tidyverse)
library(janitor)
library(randomForest)
library(pheatmap)
library(writexl)

### for main comparison, col vs cleared:
# read in data:
sample_meta <- read.csv('Data/metabolomics/Metabolon/sample_metadata_metabolon.csv', header = TRUE, sep = ",") %>%
	mutate_if(., is.character, as.factor)
data <- read.csv('Data/metabolomics/Metabolon/raw_scaled_metabolon.csv', header = TRUE, sep = ",")
compound_meta <- read.csv('Data/metabolomics/Metabolon/compound_meta_metabolon.csv', header = TRUE, sep = ",")

# run random forest:
set.seed(333)
data_forest <- data %>%
  column_to_rownames("BIOCHEMICAL") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("metaboID") %>%
  left_join(sample_meta, ., by="metaboID") %>%
  select(status:"1-methyl-5-imidazolelactate") %>% #need to fix this at some point
  clean_names() %>%
  randomForest(status ~ .,data=., proximity=T, importance=TRUE, ntree=1000)
  
#Prepare to Plot
data_forest_plot <- data_forest$importance %>%
  as.data.frame() %>%
  rownames_to_column("biochem") %>%
  left_join(., compound_meta %>%
              select(BIOCHEMICAL,SUPERPATHWAY,SUBPATHWAY, dir) %>%
              mutate(biochem=janitor::make_clean_names(BIOCHEMICAL)),by="biochem") %>%
  filter(MeanDecreaseAccuracy>0) %>%
  arrange(-MeanDecreaseAccuracy) %>%
  group_by(dir) %>%
  mutate(count=1:n()) %>%
  arrange(desc(dir), MeanDecreaseAccuracy) %>%
  mutate(plotnames=str_remove_all(BIOCHEMICAL,"\\*")) %>%
  mutate(plotnames=fct_inorder(plotnames))
  
#Plot:
ggplot(data=data_forest_plot%>%filter(count<=25), aes(y=MeanDecreaseAccuracy,x=plotnames, fill=SUPERPATHWAY))+
  geom_point(shape=21, size=2)+
  scale_fill_manual(values=beast.colors$SUPER.PATHWAY)+
  coord_flip()+
  theme_bw()+
  labs(fill="", y="Mean Decrease Accuracy",x="")+
  theme(axis.text.y=element_text(size=6),
        axis.text.x=element_text(size=6),
        axis.title.y=element_text(size=8),
        legend.position="top",
        legend.text = element_text(size=6),
        legend.key.size = unit(c(0.1,0.2),"in"),
        legend.box.margin=unit(c(0,1.5,-.15,0),"in"),
        plot.margin = unit(c(0,0.25,0,0),"in"))+
  guides(fill=guide_legend(nrow=3))
  
### if you want to get totals of compounds up/down, etc:
data_forest_plot%>%
  group_by(dir)%>%
  summarize(count=n())
#  dir             count
#  <chr>           <int>
#1 Down in cleared   150
#2 Up in cleared     107

data_forest_plot%>%
  filter(count<=25)%>%
  group_by(SUPERPATHWAY,dir)%>%
  summarise(count=n())%>%
  spread(key=dir,value=count)%>%
  ungroup()%>%
  mutate(total=rowSums(.[2:3],na.rm=T))%>%
  arrange(-total)
#    SUPERPATHWAY                        `Down in cleared` `Up in cleared` total
#  <chr>                                           <int>           <int> <dbl>
#1 "Lipid"                                             4              12    16
#2 ""                                                  4               8    12
#3 "Amino Acid"                                       10              NA    10
#4 "Xenobiotics"                                       2               3     5
#5 "Carbohydrate"                                      3              NA     3
#6 "Nucleotide"                                       NA               2     2
#7 "Cofactors and Vitamins"                            1              NA     1
#8 "Partially Characterized Molecules"                 1              NA     1

# to get order for heatmap:
top_data_forest<-as.vector((data_forest_plot %>%
	arrange(dir, MeanDecreaseAccuracy) %>%
	filter(count<=25))$BIOCHEMICAL)

###-----------------###
#### heatmap for compounds:
detach(package:plyr)

data_long <- data %>%
  gather(key="metaboID", value="ScaledImpData", -BIOCHEMICAL) %>%
  mutate(metaboID=factor(metaboID,levels=levels(sample_meta$metaboID)))%>%
  left_join(., sample_meta, by="metaboID")
#### Let's also plot a heatmap to go along with random forest:
data_by_group <-data_long %>%
  select("BIOCHEMICAL", "group", "ScaledImpData")%>%
  mutate(ScaledImpData=log10(ScaledImpData)) %>% # Log transform for plotting
  group_by(group, BIOCHEMICAL)%>%
  summarize(mean=mean(ScaledImpData)) %>% #Mean of groups
  spread(key='group', value="mean")%>%
  left_join(., compound_meta, by="BIOCHEMICAL") %>%
  filter(BIOCHEMICAL %in% top_data_forest) %>%
  mutate(BIOCHEMICAL=factor(BIOCHEMICAL, levels=rev(top_data_forest)))%>%
  arrange(BIOCHEMICAL)%>%
  column_to_rownames("BIOCHEMICAL")

pheatmap(data_by_group[c(1,3,2,4)],
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=1.5)(50),
         border_color="grey20",
         cellheight = 6.5,
         cellwidth = 10,
         fontsize = 6.5,
         #filename="RF_heatmap_sum.pdf",
         width=10,
         height=10
)

```

###-------------------------
### Figure 5D-I: SCFA targeted, fecal and cecal

### inputs:
	- HFm_scfa.fecal.seq.meta.merge.txt
	- HFm_scfa.ba.cecal.seq.meta.merge.txt

```{r}
library(dplyr)
library(tidyverse)
library(tidyr)
library(data.table)
#library(purrr)
library(readr)
library(FSA)
library("patchwork")


### filter data and set pre-FMT samples as preFMT:
# get rid of mFMT_other post-FMT, as there are only a few points (and they are not included in other metabolomic analyses)
d <- read.table("HFm_scfa.seq.meta.merge.txt", header=T, sep="\t") %>%
					filter(type == "fecal" & sampleType =="mouse") %>%
					filter(!(group =="mFMT_other" & day > 11)) %>%
					filter(day > -8) %>%
					mutate(group = case_when(
							day %in% c(-7, 0, 1, 4, 9, 11) ~ "preFMT", 			
							TRUE ~ group),
						comp = factor(group, levels=c("mFMT", "mFMT_other", "hFMT", "noFMT", "preFMT"))
						)					
# d %>% group_by(group) %>% count() %>% print(n=34)

# define colors
dcol3_preFMT=c("hFMT" = "#EACB2B", "mFMT" ="#56B29E", "preFMT" = "grey50", "noFMT" = "#F5191C")
beast.colors <- read.table("Beast.colors.R")
# now in beast.colors as dcol3_preFMT

# if graphing 1 SCFA individually:
p1 <- d %>% ggplot(aes(x=day, y=acetate, color=group, fill=group)) +
	scale_color_manual(values=beast.colors$dcol3_preFMT) + 
	scale_fill_manual(values=beast.colors$dcol3_preFMT) + 
	geom_vline(xintercept = 0, linetype = "solid", color = "red", alpha = 0.5, size = 0.5) +
	geom_vline(xintercept = 11, linetype = "solid", color = "springgreen4", alpha = 0.5, size = 0.5) +
	annotate("text", x = 0, y = 250, label = "C. diff", size = 3, color = "red", alpha=0.5) +
	annotate("text", x = 11, y = 250, label = "FMT", size = 3, color = "springgreen4") +
	#geom_ribbon(aes(ymin = Q1, ymax = Q3, group=group), alpha=0.2, color=NA) +
  	geom_jitter(width=0.25, alpha=0.4, size = 1, na.rm=TRUE) + #if you wanted points
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="line", linewidth=1,
			fun.args = list(conf.int=0.50)) +
  	stat_summary(fun.data=median_hilow, size=0.85, na.rm= TRUE, geom="point",
			fun.args = list(conf.int=0.50)) +
  	scale_y_continuous(limits=c(0,max(d$acetate))) +
	scale_x_continuous(limits=c(-7, 42), breaks=c(-7, 0, 1, 4, 9, 11, 13, 16, 20, 28, 32, 36, 41)) +
	labs(x=NULL, y="acetate (nmol)") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	theme_classic()
	
p2 <- d %>% ggplot(aes(x=day, y=propionate, color=group, fill=group)) +
	scale_color_manual(values=beast.colors$dcol3_preFMT) + 
	scale_fill_manual(values=beast.colors$dcol3_preFMT) + 
	geom_vline(xintercept = 0, linetype = "solid", color = "red", alpha = 0.5, size = 0.5) +
	geom_vline(xintercept = 11, linetype = "solid", color = "springgreen4", alpha = 0.5, size = 0.5) +
	annotate("text", x = 0, y = 250, label = "C. diff", size = 3, color = "red", alpha=0.5) +
	annotate("text", x = 11, y = 250, label = "FMT", size = 3, color = "springgreen4") +
	#geom_ribbon(aes(ymin = Q1, ymax = Q3, group=group), alpha=0.2, color=NA) +
  	geom_jitter(width=0.25, alpha=0.4, size = 1, na.rm=TRUE) + #if you wanted points
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="line", linewidth=1,
			fun.args = list(conf.int=0.50)) +
  	stat_summary(fun.data=median_hilow, size=0.85, na.rm= TRUE, geom="point",
			fun.args = list(conf.int=0.50)) +
  	scale_y_continuous(limits=c(0, max(d$propionate, na.rm=T))) +
	scale_x_continuous(limits=c(-7, 42), breaks=c(-7, 0, 1, 4, 9, 11, 13, 16, 20, 28, 32, 36, 41)) +
	labs(x=NULL, y="propionate (nmol)") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	theme_classic()

p3 <- d %>% ggplot(aes(x=day, y=butyrate, color=group, fill=group)) +
	scale_color_manual(values=beast.colors$dcol3_preFMT) + 
	scale_fill_manual(values=beast.colors$dcol3_preFMT) + 
	geom_vline(xintercept = 0, linetype = "solid", color = "red", alpha = 0.5, size = 0.5) +
	geom_vline(xintercept = 11, linetype = "solid", color = "springgreen4", alpha = 0.5, size = 0.5) +
	annotate("text", x = 0, y = 250, label = "C. diff", size = 3, color = "red", alpha=0.5) +
	annotate("text", x = 11, y = 250, label = "FMT", size = 3, color = "springgreen4") +
	#geom_ribbon(aes(ymin = Q1, ymax = Q3, group=group), alpha=0.2, color=NA) +
  	geom_jitter(width=0.25, alpha=0.4, size = 1, na.rm=TRUE) + #if you wanted points
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="line", linewidth=1,
			fun.args = list(conf.int=0.50)) +
  	stat_summary(fun.data=median_hilow, size=0.85, na.rm= TRUE, geom="point",
			fun.args = list(conf.int=0.50)) +
  	scale_y_continuous(limits=c(0, max(d$butyrate, na.rm=T))) +
	scale_x_continuous(limits=c(-7, 42), breaks=c(-7, 0, 1, 4, 9, 11, 13, 16, 20, 28, 32, 36, 41)) +
	labs(x=NULL, y="butyrate (nmol)") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	theme_classic()

# patch them together using patchwork:
p1 / p2 / p3

### stats:
subdf <- d %>% 
	mutate(status = case_when(
		group == "hFMT" ~ "hFMT", 
		group == "noFMT" ~ "noFMT",		
		group %in% c("mFMT", "mFMT_other") ~ "mFMT",
		TRUE ~ NA_character_),
		status = factor(status, levels=c("mFMT", "hFMT", "noFMT")),
		day = factor(day, levels=as.character(unique(day)))
		)
factor_levels <- unique(subdf$day)

# test for post-FMT days 19, 32, 41 (although noFMT unavailable for earlier days)
#butyrate
subdf %>% subset(., day == 13 & !is.na(butyrate) ) %>%
	kruskal.test(butyrate ~ group, data = .)		
subdf %>% subset(., day == 16 & !is.na(butyrate) ) %>%
	kruskal.test(butyrate ~ group, data = .)		
subdf %>% subset(., day == 18 & !is.na(butyrate) ) %>%
	kruskal.test(butyrate ~ group, data = .)	
subdf %>% subset(., day == 19 & !is.na(butyrate) ) %>%
	kruskal.test(butyrate ~ group, data = .)	
subdf %>% subset(., day == 32 & !is.na(butyrate) ) %>%
	kruskal.test(butyrate ~ group, data = .)	
subdf %>% subset(., day == 41 & !is.na(butyrate) ) %>%
	kruskal.test(butyrate ~ group, data = .)	

#acetate
subdf %>% subset(., day == 13 & !is.na(acetate) ) %>%
	kruskal.test(acetate ~ status, data = .)		
subdf %>% subset(., day == 16 & !is.na(acetate) ) %>%
	kruskal.test(acetate ~ status, data = .)		
subdf %>% subset(., day == 18 & !is.na(acetate) ) %>%
	kruskal.test(acetate ~ status, data = .)	
subdf %>% subset(., day == 19 & !is.na(acetate) ) %>%
	kruskal.test(acetate ~ status, data = .)	
subdf %>% subset(., day == 32 & !is.na(acetate) ) %>%
	kruskal.test(acetate ~ status, data = .)
subdf %>% subset(., day == 32 & !is.na(acetate) ) %>%
	kruskal.test(acetate ~ status, data = .)

#propionate
subdf %>% subset(., day == 13 & !is.na(propionate) ) %>%
	kruskal.test(propionate ~ status, data = .)		
subdf %>% subset(., day == 16 & !is.na(propionate) ) %>%
	kruskal.test(propionate ~ status, data = .)		
subdf %>% subset(., day == 18 & !is.na(propionate) ) %>%
	kruskal.test(propionate ~ status, data = .)	
subdf %>% subset(., day == 19 & !is.na(propionate) ) %>%
	kruskal.test(propionate ~ status, data = .)	
subdf %>% subset(., day == 32 & !is.na(propionate) ) %>%
	kruskal.test(propionate ~ status, data = .)
subdf %>% subset(., day == 41 & !is.na(propionate) ) %>%
	kruskal.test(propionate ~ status, data = .)	

### let's facet the cecal measurements here too:
####------
### let's also patch the cecal figures with this!

d <- read_tsv("HFm_scfa.ba.cecal.seq.meta.merge.txt") %>%
		filter(!group=="preFMT") %>%
		unite("group_c", day,group, sep ="_", remove=FALSE) %>%
		mutate(group_c = factor(group_c, levels=c('-7_healthy', '21_mFMT', '21_hFMT', '21_noFMT', 
										'42_mFMT', '42_hFMT', '42_noFMT')))
# what do we have?									
d %>%
  dplyr::group_by(group) %>%
  summarise(n_distinct(mouse, day))

# define colors
# specify colors for graph: also in Beast.colors.R
dcol_d21_healthy=c("hFMT" = "#EACB2B", "mFMT" ="#56B29E", "healthy" = "#9FC095", "noFMT" = "#F5191C"),
		
# graph:
d %>% ggplot(aes(x=group_c, y=Acetate_nmol, color=group_c)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_d21_healthy) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="Acetate (nmol)") +
	scale_x_discrete(breaks=c('-7_healthy', '21_mFMT', '21_hFMT', '21_noFMT', '42_mFMT', '42_hFMT', '42_noFMT'), 
			labels=c('healthy', 'mFMT', 'hFMT', 'noFMT', 'mFMT', 'hFMT', 'noFMT')) +
	scale_y_continuous(limits=c(0,max(d$Acetate_nmol, na.rm=T))) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 1.8, xend = 4.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 4.5, xend = 7.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	annotate("text", x = 2, y = 0.5, label = "day 21", size = 3, color = "grey20") +
	annotate("text", x = 5.5, y = 0.5, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey20")

# since healthy is so much higher, let's graph that separately:
# each graph can be pasted together for all SCFAs
p4<- d %>% filter(group == "healthy") %>% 
	ggplot(aes(x=group_c, y=Acetate_nmol, color=group_c)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_d21_healthy) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="acetate (nmol)") +
	scale_x_discrete(breaks=c('-7_healthy'), 
			labels=c('healthy')) +
	scale_y_continuous(limits=c(0,max(d$Acetate_nmol+20, na.rm=T))) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p5<- d %>% filter(group == "healthy") %>% 
	ggplot(aes(x=group_c, y=Propionate_nmol, color=group_c)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_d21_healthy) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="propionate (nmol)") +
	scale_x_discrete(breaks=c('-7_healthy'), 
			labels=c('healthy')) +
	scale_y_continuous(limits=c(0,max(d$Propionate_nmol+10, na.rm=T))) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p6<- d %>% filter(group == "healthy") %>% 
	ggplot(aes(x=group_c, y=Butyrate_nmol, color=group_c)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_d21_healthy) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="butyrate (nmol)") +
	scale_x_discrete(breaks=c('-7_healthy'), 
			labels=c('healthy')) +
	scale_y_continuous(limits=c(0,max(d$Butyrate_nmol+10, na.rm=T))) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# patch them together:
p4 / p5 / p6

### finally, cecal for treated groups, then patched together:
p7 <- d %>% filter(!group == "healthy") %>%
	ggplot(aes(x=group_c, y=Acetate_nmol, color=group_c)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_d21_healthy) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="acetate (nmol)") +
	scale_x_discrete(breaks=c('21_mFMT', '21_hFMT', '21_noFMT', '42_mFMT', '42_hFMT', '42_noFMT'), 
			labels=c('mFMT', 'hFMT', 'noFMT', 'mFMT', 'hFMT', 'noFMT')) +
	scale_y_continuous(limits=c(0,100)) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 0.8, xend = 3.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 3.5, xend = 6.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	annotate("text", x = 2, y = 0, label = "day 21", size = 3, color = "grey20") +
	annotate("text", x = 5, y = 0, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey20")

p8 <- d %>% filter(!group == "healthy") %>%
	ggplot(aes(x=group_c, y=Propionate_nmol, color=group_c)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_d21_healthy) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="propionate (nmol)") +
	scale_x_discrete(breaks=c('21_mFMT', '21_hFMT', '21_noFMT', '42_mFMT', '42_hFMT', '42_noFMT'), 
			labels=c('mFMT', 'hFMT', 'noFMT', 'mFMT', 'hFMT', 'noFMT')) +
	scale_y_continuous(limits=c(0,25)) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 0.8, xend = 3.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 3.5, xend = 6.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	annotate("text", x = 2, y = 0, label = "day 21", size = 3, color = "grey20") +
	annotate("text", x = 5, y = 0, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey20")

p9 <- d %>% filter(!group == "healthy") %>%
	ggplot(aes(x=group_c, y=Butyrate_nmol, color=group_c)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_d21_healthy) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="butyrate (nmol)") +
	scale_x_discrete(breaks=c('21_mFMT', '21_hFMT', '21_noFMT', '42_mFMT', '42_hFMT', '42_noFMT'), 
			labels=c('mFMT', 'hFMT', 'noFMT', 'mFMT', 'hFMT', 'noFMT')) +
	scale_y_continuous(limits=c(0,30)) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 0.8, xend = 3.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 3.5, xend = 6.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	annotate("text", x = 2, y = 0, label = "day 21", size = 3, color = "grey20") +
	annotate("text", x = 5, y = 0, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey20")

### let's try to patch ALL of these together...
# each SCFA side by side, 3 rows of each SCFA
(p1 | p4 | p7) / (p2 | p5 | p8) / (p3 | p6 | p9)

# with some editing of spacing:
layout <- "
AAAAAABCCC
DDDDDDEFFF
GGGGGGHIII
"
p1 + p4 + p7 + p2 + p5 + p8 + p3 + p6 + p9 + plot_layout(design = layout)

### stats:
d %>% subset(., group_c %in% c("21_mFMT", "21_hFMT", "21_noFMT")) %>%
	dunnTest(Acetate_nmol ~ group, data = .)
#    Comparison          Z   P.unadj     P.adj
#1  hFMT - mFMT -0.1961161 0.8445193 0.8445193
#2 hFMT - noFMT  1.0786387 0.2807488 0.5614976
#3 mFMT - noFMT  1.2747549 0.2023960 0.6071880
d %>% subset(., group_c %in% c("21_mFMT", "21_hFMT", "21_noFMT")) %>%
	dunnTest(Propionate_nmol ~ group, data = .)
#    Comparison        Z    P.unadj      P.adj
#1  hFMT - mFMT 1.274755 0.20239602 0.40479203
#2 hFMT - noFMT 2.549510 0.01078745 0.03236235
#3 mFMT - noFMT 1.274755 0.20239602 0.20239602
d %>% subset(., group_c %in% c("21_mFMT", "21_hFMT", "21_noFMT")) %>%
	dunnTest(Butyrate_nmol ~ group, data = .)
#1  hFMT - mFMT -0.5883484 0.55629846 0.55629846
#2 hFMT - noFMT  2.0592194 0.03947322 0.07894645
#3 mFMT - noFMT  2.6475678 0.00810731 0.02432193
d %>% subset(., group_c %in% c("42_mFMT", "42_hFMT", "42_noFMT")) %>%
	dunnTest(Acetate_nmol ~ group, data = .)
#    Comparison          Z    P.unadj      P.adj
#1  hFMT - mFMT -1.7764696 0.07565555 0.15131110
#2 hFMT - noFMT  0.8528029 0.39376863 0.39376863
#3 mFMT - noFMT  2.5660116 0.01028754 0.03086262
d %>% subset(., group_c %in% c("42_mFMT", "42_hFMT", "42_noFMT")) %>%
	dunnTest(Propionate_nmol ~ group, data = .)
#    Comparison         Z    P.unadj      P.adj
#1  hFMT - mFMT 0.2302831 0.81787180 0.81787180
#2 hFMT - noFMT 2.4518082 0.01421404 0.04264212
#3 mFMT - noFMT 2.0396503 0.04138517 0.08277035
d %>% subset(., group_c %in% c("42_mFMT", "42_hFMT", "42_noFMT")) %>%
	dunnTest(Butyrate_nmol ~ group, data = .)
#    Comparison         Z    P.unadj      P.adj
#1  hFMT - mFMT -2.171241 0.02991299 0.08973897
#2 hFMT - noFMT  0.000000 1.00000000 1.00000000
#3 mFMT - noFMT  2.171241 0.02991299 0.05982598


```

###-------------------------
### Figure 5J-N: BA totals and by specific BAs
	- Primary bile acids CA and TCA in main figure; mouse primary BAs aMCA and bMCA in Supplemental
	- Secondary bile acids DCA and wMCA in main

### inputs:
	- HFm_scfa.ba.cecal.seq.meta.merge.txt

```{r}
library(dplyr)
library(tidyverse)
library(tidyr)
library(data.table)
#library(purrr)
library(readr)
library(FSA)
library("patchwork")

### BA totals (primary vs. secondary):

library(dplyr)
library(tidyverse)
library(tidyr)
library(data.table)
#library(purrr)
library(readr)
library(FSA)
library(dplyr)

### plot 1 vs 2 BAs:
primary <- c("CA", "CDCA", "aMCA", "bMCA", "TCA", "TCDCA", "TaMCA_TbMCA", "HCA", "GCA", "GCDCA")
secondary <- c("DCA", "wMCA", "THDCA_TUDCA", "TLCA", "GUDCA", "TDCA", "GHDCA", "GDCA")

d <- read_tsv("HFm_scfa.ba.cecal.seq.meta.merge.txt") %>%
		mutate(group = factor(group, levels=c("healthy", "mFMT", "hFMT", "noFMT"))) %>%
			# let's get rid of undetected BAs
		select(-c(GLCA, THCA, GHCA, HDCA_UDCA)) %>%
			# calculate primary vs secondary bile acids totals:
		rowwise() %>%
		mutate(ba1 = sum(CA, CDCA, aMCA, bMCA, TCA, TCDCA, TaMCA_TbMCA, HCA, GCA, GCDCA)) %>%
		mutate(ba2 = sum(DCA, wMCA, THDCA_TUDCA, TLCA, GUDCA, TDCA, GHDCA, GDCA))

# graph means by group:
mean.bastotals<-ddply(d, c("group"), colwise(mean, is.numeric))

col.totals<-c("mistyrose", "darkturquoise")
barh<-as.data.frame(mean.bastotals[, c("ba1", "ba2")])
rownames(barh)<-mean.bastotals[,1]
bar_g<-as.matrix(t(barh))
par(mar=c(4,4,2,4))
par(xpd=T)
barplot(bar_g, main="", xlim=c(0,4), col=col.totals, ylim=c(0,20), axisnames=FALSE, cex.axis=0.8, cex.lab=0.8)
xnames=as.character(mean.bastotals[,1])
text(x =  seq(1,4,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = xnames, xpd = TRUE, cex=0.8)
mtext(side=2, text="bile acid totals \n(ug / 100 mg wet weight)", line=2, cex=0.8)
legend(3,15,legend=c("primary", "secondary"),col=col.totals,fill=col.totals,cex=0.5, title="Bile acids")


### specific BAs:

## Main figures: DCA, wMCA on top; CA, TCA, aMCA, bMCA on bottom
d <- read_tsv("HFm_scfa.ba.cecal.seq.meta.merge.txt") %>%
		unite("group_c", day,group, sep ="_", remove=FALSE) %>%
		mutate(group_c = factor(group_c, levels=c('-7_healthy', '0_preFMT', '1_preFMT', '4_preFMT', '9_preFMT',
										'21_mFMT', '21_hFMT', '21_noFMT', 
										'42_mFMT', '42_hFMT', '42_noFMT')))
										
# set colors and graphing parameters
cecal_breaks <- levels(d$group_c)
#cecal_labels <- c('day -7', "day 0", "day 1", "day 4", "day 9", 'mFMT', 'hFMT', 'noFMT', 'mFMT', 'hFMT', 'noFMT')
cecal_labels <- c('-7', "0", "1", "4", "9", '', '21', '', '', '42', '')
dcol_ba=c("hFMT" = "#EACB2B", "mFMT" ="#56B29E", "healthy" = "#9FC095", "noFMT" = "#F5191C", "preFMT" = "grey50")
	# in beast.colors$dcol_ba

# graph each BA as a separate entity, then we will patch them together
b1.1 <- d %>% ggplot(aes(x=group_c, y=DCA, color=group)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_ba) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="DCA (nmol)") +
	scale_x_discrete(breaks=cecal_breaks, 
			labels=cecal_labels) +
	scale_y_continuous(limits=c(0,max(d$DCA, na.rm=T))) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 5.6, xend = 8.4, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 8.5, xend = 11.3, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	#annotate("text", x = 7, y = 0, label = "day 21", size = 3, color = "grey20") +
	#annotate("text", x = 10, y = 0, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey20") +
	geom_vline(xintercept = 8.5, linetype = "dashed", color = "grey20")

b1.2 <- d %>% ggplot(aes(x=group_c, y=wMCA, color=group)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_ba) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="wMCA (nmol)") +
	scale_x_discrete(breaks=cecal_breaks, 
			labels=cecal_labels) +
	scale_y_continuous(limits=c(0,max(d$wMCA, na.rm=T))) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 5.6, xend = 8.4, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 8.5, xend = 11.3, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	#annotate("text", x = 7, y = 0, label = "day 21", size = 3, color = "grey20") +
	#annotate("text", x = 10, y = 0, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey20") +
	geom_vline(xintercept = 8.5, linetype = "dashed", color = "grey20")

b2.1 <- d %>% ggplot(aes(x=group_c, y=CA, color=group)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_ba) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="CA (nmol)") +
	scale_x_discrete(breaks=cecal_breaks, 
			labels=cecal_labels) +
	scale_y_continuous(limits=c(0,max(d$CA, na.rm=T))) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 5.6, xend = 8.4, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 8.5, xend = 11.3, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	#annotate("text", x = 7, y = 0, label = "day 21", size = 3, color = "grey20") +
	#annotate("text", x = 10, y = 0, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey20") +
	geom_vline(xintercept = 8.5, linetype = "dashed", color = "grey20")

b2.2 <- d %>% ggplot(aes(x=group_c, y=TCA, color=group)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_ba) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="TCA (nmol)") +
	scale_x_discrete(breaks=cecal_breaks, 
			labels=cecal_labels) +
	scale_y_continuous(limits=c(0,max(d$TCA, na.rm=T))) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 5.6, xend = 8.4, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 8.5, xend = 11.3, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	#annotate("text", x = 7, y = 0, label = "day 21", size = 3, color = "grey20") +
	#annotate("text", x = 10, y = 0, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey20") +
	geom_vline(xintercept = 8.5, linetype = "dashed", color = "grey20")
	
b2.3 <- d %>% ggplot(aes(x=group_c, y=aMCA, color=group)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_ba) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width=0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="aMCA (nmol)") +
	scale_x_discrete(breaks=cecal_breaks, 
			labels=cecal_labels) +
	scale_y_continuous(limits=c(0,max(d$aMCA, na.rm=T))) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 5.6, xend = 8.4, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 8.5, xend = 11.3, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	#annotate("text", x = 7, y = 0, label = "day 21", size = 3, color = "grey20") +
	#annotate("text", x = 10, y = 0, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey20") +
	geom_vline(xintercept = 8.5, linetype = "dashed", color = "grey20")

b2.4 <- d %>% ggplot(aes(x=group_c, y=bMCA, color=group)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_ba) + 
	stat_summary(fun = median, geom = "crossbar", size = 0.25, width=0.65 ) +
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="errorbar", width = 0.65,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="bMCA (nmol)") +
	scale_x_discrete(breaks=cecal_breaks, 
			labels=cecal_labels) +
	scale_y_continuous(limits=c(0,max(d$bMCA, na.rm=T))) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 5.6, xend = 8.4, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 8.5, xend = 11.3, y = 0.2, yend = 0.2), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	#annotate("text", x = 7, y = 0, label = "day 21", size = 3, color = "grey20") +
	#annotate("text", x = 10, y = 0, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 5.5, linetype = "dashed", color = "grey20") +
	geom_vline(xintercept = 8.5, linetype = "dashed", color = "grey20")

# patch together:
layout <- "
AAABBBCCCDDD
EEEFFFGGGHHH
"
plot_spacer() + b1.1 + b1.2 + plot_spacer() + b2.1 + b2.2 + b2.3 + b2.4 + plot_layout(design = layout)

### stats:
# DCA:
d %>% subset(., group_c %in% c("21_mFMT", "21_hFMT", "21_noFMT")) %>%
	dunnTest(DCA ~ group, data = .)
d %>% subset(., group_c %in% c("42_mFMT", "42_hFMT", "42_noFMT")) %>%
	dunnTest(DCA ~ group, data = .)
# wMCA:
d %>% subset(., group_c %in% c("21_mFMT", "21_hFMT", "21_noFMT")) %>%
	dunnTest(wMCA ~ group, data = .)
d %>% subset(., group_c %in% c("42_mFMT", "42_hFMT", "42_noFMT")) %>%
	dunnTest(wMCA ~ group, data = .)
# CA:
d %>% subset(., group_c %in% c("21_mFMT", "21_hFMT", "21_noFMT")) %>%
	dunnTest(CA ~ group, data = .)
d %>% subset(., group_c %in% c("42_mFMT", "42_hFMT", "42_noFMT")) %>%
	dunnTest(CA ~ group, data = .)
# TCA:
d %>% subset(., group_c %in% c("21_mFMT", "21_hFMT", "21_noFMT")) %>%
	dunnTest(TCA ~ group, data = .)
d %>% subset(., group_c %in% c("42_mFMT", "42_hFMT", "42_noFMT")) %>%
	dunnTest(TCA ~ group, data = .)
# aMCA:
d %>% subset(., group_c %in% c("21_mFMT", "21_hFMT", "21_noFMT")) %>%
	dunnTest(aMCA ~ group, data = .)
d %>% subset(., group_c %in% c("42_mFMT", "42_hFMT", "42_noFMT")) %>%
	dunnTest(aMCA ~ group, data = .)
# bMCA:
d %>% subset(., group_c %in% c("21_mFMT", "21_hFMT", "21_noFMT")) %>%
	dunnTest(bMCA ~ group, data = .)
d %>% subset(., group_c %in% c("42_mFMT", "42_hFMT", "42_noFMT")) %>%
	dunnTest(bMCA ~ group, data = .)





```