#### Figure 2
### 7.24.24
### Anna M. Seekatz

### Figures:
	- Fig. 2A: NMDS based on OTUs (fecal and cecal)
	- Fig. 2B: pairwise Bray-Curtis (cecal)
	- Fig. 2C: taxonimic (genus) bargraphs
	- Fig. 2D: diversity over time (fecal)
	- Fig. 2E: diversity, cecal
	- Fig. 2F: sign. different OTUs (maaslin2)
	
### Input files: Data/16S
	- HFm_summary.txt
	- HFm.otus.shared
	- HFm_summary.nmds.new.groups.txt
	- HFm_BC.2936_wide_avg.dist
	- HFm_shared.df.all.groups.input-for-dist.txt
	- allhf_runs/allhf_genfrac2p.meta.txt
	- 16S_maaslyn/16S_clear_v_colonized/significant_results.tsv
	- taxonomy_meta.txt
	- Beast.colors.R (read this into R to use with figures)


###-------------------------
### Figure 2A: NMDS (based on OTUs)

### inputs:
- Data/16S/HFm_summary.txt
- Data/16S/HFm.otus.shared

```{r}
###------------------
# use vegan (bray-curtis) to calculate and make a PCoA:
# adapted from: https://archetypalecology.wordpress.com/2018/02/21/permutational-multivariate-analysis-of-variance-permanova-in-r-preliminary/
library(tidyr)
library(reshape2)
library(tidyverse)
library(vegan)

# read in otu table; merge with selected metadata:
# example: only cecal, for the 4 groups:
# read in what samples we are using, plus get rid of samples < 2500 sequences
d <- read.table("Data/16S/HFm_summary.txt", header=T, sep="\t") %>%
		.[!is.na(.$simpson), ] %>%
		.[.$nseqs > 2500, ] %>%
		.[.$group %in% c('mFMT', 'mFMT_other', 'hFMT', 'noFMT') & .$sampleType %in% c('mouse') & .$type %in% c('cecal'), ]
otus <- read.table("Data/HFm.otus.shared", header=T, sep="\t", row.names=2) %>%
		select(., -c(label, numOtus)) %>%
		merge(d, ., by.x="seqID", by.y='row.names') %>%
		select(Otu00001:length(.))
#otus.dist <- vegdist(otus, method="bray")
otus.dist <- avgdist(otus, sample=2936)
otus.div <- adonis2(otus.dist ~ group, data = d, permutations = 999, method="bray")

dispersion <- betadisper(otus.dist, group=d$group, pairwise=TRUE)
anova(dispersion)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse

# overall tests
adonis2(dune ~ Management*A1, data = dune.env, permutations = 999, method="bray" by = NULL)


###---
# done for ALL samples of study (cecal, fecal, and inputs)
# let's read in data first
otus <- read_tsv("Data/16S/HFm.otus.shared") #%>%	# can usually combine this with the rest, but it takes too long to read in
d <- read_tsv("Data/16S/HFm_summary.txt") %>%
		filter(!is.na(simpson)) %>%
		filter(nseqs > 2500) %>%
		filter(group %in% c('mFMT', 'mFMT_other', 'hFMT', 'noFMT')) %>%		# for now, only choose these groups
		#filter(sampleType %in% c("mouse") & type %in% c('fecal'))	%>%		# if you wanted to separate out fecal from cecal
		#filter(day > 11)		# if you only wanted post-FMT days

# read in OTU table, filter out samples with low counts, etc
#days_wanted <- c(21, 42)	# if you wanted specific days 
shared_df <- otus %>%	
	select(Group, starts_with("Otu")) %>%
	pivot_longer(-Group) %>%
	#separate(Group, into=c("day", "cage", "animal"), sep="\t", remove=FALSE, convert=TRUE) %>% # can only use if the seqIDs mean this! I have too many mistakes in seqID names
	#filter(day %in% days_wanted) %>%	# filters by day; but we need to filter by meta
	group_by(Group) %>%
	#summarize(N = sum(value)) %>%	# to get number of total sequences
	#arrange(N) %>% print(n=20)	# list the sequence N
	mutate(N = sum(value)) %>%	# create a new column for N
	ungroup() %>%	# ungroup from the counts
	filter(N >=min(N)) %>% 	#filter out any sequences below your threhold; note: already did this! min(N)=2936
	select(-N) %>% #get rid of the N column (don't need it)
	pivot_wider(names_from="name", values_from="value", values_fill=0) %>%
	#inner_join(., d, by=c('Group', 'seqID')) # doesn't work
	merge(d, ., by.x="seqID", by.y="Group") %>%
	as_tibble()
# in case you have to redo this, and the order of something gets messed up:
#write_tsv(shared_df, "Data/HFm_shared.df.all.groups.input-for-dist.txt")

# separate out desired metadata:
meta <- shared_df %>%
	select(seqID, group, type, FMT_input, GROUP, TYPE, mouse, day)	#%>% 	#whatever metadata you want
	#mutate(period = if_else(day < 10, "early", "late")	# giving a period based on time constraints
	#		sex = if_else(str_detect(animal, "F"), "female", "male"))	# assigning a sex if you can based on a character in something
			
# to summarize data:
meta %>% count(group)
	
# calculate distance matrix:
# note: calculated matrix for ALL samples was used for NMDS graph
# also calculated fecal (post-FMT samples only) and cecal dataframes for adonis testing
all.otu.dist <- shared_df %>%
	#filter(TYPE == 'fecal' & day > 11) %>%	#if you want fecal samples only, post-FMT (saved as fecal.otu.dist)
	#filter(TYPE == 'cecal' & day > 11) %>%	#if you want cecal samples only, post-FMT (saved as cecal.otu.dist)
	select(seqID, starts_with("Otu")) %>%
	column_to_rownames("seqID") %>%
	select_if(colSums(.) !=0) %>% 		#eliminates some time!
	avgdist(sample=2936)
# this takes a long time, especially for 933 samples!! took ~2 hrs recommend saving:
# calculated for fecal.otu.dist (saved as "Data/HFm_fecal.only_BC.2936_wide_avg.dist") and cecal.otu.dist (saved as "Data/HFm_cecal.only_BC.2936_wide_avg.dist")
HFm.allgroups.dist <- all.otu.dist %>%
	as.matrix() %>%
	as_tibble(rownames = "samples") %>%
	#write_tsv(., "Data/HFm_BC.2936_wide_avg.dist") %>%
	pivot_longer(-samples) %>%
	filter(samples < name) #%>%
	#write_tsv(., "Data/HFm_BC.2936_long_avg.dist")
# if reading in the above wide dist object (starting from here):
#all.otu.dist <- read_tsv("Data/HFm_BC.2936_wide_avg.dist") %>%
#	column_to_rownames("samples") %>%
#	as.dist()

# run NMDS on your distance matrix
# all.otu.dist <- read_tsv("Data/HFm_BC.2936_wide_avg.dist") %>%
#	column_to_rownames("samples") %>%
#	as.dist()
set.seed(13)
all_nmds <- metaMDS(all.otu.dist) %>%
	scores() %>%
	as_tibble(rownames="seqID")

# rejoin with metadata for graphing:
meta_nmds <- inner_join(meta, all_nmds, by="seqID")
#write_tsv(meta_nmds, "Data/HFm_bray.nmds.meta.txt")	#all 933 samples

# if starting here:
#meta_nmds <- read_tsv("Data/HFm_bray.nmds.meta.txt")
all.otu.dist <- read_tsv("Data/16S/HFm_BC.2936_wide_avg.dist") %>%
	select(-samples)

# centroid calculation:
# if group_by does not work, detach plyr
#detach(package:plyr)
centroids <- meta_nmds %>%
	dplyr::group_by(GROUP) %>%
	summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))
	
# adonis test:
test <- adonis2(as.dist(all.otu.dist)~meta_nmds$group)
str(test)
p_value <- test$`Pr(>F)`[1]
adonis2(as.dist(all.otu.dist)~meta_nmds$group*meta_nmds$GROUP)		#if using + instead of *, will calculate separate interactions
#adonis2(all.otu.dist~meta$group, strata=meta$mouse)		# if you were stratifying by animal; but, we do not have early/late examples here

# dispersion:
bd <- betadisper(as.dist(all.otu.dist), meta_nmds$group)
anova(bd)
permutest(bd, pairwise=TRUE)


### let's also calculate adonis 2 across the four groups, eliminating pre-FMT and inputs; separately for fecal vs cecal
# FOR FECAL
fecal.otu.dist <- read_tsv("Data/16S/HFm_fecal.only_BC.2936_wide_avg.dist") %>%
	select(-samples)
test.fecal <- meta %>%
			filter(TYPE == 'fecal' & day > 11)
test <- adonis2(as.dist(fecal.otu.dist)~test.fecal$group)

# FOR CECAL
cecal.otu.dist <- read_tsv("Data/16S/HFm_cecal.only_BC.2936_wide_avg.dist") %>%
	select(-samples)
test.cecal <- meta %>%
			filter(TYPE == 'cecal' & day > 11)
test <- adonis2(as.dist(cecal.otu.dist)~test.cecal$group)


### graph using same color stipulations as above:
###---           
### using ALL data:
### this is what was included as Figure 2A
data <- meta_nmds					
data <- data[order(data$GROUP), ]
data$GROUP <- as.factor(data$GROUP)
data$GROUP <- as.factor(factor(data$GROUP, levels=c('pre-abx', 'post-cef', 'post-cdi', 'post-vanco', 
								'D2_A', 'D4_A', 'D5_A', 'R4_F', 'AMS001', 'AMS005', 
								'none',
								'rWT', 'yRag', 'yWT_spore', 'yWT')))
nmds.col <- c('#3B99B1', "#474747", "#A8A8A8", "#E2E2E2", 
						'#EACB2B', '#E79812', '#E8A91B', '#ED5300', '#E87700', '#E8A419', 
						'#F5191C',
						'#8BBD94', '#A6C293', '#C1C88C', '#56B29E')
# if you want to see what colors are what
#barplot(1:length(nmds.col), col=nmds.col, las=2, names=paste(levels(data$GROUP), nmds.col, sep="\n"), cex.names=0.8)

data$TYPE <- as.factor(data$TYPE)
dpch <- c(17, 19, 8) #cecal, fecal, FMTinput; colored circles
dpch2 <- c(24, 21, 8) #cecal, fecal, FMTinput
dcex <- c(1,1,2)

par(mar = c(5, 4, 4, 8), xpd = TRUE)
plot(data$NMDS1, data$NMDS2, xlim=c(min(data$NMDS1, na.rm=T), max(data$NMDS1, na.rm=T) ), ylim=c(min(data$NMDS2, na.rm=TRUE), max(data$NMDS2, na.rm=TRUE)), 
			ylab="NMDS2", xlab="NMDS1", 
			col=adjustcolor(nmds.col[data$GROUP], alpha=0.6), 
			cex=dcex[data$TYPE], pch=dpch[data$TYPE])
points(data$NMDS1[data$TYPE %in% c("cecal")], data$NMDS2[data$TYPE %in% c("cecal")], col='black', 
	pch=dpch2[data$TYPE[data$TYPE %in% c("cecal")]], cex=dcex[data$TYPE[data$TYPE %in% c("cecal")]])
legend("topright", inset = c(- 0.3, 0),                   # Create legend outside of plot
       legend = as.character(levels(data$GROUP)),
       col=adjustcolor(nmds.col, alpha=0.8),
       cex=1, pt.cex=c(1), pch=c(15))
legend("bottomright", inset = c(- 0.165, 0), 
		legend = as.character(levels(data$TYPE)), 
		col=adjustcolor(c("black"), alpha=0.8), 
		cex=1, pch=dpch2, pt.cex=1)
		


```

###-------------------------
### Figure 2B: pairwise bray-curtis distances

### inputs:
- Data/16S/HFm_summary.nmds.new.groups.txt
- Data/16S/HFm_BC.2936_long_avg.dist

```{r}
### Calculating dist differences from the above distance matrix:
# for this, will use the distance matrix generated from ALL samples
library(tidyr)
library(reshape2)
library(tidyverse)
library(vegan)

### now, onto the actual data:
# read in previously generated distance matrix and matching metadata:
all_dist <- read_tsv("Data/16S/HFm_BC.2936_long_avg.dist")	#note: this was already transformed above
d <- read_tsv("Data/16S/HFm_summary.nmds.new.groups.txt") %>%
		filter(!is.na(simpson)) %>%
		filter(nseqs > 2500) %>%
		filter(group %in% c('mFMT', 'mFMT_other', 'hFMT', 'noFMT'))

# select the metadata you want to combine together, renaming columns into a or b
meta <- d %>%
	select(seqID, group, type, GROUP, TYPE, mouse, day) %>%
	mutate(name2 = seqID, group2 = group, type2 = type, GROUP2 = GROUP, TYPE2 = TYPE, mouse2 = mouse, day2 = day) %>%
	rename(., samples = seqID)
	
		# so you don't have to rename the columns later
	
# transform distance matrix to include relevant metadata:
dist_df <- all_dist %>%
	merge(meta[, c(1:7)], ., by.x="samples", by.y="samples") %>%
	merge(meta[, c(8:length(meta))], ., by.x="name2", by.y="name") %>%		# now merge again to match 'names'
	as_tibble()
		# you now have a df where all the '2' columns match data of name2, and all the others match the samples column

# now let's make some distance comparisons that you want
# define which GROUPs we want:
wanted <- c("AMS001", "AMS005", "D2_A", "D4_A", "D5_A", "R4_F", "yWT", "rWT", "yRag", "yWT_spore", "none")
removed <- c("pre-abx", "post-cef", "post-cdi", "post-vanco")
# easy to check what you are getting rid of when you look at this:
#test %>% group_by(GROUP, TYPE, GROUP2, TYPE2) %>% count(GROUP) %>% print(n=60)


# A: match mFMT to: own input, itself (mFMT), mFMT_other, hFMT, and noFMT
dist_df %>%
	filter(!GROUP2 %in% removed | !GROUP %in% removed, 
		!GROUP2 %in% removed & !GROUP %in% removed) %>%		#get rid of other comparisons that we do not want (whether either is true, or they match each other)
	filter(group2=="mFMT" | group=="mFMT", 
		TYPE2 %in% c("cecal", "FMTinput") & TYPE %in% c("cecal", "FMTinput")) %>%					#get only mFMT comparisons, cecal for now
	#probably good for now? let's combine some columns for ease:
	mutate(	s1_group = paste(TYPE, group, sep="_"),
			s2_group = paste(TYPE2, group2, sep="_")) %>%
	#test %>% group_by(s1_group, s2_group) %>% count(s1_group) %>% print(n=50)
	#now, we should be able to set the comparisons:
	mutate(comp = case_when(
		s1_group == "cecal_mFMT" & s1_group == s2_group ~ "mFMT-mFMT", 			
		s1_group == "FMTinput_mFMT" & s2_group == "cecal_mFMT" & GROUP==GROUP2 | s2_group == "FMTinput_mFMT" & s1_group == "cecal_mFMT"  & GROUP==GROUP2 ~ "mFMT-input",
		s1_group == "cecal_hFMT" & s2_group == "cecal_mFMT" | s2_group == "cecal_hFMT" & s1_group == "cecal_mFMT" ~ "mFMT-hFMT",
		s1_group == "cecal_noFMT" & s2_group == "cecal_mFMT" | s2_group == "cecal_noFMT" & s1_group == "cecal_mFMT" ~ "mFMT-noFMT",
		s1_group == "cecal_mFMT_other" & s2_group == "cecal_mFMT" | s2_group == "cecal_mFMT_other" & s1_group == "cecal_mFMT" ~ "mFMT-mFMT_other",
		TRUE ~ NA_character_),
		comp = factor(comp, levels=c("mFMT-input", "mFMT-mFMT", "mFMT-mFMT_other", "mFMT-hFMT", "mFMT-noFMT"))
		) %>% #count(comp) # use this to check how many groupings you have
	filter(!is.na(comp)) %>% # stop here if just saving
	# if you want to continue graphing...
	ggplot(aes(x=comp, y=value)) +
	geom_jitter(width=0.25, color="grey", alpha=0.6) +
	stat_summary(fun.data=median_hilow, color="red", size=1,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="Bray-Curtis distances") +
	scale_x_discrete(breaks=c("mFMT-input", "mFMT-mFMT", "mFMT-mFMT_other", "mFMT-hFMT", "mFMT-noFMT"), 
			labels=c("mFMT-input", "mFMT-mFMT", "mFMT-mFMT_other", "mFMT-hFMT", "mFMT-noFMT")) +
	scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.2)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("FigX_bray.dist.dotplot_mFMT.png", width=4, height=4)			

# B: match hFMT to: own input, itself (hFMT), hFMT_other, mFMT, and noFMT
dist_df %>%
	filter(!GROUP2 %in% removed | !GROUP %in% removed, 
		!GROUP2 %in% removed & !GROUP %in% removed) %>%		#get rid of other comparisons that we do not want (whether either is true, or they match each other)
	filter(group2=="hFMT" | group=="hFMT", 
		TYPE2 %in% c("cecal", "FMTinput") & TYPE %in% c("cecal", "FMTinput")) %>%					#get only mFMT comparisons, cecal for now
	#probably good for now? let's combine some columns for ease:
	mutate(	s1_group = paste(TYPE, group, sep="_"),
			s2_group = paste(TYPE2, group2, sep="_")) %>%
	#test %>% group_by(s1_group, s2_group) %>% count(s1_group) %>% print(n=50)
	#now, we should be able to set the comparisons:
	mutate(comp = case_when(
		s1_group == "cecal_hFMT" & s1_group == s2_group ~ "hFMT-hFMT", 			
		s1_group == "FMTinput_hFMT" & s2_group == "cecal_hFMT" | s2_group == "FMTinput_hFMT" & s1_group == "cecal_hFMT" ~ "hFMT-input",
		s1_group == "cecal_hFMT" & s2_group == "cecal_mFMT" | s2_group == "cecal_hFMT" & s1_group == "cecal_mFMT" ~ "hFMT-mFMT",
		s1_group == "cecal_noFMT" & s2_group == "cecal_hFMT" | s2_group == "cecal_noFMT" & s1_group == "cecal_hFMT" ~ "hFMT-noFMT",
		s1_group == "cecal_mFMT_other" & s2_group == "cecal_hFMT" | s2_group == "cecal_mFMT_other" & s1_group == "cecal_hFMT" ~ "hFMT-mFMT_other",
		TRUE ~ NA_character_),
		comp = factor(comp, levels=c("hFMT-input", "hFMT-hFMT", "hFMT-mFMT_other", "hFMT-mFMT", "hFMT-noFMT"))
		) %>% #count(comp) # use this to check how many groupings you have
	filter(!is.na(comp)) %>% 
	# if you want to continue graphing...
	ggplot(aes(x=comp, y=value)) +
	geom_jitter(width=0.25, color="grey", alpha=0.6) +
	stat_summary(fun.data=median_hilow, color="red", size=1,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="Bray-Curtis distances") +
	scale_x_discrete(breaks=c("hFMT-input", "hFMT-hFMT", "hFMT-mFMT_other", "hFMT-mFMT", "hFMT-noFMT"), 
			labels=c("hFMT-input", "hFMT-hFMT", "hFMT-mFMT_other", "hFMT-mFMT", "hFMT-noFMT")) +
	scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.2)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("FigX_bray.dist.dotplot_hFMT.png", width=4, height=4)	
	
# C: compare some of the other noFMTs
dist_df %>%
	filter(!GROUP2 %in% removed | !GROUP %in% removed, 
		!GROUP2 %in% removed & !GROUP %in% removed) %>%		#get rid of other comparisons that we do not want (whether either is true, or they match each other)
	filter(group2 %in% c("noFMT", "mFMT_other") | group %in% c("noFMT", "mFMT_other"), 
		TYPE2 %in% c("cecal", "FMTinput") & TYPE %in% c("cecal", "FMTinput")) %>%					#get only mFMT comparisons, cecal for now
	#probably good for now? let's combine some columns for ease:
	mutate(	s1_group = paste(TYPE, group, sep="_"),
			s2_group = paste(TYPE2, group2, sep="_")) %>%
	#test %>% group_by(s1_group, s2_group) %>% count(s1_group) %>% print(n=50)
	#now, we should be able to set the comparisons:
	mutate(comp = case_when(
		s1_group == "cecal_noFMT" & s1_group == s2_group ~ "noFMT-noFMT", 			
		s1_group == "FMTinput_hFMT" & s2_group == "cecal_noFMT" | s2_group == "FMTinput_hFMT" & s1_group == "cecal_noFMT" ~ "noFMT-hFMT(input)",
		s1_group == "FMTinput_mFMT" & s2_group == "cecal_noFMT" | s2_group == "FMTinput_mFMT" & s1_group == "cecal_noFMT" ~ "noFMT-mFMT(input)",
		s1_group == "FMTinput_mFMT_other" & s2_group == "cecal_noFMT" | s2_group == "FMTinput_mFMT_other" & s1_group == "cecal_noFMT" ~ "noFMT-mFMT_other(input)",
		s1_group == "FMTinput_mFMT_other" & s2_group == "cecal_mFMT_other" | s2_group == "FMTinput_mFMT_other" & s1_group == "cecal_mFMT_other" ~ "mFMT_other-input",
		TRUE ~ NA_character_),
		comp = factor(comp, levels=c("noFMT-noFMT", "noFMT-hFMT(input)", "noFMT-mFMT(input)", "noFMT-mFMT_other(input)", "mFMT_other-input"))
		) %>% #count(comp) # use this to check how many groupings you have
	filter(!is.na(comp)) %>% 
	# if you want to continue graphing...
	ggplot(aes(x=comp, y=value)) +
	geom_jitter(width=0.25, color="grey", alpha=0.6) +
	stat_summary(fun.data=median_hilow, color="red", size=1,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="Bray-Curtis distances") +
	scale_x_discrete(breaks=c("noFMT-noFMT", "noFMT-hFMT(input)", "noFMT-mFMT(input)", "noFMT-mFMT_other(input)", "mFMT_other-input"), 
			labels=c("noFMT-noFMT", "noFMT-hFMT(input)", "noFMT-mFMT(input)", "noFMT-mFMT_other(input)", "mFMT_other-input")) +
	scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.2)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("FigX_bray.dist.dotplot_noFMT.png", width=4, height=4)

### let's graph these all together; first, by making each df, combining them, then graphing them how we want to:
# mFMT samples
# to check on what you are selecting: group_by(s1_group, GROUP, s2_group, GROUP2, comp) %>% count(comp) %>% print(n=60)
mFMT_comp <- dist_df %>%
	filter(!GROUP2 %in% removed | !GROUP %in% removed, 
		!GROUP2 %in% removed & !GROUP %in% removed) %>%		#get rid of other comparisons that we do not want (whether either is true, or they match each other)
	filter(group2=="mFMT" | group=="mFMT", 
		TYPE2 %in% c("cecal", "FMTinput") & TYPE %in% c("cecal", "FMTinput")) %>%					#get only mFMT comparisons, cecal for now
	#probably good for now? let's combine some columns for ease:
	mutate(	s1_group = paste(TYPE, group, sep="_"),
			s2_group = paste(TYPE2, group2, sep="_")) %>%
	#test %>% group_by(s1_group, s2_group) %>% count(s1_group) %>% print(n=50)
	#now, we should be able to set the comparisons:
	mutate(comp = case_when(
		s1_group == "cecal_mFMT" & s1_group == s2_group ~ "mFMT-mFMT", 			
		s1_group == "FMTinput_mFMT" & s2_group == "cecal_mFMT" & GROUP==GROUP2 | s2_group == "FMTinput_mFMT" & s1_group == "cecal_mFMT"  & GROUP==GROUP2 ~ "mFMT-input",
		s1_group == "cecal_hFMT" & s2_group == "cecal_mFMT" | s2_group == "cecal_hFMT" & s1_group == "cecal_mFMT" ~ "mFMT-hFMT",
		s1_group == "cecal_noFMT" & s2_group == "cecal_mFMT" | s2_group == "cecal_noFMT" & s1_group == "cecal_mFMT" ~ "mFMT-noFMT",
		s1_group == "cecal_mFMT_other" & s2_group == "cecal_mFMT" | s2_group == "cecal_mFMT_other" & s1_group == "cecal_mFMT" ~ "mFMT-mFMT_other",
		TRUE ~ NA_character_),
		comp = factor(comp, levels=c("mFMT-input", "mFMT-mFMT", "mFMT-mFMT_other", "mFMT-hFMT", "mFMT-noFMT"))
		) %>% #count(comp) # use this to check how many groupings you have
		# to check on what you are selecting: group_by(s1_group, GROUP, s2_group, GROUP2, comp) %>% count(comp) %>% print(n=60)
	filter(!is.na(comp))

#hFMT samples:
#note: removed the duplicated mFMT-hFMT row here
hFMT_comp <- dist_df %>%
	filter(!GROUP2 %in% removed | !GROUP %in% removed, 
		!GROUP2 %in% removed & !GROUP %in% removed) %>%		#get rid of other comparisons that we do not want (whether either is true, or they match each other)
	filter(group2=="hFMT" | group=="hFMT", 
		TYPE2 %in% c("cecal", "FMTinput") & TYPE %in% c("cecal", "FMTinput")) %>%					#get only mFMT comparisons, cecal for now
	#probably good for now? let's combine some columns for ease:
	mutate(	s1_group = paste(TYPE, group, sep="_"),
			s2_group = paste(TYPE2, group2, sep="_")) %>%
	#test %>% group_by(s1_group, s2_group) %>% count(s1_group) %>% print(n=50)
	#now, we should be able to set the comparisons:
	mutate(comp = case_when(
		s1_group == "cecal_hFMT" & s1_group == s2_group ~ "hFMT-hFMT", 			
		s1_group == "FMTinput_hFMT" & s2_group == "cecal_hFMT" | s2_group == "FMTinput_hFMT" & s1_group == "cecal_hFMT" ~ "hFMT-input",
		s1_group == "cecal_noFMT" & s2_group == "cecal_hFMT" | s2_group == "cecal_noFMT" & s1_group == "cecal_hFMT" ~ "hFMT-noFMT",
		s1_group == "cecal_mFMT_other" & s2_group == "cecal_hFMT" | s2_group == "cecal_mFMT_other" & s1_group == "cecal_hFMT" ~ "hFMT-mFMT_other",
		TRUE ~ NA_character_),
		comp = factor(comp, levels=c("hFMT-input", "hFMT-hFMT", "hFMT-mFMT_other", "hFMT-noFMT"))
		) %>% #count(comp) # use this to check how many groupings you have
	filter(!is.na(comp))

# noFMT and other input:
noFMT_comp <- dist_df %>%
	filter(!GROUP2 %in% removed | !GROUP %in% removed, 
		!GROUP2 %in% removed & !GROUP %in% removed) %>%		#get rid of other comparisons that we do not want (whether either is true, or they match each other)
	filter(group2 %in% c("noFMT", "mFMT_other") | group %in% c("noFMT", "mFMT_other"), 
		TYPE2 %in% c("cecal", "FMTinput") & TYPE %in% c("cecal", "FMTinput")) %>%					#get only mFMT comparisons, cecal for now
	filter(GROUP == GROUP2) %>%
	#probably good for now? let's combine some columns for ease:
	mutate(	s1_group = paste(TYPE, group, sep="_"),
			s2_group = paste(TYPE2, group2, sep="_")) %>%
	#test %>% group_by(s1_group, s2_group) %>% count(s1_group) %>% print(n=50)
	#now, we should be able to set the comparisons:
	mutate(comp = case_when(
		s1_group == "cecal_noFMT" & s1_group == s2_group ~ "noFMT-noFMT", 			
		s1_group == "FMTinput_hFMT" & s2_group == "cecal_noFMT" | s2_group == "FMTinput_hFMT" & s1_group == "cecal_noFMT" ~ "noFMT-hFMT(input)",
		s1_group == "FMTinput_mFMT" & s2_group == "cecal_noFMT" | s2_group == "FMTinput_mFMT" & s1_group == "cecal_noFMT" ~ "noFMT-mFMT(input)",
		s1_group == "FMTinput_mFMT_other" & s2_group == "cecal_noFMT" | s2_group == "FMTinput_mFMT_other" & s1_group == "cecal_noFMT" ~ "noFMT-mFMT_other(input)",
		s1_group == "FMTinput_mFMT_other" & s2_group == "cecal_mFMT_other" | s2_group == "FMTinput_mFMT_other" & s1_group == "cecal_mFMT_other" ~ "mFMT_other-input",
		s1_group == "cecal_mFMT_other" & s1_group == s2_group ~ "mFMT_other-mFMT_other",
		TRUE ~ NA_character_),
		comp = factor(comp, levels=c("noFMT-noFMT", "mFMT_other-mFMT_other", "noFMT-hFMT(input)", "noFMT-mFMT(input)", "noFMT-mFMT_other(input)", "mFMT_other-input"))
		) %>% #count(comp) # use this to check how many groupings you have
	filter(!is.na(comp))
	
# combine these together, then graph:
# if doing all levels:
#level_orders <- c("mFMT-input", "mFMT_other-input", "hFMT-input", "noFMT-mFMT(input)", "noFMT-mFMT_other(input)", "noFMT-hFMT(input)", 
#						"mFMT-mFMT", "mFMT-mFMT_other", "mFMT-hFMT", "mFMT-noFMT", 
#						"hFMT-hFMT", "hFMT-mFMT_other", "hFMT-noFMT", "noFMT-noFMT")

# if choosing specific groups
level_orders <- c("mFMT-input", "mFMT_other-input", "hFMT-input",  
						"mFMT-mFMT", "mFMT_other-mFMT_other", "hFMT-hFMT", "noFMT-noFMT",
						"mFMT-mFMT_other", "mFMT-hFMT", "mFMT-noFMT", "hFMT-mFMT_other", "hFMT-noFMT")
data <- bind_rows(list(mFMT_comp, hFMT_comp, noFMT_comp), .id = "id") %>%
	filter(comp %in% level_orders) %>%
	mutate(comp = factor(comp, levels=level_orders))
	# graph:
data %>% ggplot(aes(x=comp, y=value)) +
	geom_jitter(width=0.25, color="grey", alpha=0.6, size = 0.5) +
	stat_summary(fun.data=median_hilow, color="red", size=1,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="Bray-Curtis distances") +
	scale_x_discrete(breaks=level_orders, 
			labels=str_replace(level_orders, "-", ":")) +
	scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.2)) +
	geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey20") +
	geom_vline(xintercept = 7.5, linetype = "dashed", color = "grey20") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
	
## stats:
library(FSA)
# inputs:
data %>% subset(., comp %in% c("mFMT-input", "mFMT_other-input", "hFMT-input")) %>%
	dunnTest(value ~ comp, data = .)
#                     Comparison         Z      P.unadj        P.adj
#1 hFMT-input - mFMT_other-input  8.761986 1.918401e-18 3.836803e-18
#2       hFMT-input - mFMT-input 19.416985 5.545127e-84 1.663538e-83
#3 mFMT_other-input - mFMT-input  1.527488 1.266397e-01 1.266397e-01
data %>% subset(., comp %in% c("mFMT-mFMT", "mFMT_other-mFMT_other", "hFMT-hFMT", "noFMT-noFMT")) %>%
	dunnTest(value ~ comp, data = .)
#                           Comparison          Z      P.unadj        P.adj
#1   hFMT-hFMT - mFMT_other-mFMT_other  7.4697533 8.034532e-14 4.820719e-13
#2               hFMT-hFMT - mFMT-mFMT  0.8470408 3.969724e-01 1.000000e+00
#3   mFMT_other-mFMT_other - mFMT-mFMT -6.7737443 1.254914e-11 5.019655e-11
#4             hFMT-hFMT - noFMT-noFMT -0.2395692 8.106643e-01 8.106643e-01
#5 mFMT_other-mFMT_other - noFMT-noFMT -7.1242962 1.046139e-12 5.230694e-12
#6             mFMT-mFMT - noFMT-noFMT -0.8411676 4.002541e-01 8.005081e-01
data %>% subset(., comp %in% c("mFMT-mFMT_other", "mFMT-hFMT", "mFMT-noFMT", "hFMT-mFMT_other", "hFMT-noFMT")) %>%
	dunnTest(value ~ comp, data = .)
#                          Comparison           Z       P.unadj         P.adj
##1       hFMT-mFMT_other - hFMT-noFMT   8.8700703  7.309975e-19  2.923990e-18
#2        hFMT-mFMT_other - mFMT-hFMT  -0.7381095  4.604479e-01  4.604479e-01
##3             hFMT-noFMT - mFMT-hFMT -10.6061882  2.788962e-26  1.394481e-25
##4  hFMT-mFMT_other - mFMT-mFMT_other  19.6677457  4.074986e-86  3.259989e-85
##5       hFMT-noFMT - mFMT-mFMT_other  13.3287898  1.574306e-40  1.102015e-39
##6        mFMT-hFMT - mFMT-mFMT_other  21.4607619 3.623177e-102 3.260859e-101
##7       hFMT-mFMT_other - mFMT-noFMT  -4.7365450  2.173925e-06  6.521774e-06
##8            hFMT-noFMT - mFMT-noFMT -13.0113849  1.054142e-38  6.324853e-38
#9             mFMT-hFMT - mFMT-noFMT  -4.4222287  9.768797e-06  1.953759e-05		#significantly different
##10      mFMT-mFMT_other - mFMT-noFMT -22.4907858 5.108689e-112 5.108689e-111



```

###-------------------------
### Figure 2C: bargraphs of taxonomy

### inputs:
	- Data/16S/allhf_genfrac2p.meta.txt
	- Data/16S/allhf_runs/allhf_genfrac2p.txt --> for tax colors


```{r}

###--------------------
### C:  barcharts (cecal only)
### Look at File Prep file for how inputs were formed
library(tidyr)
library(reshape2)
library(tidyverse)

data <- read.table("Data/16S/HFm_summary.nmds.new.groups.txt", header=T, sep="\t")
samples <- as.character(data$seqID)		
phylo <- read.table(file="Data/16S/allhf_runs/allhf_genfrac2p.meta.txt", header=T) %>%
					#write_tsv(., "Data/allhf_genfrac2p.meta.txt") %>%
					.[.$seqID %in% samples, ] %>%
					#merge(data, .[, c(1, 31:length(.))], by.x="seqID", by.y="seqID")	#but doesn't work when combined?
					select(seqID, Bacteroides:other) %>%
					merge(data, ., by.x="seqID", by.y="seqID") #%>%
					#write_tsv(., "Data/HFm_genfrac2p.meta.txt") 

# let's look at what samples we have (since we have both fecal and cecal):					
phylo %>%
  group_by(TYPE, GROUP, day) %>%
  summarize(length(mouse)) %>% print(n=36)

# for genus bargraphs, let's take only cecal and FMTinput:
df <- phylo[phylo$TYPE %in% c("FMTinput", "cecal"), ]

# rename/reorder bargraph groups to reflect the groups we want to show
# names of numbered FMT_inputs:
#df[df$FMT %in% c("hFMT2"), c("TYPE", "GROUP", "FMT_input", "FMT_type", "FMT")]
	# since no hFMT 2 (??) in this dataset (that one was a C.diff-infected one), let's redo numbering scheme:
df <- df %>%
  mutate(GROUP = recode(GROUP, 'yWT' = 'mFMT1', 'yWT_spore' = 'mFMT4', 'rWT' = 'mFMT2', 'yRag' = 'mFMT3', 
  								'AMS001' = 'hFMT2', 'AMS005' = 'hFMT3', 'D2_A' = 'hFMT1', 
  								'D4_A' = 'hFMT4', 'D5_A' = 'hFMT5', 'R4_F' = 'hFMT6'))
df$GROUP[df$TYPE %in% c("cecal")] <- paste(df$GROUP[df$TYPE %in% c("cecal")], df$day[df$TYPE %in% c("cecal")], sep="_")
df$GROUP[df$GROUP %in% c('pre-abx_NA')] <- "pre-abx"
df$GROUP <- as.factor(factor(df$GROUP, levels=c('pre-abx', 'mFMT1', 'mFMT1_21', 'mFMT1_42', 'mFMT2', 'mFMT2_42', 
											'mFMT3', 'mFMT3_42', 'mFMT4', 'mFMT4_42',
											'hFMT1', 'hFMT1_21', 'hFMT1_42', 'hFMT2', 'hFMT2_21',
											'hFMT3', 'hFMT3_21', 'hFMT4', 'hFMT4_42', 
											'hFMT5', 'hFMT5_42', 'hFMT6', 'hFMT6_42',
											'none_21', 'none_42')))
# now, calculate averages for these, and graph:
library(plyr)
means <- plyr::ddply(df, c("GROUP"), colwise(mean, is.numeric))
rownames(means) <- means$GROUP
means <- select(means, Bacteroides:other)

# let's also plit these into 3 groups (for ease of separation
bar_m <- as.data.frame(means[1:10,])
bar_h <- as.data.frame(means[11:23,])
bar_no <- as.data.frame(means[24:25,])

# colors:
taxcolors <- read.table("Data/16S/allhf_runs/allhf_genfrac2p.txt", header=T, comment.char = "")
tax.col <- as.character(taxcolors$color)

# graphing both sets:
# Figure S1:
bar_g<-as.matrix(t(bar_m))
par(mar=c(5,4,2,5), mfrow=c(3,1))
par(xpd=T)
plot1 <- barplot(bar_g, main="Taxonomic comparison", ylab="Relative abundance-genera (%)", cex.names=0.8, ylim=c(0,100), col=tax.col, xlim=c(0,13), axisnames=FALSE)
text(x = plot1, y = par("usr")[3]-5, srt = 45, adj = 1, labels = unique(colnames(bar_g)), xpd = TRUE, cex=0.8)
#legend(4,110,legend=rownames(bar_g),col=tax.col,fill=tax.col,cex=0.5)	

bar_g<-as.matrix(t(bar_h))
plot2 <- barplot(bar_g, main="", ylab="Relative abundance-genera (%)", cex.names=0.8, ylim=c(0,100), col=tax.col, xlim=c(0,13), axisnames=FALSE)
text(x = plot2, y = par("usr")[3]-5, srt = 45, adj = 1, labels = unique(colnames(bar_g)), xpd = TRUE, cex=0.8)

bar_g<-as.matrix(t(bar_no))
plot3 <- barplot(bar_g, main="", ylab="Relative abundance-genera (%)", cex.names=0.8, ylim=c(0,100), col=tax.col, xlim=c(0,13), axisnames=FALSE)
text(x = plot3, y = par("usr")[3]-5, srt = 45, adj = 1, labels = unique(colnames(bar_g)), xpd = TRUE, cex=0.8)
legend(3,120,legend=rownames(bar_g)[1:32],col=tax.col[1:32],fill=tax.col[1:32],cex=0.5)
legend(7,120,legend=rownames(bar_g)[33:63],col=tax.col[33:63],fill=tax.col[33:63],cex=0.5)
legend(11,120,legend=rownames(bar_g)[64:nrow(bar_g)],col=tax.col[64:nrow(bar_g)],fill=tax.col[64:nrow(bar_g)],cex=0.5)

### let's also calculate for OVERALL treatments:
# Figure 2C
data <- read.table("Data/16S/HFm_summary.nmds.new.groups.txt", header=T, sep="\t")
samples <- as.character(data$seqID)		
phylo <- read.table(file="Data/16S/allhf_runs/allhf_genfrac2p.meta.txt", header=T) %>%
					#write_tsv(., "allhf_runs/allhf_genfrac2p.meta.txt") %>%
					.[.$seqID %in% samples, ] %>%
					#merge(data, .[, c(1, 31:length(.))], by.x="seqID", by.y="seqID")	#but doesn't work when combined?
					select(seqID, Bacteroides:other) %>%
					merge(data, ., by.x="seqID", by.y="seqID") #%>%
					#write_tsv(., "/HFm_genfrac2p.meta.txt") 

# let's look at what samples we have (since we have both fecal and cecal):					
phylo %>%
  group_by(TYPE, group, day) %>%
  summarize(length(mouse)) %>% print(n=36)

# for genus bargraphs, let's take only cecal and FMTinput:
df <- phylo[phylo$TYPE %in% c("FMTinput", "cecal"), ]

# combine with day and rerorder
df$group[df$TYPE %in% c("cecal")] <- paste(df$group[df$TYPE %in% c("cecal")], df$day[df$TYPE %in% c("cecal")], sep="_")
df$group[df$seqID %in% c('healthy_cecal_1912_1', 'healthy_cecal_1913_1', 'healthy_cecal_1913_2')] <- "pre-abx"
df$group <- as.factor(factor(df$group, levels=c('pre-abx', 'mFMT', 'mFMT_21', 'mFMT_42', 
											'mFMT_other', 'mFMT_other_42',
											'hFMT', 'hFMT_21', 'hFMT_42',
											'noFMT_21', 'noFMT_42')))
# now, calculate averages for these, and graph:
library(plyr)
means <- plyr::ddply(df, c("group"), colwise(mean, is.numeric))
rownames(means) <- means$group
means <- select(means, Bacteroides:other)

# colors:
taxcolors <- read.table("Data/16S/allhf_runs/allhf_genfrac2p.txt", header=T, comment.char = "")
tax.col <- as.character(taxcolors$color)

# graphing both sets:
bar_g<-as.matrix(t(means))
par(mar=c(5,4,2,3), mfrow=c(1,2))
par(xpd=T)
plot1 <- barplot(bar_g, main="Taxonomic comparison", ylab="Relative abundance-genera (%)", cex.names=0.8, ylim=c(0,100), col=tax.col, xlim=c(0,13), axisnames=FALSE)
text(x = plot1, y = par("usr")[3]-5, srt = 45, adj = 1, labels = unique(colnames(bar_g)), xpd = TRUE, cex=0.8)
barplot(bar_g, col = NA, border = NA, axes = FALSE, axisnames=FALSE)
legend(-3,110,legend=rownames(bar_g)[1:10],col=tax.col[1:10],fill=tax.col[1:10],cex=0.5)
legend(1,140,legend=rownames(bar_g)[11:42],col=tax.col[11:42],fill=tax.col[11:42],cex=0.5)
legend(5,140,legend=rownames(bar_g)[43:62],col=tax.col[43:62],fill=tax.col[43:62],cex=0.5)
legend(9,110,legend=rownames(bar_g)[63:nrow(bar_g)],col=tax.col[63:nrow(bar_g)],fill=tax.col[63:nrow(bar_g)],cex=0.5)


```


###-------------------------
### Figure 2D/E: diversity over time (D) and at days 21, 42 (E)

### inputs:
- Data/16S/HFm_summary.nmds.new.groups.txt


```{r}
###--------------------
### E:  diversity; overtime using fecal samples
library(tidyr)
library(reshape2)
library(tidyverse)

# fecal, overtime:
data <- read.table("Data/16S/HFm_summary.nmds.new.groups.txt", header=T, sep="\t") %>%
		filter(type == "fecal" & sampleType == "mouse") %>%
		filter(!is.na(shannon)) %>%
		dplyr::group_by(group, day) %>%
		mutate(median = median(shannon),
    			Q1 = quantile(shannon, 0.25),
    			Q3 = quantile(shannon, 0.75))
# count number of events
data %>%
  group_by(group) %>%
  summarise(n_distinct(mouse))
#1 hFMT                        49
#2 mFMT                        25
#3 mFMT_other                  18
#4 noFMT                       20

# colors and polygons:
dcol4 = c("#EACB2B", "#56B29E", "#9FC095", "#F5191C")	# now in beast.colors

# graph
data %>% ggplot(aes(x=day, y=shannon, color=group, fill=group)) +
	scale_color_manual(values=beast.colors$dcol4) + 
	scale_fill_manual(values=beast.colors$dcol4) + 
	geom_vline(xintercept = 0, linetype = "solid", color = "red", alpha = 0.5, size = 0.5) +
	geom_vline(xintercept = 11, linetype = "solid", color = "springgreen4", alpha = 0.5, size = 0.5) +
	annotate("text", x = 0, y = 115, label = "C. diff", size = 3, color = "red", alpha=0.5) +
	annotate("text", x = 11, y = 115, label = "FMT", size = 3, color = "springgreen4") +
	#geom_ribbon(aes(ymin = Q1, ymax = Q3, group=group), alpha=0.2, color=NA) +
  	geom_jitter(width=0.25, alpha=0.4, size = 1, na.rm=TRUE) + #if you wanted points
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="line", linewidth=1,
			fun.args = list(conf.int=0.50)) +
  	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="point", linewidth=1,
			fun.args = list(conf.int=0.50)) +
  	scale_y_continuous(limits=c(0,4.5)) +
	scale_x_continuous(limits=c(-7, 42), breaks=c(-7, 0, 1, 4, 9, 11, 13, 16, 20, 28, 32, 36, 41)) +
	labs(x=NULL, y="shannon diversity index") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	theme_classic()
	#geom_polygon(data = cef, aes(x = x, y = y), fill = 'grey', alpha = 0.3) + # couldn't figure this out

# stats:
subdf <- data
subdf$day <- as.factor(subdf$day)
factor_levels <- unique(subdf$day)

# Then, run KW test through oop through factor levels
for (level in factor_levels) {
  subset_data <- filter(subdf, day == level & !is.na(shannon) )  # Subset data by factor level
  
  # Check if the subsetted data contains more than one group
  if (length(unique(subset_data$group)) > 1) {
    # Run Kruskal-Wallis test
    result <- kruskal.test(shannon ~ group, data = subset_data)
  
    # Print the results
    cat("Factor level:", level, "\n")
    print(result)
    cat("\n")
  } else {
    cat("Factor level:", level, "has only one group.\n\n")
  }
}
### Note: some days do show a small significantly different weight? But, not consistent on which ones those are
#...let's check out which ones using a posthoc test
# testing days AFTER treatment
library(FSA)

subdf %>% subset(., day == 11 & !is.na(shannon) ) %>%
	dunnTest(shannon ~ group, data = .)		
#          Comparison          Z     P.unadj      P.adj
#1        hFMT - mFMT -1.5933035 0.111092117 0.33327635
#2  hFMT - mFMT_other  0.4338744 0.664379675 0.66437968
#3  mFMT - mFMT_other  1.6129770 0.106749513 0.42699805
#4       hFMT - noFMT -2.7307540 0.006318963 0.03791378
#5       mFMT - noFMT -1.1214278 0.262105796 0.52421159
#6 mFMT_other - noFMT -2.5272623 0.011495561 0.05747780
subdf %>% subset(., day == 12 & !is.na(shannon) ) %>%
	dunnTest(shannon ~ group, data = .)	
#    Comparison         Z      P.unadj        P.adj
#1  hFMT - mFMT -2.258730 0.0239001556 0.0478003113
#2 hFMT - noFMT  1.621862 0.1048330029 0.1048330029
#3 mFMT - noFMT  3.825598 0.0001304553 0.0003913658
subdf %>% subset(., day == 13 & !is.na(shannon) ) %>%
#	dunnTest(shannon ~ group, data = .)
#    Comparison         Z      P.unadj        P.adj
#1  hFMT - mFMT -2.571689 1.012036e-02 1.012036e-02
#2 hFMT - noFMT  3.599096 3.193249e-04 6.386499e-04
#3 mFMT - noFMT  5.471026 4.474374e-08 1.342312e-07
subdf %>% subset(., day == 16 & !is.na(shannon) ) %>%
	dunnTest(shannon ~ group, data = .)
#          Comparison         Z      P.unadj        P.adj
#1        hFMT - mFMT -4.103627 4.067231e-05 1.220169e-04
#2  hFMT - mFMT_other -2.622797 8.721130e-03 1.744226e-02
#3  mFMT - mFMT_other  1.198535 2.307087e-01 2.307087e-01
#4       hFMT - noFMT  4.362197 1.287628e-05 5.150512e-05
#5       mFMT - noFMT  7.178184 7.064334e-13 4.238600e-12
#6 mFMT_other - noFMT  5.896013 3.723889e-09 1.861944e-08
subdf %>% subset(., day == 19 & !is.na(shannon) ) %>%
	dunnTest(shannon ~ group, data = .)
#          Comparison         Z      P.unadj        P.adj
#1        hFMT - mFMT -2.154733 3.118275e-02 6.236549e-02
#2  hFMT - mFMT_other -4.321599 1.549022e-05 6.196088e-05
#3  mFMT - mFMT_other -1.704310 8.832312e-02 8.832312e-02
#4       hFMT - noFMT  3.164510 1.553443e-03 4.660328e-03
#5       mFMT - noFMT  4.670142 3.009918e-06 1.504959e-05
#6 mFMT_other - noFMT  6.676056 2.454582e-11 1.472749e-10
subdf %>% subset(., day == 32 & !is.na(shannon) ) %>%
	dunnTest(shannon ~ group, data = .)
#          Comparison          Z      P.unadj        P.adj
#1        hFMT - mFMT -2.2752726 2.288958e-02 4.577917e-02
#2  hFMT - mFMT_other -4.1040198 4.060328e-05 2.030164e-04
#3  mFMT - mFMT_other -0.3733818 7.088643e-01 7.088643e-01
#4       hFMT - noFMT  2.9717567 2.961013e-03 8.883038e-03
#5       mFMT - noFMT  4.0973110 4.179771e-05 1.671908e-04
#6 mFMT_other - noFMT  6.2882888 3.209844e-10 1.925906e-09
subdf %>% subset(., day == 41 & !is.na(shannon) ) %>%
	dunnTest(shannon ~ group, data = .)
#          Comparison          Z      P.unadj        P.adj
#1        hFMT - mFMT -2.5838799 9.769580e-03 2.930874e-02
#2  hFMT - mFMT_other -2.3784864 1.738388e-02 3.476776e-02
#3  mFMT - mFMT_other  0.8886532 3.741895e-01 3.741895e-01
#4       hFMT - noFMT  3.3585145 7.836262e-04 3.134505e-03
#5       mFMT - noFMT  4.7885110 1.680233e-06 8.401166e-06
#6 mFMT_other - noFMT  5.1899773 2.103198e-07 1.261919e-06


###--------------------
### E:  diversity; days 21 and 42 using cecal
# use same code as for CFU over time; just need to filter out what you want:
library(tidyr)
library(reshape2)
library(tidyverse)

### cecal, end points:
data <- read.table("Data/16S/HFm_summary.nmds.new.groups.txt", header=T, sep="\t") %>%
		filter(type == "cecal" & sampleType == "mouse") %>%
		unite("group_c", day,group, sep ="_", remove=FALSE) %>%
		mutate(group_c = factor(group_c, levels=c('21_mFMT', '21_hFMT', '21_noFMT', 
										'42_mFMT', '42_mFMT_other', '42_hFMT', '42_noFMT')))
		
# what do we have?									
data %>%
  dplyr::group_by(group) %>%
  summarise(n_distinct(mouse, day))
#1 hFMT                             41
#2 mFMT                             22
#3 mFMT_other                       15
#4 noFMT                            20

# specify colors for graph: also in Beast.colors.R
dcol2_d21 = c("#56B29E", "#EACB2B", "#F5191C", "#56B29E", "#9FC095", "#EACB2B", "#F5191C")		
	
# graph
data %>% filter(!is.na(shannon)) %>%
	ggplot(aes(x=group_c, y=shannon, color=group_c)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_d21) + 
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y="shannon diversity index") +
	scale_x_discrete(breaks=c('21_mFMT', '21_hFMT', '21_noFMT', '42_mFMT', '42_mFMT_other', '42_hFMT', '42_noFMT'), 
			labels=c('mFMT', 'hFMT', 'noFMT', 'mFMT', 'mFMT-other', 'hFMT', 'noFMT')) +
	scale_y_continuous(limits=c(0,4.5)) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 0.8, xend = 3.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 3.5, xend = 7.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	annotate("text", x = 2, y = 0.5, label = "day 21", size = 3, color = "grey20") +
	annotate("text", x = 5.5, y = 0.5, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey20")
	#ggsave("FigX_.png", width=4, height=4)

### stats:
# day 21:
library(FSA)

data %>% subset(., day == 21 & !is.na(shannon) ) %>%
	dunnTest(shannon ~ group, data = .)		
#    Comparison         Z      P.unadj        P.adj
#1  hFMT - mFMT -4.120926 3.773529e-05 7.547058e-05
#2 hFMT - noFMT  2.436509 1.482981e-02 1.482981e-02
#3 mFMT - noFMT  5.665384 1.466957e-08 4.400871e-08

data %>% subset(., day == 42 & !is.na(shannon) ) %>%
	dunnTest(shannon ~ group, data = .)		
#          Comparison          Z      P.unadj        P.adj
#1        hFMT - mFMT -2.3033716 2.125794e-02 4.251589e-02
#2  hFMT - mFMT_other -4.2521495 2.117285e-05 1.058642e-04
#3  mFMT - mFMT_other -0.5318432 5.948346e-01 5.948346e-01
#4       hFMT - noFMT  2.5809548 9.852748e-03 2.955824e-02
#5       mFMT - noFMT  3.8619097 1.125041e-04 4.500165e-04
6 mFMT_other - noFMT  6.0168115 1.778861e-09 1.067316e-08


```

###--------------------
### Figure 1F: Maaslin for significantly different OTUs:

### inputs:
- Data/16S/HFm_summary.nmds.new.groups.txt
- Data/16S/HFm_shared.df.all.groups.input-for-dist.txt
- Data/16S/16S_maaslyn/16S_clear_v_colonized/significant_results.tsv
- Data/16S/taxonomy_meta.txt


```
###--------------------
### D: using maaslyn to identify different OTUs between col vs cleared
# let's use cecal samples for these

## load packages
library(RColorBrewer)
library(gplots)
library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(Maaslin2)
library(pheatmap)
library(vegan)

# for now, using same subsampling scheme as for calculating shared OTUs
# read in meta data information (copying over from diversity info)
# cecal:
d <- read_tsv("Data/16S/HFm_summary.nmds.new.groups.txt") %>%
		filter(!is.na(simpson)) %>%
		filter(nseqs > 2500) %>%
		filter(group %in% c('mFMT', 'mFMT_other', 'hFMT', 'noFMT')) %>%		# for now, only choose these groups
		filter(sampleType %in% c("mouse") & type %in% c('cecal'))	%>%		# for now, only mouse samples that are cecal
		as.data.frame() %>%
		select(seqID, FMT, group, GROUP) %>%
		mutate(status = case_when(
			group %in% c('mFMT', 'mFMT_other') ~ "cleared",
			group %in% c('noFMT', 'hFMT') ~ "colonized"),
			status = factor(status, levels=c("cleared", "colonized")))

# fecal
d_fecal <- read_tsv("Data/16S/HFm_summary.nmds.new.groups.txt") %>%
		filter(!is.na(simpson)) %>%
		filter(nseqs > 2500) %>%
		filter(group %in% c('mFMT', 'mFMT_other', 'hFMT', 'noFMT')) %>%		# for now, only choose these groups
		filter(sampleType %in% c("mouse") & type %in% c('fecal'))	%>%		# for now, only mouse samples that are cecal
		filter(day > 11) %>%
		as.data.frame() %>%
		select(seqID, FMT, group, GROUP) %>%
		mutate(status = case_when(
			group %in% c('mFMT', 'mFMT_other') ~ "cleared",
			group %in% c('noFMT', 'hFMT') ~ "colonized"),
			status = factor(status, levels=c("cleared", "colonized")))

# how many FMT groups represented here?
d %>%
  	dplyr::group_by(FMT) %>%
 	dplyr::summarize(count=n())

# read in shared file and subsample:
shared <- read_tsv("Data/16S/HFm_shared.df.all.groups.input-for-dist.txt") %>% 
	select(seqID, starts_with("Otu")) %>%
	column_to_rownames("seqID") %>%
	select_if(colSums(.) !=0)	# separated this bc it takes a long time...

# what is the lowest number of sequences in here?
nseqs <- shared %>% 
	rownames_to_column("seqID") %>%
	pivot_longer(-seqID) %>%
	dplyr::group_by(seqID) %>%
	dplyr::summarize(N = sum(value)) %>%	# to get number of total sequences
	arrange(N) #%>% print(n=20) 		# this is 2936

# Rarefy shared file to appropriate min (subsample):
set.seed(1984)
sub_size <- min(nseqs$N)
rshared <- as.data.frame(t(shared)) %>%
	.[, colSums(.) >= sub_size]			# if you wanted to actually eliminate any samples below this threshold
# now, subsample from columns?
for (index in 1:ncol(rshared)){
  	rshared[,index] <- t(rrarefy(rshared[,index], sample=sub_size))}
	rm(index, sub_size)	

# filter out the samples you want from the shared file and set up counts and meta:
## for cecal:
counts <- rshared %>%
  	t() %>%
  	as.data.frame() %>%
  	rownames_to_column("seqID") %>%
  	left_join(d, ., by="seqID") %>%
  	select(-FMT, -GROUP, -group, -status) %>%
	column_to_rownames("seqID")

d <- d %>% column_to_rownames("seqID")

## Maaslyn comparison, cecal, status:
Maaslin2(input_data = counts, input_metadata = d, 
         output = "16S_clear_v_colonized", # maaslin will create an entire directory of results in your working directory, this is what you want it named (note: can put a path to anther directory if you want)
         normalization="NONE", transform="LOG", analysis_method ="LM", correction="BH", 
         fixed_effects = c("status")) # this is the column in meta that you want to make the comparison based on

## Maaslyn comparison, cecal, group:
Maaslin2(input_data = counts, input_metadata = d, output = "16S_groups",
         normalization="NONE", transform="LOG", analysis_method ="LM", correction="BH",
         fixed_effects = c("group"), 
         reference = 'group,mFMT', # basically the same as above but you need to call out which column to pull your reference from and what value your setting as the reference
         plot_heatmap=TRUE)

## for fecal:
counts <- rshared %>%
  	t() %>%
  	as.data.frame() %>%
  	rownames_to_column("seqID") %>%
  	left_join(d_fecal, ., by="seqID") %>%
  	select(-FMT, -GROUP, -group, -status) %>%
	column_to_rownames("seqID")

d_fecal <- d_fecal %>% column_to_rownames("seqID")

# fecal status:
Maaslin2(input_data = counts, input_metadata = d_fecal, 
         output = "16S.fecal_clear_v_colonized", 
         normalization="NONE", transform="LOG", analysis_method ="LM", correction="BH", 
         fixed_effects = c("status")) 

## fecal group:
Maaslin2(input_data = counts, input_metadata = d_fecal, output = "16S.fecal_groups",
         normalization="NONE", transform="LOG", analysis_method ="LM", correction="BH",
         fixed_effects = c("group"), 
         reference = 'group,mFMT', 
         plot_heatmap=TRUE)

### let's combine these results with a taxonomic breakdown, and plot a heatmap of abundance

## for cecal--by status:
# convert rshared to proportion -> log-transform:
# note: using sub_size from above
detach(package:plyr)
d <- d %>% rownames_to_column("seqID")

# read in additional files:
tax_meta <- read_tsv("Data/16S/taxonomy_meta.txt")

# gather for wide format
rshared_long <- rshared %>%
	rownames_to_column("Otu") %>%
	gather(key="seqID", value="otu_n", -Otu) 
	
# filter out OTUs that are very low abundance:
rshared_ranges <- rshared_long %>%
		# log-transform and calculate mean OTU abundance:
	dplyr::group_by(seqID) %>%
	mutate(RA = otu_n/sub_size*100) %>%
	#mutate(RA = log10(RA+1)) %>%
	dplyr::group_by(Otu) %>%
	mutate(otu.mean = mean(RA)) %>%
	mutate(otu.med = median(RA)) %>%
	mutate(otu.min = min(RA)) %>%
	mutate(otu.max = max(RA)) %>%
	filter(otu.mean > 0.1)
# length(unique(rshared_ranges$Otu)): this limited to 77 Otus (0.1) or 300 Otus (0.01)
# can use this list to filter out low abundance OTUs
top_otus <- unique(rshared_ranges$Otu)

# prepare to plot:
## for cecal--by status:
## this is the figure used for 2F (flipped 180)
sig_otus <- read_tsv("Data/16S/16S_maaslyn/16S_clear_v_colonized/significant_results.tsv")
rshared_norm <- rshared_long %>%
		# log-transform and calculate mean per group:
	dplyr::group_by(seqID) %>%
	mutate(RA = otu_n/sub_size*100) %>%
	mutate(RA = log10(RA+1)) %>%	# at this point, you could convert to wide form again if you wanted a table
	filter(Otu %in% top_otus) %>%
	left_join(., d, by="seqID")	%>%# join with metadata for grouping for mean within group
	filter(group %in% c("mFMT", "hFMT", "mFMT_other", "noFMT")) %>%
	dplyr::group_by(group, Otu) %>% #if this does not work, detach plyr package
  	summarize(mean=mean(RA)) %>%
  	spread(key='group', value="mean") %>%
		# join with maaslyn sig otus:
	inner_join(., sig_otus, join_by("Otu" == "feature")) %>%
		# join with taxonomy info (including colors):
	left_join(., tax_meta, join_by("Otu" == "OTU")) %>%
	mutate(., labels = paste(Otu, Genus)) %>%
	mutate(labels = gsub("Otu[0]{2}", "Otu", labels)) %>%
	mutate(labels = gsub("_unclassified", "_uncl.", labels)) %>%
	# filter(BIOCHEMICAL %in% top_data_forest) %>% # need to add filtering later
	# mutate(BIOCHEMICAL=factor(BIOCHEMICAL, levels=rev(top_data_forest)))
	arrange(Phylum, Size, Genus) %>%
	column_to_rownames("labels")
	
# plot graph:
annotation_colors = list(genera = setNames(as.character(rshared_norm$color), rshared_norm$Genus))
genera <- data.frame(genera=rshared_norm$Genus)
rownames(genera) <- row.names(rshared_norm)

pheatmap(rshared_norm[c(3,4,2,5)],
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=3)(10),
         border_color="grey20",
         cellheight = 8,
         cellwidth = 20,
         fontsize = 8,
         filename="16S_maaslyn_heatmap_sum.pdf",
         #width=14,
         #height=14,
         annotation_row = genera,
         annotation_colors = annotation_colors,
         angle_col = "45",
         annotation_legend = FALSE
         #annotation_col = rshared_norm(color)
)

# prepare to plot:
## for cecal--by group:
sig_otus <- read_tsv("Data/16S/16S_maaslyn/16S_groups/significant_results.tsv")
rshared_norm <- rshared_long %>%
		# log-transform and calculate mean per group:
	dplyr::group_by(seqID) %>%
	mutate(RA = otu_n/sub_size*100) %>%
	mutate(RA = log10(RA+1)) %>%	# at this point, you could convert to wide form again if you wanted a table
	filter(Otu %in% top_otus) %>%
	left_join(., d, by="seqID")	%>%# join with metadata for grouping for mean within group
	filter(group %in% c("mFMT", "hFMT", "mFMT_other", "noFMT")) %>%
	dplyr::group_by(group, Otu) %>% #if this does not work, detach plyr package
  	summarize(mean=mean(RA)) %>%
  	spread(key='group', value="mean") %>%
		# join with maaslyn sig otus:
	inner_join(., sig_otus, join_by("Otu" == "feature")) %>%
		# join with taxonomy info (including colors):
	left_join(., tax_meta, join_by("Otu" == "OTU")) %>%
	mutate(., labels = paste(Otu, Genus)) %>%
	mutate(labels = gsub("Otu[0]{2}", "Otu", labels)) %>%
	mutate(labels = gsub("_unclassified", "_uncl.", labels)) %>%
	# filter(BIOCHEMICAL %in% top_data_forest) %>% # need to add filtering later
	# mutate(BIOCHEMICAL=factor(BIOCHEMICAL, levels=rev(top_data_forest)))
	arrange(Phylum, Size, Genus) %>%
	dplyr::group_by(Otu) %>% 
	filter(row_number() == 1) %>% 	#note: sig Otus from a group comparison contains multiple lines per Otu (between-group info)
	column_to_rownames("labels")
	
# plot graph:
annotation_colors = list(genera = setNames(as.character(rshared_norm$color), rshared_norm$Genus))
genera <- data.frame(genera=rshared_norm$Genus)
rownames(genera) <- row.names(rshared_norm)

pheatmap(rshared_norm[c(3,4,2,5)],
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=3)(10),
         border_color="grey20",
         cellheight = 8,
         cellwidth = 20,
         fontsize = 8,
         filename="16S_maaslyn_heatmap_sum.pdf",
         #width=14,
         #height=14,
         annotation_row = genera,
         annotation_colors = annotation_colors,
         angle_col = "45",
         annotation_legend = FALSE
         #annotation_col = rshared_norm(color)
)
### I don't think this has any other differences...

### Let's check out fecal, by status:
sig_otus <- read_tsv("Data/16S/16S_maaslyn/16S.fecal_clear_v_colonized/significant_results.tsv")
rshared_norm <- rshared_long %>%
		# log-transform and calculate mean per group:
	dplyr::group_by(seqID) %>%
	mutate(RA = otu_n/sub_size*100) %>%
	mutate(RA = log10(RA+1)) %>%	# at this point, you could convert to wide form again if you wanted a table
	filter(Otu %in% top_otus) %>%
	left_join(., d, by="seqID")	%>%# join with metadata for grouping for mean within group
	filter(group %in% c("mFMT", "hFMT", "mFMT_other", "noFMT")) %>%
	dplyr::group_by(group, Otu) %>% #if this does not work, detach plyr package
  	summarize(mean=mean(RA)) %>%
  	spread(key='group', value="mean") %>%
		# join with maaslyn sig otus:
	inner_join(., sig_otus, join_by("Otu" == "feature")) %>%
		# join with taxonomy info (including colors):
	left_join(., tax_meta, join_by("Otu" == "OTU")) %>%
	mutate(., labels = paste(Otu, Genus)) %>%
	mutate(labels = gsub("Otu[0]{2}", "Otu", labels)) %>%
	mutate(labels = gsub("_unclassified", "_uncl.", labels)) %>%
	# filter(BIOCHEMICAL %in% top_data_forest) %>% # need to add filtering later
	# mutate(BIOCHEMICAL=factor(BIOCHEMICAL, levels=rev(top_data_forest)))
	arrange(Phylum, Size, Genus) %>%
	column_to_rownames("labels")
	
# plot graph:
annotation_colors = list(genera = setNames(as.character(rshared_norm$color), rshared_norm$Genus))
genera <- data.frame(genera=rshared_norm$Genus)
rownames(genera) <- row.names(rshared_norm)

pheatmap(rshared_norm[c(3,4,2,5)],
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=3)(10),
         border_color="grey20",
         cellheight = 8,
         cellwidth = 20,
         fontsize = 8,
         #filename="16S_maaslyn_heatmap_sum.pdf",
         #width=14,
         #height=14,
         annotation_row = genera,
         annotation_colors = annotation_colors,
         angle_col = "45",
         annotation_legend = FALSE
         #annotation_col = rshared_norm(color)
)
```

#### Figure S2 16S OTU engraftment
### 2.10.25
### S. Millard

### Figures:
- Fig. S2A: 16S OTU % shared (fecal, cecal)
- Fig. S2B: Relative abundance of shared OTUs recipient (fecal,cecal)
- Fig. S2C: Relative abundance of OTUs in donor that did not colonize in recipient (fecal,cecal)


```{r}
# load packages
library(gplots)
library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(glue)
library(tibble)

# read in metadata
meta<-read.table(file="/Users/home/Library/CloudStorage/Box-Box/Manuscript/ForGit/mouseCDI-SPF-hFMT/Data/16S/HFm_summary.nmds.new.groups.txt", header=T) %>%
  select(seqID, sampleID,mouse, group, FMT, FMT_type, FMT_input,Model, Experiment, day, type,sampleType) %>%
  mutate(type = ifelse(sampleType == "FMTinput", "FMTinput", type)) 

# subset cecal data 
cecal_meta <- meta %>%
  filter(type %in% c("FMTinput","cecal")) %>%
  filter(!sampleType %in% c("mouse_survey"))

# subset fecal data
fecal_meta <- meta %>%
  filter(!sampleType %in% c("human")) %>% # remove the patient fecal samples that were not used for input
  filter(type %in% c("FMTinput","fecal")) %>%
  group_by(Experiment) %>%
  mutate(day=as.numeric(day)) %>%
  mutate(new_col = ifelse(day == max(day, na.rm = TRUE), "final_fecal", "not_final")) %>% # each experiment had a different day for final fecal collection...
  ungroup() %>%
  mutate(new_col = ifelse(type %in% c("FMTinput"), "FMTinput", new_col)) %>% # changing this will give us our criteria to filter samples by
  filter(new_col %in% c("final_fecal","FMTinput")) # filter

# read in otu counts
shared <- read.table(file="/Users/home/Library/CloudStorage/Box-Box/Manuscript/Code/Data/HFm_allgroups.rarefied.shared.txt",header=T)
# subset cecal data 

cecal_df <- shared %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU) %>%
  rename(seqID = name) %>%
  filter(seqID %in% cecal_meta$seqID) %>%
  left_join(., cecal_meta, by="seqID") %>%
  mutate(pres_abs = ifelse(value != 0, 1, 0)) 

# subset fecal data 
fecal_df <- shared %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU) %>%
  rename(seqID = name) %>%
  filter(seqID %in% fecal_meta$seqID) %>%
  left_join(., fecal_meta, by="seqID") %>%
  mutate(pres_abs = ifelse(value != 0, 1, 0)) 

########################################
## Figure S2A: 16S shared OTUs
# instead of doing this for each group, lets create a function that will....
  # 1. calculate the number of unique OTUs in the input (i guess this really doesn't matter for % match but leaving it anyways)
  # 2. calculate the number of OTUs in the sample that match the input
  # 3. calculate a percentage of total OTUs from sample that were found 
process_FMT <- function(FMT_label, df) {
  #create data frame for input
  input_df <- df %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(sampleType %in% c("FMTinput")) %>%
    filter(!pres_abs == 0) %>%
    group_by(seqID) %>%
    mutate(input_unique_species_count = n_distinct(OTU)) %>%  #count unique OTUs in the input
    mutate(species_match_count = sum(OTU %in% unique(OTU)),
           percent_shared = ((species_match_count / input_unique_species_count) * 100)) %>%
    ungroup()
  #create data frame for mice
  mice_df <- df %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(!sampleType %in% c("FMTinput")) %>%
    filter(!pres_abs == 0) %>%
    group_by(seqID) %>%
    mutate(input_unique_species_count = n_distinct(OTU)) %>%  # number of unique OTUs possible in the sample
    mutate(species_match_count = sum(OTU %in% input_df$OTU),  #count matching OTUs between sample and input
           percent_shared = ((species_match_count / input_unique_species_count) * 100)) %>%  #calculate percentage of input OTUs found in the sample
    ungroup() %>%
    select(seqID, percent_shared, FMT_input, sampleType, Experiment,FMT_type) %>%
    distinct()  
  return(mice_df)
}
### run for cecal samples first
yWT_data.c <- process_FMT("yWT", cecal_df)
yWT_spore_data.c <- process_FMT("yWT_spore", cecal_df)
yRag_data.c <- process_FMT("yRag", cecal_df)
rWT_data.c <- process_FMT("rWT", cecal_df)
AMS001_data.c <- process_FMT("AMS001", cecal_df)
AMS005_data.c <- process_FMT("AMS005", cecal_df)
R4F_data.c <- process_FMT("R4_F", cecal_df)
D2A_data.c <- process_FMT("D2_A", cecal_df)
D4A_data.c <- process_FMT("D4_A", cecal_df)
D5A_data.c <- process_FMT("D5_A", cecal_df)

### now fecal 
yWT_data.f <- process_FMT("yWT", fecal_df)
yWT_spore_data.f <- process_FMT("yWT_spore", fecal_df)
yRag_data.f <- process_FMT("yRag", fecal_df)
rWT_data.f <- process_FMT("rWT", fecal_df)
AMS001_data.f <- process_FMT("AMS001", fecal_df)
AMS005_data.f <- process_FMT("AMS005", fecal_df)
R4F_data.f <- process_FMT("R4_F", fecal_df)
D2A_data.f <- process_FMT("D2_A", fecal_df)
D4A_data.f <- process_FMT("D4_A", fecal_df)
D5A_data.f <- process_FMT("D5_A", fecal_df)

# combine the data
cecal_combined <- do.call(rbind, list(yWT_data.c,yWT_spore_data.c,rWT_data.c,yRag_data.c,AMS001_data.c,AMS005_data.c,R4F_data.c,D2A_data.c,D4A_data.c,D5A_data.c)) %>%
  mutate(source="cecal") # adding label for combined plot 
fecal_combined <- do.call(rbind, list(yWT_data.f,yWT_spore_data.f,rWT_data.f,yRag_data.f,AMS001_data.f,AMS005_data.f,R4F_data.f,D2A_data.f,D4A_data.f,D5A_data.f)) %>%
  mutate(source="fecal") #adding labeld for combined plot 

## now plot them together
all_plot <- cecal_combined %>%
  rbind(.,fecal_combined) %>%
  mutate(FMT_input = factor(FMT_input, levels = c( "yWT", "rWT", "yRag","yWT_spore","D2_A","AMS001","AMS005","D4_A", "D5_A", "R4_F"))) %>%
  mutate(source =factor(source, levels=c("fecal","cecal"))) %>%
  ggplot(aes(FMT_input, percent_shared, color=FMT_input)) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(),
          legend.background = element_rect(color="black", fill=NA),
          panel.border = element_rect(fill=NA),
          strip.placement = "outside",
          strip.background = element_blank(),
          strip.text.x = element_text(size = 30)) +
    geom_boxplot(color="black",outlier.shape=NA, linewidth=.3) +
    geom_jitter(width=0.2, alpha=0.6, size = 1.5, na.rm=TRUE) +
  scale_y_continuous(breaks=seq(0, 100, by = 20),limits=c(0,100)) +
  scale_color_manual(name=NULL,
                       labels=c("mFMT1", "mFMT2", "mFMT3","mFMT4","hFMT1", "hFMT2", "hFMT3", "hFMT4", "hFMT5", "hFMT6"),
                       breaks=c( "yWT", "rWT", "yRag","yWT_spore","D2_A","AMS001","AMS005","D4_A", "D5_A", "R4_F"),
                       values=c("#56B29E","#8BBD94","#A6C293","#C1C88C","#EACB2B","#E87700","#E8A419","#E79812","#E8A91B","#ED5300")) +
    labs(x="Microbiome Input", y="% OTUs Shared") +
    facet_grid(~source, scales = "free_x")
#ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Manuscript/Figures/OTUs_shared_ALL.pdf",width=6,height=3.5,dpi=300)

# getting means for text
# mean_test <- cecal_combined %>%
#   rbind(.,fecal_combined) %>%
#   group_by(source,FMT_type) %>%
#   summarize(test = mean(as.numeric(percent_shared)), .groups = "drop")
########################################

########################################
## Figure S2B/C: Relative abundance of the OTUs that did or did not engraft
# calculate relative abundance of otus per sample
otu_rel_abund <- shared %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU) %>%
  rename(seqID = name) %>%
  filter(seqID %in% c(cecal_df$seqID, fecal_df$seqID)) %>%
  left_join(.,meta, by="seqID") %>%
  group_by(seqID) %>% # to get total # of each count for sample
  #filter(!value == 0) %>%
  mutate(rel_abund = (value / sum(value)*100)) %>% # calculate percent rel_abund
  ungroup()

# Sanity check.....make sure relative abundances add to 100 before moving on (this will be 99.something b/c filtered low abund)
otu_rel_abund %>%
  group_by(seqID) %>%
  summarize(total=sum(rel_abund)) 

# now making color df
taxonomy <- read_tsv(file="/Users/home/Library/CloudStorage/Box-Box/Manuscript/ForGit/mouseCDI-SPF-hFMT/Data/16S/allhf_runs/mothurfiles/allhumos6B.final.0.03.cons.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  mutate(Taxonomy=str_replace_all(Taxonomy, "\\(\\d+\\)",""), # remove the (100) by using regular expression d=didgit += 1 or more between two parenthesis
         Taxonomy=str_replace(Taxonomy, ";$","")) %>% # theres an extra ; at the end so need to remove that ($ means match at beginning of string whereas ^B means beginning of string)
  separate(Taxonomy, into=c("Kingdom", "Phylum","Class","Order", "Family","Genus"),
           sep=";") %>%
  filter(OTU %in% c(both_col_plot$OTU, both_plot$OTU)) %>% # NOTE!!!! this has to be done AFTER through steps below...annoying but need to get fewer color options. Just run lines up to part before comnining w/ tax_colors (need to fix this by writing df...)
  filter(OTU %in% c(otu_rel_abund$OTU)) %>%
  mutate(Phylum = replace(Phylum, Phylum == "Acidobacteria", "Other")) %>% 
  mutate(Phylum = replace(Phylum, Phylum == "Chloroflexi", "Other")) %>% 
  mutate(Phylum = replace(Phylum, Phylum == "Deferribacteres", "Other")) %>% 
  mutate(Phylum = replace(Phylum, Phylum == "Synergistetes", "Other")) %>%
  mutate(Phylum = replace(Phylum, Phylum == "Candidatus_Saccharibacteria", "Other")) %>% 
  mutate(Phylum = replace(Phylum, Phylum == "Bacteria_unclassified", "Unclassified"))

color_df <- taxonomy %>%
  select(Phylum, Genus) %>%
  mutate(Phylum = (factor(Phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria", "Firmicutes", "Bacteroidetes", "Tenericutes","Fusobacteria","Other","Unclassified")))) %>%
  arrange(Phylum) %>%
  distinct(Phylum, Genus)
# coloring at the genus level
table(color_df$Phylum)
# define colors 
verruco<-c("hotpink")
actino<-colorRampPalette(c("tan","salmon4", "brown","orangered4","darkred"))(n=15)
pro<-colorRampPalette(c("yellow2","gold","darkgoldenrod3","goldenrod2","orange2","yellow4"))(n=14)
firm<-colorRampPalette(c("midnightblue","blue","dodgerblue4","dodgerblue1","deepskyblue4","deepskyblue1","skyblue3","skyblue","steelblue4","steelblue1","royalblue4","royalblue1","slateblue4","purple3","orchid3"))(n=61)
bac<-colorRampPalette(c("darkgreen","green3","lightgreen","seagreen", "limegreen"))(n=10)
tener <- c("#7FFFD4")
fuso <- c("#FFA070")
other <- c("#D9D9D9")
unclassified<-c("black")
Color<-c(verruco, actino, pro, firm, bac, tener,fuso,other,unclassified)
color_df$gen_color <- Color

tax_colors <- taxonomy %>%
  left_join(.,color_df,by=c("Phylum","Genus")) %>%
  mutate(Phylum = (factor(Phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria", "Firmicutes", "Bacteroidetes", "Tenericutes","Fusobacteria","Other","Unclassified")))) %>%
  arrange(Phylum)

#### Figure S2B: 16S relative abundance of engrafted OTUs
# lets figure out which FMT donor OTUs colonized and what % relabund they explain in the sample
FMT_col <- function(FMT_label, df, otus) {
  # Make a df that for samples with a specified FMT_label, these are all of the possible OTUs found
  input_df <- df %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(sampleType %in% c("FMTinput")) %>%
    filter(!pres_abs == 0)
  # Now take all of the OTUs from a sample for a specified input and keep only the ones that were found in the samples
  otu_df <- otus %>%
    filter(seqID %in% df$seqID) %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(!sampleType %in% c("FMTinput")) %>%
    filter(OTU %in% input_df$OTU) %>%
    group_by(FMT_type, FMT_input,OTU) %>% #
    summarize(mean = mean(rel_abund), .groups = "drop") %>%
    #mutate(mean=mean(rel_abund)) %>% # now we can calculate per OTU, what is its mean relative abundance across this given FMT input
    ungroup()%>%
    select(OTU,FMT_type,FMT_input,mean) %>%
    distinct()
  return(otu_df)
  #return(input_df)
}
yWT_col.c <- FMT_col("yWT", cecal_df,otu_rel_abund)
yWT_spore_col.c <- FMT_col("yWT_spore", cecal_df,otu_rel_abund)
yRag_col.c <- FMT_col("yRag", cecal_df,otu_rel_abund)
rWT_col.c <- FMT_col("rWT", cecal_df,otu_rel_abund)
AMS001_col.c <- FMT_col("AMS001", cecal_df,otu_rel_abund)
AMS005_col.c <- FMT_col("AMS005", cecal_df,otu_rel_abund)
R4F_col.c <- FMT_col("R4_F", cecal_df,otu_rel_abund)
D2A_col.c <- FMT_col("D2_A", cecal_df,otu_rel_abund)
D4A_col.c <- FMT_col("D4_A", cecal_df,otu_rel_abund)
D5A_col.c <- FMT_col("D5_A", cecal_df,otu_rel_abund)

# now fecal
yWT_col.f <- FMT_col("yWT", fecal_df,otu_rel_abund)
yWT_spore_col.f <- FMT_col("yWT_spore", fecal_df,otu_rel_abund)
yRag_col.f <- FMT_col("yRag", fecal_df,otu_rel_abund)
rWT_col.f <- FMT_col("rWT", fecal_df,otu_rel_abund)
AMS001_col.f <- FMT_col("AMS001", fecal_df,otu_rel_abund)
AMS005_col.f <- FMT_col("AMS005", fecal_df,otu_rel_abund)
R4F_col.f <- FMT_col("R4_F", fecal_df,otu_rel_abund)
D2A_col.f <- FMT_col("D2_A", fecal_df,otu_rel_abund)
D4A_col.f <- FMT_col("D4_A", fecal_df,otu_rel_abund)
D5A_col.f <- FMT_col("D5_A", fecal_df,otu_rel_abund)

cecal_col <- do.call(rbind, list(yWT_col.c,yWT_spore_col.c,rWT_col.c,yRag_col.c,AMS001_col.c,AMS005_col.c,R4F_col.c,D2A_col.c,D4A_col.c,D5A_col.c)) %>%
  mutate(source="cecal") 
fecal_col <- do.call(rbind, list(yWT_col.f,yWT_spore_col.f,rWT_col.f,yRag_col.f,AMS001_col.f,AMS005_col.f,R4F_col.f,D2A_col.f,D4A_col.f,D5A_col.f)) %>%
  mutate(source="fecal")

both_col_plot <- cecal_col %>%
  rbind(.,fecal_col) %>%
  left_join(.,tax_colors, by="OTU") %>%
  arrange(desc(mean)) %>%
  mutate(Phylum = (factor(Phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria", "Firmicutes", "Bacteroidetes", "Tenericutes","Acidobacteria","Fusobacteria","Other","Unclassified")))) %>%
  arrange(Phylum) 

plot1 <- both_col_plot %>%
  mutate(Genus = factor(Genus, levels = unique(both_col_plot$Genus))) %>%
  arrange(Genus) %>%
  mutate(FMT_input = factor(FMT_input, levels = c("yWT", "rWT", "yRag","yWT_spore","D2_A","AMS001","AMS005","D4_A", "D5_A", "R4_F"))) %>%
  mutate(source = factor(source, levels = c("fecal","cecal"))) %>%
  ggplot(aes(x=FMT_input, y=mean, fill=factor(Genus))) +
  theme_classic() +
  theme(legend.position="none",
        legend.spacing.x = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank()) +
  guides(fill=guide_legend(nrow=15)) +
  geom_bar(stat="identity") + 
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual(breaks = both_col_plot$Genus, values = c(both_col_plot$gen_color)) +
  labs(x="SampleID", y= "Relative Abundance (%)") +
  facet_grid(~source, scale="free_x", space="free", switch="x") 
#ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Manuscript/Figures/16S_both_col.pdf",width=4,height=3,dpi=300)

# get means for writing
# mean_test <- cecal_col %>%
#   rbind(.,fecal_col) %>%
#   group_by(source,FMT_type,FMT_input) %>%
#   summarize(mean_sum = sum(as.numeric(mean)), .groups = "drop") %>%
#   group_by(source,FMT_type) %>%
#   summarize(test = mean(mean_sum), .groups = "drop")

#### Figure S2C: 16S relative abundance of OTUs from donor that did not engraft (in fecal or cecal samples)
# lets figure out which FMT donor OTUs did not colonize
FMT_noncol <- function(FMT_label, df, otus) {
  # Make a df that for samples with a specified FMT_label, these are all of the possible OTUs found
  mice_df <- df %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(!sampleType %in% c("FMTinput")) %>%
    filter(!pres_abs == 0)
  # Now take all of the OTUs from a specified input and remove the ones that were found in the samples leaving only OTUs that did not "colonize"
  input_df <- df %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(sampleType %in% c("FMTinput")) %>%
    filter(!pres_abs == 0) %>%
    filter(!OTU %in% mice_df$OTU) %>%
    select(OTU,sampleID) %>%
    distinct()
  # And finally lets pull the relative abundance of those OTUs from the donor sample
  otu_df <- otus %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(sampleType %in% c("FMTinput")) %>%
    filter(sampleID %in% input_df$sampleID) %>%
    filter(OTU %in% input_df$OTU) %>%
    group_by(FMT_type,FMT_input,OTU) %>% #
    summarize(mean = mean(rel_abund), .groups = "drop") %>%
    ungroup()%>%
    select(OTU,FMT_type,FMT_input,mean) %>%
    distinct()
  return(otu_df)
  #return(input_df)
}
yWT_noncol.c <- FMT_noncol("yWT", cecal_df,otu_rel_abund)
yWT_spore_noncol.c <- FMT_noncol("yWT_spore", cecal_df,otu_rel_abund)
yRag_noncol.c <- FMT_noncol("yRag", cecal_df,otu_rel_abund)
rWT_noncol.c <- FMT_noncol("rWT", cecal_df,otu_rel_abund)
AMS001_noncol.c <- FMT_noncol("AMS001", cecal_df,otu_rel_abund)
AMS005_noncol.c <- FMT_noncol("AMS005", cecal_df,otu_rel_abund)
R4F_noncol.c <- FMT_noncol("R4_F", cecal_df,otu_rel_abund)
D2A_noncol.c <- FMT_noncol("D2_A", cecal_df,otu_rel_abund)
D4A_noncol.c <- FMT_noncol("D4_A", cecal_df,otu_rel_abund)
D5A_noncol.c <- FMT_noncol("D5_A", cecal_df,otu_rel_abund)

### and fecal
yWT_noncol.f <- FMT_noncol("yWT", fecal_df,otu_rel_abund)
yWT_spore_noncol.f <- FMT_noncol("yWT_spore", fecal_df,otu_rel_abund)
yRag_noncol.f <- FMT_noncol("yRag", fecal_df,otu_rel_abund)
rWT_noncol.f <- FMT_noncol("rWT", fecal_df,otu_rel_abund)
AMS001_noncol.f <- FMT_noncol("AMS001", fecal_df,otu_rel_abund)
AMS005_noncol.f <- FMT_noncol("AMS005", fecal_df,otu_rel_abund)
R4F_noncol.f <- FMT_noncol("R4_F", fecal_df,otu_rel_abund)
D2A_noncol.f <- FMT_noncol("D2_A", fecal_df,otu_rel_abund)
D4A_noncol.f <- FMT_noncol("D4_A", fecal_df,otu_rel_abund)
D5A_noncol.f <- FMT_noncol("D5_A", fecal_df,otu_rel_abund)

# combine data and then plot
fecal_noncol <- do.call(rbind, list(yWT_noncol.f,yWT_spore_noncol.f,rWT_noncol.f,yRag_noncol.f,AMS001_noncol.f,AMS005_noncol.f,R4F_noncol.f,D2A_noncol.f,D4A_noncol.f,D5A_noncol.f)) %>%
  mutate(source="fecal")
cecal_noncol <- do.call(rbind, list(yWT_noncol.c,yWT_spore_noncol.c,rWT_noncol.c,yRag_noncol.c,AMS001_noncol.c,AMS005_noncol.c,R4F_noncol.c,D2A_noncol.c,D4A_noncol.c,D5A_noncol.c)) %>%
  mutate(source="cecal") 
  
both_plot <- fecal_noncol %>%
  rbind(.,cecal_noncol) %>%
  left_join(.,tax_colors, by="OTU") %>%
  arrange(desc(mean)) %>%
  mutate(Phylum = (factor(Phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria", "Firmicutes", "Bacteroidetes", "Tenericutes","Acidobacteria","Fusobacteria","Other","Unclassified")))) %>%
  arrange(Phylum) 

plot <- both_plot %>%
  mutate(Genus = factor(Genus, levels = unique(both_plot$Genus))) %>%
  arrange(Genus) %>%
  mutate(FMT_input = factor(FMT_input, levels = c("yWT", "rWT", "yRag","yWT_spore","D2_A","AMS001","AMS005","D4_A", "D5_A", "R4_F"))) %>%
  mutate(source = factor(source, levels = c("fecal","cecal"))) %>%
  ggplot(aes(x=FMT_input, y=mean, fill=factor(Genus))) +
  theme_classic() +
  theme(legend.position="none",
        legend.spacing.x = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks.x = element_blank()) +
  guides(fill=guide_legend(nrow=15)) +
  geom_bar(stat="identity") + 
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual(breaks = both_plot$Genus, values = c(both_plot$gen_color)) +
  labs(x="SampleID", y= "Relative Abundance (%)") +
  facet_grid(~source, scale="free_x", space="free", switch="x") 
ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Manuscript/Figures/16S_both_noncol.pdf",width=4,height=3,dpi=300)
#ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Manuscript/Figures/16S_both_noncol_legend.pdf",width=30,height=30,dpi=300)

# # calculate means for writing
# mean_test <- cecal_noncol %>%
#   rbind(.,fecal_noncol) %>%
#   group_by(source,FMT_type,FMT_input) %>%
#   summarize(mean_sum = sum(as.numeric(mean)), .groups = "drop") %>%
#   group_by(FMT_type) %>%
#   summarize(test = mean(mean_sum), .groups = "drop")
```

