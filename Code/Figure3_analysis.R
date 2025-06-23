#### Figure 3
### 8.6.24
### S. Millard

### Figures:
- Fig. 3A: species level beta diversity
- Fig. 3B: bray-curtis distances
- Fig. 3C: inverse-simpson alpha diversity
- Fig. 3D: relative abundance barchart
- Fig. 3E: Maaslin2 heatmaps


# load packages
library(RColorBrewer)
library(gplots)
library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(phyloseq)
library(glue)
library(tibble)
library(rstatix)
library(Maaslin2)
library(pheatmap)
library(ggpubr)

# set working directory 
setwd("/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript")

# read in metadata and clean it up (there is lots of unnecessary info)
meta <- read.table("./data/manuscript_meta.txt", header=T, sep='\t', quote='', na.strings= c("NA"))

# read in tax data and clean it up 
tax_counts <- read.table("./data/metagenomic_tax_relab.tsv", header=T, sep='\t', quote='', na.strings= c("NA"))

###############################################################################################################################
### Figure 3A: NMDS based on species counts
# change healthy cecal mouse sample type to cecal
beta_meta <- meta %>%
  mutate(type = replace(type, type== "mouse_survey", "cecal"))

# format count data for vegan input
beta_counts <- tax_counts %>%
  separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species","strain"), sep="\\|", fill="right",remove=FALSE) %>%
  mutate_at(vars(!matches("clade_name")),~str_remove(.,".\\_\\_")) %>%
  filter(!is.na(species)) %>%
  filter(is.na(strain)) %>%
  mutate(relative_abundance = as.numeric(relative_abundance)) %>%
  mutate(relabund = 100*(relative_abundance)) %>% # taking out of percent form
  select(species, relative_abundance, mgID) %>%
  mutate(across(where(is.numeric), round)) %>%
  pivot_wider(names_from="species", values_from="relative_abundance", values_fill=0) %>%
  column_to_rownames(var="mgID") %>%
  as.data.frame() %>%
  filter(rownames(.) %in% beta_meta$mgID) %>%
  select(-where(~sum(.) == 0))

# using vegan, calculate bray-curtis dissimilarity based on species counts
set.seed(19930514)
SPF_tax_bray <- avgdist(beta_counts, dmethod="bray", sample=79)
bray_matrix <- SPF_tax_bray%>%
  as.matrix(labels=TRUE)

### calculate nmds from bray
set.seed(2)
tax_nmds <- metaMDS(SPF_tax_bray, trymax=50)
tax_nmds # if you want to check the stress 

nmds_positions <- scores(tax_nmds)

fig3A <- nmds_positions %>%
  as_tibble(rownames="mgID") %>%
  inner_join(., beta_meta, by="mgID") %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=sample_group, shape=type)) +
  theme_classic() +
  geom_point(size = 4, alpha=0.8) + 
  theme(legend.text = element_text(),
        legend.position.inside = c(.1, .25),
        legend.key.height = unit(20,"pt"),
        legend.background = element_rect(color="black", fill=NA),
        axis.ticks.length=unit(.25, "cm"),
        panel.border = element_rect(fill=NA)) +
  scale_shape_manual(name="Sample Type",
                     breaks=c("fecal", "cecal","FMT_input"),
                     values = c(19, 17,8)) +
  scale_color_manual(name="FMT Type",
                     breaks=c("mouse_survey","pre-fmt","yWT", "D2_A","AMS001", "AMS005", "none"),
                     labels=c("Healthy","pre-fmt","yWT", "hFMT1", "hFMT2", "hFMT3","no FMT"),  
                     values=c('#3B99B1','#E2E2E2','#56B29E','#EACB2B','#E87700','#E8A419','#F5191C')) +
  ggtitle("NMDS of Bray Curtis Dissimilarity Based on Species Counts")
ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/figures/3/3A_species_bray_nmds.pdf",width=5.2,height=4,dpi=300)

## run stats
meta_stats <- meta %>%
  mutate(type = replace(type, type== "mouse_survey", "cecal")) %>%
  filter(!sample_group %in% c("pre-fmt")) %>%
  filter(!FMT_input %in% c("Healthy")) %>%
  filter(!type %in% c("FMT_input"))

stat_positions <- nmds_positions %>%
  as.data.frame() %>%
  filter(rownames(.) %in% meta_stats$mgID)
stat_positions[stat_positions < 0] <- 0.00001 # get rid of negative values

SPF_tax_bray_stats <- stat_positions %>%
  rownames_to_column("mgID") %>%
  inner_join(., meta_stats, by="mgID") %>%
  pairwiseAdonis::pairwise.adonis2(stat_positions ~ FMT_input, data = ., nperm = 1000)

# $AMS001_vs_yWT
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1   3.0951 0.53599 20.793  0.001 ***
#   Residual  18   2.6794 0.46401                  
# Total     19   5.7745 1.00000                  
# ---
# $AMS001_vs_none
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1   4.1890 0.61673 33.791  0.001 ***
#   Residual  21   2.6033 0.38327                  
# Total     22   6.7923 1.00000                  
# ---
# $AMS005_vs_yWT
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1   4.9212 0.96924 567.15  0.001 ***
#   Residual  18   0.1562 0.03076                  
# Total     19   5.0774 1.00000                  
# ---
# $AMS005_vs_D2_A
# Df SumOfSqs      R2      F Pr(>F)  
# FMT_input  1   1.2963 0.29233 8.2616  0.012 *
#   Residual  20   3.1381 0.70767                
# Total     21   4.4344 1.00000                
# ---
# $AMS005_vs_none
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1   5.6149 0.98594 1472.5  0.001 ***
#   Residual  21   0.0801 0.01406                  
# Total     22   5.6950 1.00000                  
# ---
# $yWT_vs_D2_A
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1   3.0376 0.47973 18.442  0.001 ***
#   Residual  20   3.2943 0.52027                  
# Total     21   6.3319 1.00000                  
# ---
# $yWT_vs_none
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1  1.44480 0.85946 128.42  0.001 ***
#   Residual  21  0.23626 0.14054                  
# Total     22  1.68106 1.00000                  
# ---
# $D2_A_vs_none
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1   4.5674 0.58665 32.643  0.001 ***
#   Residual  23   3.2182 0.41335                  
# Total     24   7.7856 1.00000                  
# ---
###############################################################################################################################

###############################################################################################################################
### Figure 3B: Bray-cutis distances
mat <- bray_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "RowNames")
#write.table(mat, file = "./data/tax_matrix.txt", sep = "\t", row.names = FALSE, col.names = TRUE, na = "NA")

# modify distance matrix; merge with meta data; filter out relevant comparisons
tax.dist <- mat %>%
  dplyr::rename(sampleID1 = 1) %>%
  gather("sampleID2", "distance", -sampleID1) %>%
  inner_join(., meta, by = c("sampleID1" = "mgID")) %>%
  inner_join(., meta, by = c("sampleID2" = "mgID")) %>%
  select(-type.x, -type.y) %>%
  # merge both sampleIDs with metadata
  rename_with(~ c("FMT_input1", "FMT_type1", "FMT_input2", "FMT_type2"), 
              all_of(c("FMT_input.x", "FMT_type.x", "FMT_input.y", "FMT_type.y"))
  ) %>%
  # get rid of self-comparisons (which are 0)
  filter(!sampleID1 == sampleID2) %>%	
  # define comparison
  mutate(comparison = case_when(
    FMT_type1 == "hFMT" & FMT_type2 == "mFMT" | FMT_type2 == "hFMT" & FMT_type1 == "mFMT" ~ "mFMT-hFMT",
    FMT_type1 == "noFMT" & FMT_type2 == "mFMT" | FMT_type2 == "noFMT" & FMT_type1 == "mFMT" ~ "mFMT-noFMT",
    FMT_type1 == "hFMT" & FMT_type2 == "noFMT" | FMT_type2 == "hFMT" & FMT_type1 == "noFMT" ~ "hFMT-noFMT",
    FMT_type1 == "mFMT" & FMT_type1 == FMT_type2 ~ "within-mFMT",
    FMT_type1 == "hFMT" & FMT_type1 == FMT_type2 ~ "within-hFMT",
    FMT_type1 == "noFMT" & FMT_type1 == FMT_type2 ~ "within-noFMT",
    TRUE ~ NA_character_),
    comparison = factor(comparison, levels=c("mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT", "within-mFMT", "within-hFMT", "within-noFMT"))
  ) %>% #for now, eliminating healthy...	
  filter(!is.na(comparison)) 

# graph results	
# all distances:	
tax.dist %>% ggplot(aes(x=comparison, y=distance)) +
  geom_jitter(width=0.25, color="grey", alpha=0.6, size = 0.5) +
  stat_summary(fun.data=median_hilow, color="red", size=1,
               fun.args = list(conf.int=0.50)) +
  labs(x=NULL, y="__ distances") +
  scale_x_discrete(breaks=c("mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT", "within-mFMT", "within-hFMT", "within-noFMT"), 
                   labels=c("mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT", "within-mFMT", "within-hFMT", "within-noFMT")) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("FigX_bray.dist.dotplot_mFMT.png", width=4, height=4)

# just cross-group distances:	
tax.dist %>% 
  filter(comparison %in% c("mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT")) %>%
  ggplot(aes(x=comparison, y=distance)) +
  geom_jitter(width=0.25, color="grey", alpha=0.6, size = 0.5) +
  stat_summary(fun.data=median_hilow, color="red", size=1,
               fun.args = list(conf.int=0.50)) +
  labs(x=NULL, y="__ distances") +
  scale_x_discrete(breaks=c("mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT", "within-mFMT", "within-hFMT", "within-noFMT"), 
                   labels=c("mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT", "within-mFMT", "within-hFMT", "within-noFMT")) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.2)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave("FigX_bray.dist.dotplot_mFMT.png", width=4, height=4)

## stats:
tax.dist %>% subset(., comparison %in% c("mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT")) %>%
  dunnTest(distance ~ comparison, data = .)
# Comparison          Z      P.unadj        P.adj
# 1  hFMT-noFMT - mFMT-hFMT -14.436669 3.041930e-47 6.083861e-47
# 2 hFMT-noFMT - mFMT-noFMT   8.475416 2.342403e-17 2.342403e-17
# 3  mFMT-hFMT - mFMT-noFMT  20.352373 4.423354e-92 1.327006e-91
# phrased as: "hFMT samples were significantly more dissimilar to noFMT samples than mFMT samples"

# look at means	
detach(package:plyr)
tax.dist %>% dplyr::group_by(comparison) %>% summarize(mean(distance))	
###############################################################################################################################

###############################################################################################################################
### Figure 3C: Alpha diversity
alpha <- tax_counts %>%
  filter(str_detect(clade_name,"\\|s_")) %>%
  separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species","strain"), sep="\\|",fill="right") %>%
  filter(!is.na(species)) %>%
  filter(is.na(strain)) %>%
  select(species, relative_abundance, mgID) %>%
  group_by(mgID) %>%
  summarize(sobs=specnumber(relative_abundance), # observed richness
            shannon = diversity(relative_abundance, index="shannon"),
            simpson = diversity(relative_abundance, index="simpson"),
            invisimpson = diversity(relative_abundance, index="invsimpson"),
            evenness = diversity(relative_abundance) / log(sobs), 
            n=sum(relative_abundance)) %>%
  ungroup() %>%
  merge(., meta, by="mgID") %>%
  filter(type %in% c("cecal","mouse_survey")) %>%
  mutate(FMT_type = replace(FMT_type, FMT_input== "Healthy", "Healthy")) %>%
  group_by(FMT_type) %>%
  mutate(ymin = min(invisimpson), ymax = max(invisimpson)) %>%
  ungroup()

## stats
kruskal.test(x=alpha$invisimpson, g=alpha$FMT_type)
pairwise.wilcox.test(x=alpha$invisimpson, g=alpha$FMT_type)

# using the rstatix package is the easiest way run stats and add to ggplot
kw_invisimpson <-alpha %>%
  select(invisimpson,FMT_type,mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(invisimpson ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F) %>%
  add_xy_position(x="FMT_type")

Fig3C <- alpha %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  ggplot(aes(x=FMT_type, y=invisimpson, color=FMT_type)) +
  theme_classic() +
  scale_y_continuous(breaks=c(10,20,30,40)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        #legend.text = element_text(),
        #legend.position.inside = c(0.8, 0.8),
        legend.position="none",
        #legend.background = element_rect(color="black", fill=NA),
        panel.border = element_rect(fill=NA)) +
  geom_jitter(width=0.2, alpha=0.6, size = 1, na.rm=TRUE) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.5, linewidth=.5) +
  scale_color_manual(name=NULL,
                     breaks=c("Healthy","mFMT","hFMT","noFMT"),
                     labels=c("Healthy","mFMT","hFMT","noFMT"),
                     values=c('#3B99B1','#56B29E','#EACB2B','#F5191C')) +
  labs(x="FMT Type", y="Inverse Simpson") +
  stat_pvalue_manual(filter(kw_invisimpson,p.adj<0.05), 
                     label = "p.adj.signif",
                     tip.length=0.02,
                     step.increase=0.02,
                     inherit.aes=F) +
  stat_summary(fun=median, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.1), show.legend=FALSE)
ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/figures/3/3C_FMT_input_invisimpson_stats.pdf",width=2.5,height=4,dpi=300)

###############################################################################################################################

###############################################################################################################################
### Figure 3D: Relative abundane bar-charts
# read in df to assign colors at genus and phylum levels
tax_colors <- read.table("./data/tax_colors.tsv", header=T, sep='\t', quote='', na.strings= c("NA"), comment.char = "")
tax_colors$phylum <- factor(tax_colors$phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria", "Firmicutes", "Bacteroidetes","Tenericutes","Unclassified","Other"))

# clean up taxonomic classification 
species_relab <- tax_counts %>%
  separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species","strain"), sep="\\|", fill="right",remove=FALSE) %>%
  mutate_at(vars(!matches("clade_name")),~str_remove(.,".\\_\\_")) %>%
  mutate(species = if_else(is.na(strain), species, str_glue("{species}_{strain}"))) %>%
  filter(!is.na(species)) %>%
  filter(is.na(strain)) %>%
  select(-c(kingdom,-strain)) %>%
  mutate(relative_abundance = as.numeric(relative_abundance)) %>%
  mutate(phylum = replace(phylum, phylum == "Candidatus_Saccharibacteria", "Other")) %>% 
  mutate(phylum = replace(phylum, phylum == "Synergistetes", "Other")) %>% 
  mutate(phylum = replace(phylum, phylum == "Bacteria_unclassified", "Unclassified")) 

# Metaphlann output provides counts as long format which means zero counts are not introudced, this is okay for looking at sample relative abundance but will be a problem for calculating mean so need to make df that includes 0 counts for species
zero_counts <- species_relab %>%
  select(species,relative_abundance,mgID) %>%
  mutate(relative_abundance=as.numeric(relative_abundance)) %>%
  complete(mgID, species, fill = list(relative_abundance = 0))

# now lets take the classification from the original file to add back on the relab file that includes 0 counts
class <- species_relab %>%
  select(phylum,family,species,genus) %>%
  distinct()

# and merge them back together, then re-add metadata 
species_relab_zero <- zero_counts %>%
  merge(.,class,by="species") %>% 
  merge(., meta, by="mgID") %>% 
  filter(type %in% c("cecal","mouse_survey")) %>%
  mutate(FMT_type = replace(FMT_type, FMT_input== "Healthy", "Healthy"))

# can check and make sure relative abundances add to 100 as a check 
species_relab_zero %>%
  group_by(mgID) %>%
  summarize(total=sum(relative_abundance))

# calculte mean relab by FMT_type and plot 
Fig3D <- species_relab_zero %>%
  group_by(FMT_type,phylum, genus,species) %>%
  summarize(mean_rel_abund = mean(relative_abundance), .groups = "drop") %>%
  arrange(desc(mean_rel_abund)) %>%
  ungroup() %>%
  mutate(genus = factor(genus, levels = tax_colors$genus)) %>%
  arrange(genus) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  ggplot(aes(x=FMT_type, y=mean_rel_abund, fill=factor(genus), group = factor(phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","Tenericutes","Unclassified","Other")))) + 
  theme_classic() +
  theme(legend.position="none",
        legend.spacing.x = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        axis.text.x = element_text(angle = 45),
        axis.text.y = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside") +
  guides(fill=guide_legend(nrow=25)) +
  geom_bar(position= "fill", stat="identity") + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(breaks = c(tax_colors$genus), values = c(tax_colors$gen_color)) +
  labs(x="FMT Input", y= "Relative Abundance (%)") 
ggsave(filename="./figures/3/3D_FMTinput_tax_relab.pdf",width=3.5,height=6,dpi=300)
###############################################################################################################################

###############################################################################################################################
### Figure S3: Cecal relative abundance bar-charts for each individual mouse 
S3 <- species_relab_zero %>%
  mutate(FMT_input = factor(FMT_input, levels = c("Healthy","yWT","D2_A","AMS001","AMS005","none"))) %>%
  mutate(genus = factor(genus, levels = tax_colors$genus)) %>%
  arrange(genus) %>%
  ggplot(aes(x=mgID, y=relative_abundance, fill=factor(genus), group = factor(phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","Tenericutes","Unclassified","Other")))) + 
  theme_classic() +
  theme(legend.position="none",
        legend.spacing.x = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill=guide_legend(nrow=25)) +
  geom_bar(position= "fill", stat="identity") + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(breaks = c(tax_colors$genus), values = c(tax_colors$gen_color)) +
  labs(x="SampleID", y= "Relative Abundance (%)") +
  facet_grid(~FMT_input, scale="free_x", space="free", switch="x")
ggsave(filename="./figures/3/S3_individual_tax_relab.pdf",width=8,height=4,dpi=300)
###############################################################################################################################

###############################################################################################################################
### Figure 3E: Maaslin2 heatmaps
maas_meta <- meta %>%
  filter(type %in% c("cecal")) %>%
  arrange(mgID) %>%
  column_to_rownames("mgID") %>%
  arrange(factor(Colonization, levels = c("no","yes")))

maas_counts <- tax_counts %>%
  separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species","strain"), sep="\\|", fill="right",remove=FALSE) %>%
  mutate_at(vars(!matches("clade_name")),~str_remove(.,".\\_\\_")) %>%
  filter(!is.na(species)) %>%
  filter(is.na(strain)) %>%
  mutate(relative_abundance = as.numeric(relative_abundance)) %>%
  select(species, relative_abundance, mgID) %>%
  pivot_wider(names_from="species", values_from="relative_abundance", values_fill=0) %>%
  arrange(mgID) %>%
  column_to_rownames(var="mgID") %>%
  as.data.frame() %>%
  filter(rownames(.) %in% rownames(maas_meta)) %>%
  select(-where(~sum(.) == 0))

# Compared to mFMT as baseline, what differs in abundance?
Maaslin2(input_data = maas_counts, input_metadata = maas_meta, output = "./maaslin_output/mouse_vs_others/",
         normalization="NONE", transform="LOG", analysis_method ="LM", correction="BH",
         fixed_effects = c("FMT_type"),reference = 'FMT_type,mFMT', plot_heatmap=FALSE, plot_scatter=FALSE)

# Read in results
mvo_res <- read.table("./maaslin_output/mouse_vs_others/significant_results.tsv", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  filter(qval < 0.01) %>%
  mutate(diff_express = coef) %>%
  mutate(diff_express = ifelse(diff_express > 0, "up", "down")) 

# read in df to assign colors at genus and phylum levels
tax_colors <- read.table("./data/tax_colors.tsv", header=T, sep='\t', quote='', na.strings= c("NA"), comment.char = "")
tax_colors$phylum <- factor(tax_colors$phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria", "Firmicutes", "Bacteroidetes","Tenericutes","Unclassified","Other"))

# select only the taxa identified as being significantly different 
mvo_counts <-  tax_counts %>%
  separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species","strain"), sep="\\|", fill="right",remove=FALSE) %>%
  mutate_at(vars(!matches("clade_name")),~str_remove(.,".\\_\\_")) %>%
  filter(!is.na(species)) %>%
  filter(is.na(strain)) %>%
  mutate(relative_abundance = as.numeric(relative_abundance)) %>%
  mutate(relabund = 100*(relative_abundance)) %>% # taking out of percent form
  mutate(spec_gen = paste(species, genus, sep = "+")) %>%
  select(spec_gen, relabund, mgID) %>%
  pivot_wider(names_from="spec_gen", values_from="relabund", values_fill=0) %>%
  gather(key="spec_gen", value="count", -mgID) %>%
  left_join(., meta, by="mgID") %>%
  filter(type %in% c("cecal")) %>%
  group_by(FMT_type, spec_gen) %>%
  summarize(mean=log10(mean(count+1))) %>%
  spread(key='FMT_type', value="mean") %>%
  separate(spec_gen, c("species", "genus"), sep="\\+", fill="right",remove=TRUE) %>%
  filter(species %in% mvo_res$feature) %>%
  left_join(.,tax_colors, by="genus") %>%
  filter(!grepl("^GGB", species) & !grepl("Bacteria_unclassified", genus)) %>% # will put full in supplementary
  column_to_rownames("species") %>%
  arrange(factor(phylum, levels = c("Verrucomicrobia","Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","Bacteria_unclassified"))) %>%arrange(factor(phylum, levels = c("Verrucomicrobia","Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","Bacteria_unclassified")))


mvo_an_col <- list(phylum = setNames(as.character(unique(mvo_counts$gen_color)), unique(mvo_counts$genus)))
mvo_phylum <- data.frame(phylum=mvo_counts$genus)
rownames(mvo_phylum) <- row.names(mvo_counts)

pheatmap(mvo_counts[c(3,2,4)],
         na_col = "white",
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=10)(100),
         border_color="grey20",
         cellheight = 8,
         cellwidth = 20,
         fontsize = 8,
         filename="./figures/3/3E_mouse_v_others_NAMED2.pdf",
         annotation_row = mvo_phylum,
         annotation_colors = mvo_an_col,
         angle_col = "45",
         annotation_legend = TRUE
)
dev.off()

### Figure S4: Maaslin2 heatmap for clear vs colonized comparison
Maaslin2(input_data = maas_counts, input_metadata = maas_meta, output = "./maaslin_output/colonization", 
         normalization="NONE", transform = "LOG",analysis_method ="LM", correction="BH",
         fixed_effects = c("Colonization"),plot_heatmap=FALSE, plot_scatter=FALSE) 

cvc_res <- read.table("./maaslin_output/colonization/significant_results.tsv", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  filter(qval < 0.01) %>%
  mutate(diff_express = coef) %>%
  mutate(diff_express = ifelse(diff_express > 0, "up", "down")) 

cvc_counts <-  tax_counts %>%
  separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species","strain"), sep="\\|", fill="right",remove=FALSE) %>%
  mutate_at(vars(!matches("clade_name")),~str_remove(.,".\\_\\_")) %>%
  filter(!is.na(species)) %>%
  filter(is.na(strain)) %>%
  mutate(relative_abundance = as.numeric(relative_abundance)) %>%
  mutate(relabund = 100*(relative_abundance)) %>% # taking out of percent form
  mutate(spec_gen = paste(species, genus, sep = "+")) %>%
  select(spec_gen, relabund, mgID) %>%
  pivot_wider(names_from="spec_gen", values_from="relabund", values_fill=0) %>%
  gather(key="spec_gen", value="count", -mgID) %>%
  left_join(., meta, by="mgID") %>%
  filter(type %in% c("cecal")) %>%
  group_by(FMT_type, spec_gen) %>%
  summarize(mean=log10(mean(count+1))) %>%
  spread(key='FMT_type', value="mean") %>%
  separate(spec_gen, c("species", "genus"), sep="\\+", fill="right",remove=TRUE) %>%
  filter(species %in% cvc_res$feature) %>%
  left_join(.,tax_colors, by="genus") %>%
  column_to_rownames("species") %>%
  arrange(factor(phylum, levels = c("Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","Tenericutes","Bacteria_unclassified")))

cvc_an_col <- list(phylum = setNames(as.character(unique(cvc_counts$phyl_color)), unique(cvc_counts$phylum)))
cvc_phylum <- data.frame(phylum=cvc_counts$phylum)
rownames(cvc_phylum) <- row.names(cvc_counts)

pheatmap(cvc_counts[c(3,2,4)],
         na_col = "white",
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=10)(100),
         border_color="grey20",
         cellheight = 8,
         cellwidth = 20,
         fontsize = 8,
         filename="./figures/3/S3_clear.v.col_ALL.pdf",
         annotation_row = cvc_phylum,
         annotation_colors = cvc_an_col,
         angle_col = "45",
         annotation_legend = TRUE
)
dev.off()

### Figure S5: Maaslin2 heatmap that includes all taxa (expansion of 3e)
mvo_s_counts <-  tax_counts %>%
  separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species","strain"), sep="\\|", fill="right",remove=FALSE) %>%
  mutate_at(vars(!matches("clade_name")),~str_remove(.,".\\_\\_")) %>%
  filter(!is.na(species)) %>%
  filter(is.na(strain)) %>%
  mutate(relative_abundance = as.numeric(relative_abundance)) %>%
  mutate(relabund = 100*(relative_abundance)) %>% # taking out of percent form
  mutate(spec_gen = paste(species, genus, sep = "+")) %>%
  select(spec_gen, relabund, mgID) %>%
  pivot_wider(names_from="spec_gen", values_from="relabund", values_fill=0) %>%
  gather(key="spec_gen", value="count", -mgID) %>%
  left_join(., meta, by="mgID") %>%
  filter(type %in% c("cecal")) %>%
  group_by(FMT_type, spec_gen) %>%
  summarize(mean=log10(mean(count+1))) %>%
  spread(key='FMT_type', value="mean") %>%
  separate(spec_gen, c("species", "genus"), sep="\\+", fill="right",remove=TRUE) %>%
  filter(species %in% mvo_res$feature) %>%
  left_join(.,tax_colors, by="genus") %>%
  column_to_rownames("species") %>%
  arrange(factor(phylum, levels = c("Verrucomicrobia","Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","Tenericutes","Bacteria_unclassified"))) %>%arrange(factor(phylum, levels = c("Verrucomicrobia","Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","Bacteria_unclassified")))

mvo_s_col <- list(phylum = setNames(as.character(unique(mvo_s_counts$phyl_color)), unique(mvo_s_counts$phylum)))
mvo_s_phylum <- data.frame(phylum=mvo_s_counts$phylum)
rownames(mvo_s_phylum) <- row.names(mvo_s_counts)

pheatmap(mvo_s_counts[c(3,2,4)],
         na_col = "white",
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=10)(100),
         border_color="grey20",
         cellheight = 8,
         cellwidth = 20,
         fontsize = 8,
         filename="./figures/3/S4_mouse.v.others_ALL.pdf",
         annotation_row = mvo_s_phylum,
         annotation_colors = mvo_s_col,
         angle_col = "45",
         annotation_legend = TRUE
)
dev.off()
###############################################################################################################################

###############################################################################################################################
#### Figure S6: Metagenomic species engraftment

# read in metadata and clean it up (there is lots of unnecessary info)
meta <- read.table("/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/data/manuscript_meta.txt", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  filter(type %in% c("cecal","FMT_input")) %>%
  filter(!sample_group %in% c("none")) %>%
  filter(!Experiment %in% c("GF2")) %>%
  select(mgID,sample_group, Experiment,FMT_input,type,FMT_type)

# read in tax data and clean it up 
tax_counts <- read.table("/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/data/metagenomic_tax_relab.tsv", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species","strain"), sep="\\|", fill="right",remove=FALSE) %>%
  mutate_at(vars(!matches("clade_name")),~str_remove(.,".\\_\\_")) %>%
  mutate(species = if_else(is.na(strain), species, str_glue("{species}_{strain}"))) %>%
  filter(!is.na(species)) %>%
  filter(is.na(strain)) %>%
  select(-strain) %>%
  mutate(phylum = replace(phylum, phylum == "Candidatus_Saccharibacteria", "Other")) %>% 
  mutate(phylum = replace(phylum, phylum == "Synergistetes", "Other")) %>% 
  mutate(phylum = replace(phylum, phylum == "Bacteria_unclassified", "Unclassified")) %>% 
  filter(mgID %in% meta$mgID) 

# now join them 
df2 <- tax_counts %>%
  mutate(relative_abundance = as.numeric(relative_abundance)) %>%
  mutate(pres_abs = ifelse(relative_abundance != 0, 1, 0)) %>%
  select(species, pres_abs, mgID) %>%
  left_join(.,meta,by="mgID")

#### Figure S6A: Metagenomic % species shared
### instead of doing this for each group, lets create a function that will 
# 1. calculate the number of unique species in the input (i guess this really doesn't matter for % match but leaving it anyways)
# 2. calculate the number of species in the sample that match the input
# 3. calculate a percentage of total speciess from sample that were found 
process_mg_FMT <- function(FMT_label, df) {
  #create data frame for input
  input_df <- df %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(type %in% c("FMT_input")) %>%
    group_by(mgID) %>%
    mutate(input_unique_species_count = n_distinct(species)) %>%  #count unique species in the input
    mutate(species_match_count = sum(species %in% unique(species)),
           percent_shared = ((species_match_count / input_unique_species_count) * 100)) %>%
    ungroup()
  
  #create data frame for mice
  mice_df <- df %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(!type %in% c("FMT_input")) %>%
    group_by(mgID) %>%
    mutate(input_unique_species_count = n_distinct(species)) %>%  # number of unique species in the sample
    mutate(species_match_count = sum(species %in% input_df$species),  #count matching OTUs between sample and input
           percent_shared = ((species_match_count / input_unique_species_count) * 100)) %>%  #calculate percentage of input OTUs found in the sample
    ungroup() %>%
    select(mgID, percent_shared, FMT_input, type, Experiment,FMT_type) %>%
    distinct()  
  return(mice_df)
}
yWT_mg <- process_mg_FMT("yWT", df2)
AMS001_mg <- process_mg_FMT("AMS001", df2)
AMS005_mg <- process_mg_FMT("AMS005", df2)
D2A_mg <- process_mg_FMT("D2_A", df2)

combined <- do.call(rbind, list(yWT_mg,AMS001_mg,AMS005_mg,D2A_mg))

plot <- combined %>%
  group_by(FMT_input) %>%
  mutate(ymin = min(percent_shared), ymax = max(percent_shared)) %>%
  ungroup() %>%
  mutate(FMT_input = factor(FMT_input, levels = c( "yWT","D2_A","AMS001", "AMS005"))) %>%
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
                       breaks=c("yWT","D2_A","AMS001","AMS005"),
                       values=c("#56B29E","#EACB2B","#E87700","#E8A419")) + 
    labs(x="Microbiome Input", y="% Species Shared") 
ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Manuscript/Figures/Species_shared.pdf",width=3.5,height=3,dpi=300)

# for writing:
# mean_test <- combined %>%
#   group_by(FMT_type) %>%
#   summarize(test = mean(percent_shared), .groups = "drop")

#### Figure S6B: Metagenomic relative abundance of colonized species
## now lets look at the donor species that did colonize and see what their relative abundance looks like
mg_tax_colors <- read_tsv(file="/Users/home/Library/CloudStorage/Box-Box/Manuscript/ForGit/mouseCDI-SPF-hFMT/Data/metagenomics/tax_colors.tsv")

# need add metadata to tax_counts file
tax_counts2 <- tax_counts %>%
  left_join(.,meta,by="mgID")

### now lets figure out which FMT donor specuess colonized and what % relabund they explain in the sample
FMT_mg_col <- function(FMT_label, df) {
  # Make a df that for samples with a specified FMT_label, these are all of the possible OTUs found
  input_df <- df %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(type %in% c("FMT_input"))
  # Now take all of the OTUs from a sample for a specified input and keep only the ones that were found in the samples
    mice_df <- df %>%
      filter(FMT_input %in% c(FMT_label)) %>%
      filter(!type %in% c("FMT_input")) %>%
      filter(species %in% input_df$species) %>%
      mutate(relative_abundance = as.numeric(relative_abundance)) %>%
      select(species,mgID,relative_abundance,FMT_type,FMT_input,genus) %>%
      distinct()
    return(mice_df)
  #return(input_df)
}
yWT_mg_col <- FMT_mg_col("yWT", tax_counts2)
AMS001_mg_col <- FMT_mg_col("AMS001", tax_counts2)
AMS005_mg_col <- FMT_mg_col("AMS005", tax_counts2)
D2A_mg_col <- FMT_mg_col("D2_A", tax_counts2)

all_mg_col <- do.call(rbind, list(yWT_mg_col,AMS001_mg_col,AMS005_mg_col,D2A_mg_col)) %>%
  left_join(.,mg_tax_colors, by="genus") %>%
  mutate(phylum = (factor(phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","Tenericutes","Unclassified","Other")))) %>%
  arrange(phylum) 

plot <- all_mg_col %>%
  mutate(genus = factor(genus, levels = mg_tax_colors$genus)) %>%
  arrange(genus) %>%
  mutate(FMT_input = factor(FMT_input, levels = c( "yWT","D2_A","AMS001", "AMS005"))) %>%
  ggplot(aes(x=mgID, y=relative_abundance, fill=factor(genus))) +
  theme_classic() +
  theme(legend.position="none",
        legend.spacing.x = unit(0, 'cm'),
        legend.spacing.y = unit(0, 'cm'),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill=guide_legend(nrow=15)) +
  geom_bar(stat="identity") +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual(breaks = mg_tax_colors$genus, values = c(mg_tax_colors$gen_color)) +
  labs(x="SampleID", y= "Relative Abundance (%)") +
  facet_grid(~FMT_input, scale="free_x", space="free", switch="x")
ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Manuscript/Figures/Species_col_FMT.pdf",width=5,height=3,dpi=300)

# calculate means for writing
# mean_test <- all_mg_col %>%
#   group_by(FMT_type,FMT_input,mgID) %>%
#   summarize(mean_sum = sum(as.numeric(relative_abundance)), .groups = "drop") %>%
#   group_by(FMT_type) %>%
#   summarize(test = mean(mean_sum), .groups = "drop")

#### Figure S6B: Metagenomic relative abundance of species that did not colonize
## what aabout the donor species that did not colonize?
# command that will tell us which species did not colonize
FMT_mg_noncol <- function(FMT_label, df) {
  # Make a df that for samples with a specified FMT_label, these are all of the possible species found
  mice_df <- df %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(!type %in% c("FMT_input"))
  # And finally lets pull the relative abundance of those OTUs from the donor sample
  species_df <- df %>%
    filter(FMT_input %in% c(FMT_label)) %>%
    filter(type %in% c("FMT_input")) %>%
    filter(!species %in% mice_df$species) %>%
    mutate(relative_abundance = as.numeric(relative_abundance)) %>%
    select(species,mgID,relative_abundance,FMT_input,FMT_type,genus) %>%
    distinct()
  return(species_df)
}
yWT_mg_noncol <- FMT_mg_noncol("yWT", tax_counts2)
AMS001_mg_noncol <- FMT_mg_noncol("AMS001", tax_counts2)
AMS005_mg_noncol <- FMT_mg_noncol("AMS005", tax_counts2)
D2A_mg_noncol <- FMT_mg_noncol("D2_A", tax_counts2)

all_mg_noncol <- do.call(rbind, list(yWT_mg_noncol,AMS001_mg_noncol,AMS005_mg_noncol,D2A_mg_noncol)) %>%
  left_join(.,mg_tax_colors, by="genus") %>%
  mutate(phylum = (factor(phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","Tenericutes","Unclassified","Other")))) %>%
  arrange(phylum) 

plot <- all_mg_noncol %>%
  mutate(genus = factor(genus, levels = unique(mg_tax_colors$genus))) %>%
  arrange(genus) %>%
  mutate(FMT_input = factor(FMT_input, levels = c( "yWT","D2_A","AMS001", "AMS005"))) %>%
  ggplot(aes(x=FMT_input, y=relative_abundance, fill=factor(genus))) +
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
  scale_fill_manual(breaks = mg_tax_colors$genus, values = c(mg_tax_colors$gen_color)) +
  labs(x="SampleID", y= "Relative Abundance (%)") 
ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Manuscript/Figures/Species_noncol.pdf",width=1.5,height=3,dpi=300)

# calculate means for writing
# mean_test <- all_mg_noncol %>%
#   group_by(FMT_type,FMT_input) %>%
#   summarize(mean_sum = sum(as.numeric(relative_abundance)), .groups = "drop") %>%
#   group_by(FMT_type) %>%
#   summarize(test = mean(mean_sum), .groups = "drop")
################################################################################################################################################

