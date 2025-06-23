#### Figure 4
### 8.7.24
### S. Millard

### Figures:
- Fig. 4A: KO count beta diversity
- Fig. 4B: Bray-curtis distances based on KO abundance
- Fig. 4C: Volcano plot of KOs signifficantly different in abundance between cleared and colonized
- Fig. 4D: heatmap of pathways different in abundance compared to mFMT as "baseline"


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

# read in metadata
meta <- read.table("./data/manuscript_meta.txt", header=T, sep='\t', quote='', na.strings= c("NA"))

###############################################################################################################################
### Figure 4A: NMDS based on KO abundance
# read in ko counts
ko_counts <- read.table("./data/combined_function_ko_tax.tsv", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  filter(!grepl('UN*', Gene.Family)) %>%
  mutate(Gene.Family = str_remove_all(Gene.Family, "\"")) %>%
  separate(Gene.Family, c("ko", "definition"), sep = "\\:\\s") %>%
  separate(definition, c("definition", "tax"), sep = "\\|",fill="right") %>%
  separate(definition, c("definition"), sep = "\\s\\[E", extra="drop") %>% 
  filter(is.na(tax)) %>%
  select(-c(tax, definition))

# change healthy cecal mouse sample type to cecal
nmds_meta <- meta %>%
  mutate(type = replace(type, type== "mouse_survey", "cecal"))

# format count data for vegan input
nmds_counts <- ko_counts %>%
  mutate(across(where(is.numeric), round)) %>%
  column_to_rownames(var="ko") %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% nmds_meta$mgID) %>%
  select(-where(~sum(.) == 0))

# using vegan, calculate bray-curtis dissimilarity based on KO counts
set.seed(19930514)
ko_bray <- avgdist(nmds_counts, dmethod="bray", sample=10671)
ko_bray_matrix <- ko_bray %>%
  as.matrix(labels=TRUE)

### calculate nmds from bray
set.seed(2)
ko_nmds <- metaMDS(ko_bray , trymax=50)
ko_nmds # if you want to check the stress 

nmds_positions <- scores(ko_nmds)

fig4A <- nmds_positions %>%
  as_tibble(rownames="mgID") %>%
  inner_join(., nmds_meta, by="mgID") %>%
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
  ggtitle("NMDS of Bray Curtis Dissimilarity Based on KO Counts")
#ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/figures/4/4A_ko_bray_nmds.pdf",width=5.2,height=4,dpi=300)


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

SPF_ko_bray_stats <- stat_positions %>%
  rownames_to_column("mgID") %>%
  inner_join(., meta_stats, by="mgID") %>%
  pairwiseAdonis::pairwise.adonis2(stat_positions ~ FMT_input, data = ., nperm = 1000, method="bray")

# $AMS001_vs_AMS005
# Df SumOfSqs     R2      F Pr(>F)  
# FMT_input  1  0.43128 0.2495 5.9841  0.012 *
#   Residual  18  1.29727 0.7505                
# Total     19  1.72855 1.0000                
# ---
# $AMS001_vs_yWT
# Df SumOfSqs      R2    F Pr(>F)  
# FMT_input  1   0.2013 0.10803 2.18  0.074 .
# Residual  18   1.6621 0.89197              
# Total     19   1.8634 1.00000              
# ---
# $AMS001_vs_D2_A
# Df SumOfSqs      R2    F Pr(>F)    
# FMT_input  1  0.83816 0.68051 42.6  0.001 ***
#   Residual  20  0.39350 0.31949                
# Total     21  1.23166 1.00000                
# ---
# $AMS001_vs_none
# Df SumOfSqs      R2     F Pr(>F)    
# FMT_input  1   3.5994 0.67213 43.05  0.001 ***
#   Residual  21   1.7558 0.32787                 
# Total     22   5.3552 1.00000                 
# ---
# $AMS005_vs_yWT
# Df SumOfSqs      R2      F Pr(>F)   
# FMT_input  1  0.62975 0.20635 4.6801  0.004 **
#   Residual  18  2.42206 0.79365                 
# Total     19  3.05181 1.00000                 
# ---
# $AMS005_vs_none
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1   1.5029 0.37398 12.545  0.001 ***
#   Residual  21   2.5157 0.62602                  
# Total     22   4.0186 1.00000                  
# ---
# $yWT_vs_D2_A
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1   1.2707 0.45561 16.738  0.001 ***
#   Residual  20   1.5183 0.54439                  
# Total     21   2.7890 1.00000                  
# ---
# $yWT_vs_none
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1   3.3103 0.53471 24.133  0.001 ***
#   Residual  21   2.8806 0.46529                  
# Total     22   6.1909 1.00000                  
# ---
# $D2_A_vs_none
# Df SumOfSqs      R2      F Pr(>F)    
# FMT_input  1   1.5073 0.48323 21.507  0.001 ***
#   Residual  23   1.6120 0.51677                  
# Total     24   3.1193 1.00000                  
# ---
###############################################################################################################################

###############################################################################################################################
### Figure 4B: Bray-curtis distances based on KO abundance
# for Anna
mat <- ko_bray_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "mgID")
#write.table(mat, file = "./data/ko_matrix.txt", sep = "\t", row.names = FALSE, col.names = TRUE, na = "NA")

kos.dist <- mat %>%
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
kos.dist %>% ggplot(aes(x=comparison, y=distance)) +
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
kos.dist %>% 
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
kos.dist %>% subset(., comparison %in% c("mFMT-hFMT", "mFMT-noFMT", "hFMT-noFMT")) %>%
  dunnTest(distance ~ comparison, data = .)
#               Comparison          Z     P.unadj      P.adj
#1  hFMT-noFMT - mFMT-hFMT  2.8559068 0.004291411 0.01287423
#2 hFMT-noFMT - mFMT-noFMT  1.6579011 0.097337424 0.19467485
#3  mFMT-hFMT - mFMT-noFMT -0.5680925 0.569972141 0.56997214

# look at means	
detach(package:plyr)
kos.dist %>% dplyr::group_by(comparison) %>% summarize(mean(distance))		

###############################################################################################################################

###############################################################################################################################
### Figure 4C: volcano plot of KOs different between cleared and colonized determined by Maaslin2
ko_counts <- read.table("./data/combined_function_ko_tax.tsv", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  filter(!grepl('UN*', Gene.Family)) %>%
  mutate(Gene.Family = str_remove_all(Gene.Family, "\"")) %>%
  separate(Gene.Family, c("ko", "definition"), sep = "\\:\\s") %>%
  separate(definition, c("definition", "tax"), sep = "\\|",fill="right") %>%
  separate(definition, c("definition"), sep = "\\s\\[E", extra="drop") %>% 
  filter(is.na(tax)) %>%
  select(-c(tax, definition))

# format metadata for maaslin
maas_meta <- meta %>%
  filter(type %in% c("cecal")) %>%
  arrange(mgID) %>%
  column_to_rownames("mgID")

# format count data for maaslin
maas_counts <- ko_counts %>%
  column_to_rownames(var="ko") %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% rownames(maas_meta)) %>%
  select(-where(~sum(.) == 0))

# first lets compare cleared to colonized
Maaslin2(input_data = maas_counts, input_metadata = maas_meta, output = "./function_maaslin_output/colonization",
         normalization="NONE", transform="LOG", analysis_method ="SLM", correction="BH",
         fixed_effects = c("Colonization"),plot_heatmap=FALSE, plot_scatter=FALSE) 

## Read in KEGG_db so we can classify and color by KO subcategory
volcano_db <- read.csv("./data/kegg_db_color.txt", header=T, sep='\t', na.strings= c("NA"), fill=TRUE) %>%
  separate(Definition, c("Definition"), sep = "\\[", extra="drop")

# read in results - for volcano plot, read in all results and color non-significant grey
all_res <- read.table("./function_maaslin_output/colonization/all_results.tsv", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  separate(feature, c("KO"), sep = "\\_", extra="drop") %>%
  inner_join(.,volcano_db, by="KO") %>% 
  mutate(color = ifelse(qval < 0.0001 & abs(coef) < .5, "grey", color)) %>%
  mutate(color = ifelse(qval > 0.0001, "grey", color)) %>%
  mutate(Subcategory = if_else(color == "grey", "NS", Subcategory)) 


ggplot(all_res, aes(x=coef, y=-log10(qval), color=factor(Subcategory))) +
  geom_point(size=3) +
  theme_bw() +
  theme(text = element_text(size=10), legend.position = "right") +
  ylim(-.2,27) + xlim(-3.7,3.7) +
  geom_hline(yintercept=-log10(0.0001), col="brown",linewidth=1) +
  geom_vline(xintercept = c(-.5,.5), col="navy",linewidth=1) +
  scale_color_manual(values = setNames(volcano_db$color, volcano_db$Subcategory)) +
  ggtitle("LogFC of KO Counts in Colonized vs Clear")
#ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/figures/4/4C_MAAS_ko_volcano.pdf",width=12,height=8,dpi=300)

## for results section, how many increased vs decreased 
sig_cvc_res <- all_res %>%
  filter(qval < 0.0001) %>% # using same significance parameters as volcano plot
  filter(abs(coef) > .5) %>%
  mutate(diff_express = coef) %>%
  mutate(diff_express = ifelse(diff_express > 0, "up", "down"))
up <- sig_cvc_res %>% filter(diff_express %in% c("up")) #143
down <- sig_cvc_res %>% filter(diff_express %in% c("down"))#170
###############################################################################################################################

###############################################################################################################################
### Figures S7 & S8: Heatmap of KOs found to be significantly different between cleared (mFMT) and colonized (hFMT and noFMT)
# read in significant results
col_res <- read.table("./function_maaslin_output/colonization/significant_results.tsv", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  filter(abs(coef) > .5) %>%
  mutate(diff_express = coef) %>%
  mutate(diff_express = ifelse(diff_express > 0, "up", "down")) %>%
  filter(qval < 0.0001)

# read in ko counts
ko_counts <- read.table("./data/combined_function_ko_tax.tsv", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  filter(!grepl('UN*', Gene.Family)) %>%
  mutate(Gene.Family = str_remove_all(Gene.Family, "\"")) %>%
  separate(Gene.Family, c("ko", "definition"), sep = "\\:\\s") %>%
  separate(definition, c("definition", "tax"), sep = "\\|",fill="right") %>%
  separate(definition, c("definition"), sep = "\\s\\[E", extra="drop") %>% 
  filter(is.na(tax)) %>%
  select(-c(tax, definition))

## Read in KEGG_db so we can classify and color by KO subcategory
kegg_db <- read.csv("./data/kegg_db_color.txt", header=T, sep='\t', na.strings= c("NA"), fill=TRUE) %>%
  separate(Definition, c("Definition"), sep = "\\[", extra="drop") 

## S7: KOs in the metabolism subcategory
sig_cvc_metab_counts <- ko_counts %>%
  gather(key="mgID", value="count", -ko) %>%
  rename_with(~ "KO", 1) %>%
  filter(KO %in% col_res$feature) %>%
  left_join(., meta, by="mgID") %>%
  filter(type %in% c("cecal")) %>%
  group_by(FMT_type, KO) %>%
  summarize(mean=log10(mean(count)+1)) %>%
  spread(key='FMT_type', value="mean") %>%
  left_join(.,kegg_db, by="KO") %>%
  filter(Category %in% c("Metabolism")) %>%
  arrange(Category, Subcategory) %>%
  mutate(Definition = if_else(Gene == "apgM","2,3-bisphosphoglycerate-independent phosphoglycerate mutase2", Definition)) %>%
  column_to_rownames("Definition")

cvc_metab_an_col <- list(subcat = setNames(as.character(unique(sig_cvc_metab_counts$color)), unique(sig_cvc_metab_counts$Subcategory)))
cvc_metab_subcat <- data.frame(subcat=sig_cvc_metab_counts$Subcategory)
rownames(cvc_metab_subcat) <- row.names(sig_cvc_metab_counts)

pheatmap(sig_cvc_metab_counts[c(3,2,4)],
         na_col = "white",
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=8)(120),
         border_color="grey20",
         cellheight = 8,
         cellwidth = 20,
         fontsize = 8,
         filename="./figures/4/S5_clear.v.col_METAB_ko_counts.pdf",
         #width=14,
         #height=14,
         annotation_row = cvc_metab_subcat,
         annotation_colors = cvc_metab_an_col,
         angle_col = "45",
         annotation_legend = TRUE
)

## S8: KOs in other subcategories
sig_cvc_other_counts <- ko_counts %>%
  gather(key="mgID", value="count", -ko) %>%
  rename_with(~ "KO", 1) %>%
  filter(KO %in% col_res$feature) %>%
  left_join(., meta, by="mgID") %>%
  filter(type %in% c("cecal")) %>%
  group_by(FMT_type, KO) %>%
  summarize(mean=log10(mean(count)+1)) %>%
  spread(key='FMT_type', value="mean") %>%
  inner_join(.,kegg_db, by="KO") %>%
  filter(!Category %in% c("Metabolism")) %>%
  arrange(Category, Subcategory) %>%
  mutate(Definition = if_else(Gene == "ureE","urease accessory protein ureE", Definition)) %>%
  mutate(Definition = if_else(Gene == "ureG","urease accessory protein ureG", Definition)) %>%
  mutate(Definition = if_else(Gene == "ureD,ureH","urease accessory protein ureD,ureH", Definition)) %>%
  column_to_rownames("Definition")

cvc_other_an_col <- list(subcat = setNames(as.character(unique(sig_cvc_other_counts$color)), unique(sig_cvc_other_counts$Subcategory)))
cvc_other_subcat <- data.frame(subcat=sig_cvc_other_counts$Subcategory)
rownames(cvc_other_subcat) <- row.names(sig_cvc_other_counts)

pheatmap(sig_cvc_other_counts[c(3,2,4)],
         na_col = "white",
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=8)(120),
         border_color="grey20",
         cellheight = 8,
         cellwidth = 20,
         fontsize = 8,
         filename="./figures/4/S6_clear.v.col_OTHER_ko_counts.pdf",
         #width=14,
         #height=14,
         annotation_row = cvc_other_subcat,
         annotation_colors = cvc_other_an_col,
         angle_col = "45",
         annotation_legend = TRUE
)
###############################################################################################################################

###############################################################################################################################
### Figure 4D-H: Counts of Uniref90 IDs associated with butyrate production or bile salt hydrolase
uniref_counts <- read.csv("./data/uniref_counts_only.txt", header=T, sep='\t', na.strings= c("NA"), fill=TRUE) %>%
  mutate(definition = tolower(definition))

# EC 2.8.3.8  But, butyryl-CoA:acetate CoA transferase
but_df <- uniref_counts %>%
  filter(str_detect(definition,"butyryl-coa:acetate")) %>%
  mutate(gene = "But")

#	 EC 2.8.3.9 Ato, butyryl-CoA:acetoacetate CoA transferase
ato_df <- uniref_counts %>%
  filter(str_detect(definition,"butyrate-acetoacetate coa-transferase|butyrate--acetoacetate coa-transferase|acetoacetate coa|butyryl-coa:acetoacetate")) %>%
  mutate(gene = "Ato")

#	EC 2.7.2.7 Buk, butyrate kinase
buk_df <- uniref_counts %>%
  filter(str_detect(definition,"butyrate kinase")) %>%
  mutate(gene = "Buk")

# EC:2.8.3.-	4Hbt, butyryl-CoA:4-hydroxybutyrate CoA transferase
hbt_df <- uniref_counts %>%
  filter(str_detect(definition,"4-hydroxybutyrate coa-transferase|4-hydroxybutyrate coenzyme a transferase|4-hydroxybutyrate--acetyl-coa coa transferase|
                    4-hydroxybutyrate coenzyme a transferase domain protein")) %>%
  mutate(gene = "4Hbt")

# BSH
bsh_df<- uniref_counts %>%
  mutate(definition = tolower(definition)) %>%
  filter(str_detect(definition,"choloylglycine|bile salt hydrolase")) %>%
  mutate(gene = "BSH")

merged_but_bsh <- bind_rows(but_df, ato_df, buk_df, hbt_df,bsh_df) %>%
  gather(key=mgID, value=Count, -UnirefID, -definition, -gene) %>%
  inner_join(., meta, by="mgID") %>%
  mutate(type = replace(type, type== "mouse_survey", "cecal")) %>%
  mutate(FMT_type = replace(FMT_type, FMT_input== "Healthy", "Healthy")) %>%
  filter(type %in% c("cecal")) %>%
  group_by(FMT_type, gene, mgID) %>%
  summarize(sum=sum(Count)) %>%
  ungroup() %>%
  group_by(gene, FMT_type) %>%
  mutate(ymin = min(sum), ymax = max(sum)) %>%
  ungroup() 
  
## make combined of butyrate and bsh
plot <- merged_but_bsh %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  ggplot(aes(FMT_type, sum, color=FMT_type)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position="none",
        panel.border = element_rect(fill=NA)) +
  geom_jitter(width=0.2, alpha=0.6, size = 1, na.rm=TRUE) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.5, linewidth=.5) +
  scale_color_manual(name=NULL,
                     breaks=c("Healthy","mFMT","hFMT","noFMT"),
                     labels=c("Healthy","mFMT","hFMT","noFMT"),
                     values=c('#3B99B1','#56B29E','#EACB2B','#F5191C')) +
  stat_summary(fun=median, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.1), show.legend=FALSE) +
  labs(x="FMT Type", y="BSH and Terminal Butyrate Gene Count (CPM)") +
  facet_wrap(~gene, scales = "free") # if you don't want scales to be the same
#ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/figures/4/4DtoH_combo_bsh_but_counts.pdf",width=5.5,height=4.5,dpi=300)

## stats - haven't figured out how to do this all combined so separating out into each gene
# kruskal.test(x=buk_stat$sum, g=buk_stat$FMT_type) # confirm the post-hoc test will work (need to run this on df before wilcox)
buk_stat <- merged_but_bsh %>%
  filter(gene %in% c("Buk")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)

but_stat <- merged_but_bsh %>%
  filter(gene %in% c("But")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)
Hbt_stat <- merged_but_bsh %>%
  filter(gene %in% c("4Hbt")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)
ato_stat <- merged_but_bsh %>%
  filter(gene %in% c("Ato")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)
bsh_stat <- merged_but_bsh %>%
  filter(gene %in% c("BSH")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)
###############################################################################################################################

###############################################################################################################################
### Figure 4I: heatmap of pathways different in abundance compared to mFMT as "baseline"
path_counts <- read.csv("./data/metacyc_pathabunance.txt", header=T, sep='\t', na.strings= c("NA"), fill=TRUE)

## now lets run them through Maaslin2
maas_meta <- meta %>%
  filter(type %in% c("cecal")) %>%
  arrange(mgID) %>%
  column_to_rownames("mgID") 

maas_path <- path_counts %>%
  select(-pathway) %>%
  column_to_rownames(var="description") %>%
  t() %>%
  as.data.frame() %>%
  filter(rownames(.) %in% rownames(maas_meta)) %>%
  select(-where(~sum(.) == 0))

# if we want to compare all groups to mFMT as "baseline"
Maaslin2(input_data = maas_path, input_metadata = maas_meta, output = "./data/pathway_maaslin_mvo/", 
         normalization="NONE", transform="LOG", analysis_method ="LM", correction="BH",
         fixed_effects = c("FMT_type"), reference = 'FMT_type,mFMT', plot_heatmap=FALSE, plot_scatter=FALSE) 

# read in significant results
path_res <- read.table("./data/pathway_maaslin_mvo/significant_results.tsv", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  mutate(diff_express = coef) %>%
  mutate(diff_express = ifelse(diff_express > 0, "up", "down")) %>%
  mutate(feature = str_replace_all(feature, "\\.", "-")) %>%
  filter(qval < 0.0001)

# read in database to assign colors
path_db <- read.csv("./data/metacyc_database.txt", header=T, sep='\t', na.strings= c("NA"), fill=TRUE) %>%
  filter(description %in% path_res$feature) 

sig_filt_path <- path_counts %>%
  filter(description %in% path_res$feature) %>%
  gather(key="mgID", value="count", -pathway, -description) %>%
  filter(mgID %in% rownames(maas_meta)) %>%
  left_join(., meta, by="mgID") %>%
  left_join(.,path_db, by="description") %>%
  group_by(FMT_type, category, color, pathway) %>%
  summarize(mean=log10(mean(count)+1)) %>%
  filter(category %in% c("Amine-Degradation", "Amino-Acid-Biosynthesis","Amino-Acid-Degradation","Carbohydrates-Biosynthesis","Carbohydrates-Degradation")) %>%
  spread(key='FMT_type', value="mean") %>%
  arrange(category) %>%
  column_to_rownames("pathway")

an_col <- list(subcat = setNames(as.character(unique(sig_filt_path$color)), unique(sig_filt_path$category)))
sig_subcat <- data.frame(subcat=sig_filt_path$category)
rownames(sig_subcat) <- row.names(sig_filt_path)

pheatmap(sig_filt_path[c(4,3,5)],
         na_col = "white",
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=2)(100),
         border_color="black",
         cellheight = 8,
         cellwidth = 20,
         fontsize = 8,
         filename="./figures/4/4I_mouse.vs.others_filtered.pdf",
         #width=14,
         #height=14,
         annotation_row = sig_subcat,
         annotation_colors = an_col,
         angle_col = "45",
         annotation_legend = TRUE
)
###############################################################################################################################

###############################################################################################################################
## S10: Other pathways
sig_other_path <- path_counts %>%
  filter(description %in% path_res$feature) %>%
  gather(key="mgID", value="count", -pathway, -description) %>%
  filter(mgID %in% rownames(maas_meta)) %>%
  left_join(., meta, by="mgID") %>%
  left_join(.,path_db, by="description") %>%
  group_by(FMT_type, category, color, pathway) %>%
  summarize(mean=log10(mean(count)+1)) %>%
  filter(!category %in% c("Amine-Degradation", "Amino-Acid-Biosynthesis","Amino-Acid-Degradation","Carbohydrates-Biosynthesis","Carbohydrates-Degradation")) %>%
  spread(key='FMT_type', value="mean") %>%
  arrange(category) %>%
  column_to_rownames("pathway")

an_col2 <- list(subcat = setNames(as.character(unique(sig_other_path$color)), unique(sig_other_path$category)))
sig_subcat2 <- data.frame(subcat=sig_other_path$category)
rownames(sig_subcat2) <- row.names(sig_other_path)


pheatmap(sig_other_path[c(4,3,5)],
         na_col = "white",
         cluster_cols=F,
         cluster_rows=F,
         color=colorRampPalette(c("white","grey","black"),bias=2)(100),
         border_color="black",
         cellheight = 8,
         cellwidth = 20,
         fontsize = 8,
         filename="./figures/4/S8_mouse.vs.others_OTHER.pdf",
         #width=14,
         #height=14,
         annotation_row = sig_subcat2,
         annotation_colors = an_col2,
         angle_col = "45",
         annotation_legend = TRUE
)
###############################################################################################################################


###############################################################################################################################
### Figure S9: Counts of Uniref90 IDs associated with bai operon
## read in table that contains Uniref90 IDs from Bai operon
bai_df <- read.table("./data/bai_genes.txt", header=T, sep='\t', quote='', na.strings= c("NA"))

bai <- bai_df %>%
  gather(key=mgID, value=Count, -UnirefID, -definition, -Gene) %>%
  inner_join(., meta, by="mgID") %>%
  mutate(type = replace(type, type== "mouse_survey", "cecal")) %>%
  mutate(FMT_type = replace(FMT_type, FMT_input== "Healthy", "Healthy")) %>%
  filter(type %in% c("cecal")) %>%
  group_by(FMT_type, Gene, mgID) %>%
  summarize(sum=sum(Count)) %>%
  ungroup() %>%
  group_by(Gene, FMT_type) %>%
  mutate(ymin = min(sum), ymax = max(sum)) %>%
  ungroup() %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  ggplot(aes(FMT_type, sum, color=FMT_type)) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position="none",
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        panel.border = element_rect(fill=NA)) +
  geom_jitter(width=0.2, alpha=0.6, size = 1, na.rm=TRUE) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.5, linewidth=.5) +
  facet_wrap(~Gene, scales = "free") + # if you don't want scales to be the same
  scale_color_manual(name=NULL,
                     breaks=c("Healthy","mFMT","hFMT","noFMT"),
                     labels=c("Healthy","mFMT","hFMT","noFMT"),
                     values=c('#3B99B1','#56B29E','#EACB2B','#F5191C')) +
  stat_summary(fun=median, geom="crossbar", width=0.5, 
               position=position_dodge(width=0.1), show.legend=FALSE) +
  labs(x="FMT Type", y="BAI Gene Count (CPM)")
ggsave(filename="/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/figures/4/S7_bai_operon_counts.pdf",width=4,height=4.5,dpi=300)

beta_HSDH_stat <- bai %>%
  filter(Gene %in% c("3b-HSDH")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)

alpha_HSDH_stat <- bai %>%
  filter(Gene %in% c("7a-HSDH")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)

baiA_stat <- bai %>%
  filter(Gene %in% c("baiA")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)

baiB_stat <- bai %>%
  filter(Gene %in% c("baiB")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)

baiF_stat <- bai %>%
  filter(Gene %in% c("baiF")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)

baiH_stat <- bai %>%
  filter(Gene %in% c("baiH")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)

baiK_stat <- bai %>%
  filter(Gene %in% c("baiK")) %>%
  select(sum, FMT_type, mgID) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  pairwise_wilcox_test(sum  ~ FMT_type,
                       p.adjust.method="holm",
                       ref.group = "mFMT",
                       paired=F)
###############################################################################################################################

###############################################################################################################################
### Figure 4J-P: 
path_counts <- read.csv("./data/metacyc_pathabunance.txt", header=T, sep='\t', na.strings= c("NA"), fill=TRUE)

aa_path <- path_counts %>%
  filter(description %in% path_res$feature) %>%
  left_join(.,path_db, by="description") %>%
  filter(category %in% c("Amino-Acid-Biosynthesis","Amino-Acid-Degradation"))

path_tax <- read.csv("./data/combined_pathabundance_cpm.tsv", header=T, sep='\t', na.strings= c("NA"), fill=TRUE) %>%
  separate(Pathway, c("description", "pathway"), sep = "\\:\\ ",fill="right") %>%
  separate(pathway, c("pathway", "tax"), sep = "\\|",fill="right") %>%
  filter(!is.na(tax)) %>%
  rename_with(~ str_replace(., "_Abundance.CPM", ""), -1) %>%
  gather(key="mgID", value="count", -pathway, -description, -tax) %>%
  left_join(., meta, by="mgID") %>%
  mutate(type = replace(type, type== "mouse_survey", "cecal")) %>%
  mutate(FMT_type = replace(FMT_type, FMT_input== "Healthy", "Healthy")) %>%
  filter(type %in% c("cecal"))
  
tax_class <- read.table("/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/data/metagenomic_tax_relab.tsv", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  filter(str_detect(clade_name,"\\|s_")) %>%
  separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species","strain"), sep="\\|", fill="right",remove=FALSE) %>%
  mutate_at(vars(!matches("clade_name")),~str_remove(.,".\\_\\_")) %>%
  filter(!is.na(species)) %>%
  filter(is.na(strain)) %>%
  select(phylum, genus) %>%
  distinct() %>%
  mutate(genus = str_replace_all(genus, "\\_unclassified", ""))

tax_filt <- path_tax %>%
  filter(description %in% aa_path$description) %>%
  group_by(FMT_type, pathway, tax) %>%
  summarize(mean=mean(count)) %>%
  separate(tax, c("genus", "species"), sep="\\.", fill="right") %>%
  ungroup() %>%
  mutate_at(vars(matches("genus")),~str_remove(.,".\\_\\_")) %>%
  mutate(genus = str_replace_all(genus, "\\_unclassified", "")) %>%
  mutate(genus = replace(genus, genus == "Clostridiales", "Clostridia")) %>%
  mutate(genus = replace(genus, genus == "Ruminococcaceae", "Ruminococcus")) %>%
  mutate(genus = replace(genus, genus == "Asaccharobacter", "Eggerthella")) %>%
  mutate(genus = replace(genus, genus == "Subdoligranulum", "Oscillospiraceae")) %>%
  left_join(.,tax_class, by="genus") %>%
  mutate(phylum = ifelse(is.na(phylum), "unclassified", phylum)) %>%
  mutate(pathway = gsub("superpathway of ", "", pathway)) %>%
  filter(pathway %in% c("L-histidine degradation I", "L-isoleucine biosynthesis III","L-methionine biosynthesis III","L-ornithine biosynthesis II",
                        "L-arginine biosynthesis","L-serine and glycine biosynthesis I","L-aspartate and L-asparagine biosynthesis",
                        "L-serine and glycine biosynthesis I","L-threonine biosynthesis")) %>%
  arrange(factor(phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","unclassified")), genus)

gen_color <- tax_filt %>%
  select(phylum, genus) %>%
  mutate(phylum = factor(phylum, levels = c("Verrucomicrobia","Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","unclassified")))
gen_color <- gen_color[!duplicated(gen_color$genus), ]
table(gen_color$phylum)
# define colors 
verruco<-c("hotpink")
actino<-colorRampPalette(c("tan","salmon4", "brown","orangered4","darkred"))(n=10)
pro<-colorRampPalette(c("#E4CD05","goldenrod2","#CC982B","#FFA500"))(n=5)
firm<-colorRampPalette(c("midnightblue","royalblue4","dodgerblue4","deepskyblue4","steelblue4","#1211AC","#0000CD","skyblue3","#1291AA","skyblue","lightskyblue3","deepskyblue1","dodgerblue1","steelblue1","royalblue1","#4911EA","#4911AA","#7911AA","#8911AA","darkmagenta","#7939A0","plum4","#CD69C9"))(n=41)
bac<-colorRampPalette(c("#1A6700","#6AB100", "#9AB100","#3AB100"))(n=6)
unclass <- c("#D4D4D4")
Color<-c(verruco, actino, pro, firm, bac, unclass)
gen_color$gen_color <- Color
gen_color$phylum <- factor(gen_color$phylum, levels = c("Verrucomicrobia","Actinobacteria", "Proteobacteria", "Firmicutes", "Bacteroidetes","unclassified")) 
gen_color <- gen_color[order(gen_color$phylum),]

plot <- tax_filt %>%
  mutate(pathway = factor(pathway, levels = c("L-serine and glycine biosynthesis I","L-aspartate and L-asparagine biosynthesis",
                                              "L-arginine biosynthesis","L-isoleucine biosynthesis III","L-methionine biosynthesis III",
                                              "L-ornithine biosynthesis II","L-threonine biosynthesis","L-histidine degradation I"))) %>%
  mutate(FMT_type = factor(FMT_type, levels = c("Healthy","mFMT","hFMT","noFMT"))) %>%
  ggplot(aes(FMT_type, mean, fill=factor(genus),group = factor(phylum, levels = c("Verrucomicrobia","Actinobacteria", "Proteobacteria","Firmicutes", "Bacteroidetes","unclassified")))) +
  theme_classic() +
  geom_bar(stat="identity") + 
  guides(fill = guide_legend(ncol = 2)) +
  scale_fill_manual(breaks = c(gen_color$genus), values = c(gen_color$gen_color)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        axis.title.x = element_blank(),
        legend.text = element_text(),
        legend.position="none",
        panel.border = element_rect(fill=NA),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10)) +
  facet_wrap(~pathway,strip.position = "top")
ggsave(filename="./figures/4/4JtoP_AA_pathway_TAX.pdf",width=4,height=5,dpi=300)
###############################################################################################################################
