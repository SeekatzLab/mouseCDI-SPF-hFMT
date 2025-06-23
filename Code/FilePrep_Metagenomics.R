### This is a script to prep files for downstream analysis including
# - Metaphlann output files

### 4.4.24
### S. Millard
library(gplots)
library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(glue)
library(tibble)
library(reshape2)

# set working directory 
setwd("/Users/home/Library/CloudStorage/Box-Box/Manuscript/ForGit/Data/metagenomics")

################################################ Combine metaphlann outputs ############################################
# for outputs that have one file per sample, this function will read in the file and assign mgID to all entries within that file
file_read <- function(count_data, data_path, pattern)
{
  path <- data_path  
  files <- dir(path, pattern = pattern)
  files
  count_data <- files %>%
    # read in all the files, appending the path before the filename
    map(~ read_delim(file.path(data_path, .), delim="\t", na = c("", "NA", "-----"),
                     col_names=T, escape_double=T, col_types=NULL))
  for (f in 1:length(count_data)) {
    count_data[[f]]$mgID = gsub(pattern, "", files[f])
  }
  count_data <- bind_rows(count_data)
  count_data <- as.data.frame(count_data)
}
tax.data <- file_read(count_data = tax.data, data_path = "./metaphlann/", pattern ="*_metaphlan_bugs_list.tsv")

# remove unnecessary information and write file 
meta <- read.table("/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/data/manuscript_meta.txt", header=T, sep='\t', quote='', na.strings= c("NA"))

tax.data2 <- tax.data %>%
  select(-c("NCBI_tax_id","additional_species")) %>%
  filter(mgID %in% meta$mgID)
write_tsv(tax.data2, "metagenomic_tax_relab.tsv", na = "NA")
##############################################################################################################

############################# defining colors for tax ################################################
tax_counts <- read.table("/Users/home/Library/CloudStorage/Box-Box/Working/humoFMT_manuscript/data/metagenomic_tax_relab.tsv", header=T, sep='\t', quote='', na.strings= c("NA"))

df <- tax_counts %>%
  filter(str_detect(clade_name,"\\|s_")) %>%
  separate(clade_name, c("kingdom", "phylum", "class", "order", "family", "genus", "species","strain"), sep="\\|", fill="right",remove=FALSE) %>%
  mutate_at(vars(!matches("clade_name")),~str_remove(.,".\\_\\_")) %>%
  filter(!is.na(species)) %>%
  filter(is.na(strain)) %>%
  #select(-c(kingdom)) %>%
  mutate(phylum = replace(phylum, phylum == "Candidatus_Saccharibacteria", "Other")) %>% 
  mutate(phylum = replace(phylum, phylum == "Synergistetes", "Other")) %>% 
  mutate(phylum = replace(phylum, phylum == "Bacteria_unclassified", "Unclassified")) %>% 
  merge(., meta, by="mgID") %>% #still joining with metadata because only want to define colors for tax detected in THESE samples
  filter(type %in% c("cecal","mouse_survey","FMT_input")) %>%
  select(phylum, genus)

df$phylum <- factor(df$phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria", "Firmicutes", "Bacteroidetes","Tenericutes","Unclassified","Other"))
df <- df[order(df$phylum),]
# coloring at the genus level
gen_color <- df[!duplicated(df$genus), ]
table(gen_color$phylum)
# define colors 
verruco<-c("hotpink")
actino<-colorRampPalette(c("tan","salmon4", "brown","orangered4","darkred"))(n=24)
pro<-colorRampPalette(c("yellow","yellow2","gold","darkgoldenrod3","goldenrod2","orange2"))(n=11)
firm<-colorRampPalette(c("midnightblue","royalblue4","dodgerblue4","deepskyblue4","steelblue4","#1211AC","#0000CD","skyblue3","#1291AA","skyblue","lightskyblue3","deepskyblue1","dodgerblue1","steelblue1","royalblue1","#4911EA","#4911AA","#7911AA","#8911AA","#7939A0","plum4"))(n=263)
bac<-colorRampPalette(c("#1A6700","#1A4700","#6AB100", "#9AB100","#3AB100","#9EFD80","#7EFD00"))(n=22)
tener <- c("aquamarine2")
unclassified <- colorRampPalette(c("grey40","grey","grey80","grey92","grey30"))(n=26)
other <- colorRampPalette(c("grey1", "grey10"))(n=5)
Color<-c(verruco, actino, pro, firm, bac, tener, unclassified, other)
gen_color$gen_color <- Color
gen_color$phylum <- factor(gen_color$phylum, levels = c("Verrucomicrobia", "Actinobacteria", "Proteobacteria", "Firmicutes", "Bacteroidetes", "Tenericutes","Unclassified","Other")) 
gen_color <- gen_color[order(gen_color$phylum),]

# now at the phylum level
phyl_color <- df[!duplicated(df$phylum), ]
phyl_color <- phyl_color %>% select(phylum)
table(phyl_color$phylum)
ver<-c("hotpink")
act<-c("brown")
pr<-c("orange2")
fir<-c("lightskyblue3")
ba<-c("forestgreen")
ten <- c("aquamarine2")
un <- c("grey")
oth <- c("grey10")
Color2<-c(ver, act, pr, fir, ba, ten, un, oth)
phyl_color$phyl_color <- Color2

tax_colors <- right_join(phyl_color, gen_color, by=c("phylum"))

#write_tsv(tax_colors, "/Users/home/Library/CloudStorage/Box-Box/Manuscript/ForGit/mouseCDI-SPF-hFMT/Data/metagenomics/tax_colors.tsv", na = "NA")
##############################################################################################################

########## adding colors to kegg_db ##########
kegg_db <- read.csv("kegg_db.txt", header=T, sep='\t', na.strings= c("NA"), fill=TRUE) %>%
  select(1:6) %>%
  mutate(Category = str_replace(Category, "^\\d{5}\\s", "")) %>%
  mutate(Subcategory = str_replace(Subcategory, "^\\d{5}\\s", "")) %>%
  mutate(Pathway = str_replace(Pathway, "^\\d{5}\\s", ""))
  
# lets remove some irrelevant stuff and fix a few things
kegg_subset <- kegg_db %>%
  filter(!Category %in% c("Human Diseases","Organismal Systems")) %>%
  filter(!Subcategory %in% c("Aging","Poorly characterized")) %>%
  distinct(KO, .keep_all = TRUE) %>% # remove additional KO after first instance
  mutate(Subcategory = replace(Subcategory, Subcategory == "Metabolism of other amino acids", "Amino acid metabolism")) %>%
  mutate(Subcategory = replace(Subcategory, Subcategory == "Not included in regular maps", "Other metabolism")) %>%
  mutate(Category = case_when(Subcategory == "Unclassified: metabolism" ~ "Metabolism", TRUE ~ Category)) %>%
  mutate(Subcategory = replace(Subcategory, Subcategory == "Unclassified: metabolism", "Other metabolism")) %>%
  mutate(Category = case_when(Subcategory == "Unclassified: genetic information processing" ~ "Genetic Information Processing", TRUE ~ Category)) %>%
  mutate(Category = case_when(Subcategory == "Unclassified: signaling and cellular processes" ~ "Cellular Processes", TRUE ~ Category)) %>%
  mutate(Category = case_when(Subcategory == "Protein families: metabolism" ~ "Metabolism", TRUE ~ Category)) %>%
  mutate(Subcategory = replace(Subcategory, Subcategory == "Protein families: metabolism", "Other metabolism")) %>%
  mutate(Category = case_when(Subcategory == "Protein families: signaling and cellular processes" ~ "Cellular Processes", TRUE ~ Category)) %>%
  mutate(Category = case_when(Subcategory == "Protein families: genetic information processing" ~ "Genetic Information Processing", TRUE ~ Category)) %>%
  mutate(Category = replace(Category, Category == "Brite Hierarchies", "Other")) %>%
  mutate(Subcategory = replace(Subcategory, Subcategory == "Protein families: genetic information processing", "Other genetic information processing")) %>%
  mutate(Subcategory = replace(Subcategory, Subcategory == "Unclassified: genetic information processing", "Other genetic information processing")) %>%
  mutate(Subcategory = replace(Subcategory, Subcategory == "Unclassified: signaling and cellular processes", "Other cellular processes")) %>%
  mutate(Subcategory = replace(Subcategory, Subcategory == "Protein families: signaling and cellular processes", "Other cellular processes")) %>%
  arrange(factor(Category, levels = c("Metabolism", "Cellular Processes", "Genetic Information Processing",
                                      "Environmental Information Processing","Other")))

## now picking colors
kegg_color <- kegg_subset %>%
  #distinct(Subcategory, .keep_all = TRUE) %>%
  mutate(color = Subcategory) %>%
  arrange(Subcategory) %>%
  mutate(color = replace(color, color == "Carbohydrate metabolism", "#BE53D2")) %>% 
  mutate(color = replace(color, color == "Glycan biosynthesis and metabolism", "#BE53E9")) %>% 
  mutate(color = replace(color, color == "Metabolism of cofactors and vitamins", "#FF7F50")) %>%
  mutate(color = replace(color, color == "Energy metabolism", "#2F4F4F")) %>% 
  mutate(color = replace(color, color == "Lipid metabolism", "#F9C541")) %>% 
  mutate(color = replace(color, color == "Nucleotide metabolism", "#E70005")) %>% 
  mutate(color = replace(color, color == "Amino acid metabolism", "#009AD1")) %>% 
  mutate(color = replace(color, color == "Metabolism of terpenoids and polyketides", "#ADD8E6")) %>% 
  mutate(color = replace(color, color == "Biosynthesis of other secondary metabolites", "#C1BDCD")) %>%
  mutate(color = replace(color, color == "Xenobiotics biodegradation and metabolism", "#C1BDCD")) %>%
  mutate(color = replace(color, color == "Other metabolism", "#C1BDCD")) %>% 
  mutate(color = replace(color, color == "Transport and catabolism", "#8B0A50")) %>% 
  mutate(color = replace(color, color == "Cell growth and death", "#CD1076")) %>% 
  mutate(color = replace(color, color == "Cellular community - eukaryotes", "#EE1289")) %>% 
  mutate(color = replace(color, color == "Cellular community - prokaryotes", "#EE1289")) %>% 
  mutate(color = replace(color, color == "Cell motility", "#FF69B4")) %>% 
  mutate(color = replace(color, color == "Other cellular processes", "#FFE6F5")) %>% 
  mutate(color = replace(color, color == "Transcription", "#1A4700")) %>% 
  mutate(color = replace(color, color == "Translation", "#228B22")) %>% 
  mutate(color = replace(color, color == "Folding, sorting and degradation", "#6E8B3D")) %>% 
  mutate(color = replace(color, color == "Replication and repair", "#A2CD5A")) %>%
  mutate(color = replace(color, color == "Information processing in viruses", "#C5D1BE")) %>%
  mutate(color = replace(color, color == "Other genetic information processing", "#C5D1BE")) %>%
  mutate(color = replace(color, color == "Signaling molecules and interaction", "#F5C77E")) %>%
  mutate(color = replace(color, color == "Membrane transport", "#EC9006")) %>%
  mutate(color = replace(color, color == "Signal transduction", "#FFAF7A")) %>%
  mutate(color = replace(color, color == "Viral protein families", "#E0E0E0")) %>%
  mutate(color = replace(color, color == "RNA family", "#E0E0E0"))
write_tsv(kegg_color, "kegg_db_color.txt", na = "NA")
##############################################################################################################

#################### cleaning up uniref counts ####################
uniref_counts <- read.csv("combined_genefamilies_uniref_fixed.tsv", header=T, sep='\t', na.strings= c("NA"), fill=TRUE) %>%
  rename_with(~ "gene_family", .cols = 1) %>%
  filter(!gene_family %in% c("UNMAPPED")) 

# dropping tax info tax 
uniref_counts2 <- uniref_counts %>%
  separate(gene_family, c("UnirefID", "definition"), sep = "\\:\\s", fill="right") %>%
  separate(definition, c("definition", "tax"), sep = "\\|",fill="right") %>%
  filter(!definition %in% c("NO_NAME")) %>%
  filter(is.na(tax)) %>% 
  select(-tax)
#write_tsv(uniref_counts2, "uniref_counts_only.txt", na = "NA")

# keeping tax info
uniref_counts3 <- uniref_counts %>%
  separate(gene_family, c("UnirefID", "definition"), sep = "\\:\\s", fill="right") %>%
  separate(definition, c("definition", "tax"), sep = "\\|",fill="right") %>%
  filter(!definition %in% c("NO_NAME")) %>%
  filter(!is.na(tax)) %>%
  separate(tax,into=c("genus","species"),sep = "\\.",fill="right")
#write_tsv(uniref_counts3, "uniref_counts_tax.txt", na = "NA")

bai_df <- uniref_counts %>%
  mutate(definition = tolower(definition)) %>%
  filter(str_detect(definition,"bai|bile acid|hydroxysteroid|hydroxycholanate|bile salt|cholenoic")) %>%
  mutate(Gene = definition) %>%
  select(Gene, everything())
write_tsv(bai_df, file = "bai_genes.txt", na = "NA")
##############################################################################################################

########################## Assign color to metacyc pathwyas to match kegg subcategory colors #####################################################
metacyc_db <- read.table("./data/Pathways-from-All-pathways-of-MetaCyc.txt", header=T, sep='\t', quote='', na.strings= c("NA")) %>%
  filter(description %in% path_counts$description) %>%
  mutate(color = category) %>%
  arrange(category) %>%
  mutate(color = replace(color, color == "AROMATIC-COMPOUNDS-DEGRADATION", "#0079A7")) %>% 
  mutate(color = replace(color, color == "Amine-Degradation", "#2C7FE0")) %>% 
  mutate(color = replace(color, color == "Amino-Acid-Biosynthesis", "#009AD1")) %>% 
  mutate(color = replace(color, color == "Amino-Acid-Degradation", "#00BBFA")) %>% 
  mutate(color = replace(color, color == "Antibiotic-Resistance", "#F14421")) %>% 
  mutate(color = replace(color, color == "Carbohydrates-Biosynthesis", "#BE53D2")) %>% 
  mutate(color = replace(color, color == "Carbohydrates-Degradation", "#EE82EE")) %>% 
  mutate(color = replace(color, color == "Carboxylate-Degradation", "#B3EE3A")) %>% 
  mutate(color = replace(color, color == "Cell-Structure-Biosynthesis", "#C1C633")) %>% 
  mutate(color = replace(color, color == "Cofactor-Biosynthesis", "#FF7F50")) %>% 
  mutate(color = replace(color, color == "Energy-Metabolism", "#2F4F4F")) %>% 
  mutate(color = replace(color, color == "Fatty-Acid-and-Lipid-Degradation", "#FFD700")) %>% 
  mutate(color = replace(color, color == "Lipid-Biosynthesis", "#F9C541")) %>% 
  mutate(color = replace(color, color == "Metabolic-Regulators", "#40AA84")) %>% 
  mutate(color = replace(color, color == "Noncarbon-Nutrients", "#198591")) %>% 
  mutate(color = replace(color, color == "Nucleic-Acid-Processing", "#800020")) %>% 
  mutate(color = replace(color, color == "Nucleotide-Biosynthesis", "#E70005")) %>% 
  mutate(color = replace(color, color == "Nucleotide-Degradation", "#FA180C")) %>% 
  mutate(color = replace(color, color == "Polyamine-Biosynthesis", "#B30012")) %>% 
  mutate(color = replace(color, color == "Polyprenyl-Biosynthesis", "#EA6F35")) %>% 
  mutate(color = replace(color, color == "Protein-Modification", "#F29A3B")) %>% 
  mutate(color = replace(color, color == "Secondary-Metabolite-Biosynthesis", "#5D38BD")) %>% 
  mutate(color = replace(color, color == "Tetrapyrrole-Biosynthesis", "#8E25B6")) %>%
  mutate(color = if_else(!grepl("#", color), "grey50", color))
#write_tsv(metacyc_db, "metacyc_database.txt", na = "NA")
##############################################################################################################
