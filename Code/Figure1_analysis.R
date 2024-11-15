#### Figure 1
### 7.24.24
### Anna M. Seekatz

### Figures:
	- Fig. 1B: weight change over time across groups
	- Fig. 1C: CFU over time across groups (fecal)
	- Fig. 1D: CFU at days 21, 42 (cecal) across groups
	
### Input files:
	- Data/16S/HFm_summary.txt
	- Beast.colors.R (read this into R to use with figures)

###-------------------------

```{r}
library(tidyr)
library(reshape2)
library(tidyverse)
library(readr)
library(ggplot2)


### Figure 1: design, weights, cfu

### Part 1: weights

# read in data, filter out mouse only samples, and calculate stats for weights:
data <- read.table("Data/16S/HFm_summary.txt", header=T, sep="\t") %>%
		filter(sampleType == "mouse" & type == "fecal") %>%
		filter(!is.na(percent_weight)) %>%
		dplyr::group_by(group, day) %>%
		mutate(median = median(percent_weight),
    			Q1 = quantile(percent_weight, 0.25),
    			Q3 = quantile(percent_weight, 0.75))
# count number of events
data %>%
  group_by(group) %>%
  summarise(n_distinct(mouse))
#1 hFMT                        53
#2 mFMT                        33
#3 mFMT_other                  23
#4 noFMT                       24

# colors and polygons:
dcol4 = c("#EACB2B", "#56B29E", "#9FC095", "#F5191C")	# now in beast.colors

# graph
data %>% ggplot(aes(x=day, y=percent_weight, color=group, fill=group)) +
	scale_color_manual(values=beast.colors$dcol4) + 
	scale_fill_manual(values=beast.colors$dcol4) + 
	geom_vline(xintercept = 0, linetype = "solid", color = "red", alpha = 0.5, size = 0.5) +
	geom_vline(xintercept = 11, linetype = "solid", color = "springgreen4", alpha = 0.5, size = 0.5) +
	annotate("text", x = 0, y = 115, label = "C. diff", size = 3, color = "red", alpha=0.5) +
	annotate("text", x = 11, y = 115, label = "FMT", size = 3, color = "springgreen4") +
	geom_ribbon(aes(ymin = Q1, ymax = Q3, group=group), alpha=0.2, color=NA) +
  	#geom_jitter(width=0.25, alpha=0.4, size = 1, na.rm=TRUE) + #if you wanted points
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE, geom="line", linewidth=1,
			fun.args = list(conf.int=0.50)) +
  	scale_x_continuous(limits=c(-7, 42), breaks=c(-7, 0, 1, 4, 9, 11, 13, 16, 20, 28, 32, 36, 41)) +
	labs(x=NULL, y="% weight") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	theme_classic()
	#geom_polygon(data = cef, aes(x = x, y = y), fill = 'grey', alpha = 0.3) + # couldn't figure this out
	

# stats, just in case:
subdf <- data
subdf$day <- as.factor(subdf$day)
factor_levels <- unique(subdf$day)

# Then, run KW test through oop through factor levels
for (level in factor_levels) {
  subset_data <- filter(subdf, day == level & !is.na(percent_weight) )  # Subset data by factor level
  
  # Check if the subsetted data contains more than one group
  if (length(unique(subset_data$group)) > 1) {
    # Run Kruskal-Wallis test
    result <- kruskal.test(percent_weight ~ group, data = subset_data)
  
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

subdf %>% subset(., day == 13 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# hFMT:noFMT, adj. p = 0.01082
subdf %>% subset(., day == 16 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# hFMT:mFMT, adj. p = 0.004192; hFMT:mFMT_other, adj. p = 0.03066
subdf %>% subset(., day == 17 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# hFMT:mFMT, adj p p = 0.1064
subdf %>% subset(., day == 24 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# no p adj significant
subdf %>% subset(., day == 28 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# hFMT:mFMT_other, adj p =0.01525
subdf %>% subset(., day == 30 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# hFMT:mFMT_other, adj p =0.004888
subdf %>% subset(., day == 32 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# hFMT:mFMT_other, adj p =0.001977; # mFMT:mFMT_other, adj p =0.01166
subdf %>% subset(., day == 35 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# hFMT:mFMT_other, adj p =0.004944
subdf %>% subset(., day == 36 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# noFMT:mFMT_other, adj p =0.01011
#subdf %>% subset(., day == 37 & !is.na(percent_weight) ) %>%
#	dunnTest(percent_weight ~ group, data = .)		# hFMT:mFMT_other, adj p =0.01966
subdf %>% subset(., day == 39 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# hFMT:mFMT_other, adj p =0.007573;
subdf %>% subset(., day == 41 & !is.na(percent_weight) ) %>%
	dunnTest(percent_weight ~ group, data = .)		# hFMT:mFMT_other, adj p =0.006109; # mFMT:mFMT_other, adj p =0.02074

###-------

### part 2, CFU:
data <- read.table("Data/16S/HFm_summary.txt", header=T, sep="\t") %>%
		filter(sampleType == "mouse") %>%
		mutate(CFU = if_else(CFU == 0, 100, CFU)) %>%
		mutate("logCFU" = log10(CFU)) %>%
		filter(!is.na(logCFU)) %>%
		dplyr::group_by(group, day) %>%
		mutate(median = median(logCFU),
    			Q1 = quantile(logCFU, 0.25),
    			Q3 = quantile(logCFU, 0.75))
# count number of events
data %>%
  group_by(group) %>%
  summarise(n_distinct(mouse))
#1 hFMT                        53
#2 mFMT                        33
#3 mFMT_other                  23
#4 noFMT                       21

# colors and polygons:
dcol4 = c("#EACB2B", "#56B29E", "#9FC095", "#F5191C")	# now in beast.colors

# graph:
data %>% ggplot(aes(x=day, y=logCFU, color=group, fill=group)) +
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
  	scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10), labels = expression(0, 10^2, 10^4, 10^6, 10^8, 10^10)) +
	scale_x_continuous(limits=c(-7, 42), breaks=c(-7, 0, 1, 4, 9, 11, 13, 16, 20, 28, 32, 36, 41)) +
	labs(x=NULL, y=expression("Log"["10"] * "(CFU)")) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	theme_classic()
	#geom_polygon(data = cef, aes(x = x, y = y), fill = 'grey', alpha = 0.3) + # couldn't figure this out

# stats, just in case:
subdf <- data
subdf$day <- as.factor(subdf$day)
factor_levels <- unique(subdf$day)

# Then, run KW test through oop through factor levels
for (level in factor_levels) {
  subset_data <- filter(subdf, day == level & !is.na(logCFU) )  # Subset data by factor level
  
  # Check if the subsetted data contains more than one group
  if (length(unique(subset_data$group)) > 1) {
    # Run Kruskal-Wallis test
    result <- kruskal.test(logCFU ~ group, data = subset_data)
  
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
# note: for manuscript, just indicated general significance, as it was difficult to indicate differences between ALL groups on the longitudinal graph
library(FSA)

subdf %>% subset(., day == 16 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)
subdf %>% subset(., day == 17 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)
subdf %>% subset(., day == 24 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)
subdf %>% subset(., day == 28 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)
subdf %>% subset(., day == 30 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)
subdf %>% subset(., day == 32 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)
subdf %>% subset(., day == 35 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)
subdf %>% subset(., day == 36 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)
subdf %>% subset(., day == 37 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)
subdf %>% subset(., day == 39 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)
subdf %>% subset(., day == 41 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)	



###------ 
# Figure 1C: cecal CFU
# separate by day-factor, and convert these to a different 'day'
data <- read.table("Data/16S/HFm_summary.txt", header=T, sep="\t") %>%
		filter(sampleType == "mouse") %>%
		filter(type == "cecal" & day %in% c(21, 42)) %>%
		mutate(CFU = if_else(CFU == 0, 100, CFU)) %>%
		mutate("logCFU" = log10(CFU)) %>%
		unite("group_c", day,group, sep ="_", remove=FALSE) %>%
		mutate(group_c = factor(group_c, levels=c('21_mFMT', '21_hFMT', '21_noFMT', 
										'42_mFMT', '42_mFMT_other', '42_hFMT', '42_noFMT')))
# what do we have?									
data %>%
  dplyr::group_by(group) %>%
  summarise(n_distinct(mouse, day))
#1 hFMT                        49
#2 mFMT                        30
#3 mFMT_other                  20
#4 noFMT                       21
# pretty close--did not lose too many			
# specify colors: also now in beast colors

# specify colors for graph:
dcol2_d21 = c("#56B29E", "#EACB2B", "#F5191C", "#56B29E", "#9FC095", "#EACB2B", "#F5191C")		

# graph
data %>% filter(!is.na(logCFU)) %>%
	ggplot(aes(x=group_c, y=logCFU, color=group_c)) +
	geom_jitter(width=0.25, alpha=0.6, size = 1, na.rm=TRUE) +
	scale_color_manual(values=beast.colors$dcol_d21) + 
	stat_summary(fun.data=median_hilow, size=0.65, na.rm= TRUE,
			fun.args = list(conf.int=0.50)) +
	labs(x=NULL, y=expression("Log"["10"] * "(CFU)")) +
	scale_x_discrete(breaks=c('21_mFMT', '21_hFMT', '21_noFMT', '42_mFMT', '42_mFMT_other', '42_hFMT', '42_noFMT'), 
			labels=c('mFMT', 'hFMT', 'noFMT', 'mFMT', 'mFMT-other', 'hFMT', 'noFMT')) +
		scale_y_continuous(limits=c(0,10), breaks=c(0,2,4,6,8,10), 
		labels = expression(0, 10^2, 10^4, 10^6, 10^8, 10^10)) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1) +
	#geom_vline(xintercept = c(0, 3), linetype = "dashed", color = "blue", size = 1)
	geom_segment(aes(x = 0.8, xend = 3.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	geom_segment(aes(x = 3.5, xend = 7.4, y = 1, yend = 1), color = "grey20", linetype = "solid", size = 0.25) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
	annotate("text", x = 2, y = 0.3, label = "day 21", size = 3, color = "grey20") +
	annotate("text", x = 5.5, y = 0.3, label = "day 42", size = 3, color = "grey20") +
	geom_vline(xintercept = 3.5, linetype = "dashed", color = "grey20")
#ggsave("FigX_bray.dist.dotplot_mFMT.png", width=4, height=4)

# stats:
library(FSA)

data %>% subset(., day == 21 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)		
#    Comparison          Z      P.unadj        P.adj
#1  hFMT - mFMT  5.8299933 5.542961e-09 1.662888e-08
#2 hFMT - noFMT  0.8293911 4.068831e-01 4.068831e-01
#3 mFMT - noFMT -3.6589906 2.532107e-04 5.064213e-04

subdf %>% subset(., day == 42 & !is.na(logCFU) ) %>%
	dunnTest(logCFU ~ group, data = .)		
#          Comparison          Z      P.unadj        P.adj
#1        hFMT - mFMT  3.0847272 2.037389e-03 6.112168e-03
#2  hFMT - mFMT_other  4.7396334 2.141052e-06 1.070526e-05
#3  mFMT - mFMT_other  0.1789471 8.579792e-01 8.579792e-01
#4       hFMT - noFMT -2.0464857 4.070861e-02 8.141721e-02
#5       mFMT - noFMT -4.2083619 2.572286e-05 1.028915e-04
#6 mFMT_other - noFMT -5.6964813 1.223053e-08 7.338319e-08


```
