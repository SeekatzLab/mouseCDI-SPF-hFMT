### Color palette for graphs
### A. Seekatz
### 5.5.2024

beast.colors <-list(
	# all FMT groups, by sample source:
  	FMT.allgroups.col=c("healthy" = "#3B99B1",
                "hFMT1" = "#EACB2B",
                "hFMT2" = "#E87700",
                "hFMT3" = "#E8A419",
                "mFMT" = "56B29E", 
                "mFMT_other" = "#9FC095", 
                "noFMT" = "#F5191C"),
                  
	# when only 4 groups, by treatment:               
  	dcol4=c("hFMT" = "#EACB2B",
    			"mFMT" ="#56B29E",
                "mFMT_other" = "#9FC095",
                "noFMT" = "#F5191C"),
                
	# when only 4 groups, by treatment:               
  	dcol3_preFMT=c("hFMT" = "#EACB2B",
    			"mFMT" ="#56B29E",
                "preFMT" = "grey50",
                "noFMT" = "#F5191C"),
		
	# FMT treatment groups, with day 21 and 42, in order:
	dcol_d21=c("21_mFMT" = "#56B29E",
    			"21_hFMT" ="#EACB2B",
                "21_noFMT" = "#F5191C",
                "42_mFMT" = "#56B29E",
                "42_mFMT_other" = "#9FC095",
                "42_hFMT" = "#EACB2B",
                "42_noFMT" = "#F5191C"),      

	# FMT treatment groups, with day 21 and 42, in order:
	dcol_d21_healthy=c("21_mFMT" = "#56B29E",
    			"21_hFMT" ="#EACB2B",
                "21_noFMT" = "#F5191C",
                "42_mFMT" = "#56B29E",
                "healthy" = "#9FC095",
                "42_hFMT" = "#EACB2B",
                "42_noFMT" = "#F5191C"),      
                
	# FMT treatment groups, with day 21 and 42, in order:
	dcol_ba=c("hFMT" = "#EACB2B",
    			"mFMT" ="#56B29E",
                "preFMT" = "grey50",
                "noFMT" = "#F5191C", 
                "healthy" = "grey50"),      
		
	# 16S NMDS colors:
	dcol_nmds=c("pre-abx" = "#3B99B1",
    			"post-cef" ="#474747",
                "post-cdi" = "#A8A8A8",
                "post-vanco" = "#E2E2E2",
                "D2_A" = "#EACB2B",
                "D4_A" = "#E79812",
                "D5_A" = "#E8A91B", 
                "R4_F" = "#ED5300", 
                "AMS001" = "#E87700",
                "AMS005" = "#E8A419", 
                "none" = "#F5191C", 
                "rWT" = "#8BBD94", 
                "yRag" = "#A6C293", 
                "yWT_spore" = "#C1C88C", 
                "yWT" = "#56B29E"),  
	
	# for metabolomics categories:
	SUPER.PATHWAY=c("Amino Acid" = "deepskyblue2",
                  "Carbohydrate" = "darkorchid1",
                  "Cofactors and Vitamins"  = "coral",
                  "Energy" = "darkslategrey",
                  "Lipid" = "gold",
                  "Nucleotide" = "firebrick2",
                  "Peptide" = "lightblue",
                  "Xenobiotics" = "olivedrab2")
  
)