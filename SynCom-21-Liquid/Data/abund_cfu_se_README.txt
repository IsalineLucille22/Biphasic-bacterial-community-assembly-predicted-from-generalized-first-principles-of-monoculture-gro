README
abund_se.csv
From: "Enhanced Toluene Removal by Pseudomonas veronii Through Phosphorus Metabolic Control" by Bailey et al.

Relative abundance data for SynCom and P. veronii ptxD+ in liquid SE derived from the 16S rRNA gene amplicon sequencing. All samples containing SynCom background community and P. veronii have 600 ug toluene/mL. SynCom background community only has no toluene.Additional conditions are specified (inoculant, phosphite). Library preparation is detailed under "Community 16S rRNA Gene Amplicon Sequencing" in Materials and Methods.

Variables:

-sample_id: Identifier for soil microcosm conditions and the time point. The conditions present are detailed in neighboring columns. 
-timepoint: Time past inoculation samples were taken, in hours.
-inoculant: Indicates whether P. veronii ptxD+ was ("y") or was not ("n") inoculated into the liquid SE.
-phosphite: Indicates whether 1.3 mM phosphite was ("y") or was not ("n") added to liquid SE during inoculation.
-total_grep: Total counts detected in the sample (of any strain) using the Bash “grep” command for regions of V3/V4 that were unique to each strain.
-species: Syncom member being counted.
-abs_abund:  Read counts for each organism normalized to the number of 16S rRNA gene operons found in the given strain.
-rel_abund: Value of "abs_abund" normalized to the total read counts for that sample. This value defines the "relative abundance" as described in the text.
-Replicate: Replicate identifier, first portion of the label in "sample_id".
-cfu_ml: log(CFU/mL) determined for the sample at the given timepoint by CFU counting, calculated as the average of 3-4 technical repliactes. 
-cell_abund: rel_abund * cfu_mL. This value defines the "absolute abundance" as described in the text.