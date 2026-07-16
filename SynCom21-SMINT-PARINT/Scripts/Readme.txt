README

This folder has the following scripts (with a short description).

1) AsymmetryIndex.m

Calculates the asymmetry of the byproduct utilization (main Figure 5D and 5E). Computes the asymmetry index of a matrix using the Frobenius norm, 0 = completely symmetric, 2 = completely skewed or symmetric.

2) Average_Mat.m

Calculates the average of the same interaction matrices coming from different fitted sets. 

3) BrayCurtisDistance.m

Calculates the Bray-Curtis dissimilarity of the relative species abundances in the communities, between experimental and/or simulated datasets.

4) cliffsDelta.m
 
Calculates the Cliff’s Delta measure between observed and simulated final biomass.

5) CompSubcomProps.m

Compares abundances of subcommunities. Does this for observed and simulated sub-communities. Refers to Supplementary figure 8.

6) CorrectionT0.py

Infers the missing initial species abundance using a linear fitting on the log abundance from timepoint 12 hours,the first experimentally measured abundance.

7) CorrectionT0Monod.py

Infers the missing initial species abundance assuming a Monod growth. First measured abundance at time point t = 12 hours

8) distinguishable_colors.m

Produces the colors to differentiate the data of the individual soil isolates in the community plots

9) fun_CF_Death_Lyso.m

Function to simulate the SMINT-PARINT model. Calculates the outcome of the combined matrices with SMINT-PARINT parameters.

10) fun_MH_Candidate_Rob_Death.m

SAO-algorithm to search parameter space iteratively. Computes the difference of energy between the current parameter states and the new parameter candidates. Keep new candidate values if rand_val < ratio.  
Generates a new candidate in the space S_consumer*S_consumed
 
11) fun_MH_Candidate_Rob.m

SAO-algorithm to search parameter space iteratively. Computes the difference of energy between the current parameter states and the new candidate values. Keeps the new candidate if rand_val < ratio.  
Generates a new candidate in the space S_consumer

12) Increased_Mat.m

Function used to increase the size of an existing matrix by the addition of a 0 row and/or column. For instance, when we have a SxS matrix and we want to initialize a S+1xS+1 matrix from this SxS matrix by adding a new row and column (n_row = index new row, n_col = index new column).

13) MH_Algorithm_21SynCom_CF_Death_PARFOR_Lysobacter_Alone.m

Main script to start the SAO-algorithm for parameter searches. It only fits SMINT-parameters for Lysobacter keeping the interactions between the other SynCom20 members fixed. (This reference to the dataset_3 generated from the script MH_Algorithm_CF_Death_PARFOR.m in SynCom-20-SMINT is in line 65. This index can be changed, or averaged over the entire fitted sets by using the function Average_Mat.m).
The script processes in 4 SAO loops to fit Cross-feeding, Self-inhibition, Resource consumption and Predation. At the end of the last loop it starts again at the first loop until stop criteria is reached (Number of iterations p < N).

14) permanova_stats.m

Calculates the PERMANOVA statistic (F) and R2. Values obtained as in figure 3 panel B.

15) permanova1_exact.m

Calculates a permutation test based on a PERMANOVA statistic. Outputs as in figure 3 panel B.

16) permanova1.m

Computes the one-factor PERMANOVA pseudo-F statistic, the coefficient of determination R2 and the associated permutation-based p-value from the Bray-Curtis dissimilarity.

17) Permutation_Test_Comm.m

Performs a permutation test based on statistics computed from the Bray-Curtis dissimilarity index between two communities (for instance observed and simulated). 

18) Permutation_Test.m

Performs a permutation test based on Cliff’s delta between two groups (for example, observed and simulated community measurements).

19) Permutt_cosin.m

Tests the similarities and P-value for the comparison of the theta cosine trajectories of simulated versus empirical data (values reported in the Table 1, or Old_figures and not-used => Suppldata-MDS-plots)

20) Plot_syncom_from_multi_simulation_data_CF_SynCom21.m

Plots byproduct formation by species and species biomass formation on byproducts in a consumer-producer matrix, as in main Figure 4A-F, Figure 5B-C. Compares SMINT-PARINT model to a model without SMINT-PARINT interaction using fitted parameter set for CF, death, primary resource consumption and predation. Data from S20_S21_abs_abund_cfu_Senka.xlsx excel file, sheet 7 SMINT_PARINT (soil microcosm no MT)

21) Quantile_Transformer.m

Quantile transformation of the observed biomass. Used in the SAO functions fun_MH_Candidate_Rob_Death.m and fun_MH_Candidate_Rob.m

22) Shannon_Simpson_Indices.m

Calculates the Shannon and Simpson diversity indices on the relative species abundances.

23) Sort_diff.m

Sorts species by their biomass differences between observations and simulations

24) Testing_competition_SnyCom21_individual_CF.m

Plots the inividual observed biomass in soil using the Excel file MergedData_copy.xlsx. First section used to generate Figure 2C. Last section compares observed community with 20 to observed community with 21 members Figure 4B.

25) Transferred_Plot_SynCom.m

Generates a soil microcosm community with 21 members with weekly transfer into a fresh environment. Simulations use SMINT-PARINT interactions with the fitted set or without SMINT-PARINT interactions (only resource competition). Supplementary Figure  9.

26) TrajectoryDist.m

Calculates the theta cosine trajectory. Calls the function Permutt_cosin.m to perform a statistical test in order to assess trajectories similarities.


27) Weighted_Struct.m

Returns a structure containing, for each fields, a weigthed average of the interaction matrices coming from different SAO-optimizations. 

28) Yield_per_Byproducts_Time.m. 

Calculates biomass formation from byproducts per strain and time point. Creation of the heatmap of byproduct exchanges. Figure 5E. Data coming from the S20_S21_abs_abund_cfu_Senka.xlsx excel file, sheet 7 SMINT_PARINT (soil microcosm no MT)


29) Yield_per_Predation_Time.m

Calculates biomass formation from predation per strain and time point. Creation of the heatmap of predation. Figure 4G. Data coming from the S20_S21_abs_abund_cfu_Senka.xlsx excel file, sheet 7 SMINT_PARINT (soil microcosm no MT)