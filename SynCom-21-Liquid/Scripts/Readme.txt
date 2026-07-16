README

This folder has the following scripts (short description is provided). The community growth model can be initiated in the script 'Liquid_Plot_syncom_from_multi_simulation_data_CF_SynCom21.m'

The SAO algorithm for parameter simulation is initiated in the script 'Liquid_MH_Algorithm_CF_Death_PARFOR.m'

1) AsymmetryIndex.m

Calculates the asymmetry of the byproduct utilization (as in main Figure 5D and 5E). Computes the asymmetry index of a matrix using the Frobenius norm, 0 = completely symmetric, 2 = completely skewed or symmetric.

2) BrayCurtisDistance.m

Calculates the Bray-Curtis dissimilarity of the relative species abundances in the communities, between experimental and/or simulated datasets.

3) cliffsDelta.m
 
Calculates the Cliff’s Delta measure between observed and simulated final biomass.

4) distinguishable_colors.m

Produces the colors to differentiate the data of the individual soil isolates in the community plots

5) fun_CF_Death_Lyso.m

Function to simulate the SMINT-PARINT model. Calculates the outcome of the combined matrices with SMINT-PARINT parameters.

6) fun_MH_Candidate_Rob_Death.m

SAO-algorithm to search parameter space iteratively. Computes the difference of energy between the current parameter states and the new parameter candidates. Keeps new candidate values if rand_val < ratio.  
Generates a new candidate in the space S_consumer*S_consumed
 
7) fun_MH_Candidate_Rob.m

SAO-algorithm to search parameter space iteratively. Computes the difference of energy between the current parameter states and the new candidate values. Keeps the new candidate if rand_val < ratio.  
Generates a new candidate in the space S_consumer

8) Increased_Mat.m

Function used to increase the size of an existing matrix by the addition of a 0 row and/or column. For instance, when we have a SxS matrix and we want to initialize a S+1xS+1 matrix from this SxS matrix by adding a new row and column (n_row = index new row, n_col = index new column).

9) Liquid_MH_Algorithm_CF_Death_PARFOR.m

Main script to start the SAO-algorithm for parameter searches. The script processes in three SAO loops to fit Cross-feeding, Self-inhibition and Resource consumption. The predation rate of Lysobacter on its four prey species is fixed at 1e-03. At the end of the last loop it starts again at the first loop until stop criteria is reached (Number of iterations p < N).
10) Liquid_Plot_syncom_from_multi_simulation_data_CF_SynCom21.m

Plots byproduct formation by species and species biomass formation on byproducts in a consumer-producer matrix, as in supplementary Figure 11A-C. Compare SMINT-PARINT model to a model without SMINT-PARINT interaction using fitted parameter set for CF, death, primary resource consumption (simulated as one resource dues to the liquid environment) and predation. Data from Excel file abund_cfu_se.xlsx, sheet 3 SMINT_PARINT (soil community in liquid)

11) Liquid_Yield_per_Byproducts_Time.mCalculates biomass formation from byproducts per strain and time point. Produces a heatmap of byproduct exchanges (Supplementary Figure 11D). Input data from Excel 'abund_cfu_se.xlsx', sheet 3 SMINT_PARINT (soil community in liquid)

12) permanova_stats.m

Calculates the PERMANOVA statistic (F) and R2.

13) permanova1_exact.m

Calculates a permutation test based on a PERMANOVA statistic.

14) Permutt_cosin.m

Tests the similarities and P-value for the comparison of the theta cosine trajectories of simulated versus empirical data (values reported in the Table 1, or Old_figures and not-used => Suppldata-MDS-plots)

15) Permutation_Test_Comm.m

Performs a permutation test based on statistics computed from the Bray-Curtis dissimilarity index between two communities (for instance observed and simulated). 

16) Permutation_Test.m

Performs a permutation test based on Cliff’s delta between two groups (for example, observed and simulated community measurements).

17) Quantile_Transformer.m

Quantile transformation of the observed biomass. Used in the SAO functions fun_MH_Candidate_Rob_Death.m and fun_MH_Candidate_Rob.m

18) Shannon_Simpson_Indices.m

Calculates the Shannon and Simpson diversity indices on the relative species abundances.

19) Sort_diff.m

Sorts species by their biomass differences between observations and simulations


20) TrajectoryDist.m

Calculates the theta cosine trajectory. Calls the function Permutt_cosin.m to perform a statistical test in order to assess trajectories similarities.