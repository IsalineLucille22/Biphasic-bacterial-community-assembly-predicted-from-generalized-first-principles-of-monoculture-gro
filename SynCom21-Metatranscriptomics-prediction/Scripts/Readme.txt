README

This folder has the following scripts (with a short description).

1) AsymmetryIndex.m

Calculates the asymmetry of the byproduct utilization (main Figure 5D and 5E). Computes the asymmetry index of a matrix using the Frobenius norm, 0 = completely symmetric, 2 = completely skewed or symmetric.

2) Average_Mat.m

Calculates the average of the same interaction matrices coming from different fitted sets. 

3) BrayCurtisDistance.m

Calculates the Bray-Curtis dissimilarity of the relative species abundances in the communities, between experimental and/or simulated datasets.


4) cliffsDelta.m
 
Calculates the Cliff’s Delta measure between observed and simulated final biomass

5) Development_Sim_48hoursMT.m

Runs the SMINT-PARINT model to predict biomass development of the 21-member SynCom in the soil microcosms during the first 48 hours starting from the initial species abundances of the SynCom21-MT experiment. Used to generate supplementary figure 2F

6) Development_Sim_500hoursMT.m

Runs the SMINT-PARINT model to predict biomass development of the 21-member SynCom in the soil microcosms during the first 500 hours. Used to generate figure 5A.

7) distinguishable_colors.m

Produces the colors to differentiate the data of the individual soil isolates in the community plots

8) fun_CF_Death_Lyso.m

Function to simulate the SMINT-PARINT model. Calculates the outcome of the combined matrices with SMINT-PARINT parameters.

9) Increased_Mat.m

Function used to increase the size of an existing matrix by the addition of a 0 row and/or column. For instance, when we have a SxS matrix and we want to initialize a S+1xS+1 matrix from this SxS matrix by adding a new row and column (n_row = index new row, n_col = index new column).

10) permanova_stats.m

Calculates the PERMANOVA statistic (F) and R2. Values obtained in the figure 3 panel B.

11) permanova1_exact.m

Calculates a permutation test based on a PERMANOVA statistic. Values obtained in the figure 3 panel B.

12) Permutation_Test_Comm.m

Performs a permutation test based on statistics computed from the Bray-Curtis dissimilarity index between two communities (for instance observed and simulated). 

13) Permutation_Test.m

Performs a permutation test based on Cliff’s delta between two groups (for example, observed and simulated community measurements).

14) Permutt_cosin.m

Tests the similarities and P-value for the comparison of the theta cosine trajectories of simulated versus empirical data (values reported in the Table 1, or Old_figures and not-used => Suppldata-MDS-plots)

15) Plot_syncom_from_multi_simulation_data_CF_SynCom21MT.m

Plots byproduct formation by species and species biomass formation on byproducts in a consumer-producer matrix, as in main Figure 6B. Compare SMINT-PARINT model to a model without SMINT-PARINT interaction using fitted parameter set for CF, death, primary resource consumption and predation. Data from SSC21_genera_relative-abundances.xlsx excel file, sheet 4 SMINT_PARINT (soil microcosm MT)

16) Shannon_Simpson_Indices.m

Calculates the Shannon and Simpson diversity indices on the relative species abundances.

17) Sort_diff.m

Sorts species by their biomass differences between observations and simulations


18) TrajectoryDist.m

Calculates the theta cosine trajectory. Calls the function Permutt_cosin.m to perform a statistical test in order to assess trajectories similarities.

19) Weighted_Struct.m

Returns a structure containing, for each fields, a weighted average of the interaction matrices coming from different SAO-optimizations. 

20) Yield_per_Byproducts_TimeMT.m

Calculates biomass formation from byproducts per strain and time point. Creation of the heatmap of byproduct exchanges. Figure 6C. Data coming from the SSC21_genera_relative-abundances.xlsx excel file, sheet 4 SMINT_PARINT (soil microcosm MT)
