README

This folder has the following scripts (short description is provided). The community growth model can be initiated in the script 'Plot_syncom_from_multi_simulation_data_CF_SynCom20.m'

The SAO algorithm for parameter simulation is initiated in the script 'MH_Algorithm_CF_Death_PARFOR.m'

1) AsymmetryIndex.m

Calculates the asymmetry of the byproduct utilization (main Figure 5D and 5E). Computes the asymmetry index of a matrix using the Frobenius norm, 0 = completely symmetric, 2 = completely skewed or symmetric.

2) Average_Mat.m

Calculates the average of interaction matrices coming from different SAO-optimizations. 

3) BrayCurtisDistance.m

Calculates the Bray-Curtis dissimilarity of the relative species abundances in the communities, between experimental and/or simulated datasets.

4) cliffsDelta.m
 
Calculates the Cliff’s Delta measure between observed and simulated final biomass.

5) Com_Scatter.m

Plots the absolute species abundances in the simulated versus the experimental communities.

6) Compound_headmap.m

Creates the resource heatmap of the individual soil isolates and the log-normal fitted distributions of the growth rates per substrate for each isolate (main Figure 2A and B).

7) distinguishable_colors.m

Produces the colors to differentiate the data of the individual soil isolates in the community plots (e.g., stacked bar as in Fig. 3A)

8) fun_CF_Death.m

Function to simulate the model without predation. Calculates the outcome of the combined matrices with SMINT parameters only. 

9) fun_MH_Candidate_Rob_Death.m

SAO-algorithm to iteratively search parameter space. Computes the difference of energy between the current parameter states and the new candidate values. Keeps new candidates if rand_val < ratio.  
Generates new candidates in the space S_consumer*S_consumed.
 
10) fun_MH_Candidate_Rob.m

SAO-algorithm to search parameter space iteratively. Computes the difference of energy between the current parameter states and the new candidate values. Keeps new candidates if rand_val < ratio.  
Generates candidate in the space S_consumer.

11) ImpactInitialProps.m

Specific script to test the effect of differences in initial abundances on the development of the community (as in Supplementary figure 10)

12) Increased_Mat.m

Function used to increase the size of an existing matrix by the addition of a 0 row and/or column. For instance, when we have a SxS matrix and we want to initialize a S+1xS+1 matrix from this SxS matrix by adding a new row and column (n_row = index new row, n_col = index new column).

13) MH_Algorithm_CF_Death_PARFOR.m

Main script to start the SAO-algorithm for parameter searches. The script processes in 3 SAO loops to fit Cross-feeding, Self-inhibition and Resource consumption. At the end of the last loop it starts again at the first loop until stop criteria is reached (Number of iterations p < N).

14) permanova_stats.m

Calculates the PERMANOVA statistic (F) and R2. Values obtained in figure 3 panel B.

15) permanova1_exact.m

Calculates a permutation test based on a PERMANOVA statistic. Values obtained in the figure 3 panel B.

16) Permutation_Test_Comm.m

Performs a permutation test based on statistics computed from the Bray-Curtis dissimilarity index between two communities (for instance observed and simulated). 

17) Permutation_Test.m

Performs a permutation test based on Cliff’s delta between two groups (for example, observed and simulated community measurements).

18) Permutt_cosin.m

Tests the similarities and P-value for the comparison of the theta cosine trajectories of simulated versus empirical data (values reported in the Table 1, or Old_figures and not-used => Suppldata-MDS-plots)

19) Plot_observed_individual_biomass_in_soil.m

Plots the observed individual biomass in soil (as in main Figure 2C).

20) Plot_syncom_from_multi_simulation_data_CF_SynCom20.m

Plots byproduct formation by species and species biomass formation on byproducts in a consumer-producer matrix, as in main Figure 3A-E, Figure 5A. Compares SMINT model to a model without SMINT interactions using fitted parameter sets for CF, death and primary resource consumption. 

21) Plot_syncom_from_multi_simulation_data_CF_SynCom20Random.m

Plots byproduct formation by species and biomass formation on byproducts in a consumer-producer matrix, from the randomly sampled SMINT parameters, as in supplementary figure 1. Simulates a SMINT model using random interactions for CF.

22) Quantile_Transformer.m

Quantile transformation of the observed biomass. Used in the SAO functions fun_MH_Candidate_Rob_Death.m and fun_MH_Candidate_Rob.m

23) Shannon_Simpson_Indices.m

Calculates the Shannon and Simpson diversity indices on the relative species abundances.

24) Sort_diff.m

Sorts species by their biomass differences between observations and simulations

25) TrajectoryDist.m

Calculates the theta cosine trajectory. Calls the function Permutt_cosin.m to perform a statistical test in order to assess trajectories similarities.

26) Yield_per_Byproducts_Time_SynCom20.m

Calculates biomass formation from byproducts per strain and time point. Creation of the heatmap of byproduct exchanges. Figure 5D.