README

This folder contains the following input and output files:

The experimental datasets used for the training and validation of the SMINT-parameters by the simulated-annealing optimization (SAO) algorithm
1) MergedData.xlsx sheet 8 yields and sorted names
2) MergedData.xlsx sheet 9 maximum growth rates and lag times
3) S20_S21_abs_abund_cfu_Senka.xlsx sheet 7 observed abundances
4) S20_S21_abs_abund_cfu_Senka.xlsx sheet 12 observed abundances with corrected initial abundances for low abundance species
5) S20_S21_abs_abund_cfu_Senka.xlsx sheet 13 observed abundances for week transfer experiment
6) SSC21_genera_relative-abundances.xlsx abundances for the SnyCom21MT community (for testing)
7) SynComNX_data_SenkaCB_IG.xlsx Subcommunities abundances (Especially sheet 2:No Mucilaginibacter, sheet 3: No Pseudomonas1, sheet 4: No Lysobacter)

Matlab files:

8) data_to_save_Inter_v2.mat the output of the SAO for the SynCom20 saved in the .mat file
9) data_to_save_SynCom21_Soil_only_Lyso_Parfor_Predation.mat the output of the SAO for the SynCom21 saved in the .mat file without pre-defined preys for Lysobacter
10)Phil_SynCom21_21wk_Paper_Transfer_Predation_Sim_y0mean_y_0_T.mat MATLAB table containing the mean abundances of each species after each transfer, averaged across multiple simulations.
11) struct_tot_SynCom21_Soil_only_Lyso_Parfor.mat the output of the SAO for the SynCom21 saved in the .mat file with pre-defined preys for Lysobacter
12) Zenodo_Senka_Average_InteractionStacked_abundances.mat output of the script Plot_syncom_from_multi_simulation_data_CF_SynCom21.m containing stacked abundances
13)Zenodo_Senka_Average_Interactionz.mat output of the script Plot_syncom_from_multi_simulation_data_CF_SynCom21.m containing day-21 abundances


Text file:

14) Flow_Cytometry_README.txt README file made by Senka to describe SSC21_genera_relative-abundances.xlsx files 
