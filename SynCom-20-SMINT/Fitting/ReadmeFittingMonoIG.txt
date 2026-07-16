README

This folder contains the following scripts (a short description of each is provided). These R scripts are used to fit the growth kinetics of the soil community species grown in monoculture using a Monod growth model. The fitted parameter values are saved in ../Data/MergedData.xlsx, in the sheet 'Param. Val. Logistic Monod'. The collection of fitted growth kinetics parameters is then read by FitDist.R, which estimates the log-normal distributions that best describe the fitted parameter values.

1)FitDist.R

Reads the fitted parameters from the different monoculture experiments that were produced by  one of the scripts 'SenkaFitting(..).R' (see below). Determines the parameters of the log-normal distribution that best fits the histogram of the growth kinetics.

2)MicrobacteriumFitMultipleCarbons.R

Reads the input file in Data/Transfer of data for SynCom model/ 'Microbacterium_r2a50xdil&ecoplate.xlsx' containing the monoculture growth of Microbacterium. Fits a Monod growth model to the observed growth curves to estimate the growth kinetics parameters (lag times, yields, and maximum growth rates).


3)SenkaFitting.R

Reads the input file in Data/Transfer of data for SynCom model/ 'Senka_OD600_succinategrowth_pooleddata.xlsx' containing monoculture growth curves on succinate by optical density. Fits a Monod growth model to the observed growth curves to estimate the growth kinetics parameters (lag times, yields, and maximum growth rates).

4)SenkaFitting2.R

Reads the input file in Data/Transfer of data for SynCom model/  'PTYG_OD600_21strain_Senka&Tania.xlsx' with monoculture growth curves on PTYG. Fits a Monod growth model to the observed growth curves to estimate the growth kinetics parameters (lag times, yields, and maximum growth rates).

5)SenkaFitting3

Reads the input file in Data/Transfer of data for SynCom model/ 'Senka_Biolog_PlateReader_allstrainscombined.xlsx' with all other monoculture growth from Biolog plates. Fits a Monod growth model to the observed growth curves to estimate the growth kinetics parameters (lag times, yields, and maximum growth rates).