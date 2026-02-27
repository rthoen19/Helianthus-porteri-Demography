# Helianthus-porteri-Demography
Scripts used in the analysis for "Negative density dependance and the soil seed bank buffer an endemic plant from climate-induced population declines"
Thoen and DeMarche (2026)

# Contents

- FitVRMods_RunBoot_RunSim.R - The script which 1) fits all vital rate models, 2) runs the matrix model, 3) performs bootstraps, and 4) simulates individual-based models. This script requires all other files to run completely.
- FunctinsForPorteriDemography.R - contains all functions required for the analyses. Includes vital rate kernels, functions to rescale predictor variables, the matrix model function, and lifespan function.
- prepare_BaseIBMDataForSimulations.R - this script builds the dataframe used to simulate and store data in the indivdual-based model simulations. It outputs a large list.  
- VitalRateDataSets.Rdata - contains most the datasets required for analyses. Specifically, data for all vital rate models and the base dataframes to run the matrix model and bootstrap the matrix model.
- climateDrawData_For1000IBMs.Rdata - includes random draws of climate data which are used to predict vital rates in the IBM simulations. It is required to build the base IBM list.
- StableStabes_ForIBM.Rdata - has a file with the stable stage distribution across all starting seedbank and density scenarios used in the individual based model. Used to build the base IBM list. 

# Notes

- all code were initially run in R version 4.0.2
