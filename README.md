# Helianthus-porteri-Demography
Scripts used in the analysis for "Negative density dependance and the soil seed bank buffer an endemic plant from climate-induced population declines"
Thoen and DeMarche (2026)

# Contents

- FitVRMods_RunBoot_RunSim - The script which 1) fits all vital rate models, 2) runs the matrix model, 3) performs bootstraps, and 4) simulates individual-based models. This script requires  VitalRateDataSets.Rdata, IBMListForSimulations.Rdata, and FunctinsForPorteriDemography.R.
- FunctinsForPorteriDemography.R - contains all functions required for the analyses. Includes vital rate kernels, functions to rescale predictor variables, the matrix model function, and lifespan function.
- VitalRateDataSets.Rdata - contains most the datasets required for analyses. Specifically, data for all vital rate models and the base dataframes to run the matrix model and bootstrap the matrix model.
- IBMListForSimulations.Rdata - contains the list which prompts and stores data from individual-based model simulations.

# Notes

- all code were run in R version 4.0.2
