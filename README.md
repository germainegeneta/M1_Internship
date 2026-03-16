# M1_Internship
R scripts for cross-correlation mapping of mosquito population responses and climate variables.

0.Preprocessing Aedes datasets - cleans and joins VectAbundance (Da Re et al., 2024) dataset and Nice ovitrap dataset (EID) and assigns biogeographical EEA-based regions. Outputs are raw joined data tables, data table excluding NA trap counts, data table excluding nonzero trap counts, and data table excluding below 55 trap counts as set by Da Re et al. for their proposed period-over-threshold (POT) phenological index.

0.1 Data visualization - abundance timeseries only - plots and aggregates abundance timeseries into relative weekly abundance, relative yearly abundance, and relative abundance per timeseries splits. Outputs are plots

0.2 Data visualization - abundance and climate timeseries only - plots, aggregates, and overlays climate and abundance timeseries into relative weekly abundance, relative yearly abundance, and relative abundance per timeseries splits. Outputs are plots

2.OpenMeteo data retrieval - downloads historical climate data from OpenMeteo via r/openmeteo. Outputs are data tables containing climate measures for the specified time and location.

3.Cross-correlation mapping using ecoXCorr - fits models and visualizes best time lag windows using CCMs via r/ecoXCorr. Outputs are data tables for fit models, regional-scale CCMs, timeseries split-scale CCMs, and additional plots for best time lag start, end, and length.

3.1Cross-correlation mapping using ecoXCorr - above-threshold datasets - same code as #3 but for above-threshold abundance datasets.

Transforming to DwC format - harmonizes entomological datasets to comply with DarwinCore (DwC) standards. Outputs are data tables following DwC sampling-event, occurrence, and metadata dataset requirements.

Cross-correlation mapping using climwin - generates CCMs via r/climwin. Outputs are site-scale and regional-scale maps computed using beta coefficients and deltaAICc. 
