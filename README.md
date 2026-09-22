# Large-scale bird co-occurrence stability over time
This repository contains the R-code and some data for data processing and analysis and figures in
our manuscript _Large-scale bird co-occurrence stability over time_. For all data, please refer to the Zenodo repository.


## Summary
We analyzed co-occurrences among bird species pairs using four
large-scale datasets covering Czechia, Europe, New York State, and New Zealand.
Each dataset covers approx. 30 years, offering a unique temporal view on
co-occurrences. Using increasingly realistic counterfactual scenarios, we
compare the observed changes in spatial associations to changes expected due to
basic geometric properties, simple environmental sorting and climate- and
landuse change, and we analyse potential traits as drivers of differences
between species pairs.

![](figures2/fig1.png "Data and Processing")


## Desription of content & How to use

If you decide to use code or data from this repo, please simply cite the
publication [doi will be inserted once published].

You can download the entire repository and use the derived datasets included
here to explore data, analysis, results, and figures. The computationally
intensive counterfactual and z-score steps are run separately on a server using
the provided shell scripts. You will want to run these parts in parallel to
reduce runtime. For the initial occupancy-probability data and other external
input data, a second repository is necessary. These are provided separately both
for peer review and eventually publication. Link to this data repository is
provided in the manuscript file and will be added here once it goes public.

### Data and files included

| File name | Description | Reads | Produces |
|---|---|---|---|
| `00_climate_download.R` | Downloads and prepares annual climate data and extracts environmental information for the atlas grids. | Elevation raster; atlas grid files; external climate data. | Annual climate rasters; T1/T2 climate rasters and differences; environmental extraction files. |
| `00_counterfactuals_prep.R` | Prepares occupancy and spatial data for the counterfactual analyses. | `data/occ_<atlas>/*`; atlas grid files; elevation raster. | `data/processed_occupancy_for_CFs_2_<atlas>.rds`; `data/one_scale_atlas_with_elevation<atlas>.gpkg`. |
| `00_climate_models.R` | Fits the climate-based counterfactual models for each species and dataset. | Processed occupancy; annual climate rasters; land-use data; atlas grids; uncertainty simulations. | `data/environment/climate/env_input_<atlas>.rds`; climate-model counterfactual RDS files. |
| `01_counterfactuals_server.R` | Generates the stochastic and environmental counterfactual simulations. | Processed occupancy; atlas/grid data; uncertainty simulations; climate-model counterfactuals. | Uncertainty simulations and standard, elevation, environmental, and climate counterfactual simulations. |
| `02_counterfactuals_pairwise_analysis_server.R` | Calculates pairwise co-occurrence/overlap statistics for the counterfactual simulations. | Processed occupancy; uncertainty and counterfactual simulation files. | Chunked overlap and pairwise-result RDS files. |
| `03_z_scores_server.R` | Calculates z-scores comparing observed and simulated pairwise associations. | Uncertainty simulations and pairwise counterfactual results. | `data/z_scores.rds`. |
| `04_all_analysis.R` | Combines results and performs the final statistical, transition, and trait analyses. | Counterfactual results; uncertainty results; overlap files; z-scores; trait data; processed occupancy. | `output/corT1T2.rds`; `data/final_analysis_input2.rds`; `data/qualitative_transition_analysis.rds`; `output/transition_stats_last.csv`; phylogenetic-signal results. |
| `05_make_figures.R` | Produces the figures and figure-specific summaries from the analysis outputs. | Final analysis data; transition analysis; z-scores; processed occupancy; climate metadata; simulation and spatial data. | Files in `figures2/`, including the main figures and supporting plots. |
| `ssh_run_CF.sh` | Server wrapper for running the counterfactual simulations. | `01_counterfactuals_server.R`. | Server run/log output. |
| `ssh_run_CF_analysis.sh` | Server wrapper for running the pairwise counterfactual analysis. | `02_counterfactuals_pairwise_analysis_server.R`. | Server run/log output. |
| `ssh_z_scores.sh` | Server wrapper for calculating z-scores. | `03_z_scores_server.R`. | Server run/log output. |

### Software

All analyses were run in R version 4.5.2. The R packages actually used by the
active code are: `ape`, `broom`, `broom.mixed`, `caret`, `circlize`, `clootl`,
`cluster`, `cowplot`, `data.table`, `exactextractr`, `gbm`, `geodata`,
`ggeffects`, `ggplot2`, `ggpubr`, `ggpmisc`, `ggsci`, `ggtext`, `gt`, `knitr`,
`lmerMultiMember`, `MASS`, `parallel`, `performance`, `phytools`,
`RhpcBLASctl`, `scico`, `scattermore`, `sf`, `sp`, `terra`, and `vegan`.
Tested on Linux Mint 22.3, 32.5 GB RAM.


## Files read and produced by the scripts

| File | Description |
|---|---|
| `data/occ_<atlas>/*` | Species occurrence/occupancy data for each atlas dataset. |
| `data/all_scales_atlas_<atlas>.gpkg` | Atlas spatial data containing the sampling grid and spatial units at different scales. |
| `data/environment/elevation/wc2.1_30s/wc2.1_30s_elev.tif` | Elevation raster used in the spatial preparation and counterfactual analyses. |
| `data/environment/climate/annual/MAT_<year>_<atlas>.tif` | Annual mean annual temperature raster. |
| `data/environment/climate/annual/TAP_<year>_<atlas>.tif` | Annual total annual precipitation raster. |
| `data/environment/climate/annual/SEA_T_<year>_<atlas>.tif` | Annual temperature seasonality raster. |
| `data/environment/climate/annual/SEA_P_<year>_<atlas>.tif` | Annual precipitation seasonality raster. |
| `<atlas>_<variable>_T1.tif` | Climate-variable raster representing conditions at time 1. |
| `<atlas>_<variable>_T2.tif` | Climate-variable raster representing conditions at time 2. |
| `<atlas>_<variable>_diff.tif` | Difference between the climate variable at time 2 and time 1. |
| `data/landuse/*` | Land-use/environmental predictor data used in the environmental counterfactual models. |
| `data/uncertainty/simulated_<species>_<atlas>.rds` | Species-level stochastic simulations used to represent uncertainty. |
| `data/processed_occupancy_for_CFs_2_<atlas>.rds` | Processed occupancy and spatial information prepared for the counterfactual analyses. |
| `data/one_scale_atlas_with_elevation<atlas>.gpkg` | Analysis-scale atlas grid with elevation information attached. |
| `data/environment/climate/env_input_<atlas>.rds` | Prepared climate/environmental predictor data used by the climate-based counterfactual models. |
| `data/counterfactuals2/c+l_sim_T2_species_<species>_<atlas>.rds` | Climate-based counterfactual model output for an individual species and atlas. |
| `data/counterfactuals2/simulated_<species>_<atlas>.rds` | Standard counterfactual simulations. |
| `data/counterfactuals2/simulated_elevation_<species>_<atlas>.rds` | Elevation-based counterfactual simulations. |
| `data/counterfactuals2/simulated_env_<species>_<atlas>.rds` | Environmental counterfactual simulations. |
| `data/<atlas>_overlap_chunk_<chunk>.rds` | Pairwise overlap results calculated in chunks for each atlas. |
| `data/uncertainty/results/results_<atlas>_<method>_chunk_<chunk>.rds` | Pairwise results from the uncertainty simulations. |
| `data/counterfactuals2/results/results_<atlas>_<method>_chunk_<chunk>.rds` | Pairwise results from the counterfactual simulations. |
| `data/z_scores.rds` | Z-scores comparing observed pairwise associations with simulated expectations. |
| `data/trait_analysis_input2.rds` | Trait data prepared for the trait-based analyses. |
| `data/AVONET1_BirdLife.csv` | Bird trait data used in the trait analyses. |
| `data/final_analysis_input2.rds` | Combined final analysis dataset used for the main statistical analyses and figures. |
| `data/qualitative_transition_analysis.rds` | Data summarizing qualitative transitions in co-occurrence associations. |
| `output/corT1T2.rds` | Results describing the relationship between co-occurrence patterns at time 1 and time 2. |
| `output/transition_stats_last.csv` | Summary statistics for the transition analysis. |
| `output/model_phy_sig_output_<atlas>.rds` | Results of the phylogenetic-signal analyses. |
| `data/clim_file_names.csv` | Climate-file metadata used to link climate rasters to the analyses and figures. |
| `figures2/*` | Figure outputs generated by `05_make_figures.R`. |

# Analysis Pipeline

## `00_climate_download.R`

Downloads and prepares the climate and environmental raster data required for the analyses. It extracts elevation data and prepares climate layers for the two time periods and their differences.

## `00_climate_models.R`

Fits climate/environment-based species distribution models using the first time period and predicts species occurrences in the second time period. These predictions provide the climate-based counterfactual used in subsequent analyses.

## `00_counterfactuals_prep.R`

Processes and quality-controls the species occurrence data and prepares the spatial input data for the counterfactual simulations. It also creates the spatial layers required by the elevation- and environment-informed counterfactuals.

## `01_counterfactuals_server.R`

Generates the spatial counterfactual communities and their associated uncertainty using the selected counterfactual method. This computationally intensive script is designed to be run on the server with multiple cores.

## `02_counterfactuals_pairwise_analysis_server.R`

Calculates pairwise species associations from the observed and counterfactual communities. It compares associations across time and across counterfactual scenarios and saves the resulting pairwise statistics.

## `03_z_scores_server.R`

Quantifies qualitative changes in species associations by comparing observed association changes with randomized expectations. The computationally intensive permutation procedure is run separately on the server using many cores.

## `04_all_analysis.R`

Combines the results from the preceding analyses to quantify changes in community associations and their relationships with species traits and other community properties. It produces the statistical results and processed datasets used to generate the figures.

## `05_make_figures.R`

Generates all main and supplementary figures from the saved analysis results. No substantive analyses are performed here beyond calculations directly required for visualization.
