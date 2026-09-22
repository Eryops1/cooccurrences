---
output:
  html_document: default
  pdf_document: default
---
# Large-scale bird co-occurrence stability over time
This repository contains the R-code and processed data for data processing and analysis and figures in
our manuscript _Large-scale bird co-occurrence stability over time_. 


## Summary
Focusing on birds, we analyzed co-occurrences among species pairs using four
large-scale datasets covering Czechia, Europe, New York State, and New Zealand.
Each dataset covers approx. 30 years, offering a unique temporal view on
co-occurrences. Using increasingly realistic counterfactual scenarios, we
compare the observed changes in spatial associations to changes expected due to
basic geometric properties, simple environmental sorting and climate- and
landuse change, and we analyse potential traits as drivers of differences
between pairs.

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

| File name | Description |
|:---|:---|:---|:---|
| **SCRIPTS** | **Description** | **Reads** | **Produces** |
|:---|:---|:---|:---|
| `00_climate_download.R` | Downloads and prepares climate and environmental raster data, including annual climate layers and elevation extracts. | `data/all_scales_atlas_<atlas>.gpkg`; external climate and elevation data | `data/climate/annual/*`; climate T1/T2/difference rasters; `data/counterfactuals/env/elevation_extract_atlas<atlas>.rds` |
| `00_climate_models.R` | Extracts climate and land-use predictors and fits climate/environment-based species distribution models used to generate the climate counterfactual. | `data/processed_occupancy_for_CFs_2_<atlas>.rds`; climate rasters; land-use rasters; atlas grids | `data/climate/env_input_<atlas>.rds`; `data/counterfactuals2/c+l_sim_T2_species_<species>_<atlas>.rds` |
| `00_counterfactuals_prep.R` | Processes the occupancy-probability input data and prepares the spatial data used by the counterfactual simulations, including elevation-enhanced grids. | `data/occ_<atlas>/*`; elevation extract; atlas grids | `data/processed_occupancy_for_CFs_2_<atlas>.rds`; `data/one_scale_atlas_with_elevation<atlas>.gpkg` |
| `01_counterfactuals_server_with_uncertainty.R` | Generates the uncertainty simulations and spatial counterfactual communities; designed to run on the server with multiple cores. | Processed occupancy; atlas/elevation data; uncertainty and climate/environment inputs | `data/uncertainty/simulated_<species>_<atlas>.rds`; `data/counterfactuals2/simulated_<species>_<atlas>.rds`; `simulated_elevation_*`; `simulated_env_*` |
| `02_counterfactuals_pairwise_analysis.R` | Calculates pairwise association statistics for the observed and counterfactual communities. | Counterfactual simulation files; processed occupancy; spatial-overlap inputs | `data/counterfactuals2/results/results_<atlas>_<method>_chunk_<chunk>.rds`; `data/uncertainty/results/results_<atlas>_<method>_chunk_<chunk>.rds`; overlap chunk files |
| `03_z_scores_server.R` | Calculates z-scores for qualitative changes in species associations using permutation-based null distributions; designed to run on the server with many cores. | Pairwise analysis results from script 02 | `data/z_scores.rds` |
| `04_all_analysis.R` | Combines the counterfactual and z-score results and performs the main association, transition, trait, and phylogenetic analyses. | Counterfactual and uncertainty results; overlap results; `data/z_scores.rds`; trait input data; AVONET data | `data/final_analysis_input2.rds`; `data/qualitative_transition_analysis.rds`; `corT1T2.rds`; `output/transition_stats_last.csv`; `output/model_phy_sig_output_<atlas>.rds` |
| `05_make_figures.R` | Generates the main and supplementary figures from the saved analysis results. | Saved analysis results; atlas spatial data; geographic boundary files | Figure files in `figures2/` |
| `99_functions.R` | Shared functions sourced by the analysis scripts. | — | — |
| `ssh_run_CF.sh` | Server wrapper for `01_counterfactuals_server_with_uncertainty.R`. | — | Server log files |
| `ssh_run_cf_analysis.sh` | Server wrapper for `02_counterfactuals_pairwise_analysis.R`. | — | Server log files |
| `ssh_z_scores.sh` | Server wrapper for `03_z_scores_server.R`. | — | Server log files |
| **INPUT DATA** | |
| `data/occ_<atlas>/` | Occupancy-probability RDS files used by `00_counterfactuals_prep.R`; these are provided separately. |
| `data/all_scales_atlas_<atlas>.gpkg` | Spatial atlas grids used by the preparation, climate-model, counterfactual, and figure scripts for atlases 5, 6, 17, and 26. |
| `data/counterfactuals/env/elevation/wc2.1_30s/wc2.1_30s_elev.tif` | Elevation raster used to extract elevation for the spatial grids. |
| `data/climate/annual/` | Annual climate raster layers used by the climate models. |
| `data/landuse/` | Land-use raster layers used by the climate/environment models. |
| `data/trait_analysis_input2.rds` | Trait-analysis input data used by `04_all_analysis.R`. |
| `data/AVONET1_BirdLife.csv` | AVONET bird trait data used to attach species traits to the pairwise results. |
| `data/europe.gpkg` | Europe outline used for mapping. |
| `data/czechia.gpkg` | Czechia outline used for mapping. |
| `data/nz.gpkg` | New Zealand outline used for mapping. |
| `data/new_york_state.gpkg` | New York State outline used for mapping. |
| **GENERATED DATA** | |
| `data/counterfactuals/env/elevation_extract_atlas<atlas>.rds` | Elevation extracted for the atlas grids. |
| `data/processed_occupancy_for_CFs_2_<atlas>.rds` | Processed occupancy data used as input for the counterfactual analyses. |
| `data/one_scale_atlas_with_elevation<atlas>.gpkg` | Atlas grids with extracted elevation values. |
| `data/climate/env_input_<atlas>.rds` | Climate and land-use predictors extracted for the atlas grids. |
| `data/counterfactuals2/c+l_sim_T2_species_<species>_<atlas>.rds` | Climate/environment-based species predictions for the second time period. |
| `data/uncertainty/simulated_<species>_<atlas>.rds` | Species-level uncertainty simulations. |
| `data/counterfactuals2/simulated_<species>_<atlas>.rds` | Standard spatial counterfactual simulations. |
| `data/counterfactuals2/simulated_elevation_<species>_<atlas>.rds` | Elevation-informed counterfactual simulations. |
| `data/counterfactuals2/simulated_env_<species>_<atlas>.rds` | Environment-informed counterfactual simulations. |
| `data/<atlas>_overlap_chunk_<chunk>.rds` | Pairwise spatial-overlap results. |
| `data/counterfactuals2/results/results_<atlas>_<method>_chunk_<chunk>.rds` | Pairwise counterfactual analysis results. |
| `data/uncertainty/results/results_<atlas>_<method>_chunk_<chunk>.rds` | Pairwise uncertainty-analysis results. |
| `data/z_scores.rds` | Permutation-based z-score results from the server analysis. |
| `corT1T2.rds` | Pairwise uncertainty results retained for the T1–T2 association analysis. |
| `data/final_analysis_input2.rds` | Main processed pairwise analysis data used by subsequent analyses and figures. |
| `data/qualitative_transition_analysis.rds` | Processed results from the qualitative association-transition analysis. |
| `output/transition_stats_last.csv` | Summary statistics for qualitative association transitions. |
| `output/model_phy_sig_output_<atlas>.rds` | Trait-model and phylogenetic-signal results for each atlas. |
| **CLIMATE DATA GENERATED BY `00_climate_download.R`** | |
| `data/climate/annual/MAT_<year>_<atlas>.tif` | Annual mean annual temperature layers. |
| `data/climate/annual/TAP_<year>_<atlas>.tif` | Annual total annual precipitation layers. |
| `data/climate/annual/SEA_T_<year>_<atlas>.tif` | Annual temperature seasonality layers. |
| `data/climate/annual/SEA_P_<year>_<atlas>.tif` | Annual precipitation seasonality layers. |
| `<atlas>_<variable>_T1.tif`, `<atlas>_<variable>_T2.tif`, `<atlas>_<variable>_diff.tif` | Time-period summaries and differences of the climate variables used in the analyses. |

### Software

All analyses were run in R version 4.5.2. The R packages actually used by the
active code are: `ape`, `broom`, `broom.mixed`, `caret`, `circlize`, `clootl`,
`cluster`, `cowplot`, `data.table`, `exactextractr`, `gbm`, `geodata`,
`ggeffects`, `ggplot2`, `ggpubr`, `ggpmisc`, `ggsci`, `ggtext`, `gt`, `knitr`,
`lmerMultiMember`, `MASS`, `parallel`, `performance`, `phytools`,
`RhpcBLASctl`, `scico`, `scattermore`, `sf`, `sp`, `terra`, and `vegan`.
We use groundhog for package version control.
Tested on Linux Mint 22.3, 32.5 GB RAM.


# Analysis Pipeline

## `00_climate_download.R`

Downloads and prepares the climate and environmental raster data required for the analyses. It extracts elevation data and prepares climate layers for the two time periods and their differences.

## `00_climate_models.R`

Fits climate/environment-based species distribution models using the first time period and predicts species occurrences in the second time period. These predictions provide the climate-based counterfactual used in subsequent analyses.

## `00_counterfactuals_prep.R`

Processes and quality-controls the species occurrence data and prepares the spatial input data for the counterfactual simulations. It also creates the spatial layers required by the elevation- and environment-informed counterfactuals.

## `01_counterfactuals_server_with_uncertainty.R`

Generates the spatial counterfactual communities and their associated uncertainty using the selected counterfactual method. This computationally intensive script is designed to be run on the server with multiple cores.

## `02_counterfactuals_pairwise_analysis.R`

Calculates pairwise species associations from the observed and counterfactual communities. It compares associations across time and across counterfactual scenarios and saves the resulting pairwise statistics.

## `03_z_scores_server.R`

Quantifies qualitative changes in species associations by comparing observed association changes with randomized expectations. The computationally intensive permutation procedure is run separately on the server using many cores.

## `04_all_analysis.R`

Combines the results from the preceding analyses to quantify changes in community associations and their relationships with species traits and other community properties. It produces the statistical results and processed datasets used to generate the figures.

## `05_make_figures.R`

Generates all main and supplementary figures from the saved analysis results. No substantive analyses are performed here beyond calculations directly required for visualization.
