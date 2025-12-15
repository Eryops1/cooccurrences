# Large-scale avian co-occurrence stability over time
This repository contains the R-code and processed data for data processing and analysis and figures in
our manuscript _Large-scale avian co-occurrence stability over time_. 


## Abstract (shortened)
Focusing on birds, we analyzed co-occurrences among species pairs using four
large-scale datasets covering Czechia, Europe, New York State, and New Zealand.
Each dataset covers approx. 30 years, offering a unique temporal view on
co-occurrences. 

Overall, we show that most species maintain their co-occurrences over
decades, suggesting that co-occurring species show comparable reactions to
environmental changes. Exploring how this stability in co-occurrence patterns
aligns with previously observed large-scale biotic change could be an exciting
direction for future work.


![](figures/fig1_V2.png "Data and Processing")


## How to use this repository

You can simply download the entire repository and use the derived datasets
included here to explore data, analysis, results, and figures - those are
scripts 03 and 04. For scripts 01 and 02, a second repository with the occupancy
probability data is necessary. These are provided separately both for peer
review and eventually publication. Link to this data repository is provided in
the manuscript file.

Following is a list of files included in the repository:

### Data and files included

|file name                          |  description  | 
|-----------------------------------:|:--------------|
|SCRIPTS                            |               | 
|01_get_co-occurrence_estimates.R   | Loads occupancy probability data, calculates co-occurrence metrics, save processed data for analysis |
|02_get_dissimilarities.R           | Calculate Mantel tests for influence of phylogeny and traits |
|03_get_changes_in_co_occurrences.R | Calculate changes between sampling periods, save the processed data for analysis |
|04_make_figures.R                  | Make figures, get stats for differences between e.g. taxonpomic groups |
|05_figure5.R                       | Do the math for schematic figure 5|
|DATA                               |               |
| all_scales_atlas_17.gpkg          | spatial grids for New Zealand atlas |
| all_scales_atlas_26.gpkg          | spatial grids for Europe atlas |
| all_scales_atlas_5.gpkg           | spatial grids for Czech atlas |
| all_scales_atlas_6.gpkg           | spatial grids for New York atlas |
| czechia.gpkg                      | Czechia outline shapefile | 
| europe.gpkg                       | Europe outline shapefile | 
| new_york_state.gpkg               | New York State outline shapefile | 
| nz.gpkg                           | New Zealand outline shapefile |
| braycurtis.rds                    | Bray-Curtis dissimilarities species level | 
| mantel_results.rds                | atlas level mantel test results | 
| AVONET1_BirdLife.csv              | Avonet Data, available at https://doi.org/10.6084/m9.figshare.16586228 |
| processed_spass.rds               | changes in spatial association (spass) data | 
| range_change.rds                  | range size data on the species level |
| sum_mean_psi_for_maps.rds         | Occupancy sums per grid cell for figure 1 | 
| more_lists_atlas=17_SES_cor_scales=1_2025-11-28.rds | list object containing Spearman correlation, C-score, and null distributions for both measures for each pair, sampling period, and atlas dataset. For example, dat[[1]][[1]] would contain the data for the first sampling period, first species pair |
| more_lists_atlas=5_SES_cor_scales=1_2025-11-28.rds | see above  |
| more_lists_atlas=6_SES_cor_scales=2_2025-11-28.rds | see above  |
| more_lists_atlas=26_SES_cor_scales=1_2025-11-28_chunk_1.rds | see above, but divided in 4 chunks to avoid large data files |
| more_lists_atlas=26_SES_cor_scales=1_2025-11-28_chunk_2.rds | see above |
| more_lists_atlas=26_SES_cor_scales=1_2025-11-28_chunk_3.rds | see above |
| more_lists_atlas=26_SES_cor_scales=1_2025-11-28_chunk_4.rds | see above |
|OUTPUT                             |               |
| | csv files with stats from permutation and mantel tests |



