# Large-scale avian co-occurrence stability over time
This repository contains the R-code and processed data for data processing and analysis and figures in
our manuscript _Large-scale avian co-occurrence stability over time_. 


## Abstract 
Past decades have seen changes in species abundances, distributions, and species
composition, alongside environmental changes. Given this omnipresent sense of
change, we expect most facets of biodiversity to also be dynamic. One such facet
is co-occurrence among species, but we know little about its temporal change,
even though it may indicate disruptions of established species interactions,
dispersal barriers, or shared niches.

Studies of biodiversity dynamics in recent decades have shown changes in
community composition, yet the exact nature of how species assemblages change
remains unclear. Given the commonly observed temporal community composition
turnover and environmental changes, we expected to see changes in the
co-occurrences among species.

Focusing on birds, we analyzed co-occurrences among species pairs using four
large-scale datasets covering Czechia, Europe, New York State, and New Zealand.
Each dataset covers approx. 30 years, offering a unique temporal view on
co-occurrences. 

Surprisingly, we found that bird co-occurrence patterns remained
remarkably stable through time. Changes in co-occurrences showed no tendency
toward positive or negative values and did not correspond with phylogenetic or
functional distances among species pairs, nor with shifts in their range sizes.
In addition, the composition of co-occurrence values within entire groups of
co-occurring species showed little change over time. There was a mild shift
towards more co-occurrence in aquatic species, and less in forest-associated
species. 

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

file name                          |  description | included?
-----------------------------------|--------------|-----------
01_get_co-occurrence_estimates.R   |     bla     | no
02_get_dissimilarities.R
03_get_changes_in_co_occurrences.R
04_make_figures.R
05_figure5.R
