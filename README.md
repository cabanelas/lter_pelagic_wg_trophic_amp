# Trophic Amplification 
## Pelagic Community Structure: Interannual Variability and Long-Term Change in Pelagic Community Structure Across a Latitudinal Gradient
### LTER Pelagic Synthesis Working Group

The Pelagic Community Structure Synthesis Working Group investigates interannual variability and long-term changes in pelagic community structure across a latitudinal gradient. Recent synthesis efforts have identified both similarities and differences in how pelagic marine ecosystems respond to cyclical and long-term environmental changes. This working group aims to test conceptual models describing how pelagic communities respond to stochastic and long-term environmental shifts by leveraging comparative data from four Long-Term Ecological Research (LTER) sites:
* CCE – California Current Ecosystem
* NES – Northeast U.S. Shelf
* NGA – Northern Gulf of Alaska
* PAL – Palmer LTER (Antarctic Peninsula)

Science at the four pelagic LTER sites has matured to the point that we can bring together detailed, comparable, multi‐decadal time series data. Our focus is on understanding how species and communities are structured and how they vary seasonally and interannually. 

**Key Research Questions**
1. Mechanistic Explanations for Long-Term Population Variability
* Linear Tracking Window Hypothesis: Populations track stochastic environmental forcing most effectively when their generation time matches the characteristic timescale of the environmental signal (Hsieh & Ohman 2006). This hypothesis has been supported by studies of phytoplankton, copepods, krill, and fish in the California Current Ecosystem (CCE).
* Double-Integration Hypothesis: Marine population responses to white-noise atmospheric forcing may exhibit strong transitions and prolonged apparent state changes due to cumulative integrations (Di Lorenzo & Ohman 2013). Evidence from euphausiid species in the CCE suggests that responses depend on whether species are influenced by upwelling-driven prey availability or large-scale advection and circulation buffering.
2. Testing Global Model Predictions with Field Data
* Trophic Amplification Hypothesis: Climate-driven ocean biomass declines disproportionately impact higher trophic levels and lower latitudes (Kwiatkowski et al. 2019, Lotze et al. 2019, Petrik et al. 2020). This pattern has been consistent across multiple coupled marine ecosystem models, but field-based validation is lacking. The long-term, multi-trophic level data from the four LTER sites provides an opportunity to assess whether global model predictions hold when confronted with empirical data, and whether simplistic trophic compartmentation in models affects these predictions.
3. Plankton Size Structure and Community Organization
* Normalized Biomass Size Spectrum (NBSS) Hypothesis: Size and size‐structure are defining traits of marine organisms. Analysis of size spectra provides insights into ecosystem resource availability, primary production, and biomass transfer efficiency within food webs (San Martin et al. 2006; Taniguchi et al. 2014). Testing how size spectra and biomass distributions vary across LTER sites, assessing their role in structuring food webs, and evaluating their response to environmental variability over seasonal and interannual timescales.

## Repository Scope
This repository contains scripts specifically related to the trophic amplification hypothesis at the Northeast U.S. Shelf (NES) site. The included script processes NOAA EcoMon zooplankton and larval fish data collected via bongo tows and performs calculations to assess trophic amplification trends over time.

Key Script: trophic_amplification_NES.R
This script:
* Cleans and tidies the EcoMon dataset.
* Computes metrics related to trophic amplification using the following approach:
1. Log transformation: log10(x + (min/2)) at each station, where x is the sum of trophic levels.
2. Aggregation: Averages values across stations for a given year and season.
3. Smoothing: Applies a 5 year running mean.
4. SD: Computes the standard deviation of the resulting time series.

Here, we use volume of zooplankton and sum of 4 forage fish species. 

