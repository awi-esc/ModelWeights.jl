# Todos

## Data

- Note: AMOC observational data is stored here: ``/albedo/work/projects/p_pool_clim_data/ERA5/amoc/Johns_2023_mht_data_2020_ERA5/mocha_mht_data_ERA5_v2020.nc``.

- Observational data: for reducing the uncertainties of the observational data, it is recommended to use several different observational datasets (also done in Brunner et al.)

- sea-surface temperature (tos), use mask?

- check the error with Branch_time_in_parent (now files are just excluded)

### Plots

- Fix tos-plot something went wrong here!
- Finish AMOC plot


## Weights

- implement different definition of weighted quantiles
- implement distances for sea-ice extent

## Program

- add sanity checks correct configs
- implement proper logging instead of logging to console
- add tests where missing
- in plot-uitls check Lookups for negative values when lookup values had been changed!

## CLIMWIP

- add used diagnostics to metadata, now this is not saved, for every variable
different diagnostics might be used
- check units in performance weights plots
- plot weighted temp. map check why I get more missing data than they
- recipe fails with forked-esmvaltool version.. for mask_sea preprocessor, might be solved due to updated packages 
