# Todos

## Data

- AMOC observational data is stored here: ``/albedo/work/projects/p_pool_clim_data/ERA5/amoc/Johns_2023_mht_data_2020_ERA5/mocha_mht_data_ERA5_v2020.nc``. The depth-levels are in "z". 
"maxmoc" corresponds to "amoc" in climate models. Do we further need "moc"? 

- Observational data: for reducing the uncertainties of the observational data, it is recommended to use several different observational datasets (also done in Brunner et al.)

- sea-surface temperature (tos), use mask?

- check the error with Branch_time_in_parent (now files are just excluded)

- Start with the Paleodata

### Plots

- Make proper functions for plotting  ensemble members and weigths! 
- Fix tos-plot something went wrong here!
- Finish AMOC plot


## Weights

- implement weighted quantiles!!

- We should also use more than just the climatological average as diagnostic?! In Brunner et al. a set of different diagnostics is used. This should be fairly easy to add.

- implement distances for sea-ice extent

## Program

- save the computed weights to target folder from config
- add sanity checks correct configs
- implement proper logging instead of logging to console
- add tests where missing


## CLIMWIP

- add used diagnostics to metadata, now this is not saved, for every variable
different diagnostics might be used 
- check units in performance weights plots
- add plots weighted temp. graph/map
- recipe fails with forked-esmvaltool version.. for mask_sea preprocessor 