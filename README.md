# Cumulative_Emissions_Attribution
Code for "Extreme heat and rainfall risk attributed to cumulative CO2 emissions from countries and firms"

This analysis was performed on a high-performance computing system. The underlying data is not available here due to large file sizes, but each dataset is available for download at public locations. The Multi-Model Large Ensemble Archive is available at: https://www.cesm.ucar.edu/community-projects/mmlea. Global Carbon Budget data are available at: https://globalcarbonbudget.org/. The Carbon Majors database is available at: https://carbonmajors.org. CEDS data are available at: https://github.com/JGCRI/CEDS. ERA5 data are available at: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels. CPC data are available at: https://psl.noaa.gov/data/gridded/data.cpc.globalprecip.html.

Specific scripts:
- `MMLEA_Fit_GEV.R` does the actual GEV fitting.
- `Obs_Mask.ipynb` assesses whether the parameters fitted to the observations fall within the model distributions.
- `Figures_Public.ipynb` does the final analysis and makes the main text figures found in the manuscript. 
