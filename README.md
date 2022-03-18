# Snow-Persistence-High-Low-Flow
Data and Code for "Snow Persistence Explains Stream High Flow and Low Flow Signatures with Differing Relationships by Aridity and Climatic Seasonality" (Under Review)

The code, depending on IDE and how the code is opened, may need to have relative directories (as is currently set) hard coded to your system's local directory(ies). 

## Contents

**Data**: The dataset used for the paper (water_climate_data.csv)

**Code**: The code used to generate the plots (plotter_sites_signatures_climatic_metrics.R) and the code to analyze the data (water_climate_analysis.R)

Within water_cliamte_analysis.R, there are some commentted out lines of code to manually examine objects/models (e.g., linear model objects), if such detailed is desired. 

**Variables/Columns Key**:

* mean_sp: Mean Snow Persistence

* mean_ai: Mean Aridity Index

* mean_si: Mean Seasonality Index

* bfi: Baseflow Index

* q5Frac: Normalized Q5

* low_fdc: Low Slope of the Flow Duration Curve

* q95Frac: Normalized Q95

* high_fdc: High Slope of the Flow Duration Curve

