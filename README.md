# AlaskaFireProjection
Contains code and not-publicly-available input data for research projecting fire and management in Alaska through 2100

Start with ignitions.R, which cleans ALFD dataset, makes a grid for the study area, and assigns each fire to a grid cell. 
You will need to download the fire points data and FMZ boundaries; other input data is on github.
Fire points - https://www.frames.gov/catalog/10465
FMZ boundaries - https://hub.arcgis.com/datasets/ffc987188fd64e6b9950fdfac3cdfe0f/explore

Next, run match_historical.R, which is the parameter selection model. 
It generates a parameter for each grid cell so that simulations will match historical burned area averages.

Third, run mgmtxgridcell.R, which calculates the management spending weight for each grid cell based on 2020 FMZ boundaries.

Finally, management.R is the main projection model, which simulates fire and management regimes through 2100. 
Different scenarios can be selected using the Booleans near the beginning of the code. 
