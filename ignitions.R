#this script cleans fire points data and assigns them to grid cells

library(tidyverse)  
library(sf)         
library(raster)     
library(fasterize)  
library(rgdal)  
library(ggplot2)
library(doSNOW)
library(parallel)
library(exactextractr)
library(terra)


#set desired output pathway
out_path <- "C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly"
dir.create(out_path, recursive = T)

cores <- detectCores() - 2
cl <-makeCluster(cores)
registerDoSNOW(cl)

#read in the fire points for alaska
#this is from the Alaska Large Fire Database, https://www.frames.gov/catalog/10465
firepoints <- read_sf("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/ALFD_Points.shp")
  #add a buffer around each point, solved some errors for polygons
    #firepoints = sf::st_buffer(firepoints, dist = 0)

#set parameters and variable names
in_shape = firepoints
start_year = 1970      #likely need to revisit this
end_year = 2019
grid_size = 100000     #length of one side of grid cell, in METERS

#need to convert FIRESEASON to as.numeric
firepoints$FIRESEASON <- as.numeric(firepoints$FIRESEASON)

 
#discard false alarms (1116 of 31909 observations dropped)
in_shape <- in_shape %>% filter(FALSEALARM == "N") 
  
#discard prescribed burns (229 of 30793 observations dropped)
in_shape <- in_shape %>% filter(PRESCRIBED == "N") 
  
#discard fires of size 0 (1733 of 30564 observations dropped)
in_shape <- in_shape %>% filter(ESTIMATEDT != 0)
  
#keep fire points' location, year (FIRESEASON), and fire size *in acres* (ESTIMATEDT) 
#drop fires outside of the year range 
ak <- in_shape %>% dplyr::select(FIRESEASON, ESTIMATEDT) %>%
  filter(FIRESEASON >= as.numeric(start_year) & FIRESEASON <= as.numeric(end_year)) 
  
#convert shapefile to alaskan albers so units are in meters
ak <- st_transform(ak, '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  
#create a grid over the state
full_grid <- st_make_grid(ak, square = T, cellsize = c(grid_size, grid_size)) %>% 
  st_sf() 
  
#calculate area of grid in km2
full_grid <- full_grid %>% mutate(grid_km2 = (as.numeric(st_area(full_grid) * 1e-6))) 
     #the 1e-6 converts from square meters to square kilometers 
  
#assign each grid cell an ID number 
full_grid <- full_grid %>% dplyr::mutate(ID = row_number())
  
#combine the fire points and the grid
ak_intersect <- st_intersection(ak, full_grid) 
  
#add a fire area column in km2 (called firekm) - was acres in original data
#conversion factor is 1 acre = 0.00404686 km2
ak_intersect <- ak_intersect %>% mutate(firekm = ESTIMATEDT*0.00404686)
  
#count ignitions per grid cell
ak_intersect <- ak_intersect %>% group_by(ID) %>% mutate(count = n())
  
#calculate total burned area per grid cell (all years)
ak_intersect <- ak_intersect %>% group_by(ID) %>% mutate(totalkm = sum(firekm, na.rm = F))
  
#calculate annual burned area per grid cell
ak_intersect <- ak_intersect %>% mutate(ann_km = totalkm/(as.numeric(end_year)-as.numeric(start_year)))
  
#get median ignitions per grid cell per year
#first, squash ak_intersect so there's one observation per grid cell per year
#so, the number of observations per grid cell gives the number of years there were ignitions
test <- ak_intersect %>% group_by(ID, FIRESEASON) %>% summarize(ig_by_year = n()) 
final_list <- list()
test <- test %>% mutate(FIRESEASON=as.integer(FIRESEASON))
#make a list of all fire seasons in the date range 
all_years <- data.frame('FIRESEASON'= start_year:end_year)
  
foreach(id = unique(test$ID), .packages=c("sf", "tidyverse")) %do% {
#foreach(id = unique(filter(test, ID ==24)$ID), .packages=c("sf", "tidyverse")) %do% {
    #use %dopar% for parallel, %do% for not parallel
    #isolate the rows from grid cell in question
  print(paste0('This id = ', id))     #this doesn't do anything in parallel
  #pull out grid cell of interest 
  sub_df <- test %>% dplyr::filter(ID == id)  
    #expand to include row for every fire season
  sub_df <- sub_df %>% full_join(all_years, by=c('FIRESEASON'))
    #each year now has the number of ignitions, which is sometimes NA
  #replace NA's with 0's
  sub_df <- sub_df %>% mutate(ig_by_year = ifelse(is.na(ig_by_year), 0, ig_by_year))
  #similarly, make sure ID always contains the cell ID number 
  sub_df <- sub_df %>% mutate(ID = ifelse(is.na(ID), id, ID))
  sub_df <- sub_df %>% summarise(ID=max(ID), median_ig = median(ig_by_year))
  final_list <- bind_rows(final_list, sub_df)
}
  
final_list <- st_as_sf(final_list)

nrow(ak_intersect)
ak_intersect <- ak_intersect %>% st_join(final_list, by = c("ID"))
nrow(ak_intersect)
  
ak_intersect <- ak_intersect %>% mutate(ID = ID.x)
ak_intersect = subset(ak_intersect, select = -c(ID.y, ID.x))

#select only grid cells in full_grid which have a burned area within them
ak_intersect <- ak_intersect %>% filter(count >= 1)
  
#currently, ak_intersect has every ignition point, sorted into grid cell ID categories
  
#shrink down by grid cell ID --> one observation per cell, record average fire size
ig_count <- ak_intersect %>% group_by(ID) %>% 
  dplyr::summarize(avg_firekm = mean(firekm, na.rm = T), count=max(count), 
                     med = max(median_ig), ann_km=max(ann_km))
  
#put ignition count, avg fire size into grid cell dataframe
full_grid <- full_grid %>% st_join(ig_count, by=c("ID"), left=T)
full_grid <- full_grid %>% mutate(ID = ID.x)
full_grid = subset(full_grid, select = -c(ID.y, ID.x))
  
#get annual number of fires
full_grid <- full_grid %>% mutate(annual_count = count/(as.numeric(end_year)-as.numeric(start_year)))
  
#get rid of grid cells that don't have any burnable area
#read in the land cover file
#land cases are as follows: 
# 1: tundra, 2: cavm barren, 3: cavmvwetlands, 4:water, 5:open boreal forest, 6: boreal wetlands, 7: grasslands, 8:barren/ice, 9: crops, 10: urban, 11: temperate, 12: tundra water, 13: closed boreal forest
lc <- raster('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/BorMask_NDVI_VPD_V5.tif')
  
#save the resolution for later
res <- res(lc)[1] #meters
# 
# #all the values to remove, barren, water etc
bad_vals <- c(2, 4, 8, 10, 11, 12)
# lc[lc %in% bad_vals] = NA
  
#turn bad to 0
lc[lc %in% bad_vals]  = 0
  
#turn bad to 1
lc[lc != 0] = 1
  
#make crs of burned_fraction the same as lc
full_grid <- st_transform(full_grid, crs(lc))
  
#faster extract
ex <-  as_tibble(terra::extract(rast(lc), vect(full_grid), fun = 'sum', list = F, na.rm = T, method = 'simple'))
  
#change ID names
ex$ID = full_grid$ID
names(ex) <- c('ID', 'Count')
  
#remove 0 values (this removes six grid cells)
ex <- ex %>% dplyr::filter(Count >0)
  
#add the area based on the land cover source cell size
ex <- ex %>% mutate(lc_area = (res * res) * 1e-6)
  
#get total area
ex <- ex %>% mutate(final_area = Count * lc_area)
  
# #join back to original burned fraction object, just so we can have the 'grid_m2' field
ex <- left_join(full_grid, as_tibble(ex), by = 'ID') %>% drop_na()
  
#calculate percent of grid cell which is appropriate land cover
ex <- ex %>% mutate(bor_percent = final_area / grid_km2)
  
#scale the burned fraction by the percent boreal pixels
ex <- ex %>% mutate(bf_scaled = (ann_km) / (grid_km2 * bor_percent))
  
#calcuate the fire return interval
ex <- ex %>% mutate(FRI = 1 / bf_scaled)
  
#cap FRI at 1000 if it was larger
ex <- ex %>% mutate(FRI = ifelse(FRI > 1000, 1000, FRI))
  
#full_grid is the dataframe I want
full_grid <- ex
  
#drop bad grid cells in ak_intersect too
ak_intersect <- left_join(ak_intersect, as_tibble(ex), by = 'ID') %>% drop_na()
  
  
full_grid = st_transform(full_grid, '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
ak_intersect = st_transform(ak_intersect, '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
  
#delete output if it exists
if (file.exists(file.path(out_path, 'full_grid.shp'))){file.remove(file.path(out_path, 'full_grid.shp'))}
if (file.exists(file.path(out_path, 'ak_intersect.shp'))){file.remove(file.path(out_path, 'ak_intersect.shp'))}
  
#write out file
write_sf(full_grid, file.path(out_path, 'full_grid.shp'))
write_sf(ak_intersect, file.path(out_path, 'ak_intersect.shp'))
  
  
  
  
  

#try some histograms 
#high-ignition cells include 75, 129, 130, 56, 116
#medium-ignition cells include 100, 125, 115, 140, 153
#low-ignition cells include 41, 24, 36, 135, 168, 84

#histogram of sizes of all fires:
ggplot(as.data.frame(ak_intersect), aes(x = firekm)) + 
  geom_histogram(binwidth = 30) + geom_density()

#histogram of sizes for a specific grid cell:
ggplot(subset(as.data.frame(ak_intersect), ID==110), aes(x = firekm)) + 
  geom_histogram(binwidth = 30) + geom_density()
#all histograms look basically the same - very high left density, long right tail

#plot of ignitions by year
ggplot(subset(as.data.frame(ak_intersect), ID==125) %>% group_by(FIRESEASON) %>% summarize(annual_ignitions = n()), aes(x=annual_ignitions)) +
  geom_histogram(binwidth = 2) + geom_density()                 
#histograms for high-ignition cells show a normal-ish pattern, but most cells have a noticeable right skew
#in the next script, i will use MEDIAN ignitions to avoid overestimating the expected number of ignitions





