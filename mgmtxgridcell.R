#this script calculates management spending weight for each grid cell 

library(tidyverse)  
library(sf)         
library(raster)     
library(fasterize)  
library(rgdal)  
library(ggplot2)
library(lwgeom)
library(dplyr)
library(doSNOW)
library(parallel)


#set desired output pathway
out_path <- "C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly"
dir.create(out_path, recursive = T)

#load dataframes
#one observation per grid cell
full_grid <- read_sf("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/full_grid.shp")
#all historical fires
ak_intersect <- read_sf("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/ak_intersect.shp")
#2020 FMZs
zones <- read_sf("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/managementzoneshope.shp")

#combine the FMZs and the grid
#change projection to meters, to calculate area 
zones = st_transform(zones, '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
full_grid = st_transform(full_grid, '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
#deal with weird vertices due to zones being complicated polygons 
zones <- zones %>% st_make_valid()

#combine FMZs and grid cells 
abc <- st_intersection(zones, full_grid) 

#calculate area in meters
abc$area <- st_area(abc)

#drop some unused columns
abc = subset(abc, select = -c(Shape_Leng, Shape_Area, DATE_, GIS_ACRES) )

#combine into one observation per zone per grid cell
abc <- abc %>% group_by(ID, PROT) %>% summarise(br_prcn=max(br_prcn), FRI=max(FRI), ann_km=max(ann_km), zone_area=sum(area))

#get rid of undefined area
abc <- abc %>% filter(PROT != "U")

#generate total area so we can get a fraction
abc <- abc %>% group_by(ID) %>% mutate(tot_area = sum(zone_area))

#generate fraction of area in each zone
abc <- abc %>% mutate(zone_frac = zone_area/tot_area)

#generate a management multiplication factor for each grid cell
#the factors are
  #critical -- 6.53
  #full -- 3.70
  #modified -- 1.26
  #limited -- 0.41
#drop some unused columns
abc = subset(abc, select = -c(FRI, tot_area, ann_km) )
abc <- abc %>% mutate(mgmtfactor=0)
abc <- as.data.frame(abc)
abc = subset(abc, select = -c(geometry) )
foreach(id = unique(abc$ID), .packages=c("sf", "tidyverse")) %do% {
  mgmtfac = 
    6.53*sum(as.numeric(abc$zone_frac[(abc$PROT=="C")&(abc$ID==id)])) + 
    3.70*sum(as.numeric(abc$zone_frac[(abc$PROT=="F")&(abc$ID==id)])) +
    1.26*sum(as.numeric(abc$zone_frac[(abc$PROT=="M")&(abc$ID==id)])) + 
    0.41*sum(as.numeric(abc$zone_frac[(abc$PROT=="L")&(abc$ID==id)]))
  abc$mgmtfactor[abc$ID==id] <- mgmtfac
}

#squash abc so we just have mgmtfac for each grid cell, and attach it to full_grid
abc <- abc %>% group_by(ID) %>% summarise(mgmtfactor=max(mgmtfactor))
full_grid = left_join(full_grid, abc, by = 'ID')
full_grid <- full_grid[, !names(full_grid) %in% c("mgmtfctr_y", "mgmtfctr_x")]

#save full_grid file to keep management multiplier 
if (file.exists('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/full_grid.shp')){file.remove('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/full_grid.shp')}
write_sf(full_grid, file.path('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/full_grid.shp'))

