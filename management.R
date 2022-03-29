#this script runs projections with increased management activity 

library(tidyverse)  
library(sf)     
library(doSNOW)
library(raster)     
library(fasterize)  
library(rgdal)  
library(ggplot2)
library(ggridges)
library(RColorBrewer)  
library(ggpubr)
library(dplyr)

options(warn = -1)

#set desired output pathway
out_path <- "C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/projection"
dir.create(out_path, recursive = T)

#load dataframes 
exponents <- read_sf("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/exponents.shp")
prob <- read_csv('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/ProbAdjust.csv')
full_grid <- read_sf('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/full_grid.shp')
ak_intersect <- read_sf("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/ak_intersect.shp")
cc_adjust <- read_sf("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/CCadjustmentFactor.csv")
SCCoutput <- read_csv('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/SCCoutput.csv')

#fix projection to be all the same
exponents <- st_transform(exponents, crs = st_crs(full_grid))

#grab burnable fraction for each grid cell to get FRI later - it is the br_prcn column in full_grid
#grab the boreal_percent (br_prcn) and ID columns from full_grid data frame
exponents <- left_join(exponents, as.data.frame(full_grid[,c("br_prcn", "ID")]), by = 'ID')
#get rid of new repeat geometry column and fix the remaining geometry column's name 
exponents = subset(exponents, select = -c(geometry.y))
names(exponents)[names(exponents)=="geometry.x"] <- "geometry"
st_geometry(exponents) <- "geometry"
#each grid cell is 10000 km^2 in area

#load fire size cdf
firesorted <- sort(ak_intersect$firekm)
e_cdf <- 1:length(firesorted)/length(firesorted)

#set some initial parameters and dataframes
start_year <- 2020
end_year <- 2100
histyearfirst = 1970      
histyearlast = 2019
histburned = sum(ak_intersect$firekm[which(ak_intersect$FIRESEA<1990)])/20
#average annual burned area 1970-1989 is our  baseline burned area average,  2061.797 km2/yr
modernburned = sum(ak_intersect$firekm[which(ak_intersect$FIRESEA>2009)])/10
#average annual burned area 2010-2019 is our modern baseline burned area average, 5114.747 km2/yr

#these are for choosing ignitions by draw
test <- ak_intersect %>% group_by(ID, FIRESEA) %>% summarize(ig_by_year = n()) 
test <- test %>% mutate(FIRESEA=as.integer(FIRESEA))
all_years <- data.frame('FIRESEA'= histyearfirst:histyearlast)

#mgmtfactor in full_grid dataframe gives the management multiplier for each grid cell 
#need to determine saturation point in each grid cell, when mgmt increases s.t. fire size -> 0
#past the saturation point, continuing to increase mgmt results in no further burned area changes in that grid cell 
full_grid <- full_grid %>% mutate(saturation = ((log(0.05)/log(1-0.002063))/(100*mgmtfactor))+1)

#options to turn on/off
full_mgmt_effectiveness = 1
half_mgmt_effectiveness = 0
full_scc = 1
half_scc = 0
mgmt_x2 = 1
mgmt_x1 = 0
mgmt_historical_baseline = 0
mgmt_modern_baseline = 0

costfunc<-function(mgmtlvl){
  #adjust fire sizes: 1% increase in spending -> 0.2063% decrease in burned area (Phillips et al forthcoming)
  #weighted by mgmtfactor, which is based on coverage of different FMZs
  foreach(id = unique(manageproj$ID[(manageproj$year==year) & (manageproj$iteration==i)]), .packages=c("sf", "tidyverse")) %do% {
    manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)] <- 
      manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)& (manageproj$ID==id)]*
      (1-0.002063)^((mgmtlvl-1)*100*full_grid$mgmtfactor[full_grid$ID==id])*full_mgmt_effectiveness +      
      manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)& (manageproj$ID==id)]*
      (1-(0.002063*0.5))^((mgmtlvl-1)*100*full_grid$mgmtfactor[full_grid$ID==id])*half_mgmt_effectiveness
    if(manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)]<=0) {
      manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)]<-0.0001 
    }#end if, which makes sure burned area can't go negative
  }#end foreach
  #cost is in 2015 dollars 
  #add additional mgmt to cost
  #add excess emissions to cost, if fire_sum for the year greater than historical baseline 
  #3325gC/m2 -> 12191.67 metric tons CO2/km2
  #while most combustion emissions become atmospheric CO2, not all do, following Akagi et al (2011).
  #conversion factor is 0.84
  #convert scc from $2010/ton CO2 to $2015/ton CO2 using (237.002/218.076)
  cost <-
    ((sum(manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i)])-histburned)*12191.67*0.84*
       SCCoutput$SCC[SCCoutput$Year==2019+year]*(237.002/218.076)*full_scc +          
       (sum(manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i)])-histburned)*12191.67*0.84*
       SCCoutput$SCC[SCCoutput$Year==2019+year]*(237.002/218.076)*0.5*half_scc  
    )#*(sum(manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i)]) > histburned) 
  
  foreach(id = unique(manageproj$ID[(manageproj$year==year) & (manageproj$iteration==i)]), .packages=c("sf", "tidyverse")) %do% {
    if (mgmtlvl < full_grid$saturation[full_grid$ID==id]) {
      cost <- cost + (full_grid$mgmtfactor[full_grid$ID==id]/sum(full_grid$mgmtfactor))*(mgmtlvl)*99386014.75 
    }
    if (mgmtlvl >= full_grid$saturation[full_grid$ID==id]) {
      cost <- cost + (full_grid$mgmtfactor[full_grid$ID==id]/sum(full_grid$mgmtfactor))*
        (full_grid$saturation[full_grid$ID==id])*99386014.75 
    }
  }
  
  return(cost)
}


for(i in 76:100) {
  #for(i in 1:1) {
  itercount = i
  print(paste0('This iteration = ', i))
  
  #generate an empty dataframe to store results
  manageproj <- data.frame(matrix(ncol = 12, nrow = 0))
  colnames(manageproj) <- c('year', 'ID', 'br_pcrn', 'fire_sum', 'FRI', 'geometry', 'iteration', 'mgmtlvl', 'mgmtbound', 'cost', 'optfirekm', 'socialcost')
  
  
  for(year in 1:80) {
    #for(year in 1:1) {
    
    print(paste0('This year = ', year))
    
    #foreach (id = unique(filter(exponents, ID >190)$ID), .packages=c("sf", "tidyverse")) %do% {
    foreach(id = unique(exponents$ID), .packages=c("sf", "tidyverse")) %do% {
      
      #print(paste0('This ID = ', id))
      
      #decide number of ignitions --
      #make a dataframe to use as the ignitions distribution
      #pull out grid cell of interest from test
      ig_dist <- test %>% dplyr::filter(ID == id) 
      #expand to include row for every fire season
      ig_dist <- ig_dist %>% full_join(all_years, by=c('FIRESEA'))
      #each year now has the number of ignitions, which is sometimes NA
      #replace NA's with 0's
      ig_dist <- ig_dist %>% mutate(ig_by_year = ifelse(is.na(ig_by_year), 0, ig_by_year))
      #similarly, make sure ID always contains the cell ID number 
      ig_dist <- ig_dist %>% mutate(ID = ifelse(is.na(ID), id, ID))
      #make a year column, which starts at 1 rather than 1960 or w/e histyearfirst is
      ig_dist <- ig_dist %>% mutate(year = ig_dist$FIRESEA + 1 - histyearfirst)
      
      #using random draw, choose number of ignitions for this grid cell in this year
      #use the same seed for all grid cells in each year
      #basically each projection year matches ignitions from a historic year across the whole state
      set.seed(i*year)
      ignitions <- ig_dist$ig_by_year[ig_dist$year==min(ig_dist$year[which(ig_dist$year >= runif(1)*50[1])])]
      
      #adjust for climate factor, which increases ignitions
      ignitions<-as.numeric(cc_adjust$CCFactor[which(cc_adjust$number==year)]) * ignitions
      ignitions<-round(ignitions)
      
      #assuming there are ignitions in this grid cell in this year...
      if (ignitions>0) {
        
        #i'm going to make a lil subdf to store this cell in this year, temporarily
        sub_df <- exponents %>% dplyr::filter(ID == id)
        
        #expand the current sub_df by the number of ignitions 
        sub_df <- sub_df[rep(seq_len(nrow(sub_df)), each = ignitions), ]
        
        #add columns for year, firekm, FRI
        sub_df <- sub_df %>% mutate(year=year, firekm=0, FRI=0)
        
        #populate fire sizes for each ignition, by random draw 
        #for the first 30 years we use historical FRI
        if(year<31) {
          #calculate fire size
          sub_df$firekm <- replicate(ignitions, firesorted[which(e_cdf >= runif(1)^sub_df$exponent[1])[1]]*1.76)
          #adjust for low FRI
          if (full_grid[full_grid$ID==id,]$FRI < 84.0) {
            sub_df$firekm <- sub_df$firekm*prob$ProbAdjustment[which(prob$FRI >= full_grid[full_grid$ID==id,]$FRI)[1]]
          }
          #populate FRI (just historical for years 1-30)
          sub_df$FRI <- full_grid[full_grid$ID==id,]$FRI
        }
        
        if(year>30) {
          #calculate fire size
          sub_df$firekm <- replicate(ignitions, firesorted[which(e_cdf >= runif(1)^sub_df$exponent[1])[1]]*1.76)
          #the *1.76 adjusts the baseline from 1996 to 2020
          
          #adjust for low FRI if necessary 
          #grab FRI from nearest existing year
          #first if statement makes sure there has been at least one ignition by Year-1
          if (any(manageproj$year<=(year-1) & manageproj$ID==id & manageproj$iteration==i)) {
            #figure out the closest nonzero year, to use for subdf year reference 
            cutoff <- unique(manageproj$year[(manageproj$ID==id) & (manageproj$iteration==i)])
            cutoff <- max(cutoff[cutoff<=(year-1)])
            #if the FRI for the most recent year with burning is small, adjust
            if (manageproj$FRI[(manageproj$ID==id) & (manageproj$iteration==i) & (manageproj$year==cutoff)][1] < 84.0) {
              sub_df$firekm <- sub_df$firekm*prob$ProbAdjustment[which(prob$FRI >=
                                                                         manageproj$FRI[(manageproj$ID==id) & (manageproj$iteration==i) & (manageproj$year==cutoff)][1])[1]]
            }
            #if this is the first ignition of the run, FRI should be high and no adjustment will be needed
            
            #calculate FRI. we want to use optfirekm not fire_sum, since we want the chosen burned area, not the original
            #if there was at least one fire before 30 years ago
            if (!is.na( sum(manageproj$fire_sum[(manageproj$year<=(year-30)) & (manageproj$ID==id) & (manageproj$iteration==i)], na.rm=TRUE))) {
              sub_df$FRI <- 1 / (((sum(manageproj$optfirekm[(manageproj$year<=year) & (manageproj$ID==id) & (manageproj$iteration==i)], na.rm=TRUE)-
                                     sum(manageproj$optfirekm[(manageproj$year<=(year-30)) & (manageproj$ID==id) & (manageproj$iteration==i)], na.rm=TRUE))/30)/
                                   (10000*full_grid[full_grid$ID==id,]$br_prcn))
              if (sub_df$FRI[1]>1000) {sub_df$FRI <- 1000}
            }
            #if there were no fires before 30 years ago
            if (is.na( sum(manageproj$optfirekm[(manageproj$year<=(year-30)) & (manageproj$ID==id) & (manageproj$iteration==i)], na.rm=TRUE))) {
              sub_df$FRI <- 1 / (((sum(manageproj$optfirekm[(manageproj$year<=year) & (manageproj$ID==id) & (manageproj$iteration==i)], na.rm=TRUE))/30)/
                                   (10000*full_grid[full_grid$ID==id,]$br_prcn))
              if (sub_df$FRI[1]>1000) {sub_df$FRI <- 1000}
            }
          }
          #if this is the first fire, set FRI=1000
          if (!any(manageproj$year<=(year-1) & manageproj$ID==id & manageproj$iteration==i)) {sub_df$FRI <- 1000}
        } #this closes the if(year>30) condition, calculating fire size based on previous FRI
        
        #store sub_df as csv
        #if (file.exists(file.path(out_path, paste0("i", id, "y", year,'.csv')))){file.remove(file.path(out_path, paste0("i", id, "y", year,'.csv')))}
        #write_csv(sub_df, file.path(out_path, paste0("i", id, "y", year,'.csv')))    #using write_csv here preserves the geometry column
        
        #shrink sub_df to have one observation per year, and keep track of the iteration
        sub_df <- sub_df  %>% summarise(year=max(year), ID=max(ID), br_prcn=max(br_prcn), 
                                        fire_sum=sum(firekm), FRI=mean(FRI))
        #add mgmt data columns, to be filled in after optimization
        sub_df <- sub_df %>% mutate(iteration=as.numeric(itercount), mgmtlvl=1, mgmtbound=0, cost=0, optfirekm=0, socialcost=0)
        
        #append sub_df results to fireproj dataframe 
        manageproj <- rbind(manageproj, sub_df)
        
      }#this closes the if condition, of "if there are any fires this year+grid cell"
    }#this closes the foreach loop, going through grid cells
    
    #adjust management cost to minimize total costs (excess mgmt + excess emissions)
    cost<-0
    mgmtlvl<-1
    mgmtbound<-1
    
    #first step -- figure out upper mgmt bound, mgmt to hit historical baseline
    ExitFlag <- 0
    
    if (mgmt_historical_baseline==1) {
      #choose ecological bound so that 10 year average does not exceed 5% of historical burned area 
      #for the first 9 simulation years, use a shorter moving average of only the existing years
      #if current burned area more than 5% larger than historical, increase mgmt
      if (year<10) {
        if (((histburned*year)-
             (sum(manageproj$optfirekm[(manageproj$year<year) (manageproj$iteration==i)]) + 
              sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)]) 
             ))/(histburned*year) < (-0.05)) {
          mgmtbound <- mgmtbound +.01
        }
        if (((histburned*year)-
             (sum(manageproj$optfirekm[(manageproj$year<year) & (manageproj$iteration==i)]) + 
              sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)]) 
             ))/(histburned*year) >= (-0.05)) {
          ExitFlag<-1
        }}
      if (year>=10) {
        if (((histburned*10)-
             (sum(manageproj$optfirekm[(manageproj$year>(year-10)) & (manageproj$year!=year) & (manageproj$iteration==i)]) +
              sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)])))/(histburned*10) < (-0.05)) {
          mgmtbound <- mgmtbound +.01
        }
        if (((histburned*10)-
             (sum(manageproj$optfirekm[(manageproj$year>(year-10)) & (manageproj$year!=year) & (manageproj$iteration==i)]) +
              sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)])))/(histburned*10) >= (-0.05)) {
          ExitFlag<-1
        }}
      
      while (ExitFlag==0) {
        #adjust fire sizes: 1% increase in spending -> 0.2063% decrease in burned area (Phillips et al forthcoming)
        #weighted by mgmtfactor, which is based on coverage of different FMZs
        subdf<-0
        foreach(id = unique(manageproj$ID[(manageproj$year==year) & (manageproj$iteration==i)]), .packages=c("sf", "tidyverse")) %do% {
          #print(paste0('This ID = ', id))
          #print(paste0('This mult = ', (1-0.002063)^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])))
          
          if((1-0.002063)^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])>0.05){
            subdf <- subdf +
              manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)][1]*
              (1-0.002063)^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])
            #print(paste0('This subdf = ', subdf))
          }}
        #if current burned area close to historical, or smaller than historical, change exit flag
        if (year<10) {
          if (((histburned*year)-
               (sum(manageproj$optfirekm[(manageproj$year<year) & (manageproj$iteration==i)]) + subdf))/(histburned*year) < (-0.05)) {
            mgmtbound <- mgmtbound +.01
          }
          if (((histburned*year)-
               (sum(manageproj$optfirekm[(manageproj$year<year) & (manageproj$iteration==i)]) + subdf))/(histburned*year) >= (-0.05)) {
            ExitFlag<-1
          }}
        if (year>=10) {
          if (((histburned*10)-
               (sum(manageproj$optfirekm[(manageproj$year>(year-10)) & (manageproj$year!=year) & (manageproj$iteration==i)]) + 
                subdf))/(histburned*10) < (-0.05)) {
            mgmtbound <- mgmtbound +.01
          }
          if (((histburned*10)-
               (sum(manageproj$optfirekm[(manageproj$year>(year-10)) & (manageproj$year!=year) & (manageproj$iteration==i)]) + 
                subdf))/(histburned*10) >= (-0.05)) {
            ExitFlag<-1
          }}
        if (mgmtbound >= 12) {ExitFlag<-1}
      }#this closes the while loop with the exit flag, adjusting management 
    }  
    
    if (mgmt_modern_baseline==1) {
      #choose ecological bound so that 10 year average does not exceed 5% of modern burned area 
      #for the first 9 simulation years, use average of shorter period 
      #if current burned area more than 5% larger than modern, increase mgmt
      if (year<10) {
        if (((modernburned*year)-
             (sum(manageproj$optfirekm[(manageproj$year<year) (manageproj$iteration==i)]) + 
              sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)]) 
             ))/(modernburned*year) < (-0.05)) {
          mgmtbound <- mgmtbound +.01
        }
        if (((modernburned*year)-
             (sum(manageproj$optfirekm[(manageproj$year<year) & (manageproj$iteration==i)]) + 
              sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)]) 
             ))/(modernburned*year) >= (-0.05)) {
          ExitFlag<-1
        }}
      if (year>=10) {
        if (((modernburned*10)-
             (sum(manageproj$optfirekm[(manageproj$year>(year-10)) & (manageproj$year!=year) & (manageproj$iteration==i)]) +
              sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)])))/(modernburned*10) < (-0.05)) {
          mgmtbound <- mgmtbound +.01
        }
        if (((modernburned*10)-
             (sum(manageproj$optfirekm[(manageproj$year>(year-10)) & (manageproj$year!=year) & (manageproj$iteration==i)]) +
              sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)])))/(modernburned*10) >= (-0.05)) {
          ExitFlag<-1
        }}
      
      while (ExitFlag==0) {
        #adjust fire sizes: 1% increase in spending -> 0.2063% decrease in burned area (Phillips et al forthcoming)
        #weighted by mgmtfactor, which is based on coverage of different FMZs
        subdf<-0
        foreach(id = unique(manageproj$ID[(manageproj$year==year) & (manageproj$iteration==i)]), .packages=c("sf", "tidyverse")) %do% {
          #print(paste0('This ID = ', id))
          #print(paste0('This mult = ', (1-0.002063)^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])))
          
          if((1-0.002063)^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])>0.05){
            subdf <- subdf +
              manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)][1]*
              (1-0.002063)^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])
            #print(paste0('This subdf = ', subdf))
          }}
        #if current burned area close to modern, or smaller than modern, change exit flag
        if (year<10) {
          if (((modernburned*year)-
               (sum(manageproj$optfirekm[(manageproj$year<year) & (manageproj$iteration==i)]) + subdf))/(modernburned*year) < (-0.05)) {
            mgmtbound <- mgmtbound +.01
          }
          if (((modernburned*year)-
               (sum(manageproj$optfirekm[(manageproj$year<year) & (manageproj$iteration==i)]) + subdf))/(modernburned*year) >= (-0.05)) {
            ExitFlag<-1
          }}
        if (year>=10) {
          if (((modernburned*10)-
               (sum(manageproj$optfirekm[(manageproj$year>(year-10)) & (manageproj$year!=year) & (manageproj$iteration==i)]) + 
                subdf))/(modernburned*10) < (-0.05)) {
            mgmtbound <- mgmtbound +.01
          }
          if (((modernburned*10)-
               (sum(manageproj$optfirekm[(manageproj$year>(year-10)) & (manageproj$year!=year) & (manageproj$iteration==i)]) + 
                subdf))/(modernburned*10) >= (-0.05)) {
            ExitFlag<-1
          }}
        if (mgmtbound >= 12) {ExitFlag<-1}
      }#this closes the while loop with the exit flag, adjusting management 
    }  
    
    if (sum(mgmt_historical_baseline + mgmt_modern_baseline)==0) {
      #if 10 year average of burned area more than 5% larger than historical, increase mgmt
      if (year<10) {
        if(((histburned*year)-sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)])-
            sum(manageproj$optfirekm[(manageproj$iteration==i) & (manageproj$year<year)]))/(histburned*year)<(-0.05)) {
          mgmtbound <- mgmtbound +.01
        }
        if(((histburned*year)-sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)])-
            sum(manageproj$optfirekm[(manageproj$iteration==i) & (manageproj$year<year)]))/(histburned*year)>=(-0.05)) {
          ExitFlag<-1
        }
      }
      if (year>=10) {
        if(((histburned*10)-sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)])-
            sum(manageproj$optfirekm[(manageproj$iteration==i) & (manageproj$year<year) & (manageproj$year>(year-10))]))/(histburned*10)<(-0.05)) {
          mgmtbound <- mgmtbound +.01
        }
        if(((histburned*10)-sum(manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i)])-
            sum(manageproj$optfirekm[(manageproj$iteration==i) & (manageproj$year<year) & (manageproj$year>(year-10))]))/(histburned*10)>=(-0.05)) {
          ExitFlag<-1
        } 
      }
      while (ExitFlag==0) {
        #adjust fire sizes: 1% increase in spending -> 0.2063% decrease in burned area (Phillips et al forthcoming)
        #weighted by mgmtfactor, which is based on coverage of different FMZs
        subdf<-0
        foreach(id = unique(manageproj$ID[(manageproj$year==year) & (manageproj$iteration==i)]), .packages=c("sf", "tidyverse")) %do% {
          #print(paste0('This ID = ', id))
          #print(paste0('This mult = ', (1-0.002063)^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])))
          
          if (full_mgmt_effectiveness==1) {
            if((1-0.002063)^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])>0.05) {        
              subdf <- subdf +
                manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)][1]*
                (1-0.002063)^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])                 
            } 
          }
          
          if (half_mgmt_effectiveness==1){
            if((1-(0.002063*0.5))^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])>0.05){    
              subdf <- subdf +
                manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)][1]*
                (1-(0.002063*0.5))^((mgmtbound-1)*100*full_grid$mgmtfactor[full_grid$ID==id])              
            } #closes if
          }
        } #closes foreach
        #add previous years' burned area
        if (year<10) {
          subdf <- subdf + sum(manageproj$optfirekm[(manageproj$iteration==i) & (manageproj$year<year)])
        }
        if (year>=10) {
          subdf <- subdf + sum(manageproj$optfirekm[(manageproj$iteration==i) & (manageproj$year<year) & (manageproj$year>(year-10))])
        }
        #if current burned area close to historical, or smaller than historical, change exit flag
        if(((histburned*10)-subdf)/(histburned*10)>=(-0.05)) {
          ExitFlag <- 1
        }
        if(((histburned*10)-subdf)/(histburned*10)<(-0.05)) {
          mgmtbound <- mgmtbound +.01
        }
        if (mgmtbound >= 12) {ExitFlag<-1}
      }#this closes the while loop with the exit flag, adjusting management 
    }
    
    #optimize 
    if(sum(mgmt_x1+mgmt_x2+mgmt_historical_baseline+mgmt_modern_baseline)==0){
      if(mgmtbound>1){
        xmin<-optimize(costfunc, interval = c(1,mgmtbound), tol = 0.001)
        #store optimal value, bound value, cost
        manageproj$mgmtlvl[(manageproj$year==year) & (manageproj$iteration==i)] <- xmin$minimum
        manageproj$cost[(manageproj$year==year) & (manageproj$iteration==i)] <- xmin$objective
        manageproj$mgmtbound[(manageproj$year==year) & (manageproj$iteration==i)] <- mgmtbound
      }
      if(mgmtbound==1){
        #store optimal value, bound value, cost
        manageproj$mgmtlvl[(manageproj$year==year) & (manageproj$iteration==i)] <- 1
        manageproj$cost[(manageproj$year==year) & (manageproj$iteration==i)] <- costfunc(mgmtlvl = 1)
        manageproj$mgmtbound[(manageproj$year==year) & (manageproj$iteration==i)] <- 1  
      } 
    }  
    
    #set mgmt to a vertain level
    #store optimal value, bound value, cost
    if(mgmt_x1==1){
      mgmtlvl<-1
      manageproj$mgmtlvl[(manageproj$year==year) & (manageproj$iteration==i)] <- mgmtlvl
      manageproj$mgmtbound[(manageproj$year==year) & (manageproj$iteration==i)] <- mgmtbound
      manageproj$cost[(manageproj$year==year) & (manageproj$iteration==i)] <- costfunc(mgmtlvl = 1)
    }
    if(mgmt_x2==1){
      mgmtlvl<-2
      manageproj$mgmtlvl[(manageproj$year==year) & (manageproj$iteration==i)] <- mgmtlvl
      manageproj$mgmtbound[(manageproj$year==year) & (manageproj$iteration==i)] <- mgmtbound
      manageproj$cost[(manageproj$year==year) & (manageproj$iteration==i)] <- costfunc(mgmtlvl = 2)
    }  
    if (mgmt_historical_baseline==1) {
      mgmtlvl<-mgmtbound
      manageproj$mgmtlvl[(manageproj$year==year) & (manageproj$iteration==i)] <- mgmtlvl
      manageproj$mgmtbound[(manageproj$year==year) & (manageproj$iteration==i)] <- mgmtbound
      manageproj$cost[(manageproj$year==year) & (manageproj$iteration==i)] <- costfunc(mgmtlvl = mgmtbound)
    }
    
    if (mgmt_modern_baseline==1) {
      mgmtlvl<-mgmtbound
      manageproj$mgmtlvl[(manageproj$year==year) & (manageproj$iteration==i)] <- mgmtlvl
      manageproj$mgmtbound[(manageproj$year==year) & (manageproj$iteration==i)] <- mgmtbound
      manageproj$cost[(manageproj$year==year) & (manageproj$iteration==i)] <- costfunc(mgmtlvl = mgmtbound)
    }
    
    #update optfirekm to store new firekm in each grid cell under optimized mgmt
    foreach(id = unique(manageproj$ID[(manageproj$year==year) & (manageproj$iteration==i)]), .packages=c("sf", "tidyverse")) %do% {
      manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)] <- 
        manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)]*
        (1-0.002063)^((manageproj$mgmtlvl[(manageproj$year==year) & (manageproj$iteration==i)]-1)*100*full_grid$mgmtfactor[full_grid$ID==id])*full_mgmt_effectiveness+                 
        manageproj$fire_sum[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)]*
        (1-(0.002063*0.5))^((manageproj$mgmtlvl[(manageproj$year==year) & (manageproj$iteration==i)]-1)*100*full_grid$mgmtfactor[full_grid$ID==id])*half_mgmt_effectiveness             
      if(manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)]<=0) {
        manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i) & (manageproj$ID==id)]<-0.0001 
        #don't want to divide by zero for FRI
      }
    }
    
    #update socialcost to reflect emissions x scc
    manageproj$socialcost[(manageproj$year==year) & (manageproj$iteration==i)] <- 
      ((sum(manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i)])-histburned)*12191.67*0.84*
         SCCoutput$SCC[SCCoutput$Year==2019+year]*(237.002/218.076)*full_scc +          
         (sum(manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i)])-histburned)*12191.67*0.84*
         SCCoutput$SCC[SCCoutput$Year==2019+year]*(237.002/218.076)*0.5*half_scc  
      )#*(sum(manageproj$optfirekm[(manageproj$year==year) & (manageproj$iteration==i)]) > histburned) 
    #the last term is logical, social costs only !=0 if burned area greater than historical 
    
    
  }#this closes the for loop, going through years
  
  #save fireproj file
  #save fireproj file, with different options for each scenario
  #scenario 1: optimized mgmt (full effectiveness and full SCC)
  if (full_mgmt_effectiveness==1 & full_scc==1 & sum(mgmt_historical_baseline+mgmt_modern_baseline+mgmt_x1+mgmt_x2+ccadjusted)==0) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_bound_each_year", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_bound_each_year", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_bound_each_year", paste0(i, '.csv'))))
  }
  #scenario 2: optimized management, management half as effective, full SCC
  if (half_mgmt_effectiveness==1 & full_scc==1 & sum(mgmt_historical_baseline+mgmt_modern_baseline+mgmt_x1+mgmt_x2+ccadjusted)==0) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_half_effectiveness", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_half_effectiveness", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_half_effectiveness", paste0(i, '.csv'))))
  }
  #scenario 3: optimized management, half SCC
  if (full_mgmt_effectiveness==1 & half_scc==1 & sum(mgmt_historical_baseline+mgmt_modern_baseline+mgmt_x1+mgmt_x2+ccadjusted)==0) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_half_SCC", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_half_SCC", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_half_SCC", paste0(i, '.csv'))))
  }
  #scenario 4: maintain historical baseline burning 
  if (mgmt_historical_baseline==1) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_historical_levels", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_historical_levels", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_historical_levels", paste0(i, '.csv'))))
  }
  #scenario 5: maintain modern baseline burning 
  if (mgmt_modern_baseline==1) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_modern_levels", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_modern_levels", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_modern_levels", paste0(i, '.csv'))))
  }
  #scenario 6: double management 
  if (mgmt_x2==1) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_x2", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_x2", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_x2", paste0(i, '.csv'))))
  }
  #scenario 7: maintain current levels of management 
  if (mgmt_x1==1) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_no_mgmt_change", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_no_mgmt_change", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_no_mgmt_change", paste0(i, '.csv'))))
  }
  #scenario 8: minimum burned area climate scenario
  if (ccmin==1) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_ccmin", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_ccmin", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_ccmin", paste0(i, '.csv'))))
  }
  #scenario 9: 25th pct burned area climate scenario
  if (cc25==1) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_cc25", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_cc25", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_cc25", paste0(i, '.csv'))))
  }
  #scenario 10: median burned area climate scenario
  if (ccmed==1) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_ccmed", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_ccmed", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_ccmed", paste0(i, '.csv'))))
    #scenario 11: 75th pct burned area climate scenario
  }
  if (cc75==1) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_cc75", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_cc75", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_cc75", paste0(i, '.csv'))))
  }
  #scenario 12: maximum burned area climate scenario
  if (ccmax==1) {
    if (file.exists(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_ccmax", paste0(i, '.csv'))))
    {file.remove(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_ccmax", paste0(i, '.csv')))}
    write_csv(manageproj, file.path(file.path("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/proj_mgmt_ccmax", paste0(i, '.csv'))))
  }
}#this closes the iteration loop, going through the whole thing 


ggplot(manageproj, aes(x=mgmtbound, y=mgmtlvl)) +
  geom_point() 
  #+ geom_smooth(method = lm) 

#save all runs together
optmgmt <- do.call(rbind, lapply(paste0(out_path,"/",list.files(path = out_path)), read.csv))
if (file.exists('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/optmgmt.csv')){file.remove('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/optmgmt.csv')}
write_csv(optmgmt, file.path('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/optmgmt.csv'))




















  






