#this script calculates a parameter for each grid cell so that simulations match historical average burned area

library(tidyverse)  
library(sf)     
library(doSNOW)
library(parallel)
library(raster)     
library(fasterize)  
library(rgdal)  
library(ggplot2)


#set desired output pathway
out_path <- "C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/matching"
dir.create(out_path, recursive = T)


#set things up to go parallel
cores <- detectCores()-2
cl <-makeCluster(cores)
registerDoSNOW(cl)


#load dataframes from ignitions.R --
  #ak_intersect is a (somewhat cleaned) full data set
  #full_grid has one observation per grid cell 
  #exp_initial_guess gives a guess for each grid cell's exponent from previous model output
full_grid <- read_sf("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/full_grid.shp")
ak_intersect <- read_sf("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/ak_intersect.shp")
exp_initial_guess <- read_csv("C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/exp_initial_guess.csv")
prob <- read_csv('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/ProbAdjust.csv')

#generate a distribution of fire size for all Alaska
ggplot(as.data.frame(ak_intersect), aes(x = firekm)) + xlim(0, 500) +
  ylim(0,300) + geom_histogram(color = '#E69F00') + geom_density(color = '#56B4E9')
#km_cdf <- ecdf(ak_intersect$firekm)
  #ecdf = empirical cumulative distribution function
  plot.ecdf(ak_intersect$firekm, col.hor = "#56B4E9", 
            ylab = "Probability", xlab = "Square Km")

  
#generate a CDF for all fires in the dataset
#syntax from https://stackoverflow.com/questions/38537311/how-to-find-quantiles-of-an-empirical-cumulative-density-function-ecdf
firesorted <- sort(ak_intersect$firekm)
e_cdf <- 1:length(firesorted)/length(firesorted)
plot(firesorted, e_cdf, type = "l", xlim = c(0, 1000), ylab = "Cumulative Probability",
     xlab = "Fire Size, Square Kilometers")
      #, main = "CDF of All Historical Fire Sizes in Alaska"
#this line takes a probability cutoff and gives the next largest fire
firesorted[which(e_cdf >= 0.9999627)[1]]

#same plot, but with different size labels and fire size in thousands of hectares
firesorted_thouHA <- firesorted/10
cdf_df <- data.frame(firesorted_thouHA, e_cdf, firesorted)
cdf_df <- cdf_df %>% mutate(filling=1)
ggplot(data = cdf_df, aes(x=firesorted, y=e_cdf)) + 
  geom_line(color="#440154FF", size=1.5) +
  labs(x = "Fire Size (ha)",
       y = "Cumulative Probability") +
  scale_x_continuous(trans="log10", breaks = c(1, 10, 100, 1000)) +
  theme_bw()
  
#look at each grid cell in each year
#each year, there will be the median number of ignitions for that grid cell
#each ignition will cause a fire of random size, drawn from cdf
#modify where draw is from to match burned area
  
#make empty df - will be populated as we run the model  
#headers <- c("ID","exponent", "iterations")
#final_df <- as.data.frame(matrix(,ncol=3,nrow=0))
#names(final_df)<-headers

#set some initial parameters and dataframe
iter_years <- 10000
test <- ak_intersect %>% group_by(ID, FIRESEA) %>% summarize(ig_by_year = n()) 
test <- test %>% mutate(FIRESEA=as.integer(FIRESEA))
histyearfirst = 1970      
histyearlast = 2019
all_years <- data.frame('FIRESEA'= histyearfirst:histyearlast)

#loop through each grid cell, or choose one to practice on
#foreach (id = unique(full_grid$ID), .packages=c("sf", "tidyverse")) %do% {
foreach (id = unique(filter(full_grid, ID ==47)$ID), .packages=c("sf", "tidyverse")) %do% {
  
  print(paste0('This ID = ', id))
  
  #a dataframe with just a single id in it
  sub_df <- full_grid %>% dplyr::filter(ID == id)
  
  #add columns for exponent
  sub_df <- sub_df %>% mutate(exponent=filter(exp_initial_guess, ID==id)$exponent)
  
  
  #decide annual ignitions
    #just doing median 
      #if (sub_df$med==0) {sub_df$ignitions <- 1}
      #if (sub_df$med>=1){sub_df$ignitions <-as.integer(sub_df$med)}
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
    
    #make a list iter_years long
    ig_per_year <- list()
    ig_per_year <- ig_dist$ig_by_year[ig_dist$year==min(ig_dist$year[which(ig_dist$year >= runif(1)*50[1])])]
    #each row gives number of ignitions in that year
    #multiplying by 50 because there are 50 years in the historical sample at this point 
    #if historical data is narrowed, 50 will need to change !!
    while (length(ig_per_year) < iter_years) {
      ig_per_year[length(ig_per_year)+1] <- 
        ig_dist$ig_by_year[ig_dist$year==min(ig_dist$year[which(ig_dist$year >= runif(1)*50[1])])]
    }
  
  
  #expand the current sub_df by the total number of ignitions 
        #sub_df <- sub_df[rep(seq_len(nrow(sub_df)), each = sub_df$ignitions*iter_years), ]
  sub_df <- sub_df[rep(seq_len(nrow(sub_df)), each = sum(ig_per_year)), ]
  
  #add year and firekm columns
  sub_df <- sub_df %>% mutate(year=0, firekm=0)
  #add a column to identify each year
  #the number year observations is determined by the corresponding ignition number in the x spot in the ig_per_year list
        #sub_df$year <- rep(1:iter_years, max(sub_df$ignitions))
  sub_df$year <- rep(1:iter_years, ig_per_year)
  
  #populate fire sizes for each ignition, by random draw, and adjust if FRI is low 
  sub_df$firekm <- replicate(sum(ig_per_year), firesorted[which(e_cdf >= runif(1)^sub_df$exponent[1])[1]])
  if (full_grid[full_grid$ID==id,]$FRI < 84.02) {
    sub_df$firekm <-  sub_df$firekm*prob$ProbAdjustment[which(prob$FRI >= full_grid[full_grid$ID==id,]$FRI)[1]]
  }
  
  #calculate annual burned area 
  sub_df <- sub_df %>% group_by(year) %>% mutate(sumkm = sum(firekm))
  
  #calculate average annual burned area
  #squash down by year
  sub_df <- sub_df %>% group_by(year) %>% summarize(ID=max(ID), sumkm=max(sumkm), year=max(year), 
                                                     exponent=max(exponent))
  #squash down completely, keeping only average annual burned area and grid cell ID 
  sub_df <- sub_df %>% summarise(ID = max(ID), avg_km = sum(sumkm)/iter_years,
                                 exponent=max(exponent))
  
  #extract historical average annual burned area for comparison to simulated 
  sub_df <- sub_df %>% mutate(hist_km = as.data.frame(full_grid)[full_grid$ID==id,]$ann_km)
  
  #set the exit flag, which means when the model should stop or run criteria, if exit flag is zero keep the model running, if it is 1 it will stop
  sub_df <- sub_df %>% mutate(ExitFlag = 0)
  
  #count number of iterations
  sub_df <- sub_df %>% mutate(IterNum = 1)
  
  #empty list to store hist_km - avg_km
  sub_df <- sub_df %>% mutate(current_km_diff = 1)
  
  
  #keep looping until ExitFlag=1
  while (sub_df$ExitFlag==0) {
    
    #store the current hist_km - avg_km if IterNum == 1
    if (as.numeric(sub_df$IterNum) == 1){sub_df$current_km_diff= (sub_df$hist_km - sub_df$avg_km)}
    
    #increase iteration number by 1
    sub_df$IterNum <- sub_df$IterNum + 1
    
    #if the calculated km is close to real km then stop iteration by changing the exit to 1
    #if it is 0 loop will iterate again
    #for now, "close" means within 10% of historical
    sub_df$ExitFlag <- ifelse(abs((sub_df$hist_km - sub_df$avg_km)/sub_df$avg_km) <= 0.1, 1, 0)
    
    #assuming convergence not achieved, adjust exponent on random draw
    if ((sub_df$hist_km - sub_df$avg_km)>0) {sub_df$exponent <- sub_df$exponent-.025}
        #simulated km too small --> make exponent smaller
    if ((sub_df$hist_km - sub_df$avg_km)<0) {sub_df$exponent <- sub_df$exponent+.025}
        #simulated km too large --> make exponent larger
    
    #rerun the simulation with new adjustment parameter to see if it matches
    #expand by years x ignitions
    sub_df <- sub_df[rep(seq_len(nrow(sub_df)), each = sum(ig_per_year)), ]
    
    #add a column to identify each year
    sub_df <- sub_df %>% mutate(year=0, firekm=0)
    sub_df$year <- rep(1:iter_years, ig_per_year)
    
    #populate fire sizes for each ignition, by random draw
    sub_df$firekm <- replicate(sum(ig_per_year), firesorted[which(e_cdf >= runif(1)^sub_df$exponent[1])[1]])
    if (full_grid[full_grid$ID==id,]$FRI < 84.02) {
      sub_df$firekm <-  sub_df$firekm*prob$ProbAdjustment[which(prob$FRI >= full_grid[full_grid$ID==id,]$FRI)[1]]
    }
    
    #calculate annual burned area 
    sub_df <- sub_df %>% group_by(year) %>% mutate(sumkm = sum(firekm))
    
    #calculate average annual burned area
    #squash down by year
    sub_df <- sub_df %>% group_by(year) %>% summarize(ID=max(ID), sumkm=max(sumkm), year=max(year), hist_km = max(hist_km),
                                                       exponent=max(exponent), IterNum = max(IterNum))
    #squash down completely, keeping only average annual burned area and grid cell ID 
    #average annual burned area is the sum of all burned area, divided by total years. 
        #not mean(sumkm), because some years have no fires
    sub_df <- sub_df %>% summarise(ID = max(ID), avg_km = sum(sumkm)/iter_years,  hist_km = max(hist_km),
                                    exponent=max(exponent), IterNum = max(IterNum))
    
    #check exit flag with new simulated mean burned area
    sub_df$ExitFlag <- ifelse(abs((sub_df$hist_km - sub_df$avg_km)/sub_df$avg_km) <= 0.1, 1, 0)
    
    #if IterNum > 600 make the loop exit (upper limit on how large exponent can get)
    if (sub_df$IterNum > 600){sub_df$ExitFlag <- 1}
    
    #if about to go negative, make the loop exit (lower limit)
    if (sub_df$exponent <= 0.03){sub_df$ExitFlag <- 1}
  }
  
  #store sub_df as csv
  if (file.exists(file.path(out_path, paste0(id,'.csv')))){file.remove(file.path(out_path, paste0(id,'.csv')))}
  write_csv(sub_df, file.path(out_path, paste0(id,'.csv')))    #using write_csv here preserves the geometry column 
  #store final parameter and ID in final_df
  #final_df<-rbind(final_df, c(id, exponent, IterNum))
  #colnames(final_df) <- c('ID', 'exponent', 'iterations')
}

#bind all csv's together to one final df
exponents <- do.call(rbind, lapply(paste0(out_path,"/",list.files(path = out_path)), read.csv))

#make the exponents dataframe into a spatial dataframe
#keep only the ID and the spatial column from full_grid
shape = subset(full_grid, select= c(ID))
#join shape (spatial base) with exponent (variables we want)
joined = left_join(shape, exponents, by = 'ID')
joined = subset(joined, select = -c(geometry.y))
joined = st_transform(joined, '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')


#save exponents file, as shp
if (file.exists('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/exponents.shp')){file.remove('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/exponents.shp')}
write_sf(joined, file.path('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/exponents.shp'))
if (file.exists('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/exponents.csv')){file.remove('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/exponents.csv')}
write_csv(joined, file.path('C:/Users/Molly/OneDrive - Tufts/Job Stuff/WHRC/Future_Fire_Projections/molly/molly/exponents.csv'))





#linear regression of burned area on time
  #since the parameter is chosen using data from the entire sample, it actually would give results for the sample average/middle
  #which is 1990. For the future projections, we need to adjust using a linear time trend to project from a 2020 baseline
#collapse ak_intersect so there's just total burned area and fire season (year)
ak_intersect <- ak_intersect %>% mutate(lnfirekm = log(firekm))
ak_intersect <- ak_intersect %>% mutate(lnfireha = log(firekm*100))
ak_intersect$FIRESEA<-as.numeric(ak_intersect$FIRESEA)  #don't do this for boxplots
burnkmtime.lm <- lm(lnfirekm ~ FIRESEA, data = ak_intersect) 
summary(burnkmtime.lm)
pctfirekmincrease <- (exp(-0.017858)-1)*100       #small transformation to get exact interpretation 
pctfirekmincrease
  #the average fire size DECREASES by 1.77%/year, which seems to be driven by the lower bound on
    #recorded fire size being higher in 1960-1973, 1985, and 1987-1989 

burnedannual <- ak_intersect %>% group_by(FIRESEA) %>% summarise(mediankm = median(firekm), firekm=sum(firekm), count=n())
burnedannual <- burnedannual %>% mutate(meankm = firekm/count, lnmed = log(mediankm),
                                        lnfirekm = log(firekm), lncount = log(count), lnmean = log(meankm))

ignitionstime.lm <- lm(count ~ FIRESEA, data = burnedannual)
summary(ignitionstime.lm)
  #the number of annual ignitions increases by 6.489/year, intersect = -12508.111
  #ignitions increase by about 50% from 1990 to 2019
  #this means that more recent years are weighted more heavily in a random draw of fire size (but not number of ignitions)
lnignitions.lm <- lm(lncount ~ FIRESEA, data = burnedannual)
summary(lnignitions.lm)
pctcountincrease <- (exp(0.023836)-1)*100       #small transformation to get exact interpretation 
pctcountincrease
  #annual ignition increase of 2.41% on average year-to-year

#figure out where baseline is - average FIRESEA, weighted by count
burnedannual <- burnedannual %>% mutate(weight = count/sum(count))
burnedannual <- burnedannual %>% mutate(weighted_year = weight*as.numeric(FIRESEA))
sum(burnedannual$weighted_year)  
#baseline year is 1996.334

annburnedareabytime.lm <- lm(firekm ~ FIRESEA, data = burnedannual)
summary(annburnedareabytime.lm)
annualpctincrease <- (exp(0.04698)-1)*100       #small transformation to get exact interpretation 
annualpctincrease
  #annual burned area increases by 4.81% on average year-to-year

#look at each half of period separately --1970-1989, 1990-2019
###################
burnkmtime.lm <- lm(lnfirekm ~ FIRESEA, data = ak_intersect[ak_intersect$FIRESEA<1990,]) 
summary(burnkmtime.lm)
pctfirekmincrease <- (exp(-0.004302)-1)*100     
pctfirekmincrease
#in first half of period, no significant change in fire size
burnkmtime.lm <- lm(lnfirekm ~ FIRESEA, data = ak_intersect[ak_intersect$FIRESEA>1989,]) 
summary(burnkmtime.lm)
pctfirekmincrease <- (exp(0.033149)-1)*100     
pctfirekmincrease
  #in second half of period, 3.37% increase in average fire size 

lnignitions.lm <- lm(lncount ~ FIRESEA, data = burnedannual[burnedannual$FIRESEA<1990,])
summary(lnignitions.lm)
#in the first half of the period, there is no statistically significant trend in number of ignitions 
lnignitions.lm <- lm(lncount ~ FIRESEA, data = burnedannual[burnedannual$FIRESEA>1989,])
summary(lnignitions.lm)
  #in the second half of the period, there is no statistically significant trend in number of ignitions 
lnignitions.lm <- lm(lncount ~ FIRESEA, data = burnedannual)
summary(lnignitions.lm)



















            