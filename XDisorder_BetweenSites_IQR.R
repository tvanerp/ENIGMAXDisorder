# Cross Dx Outlier Script

###### SETTING UP ######
library(ggplot2)
library(dplyr)
library(psych)
library(tidyverse)
library(ggpubr)

# set working directory where the data is & where you want plots to end up
wrkdir <- "~/Desktop/XDisorder/"
iqr_range <- 2.5 # this is the IQR range you would like to use

options(pillar.sigfig=20) # sets dplyr to 20 significant figures to increase accuracy
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")
dir.create(paste0(wrkdir, "IQR_OutlierSites"))
  dir.create(paste0(wrkdir, "IQR_OutlierSites/1.5IQR"))
  dir.create(paste0(wrkdir, "IQR_OutlierSites/2.5IQR"))
  
if(iqr_range == 1.5){
  wrkdir_iqr <- paste0(wrkdir, "IQR_OutlierSites/1.5IQR/")
}else if(iqr_range == 2.5){
  wrkdir_iqr <- paste0(wrkdir, "IQR_OutlierSites/2.5IQR/")
}


all3 <- read.csv(paste0(wrkdir, "ALL2_13July2020.csv")) # read in the latest data sheet
all3$Site_new <- as.factor(all3$Site_new) # change everything to factors for plotting later
all3$Dx <- as.factor(all3$Dx)
all3$Site <- as.factor(all3$Site)

####### FUNCTIONS ######

# this function allows us to calculated the weighted sd later in the script
# adapted from a script found on GitHub
weighted.sd <- function(x, wt, na.rm = TRUE) {
  if (na.rm) {
    ind <- is.na(x) | is.na(wt)
    x <- x[!ind]
    wt <- wt[!ind]
  }
  wt <- wt / sum(wt)
  wm <- weighted.mean(x, wt)
  sqrt(sum(wt * (x - wm) ^ 2))
}

####### CREATE WEIGHTS FOR WEIGHTED MEAN/SD #####
weights_grouped <- all3 %>% count(Site, Dx)
all3$Weights <- NA

for(i in 1:length(weights_grouped$n)){
  index <- which((all3$Site == weights_grouped$Site[i]) & (all3$Dx == weights_grouped$Dx[i]))
  all3[index,which(colnames(all3) == "Weights")] <- weights_grouped$n[i]
}

####### WITHIN SITE OUTLIER ######

# first, get rid of individual outliers within the region of interest & plot
start <- which(colnames(all3)=="LLatVent")
end <- which(colnames(all3)=="FullSurfArea")

# As a note: this is going by 1.5*IQR WITHIN sites, through EACH of the regions of interest.
# first taking these out to replicate the plots and not to skew the mean ridiculously
regions <- start:end
#divide by site
df_list <- split(all3, as.factor(all3$Site))
#loop through regions 
for (j in df_list){
  outliers <-list()
  #run for loop for each region, get a vector of row numbers that are outliers
  for (i in regions){
    outliers[[(i-(start-1))]] <- boxplot(j[,i], plot=FALSE, range=1.5)$out
    #for each outlier, get the row number in the variable excl and then set that row to NA
    for (h in 1:length(outliers[[(i-(start-1))]])){
      excl <- which(all3[,regions[(i-(start-1))]]==outliers[[(i-(start-1))]][h])
      all3[excl,regions[(i-(start-1))]] <- NA
    }
  }
}

###### OUTLIER SITES ######
# setting up variables needed for the next for loop
site_out_index <- NA
count <- 1
# site_out_df <- as.data.frame(matrix(NA, ncol = 3)) # data frame that holds all the sites that are outliers
# colnames(site_out_df) <- c("Site", "Dx", "ROI")
# 
site_out_df <- as.data.frame(matrix(NA, ncol = 2)) # data frame that holds all the sites that are outliers
colnames(site_out_df) <- c("Site", "ROI")

for(i in start:end){
  
  current <- colnames(all3[i]) # current ROI we are working with
  
  #site_mean <- as.data.frame(all3 %>% select(2, i, 238, 239, 270) %>% na.omit() %>%
   #                 group_by(Site, Dx) %>% summarise_at(current, funs(weighted.mean(., w = Weights, na.rm = TRUE)))) # take the mean by site and dx for ROI

  # this is the weighted mean
  
  
  # site_mean <- as.data.frame(all3 %>% select(2,i,238,239,270) %>% na.omit() %>% group_by(Site, Dx) %>% summarise_at(current, mean))
  # this is the regular mean by site
  
  site_mean <-  site_mean <- as.data.frame(all3 %>% select(2,i,238,239,270) %>% na.omit() %>% group_by(Site) %>% summarise_at(current, mean))
  # JUST by group

  
  site_outliers <- boxplot(site_mean[,current], plot = FALSE, range = iqr_range)$out # determine outliers
  
  if(length(site_outliers) == 0){ # if there were no outliers were detected
  
    # histogram of values that have no outliers
    histzero <- ggplot(all3, aes(x=all3[,i])) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666")+
      xlab(paste0(current))
    
    filenamezero <- paste0(wrkdir_iqr, colnames(all3)[i], '_hist_nooutliers.png')
    ggsave(filename = filenamezero, plot = histzero, width = 16, height = 9)
    
    # box plot of no outlier sites
    plotzero <- ggplot(all3, aes(x=Site, y=all3[,i], fill=Dx)) + 
      scale_fill_manual(values=cbPalette)+
      geom_boxplot(outlier.colour="black") +
      scale_x_discrete(breaks = levels(all3$Site)[c(T,F,F)]) +
      xlab("Site")+
      ylab(paste0("Units for ", current))

      filenamezero2 <- paste0(wrkdir_iqr, colnames(all3)[i], '_box_nooutliers.png')
      ggsave(filename = filenamezero2, plot = plotzero, width = 16, height = 9)
    
      }else{

  
  site_out_index <- NA # what index in the data frame has the outliers?
  site_mean_noout <- site_mean
  
  for(j in 1:length(site_outliers)){ # what index in the site mean data frame has the outliers?
    site_out_index[j] <- which(site_mean[,current] == site_outliers[[j]])
    
  }

site_mean_noout[site_out_index,current] <- NA # mark the outliers as NA

  # append to the current dataframe to keep a running list by ROI of all the site / dx outliers
  # site_out_df <- site_out_df %>% add_row(Site = site_mean$Site[site_out_index], Dx = site_mean$Dx[site_out_index], ROI = current)
  site_out_df <- site_out_df %>% add_row(Site = site_mean$Site[site_out_index], ROI = current)

  site_out_ROI <- site_out_df %>% filter(ROI == current) # pull the outlier sites just for the current ROI
  
  # take out the newly defined outliers to plot: ALL/raw values, not the means
  
  takeout <- NA
  all_noout <- all3
  
  for(k in 1:dim(site_out_ROI)[1]){
    # takeout <- which(all_noout$Site == site_out_ROI[k,1] & all_noout$Dx == site_out_ROI[k,2])
    takeout <- which(all_noout$Site == site_out_ROI[k,1])
    all_noout[takeout,i] <- NA 
  }
  
  # print what it looked like before: boxplot
  range_axis <- range(all3[,i], na.rm = TRUE)
  
  plot <- ggplot(all3, aes(x=Site, y=all3[,i], fill=Dx)) + 
    scale_fill_manual(values=cbPalette)+
    geom_boxplot(outlier.colour="black") +
    scale_x_discrete(breaks = levels(all3$Site)[c(T,F,F)]) +
    xlab("Site")+
    ylab(paste0("Units for ", current))+
    lims(y = range_axis)
  
  # print what it looks like now: boxplot
  plot2 <- ggplot(all_noout, aes(x = Site, y = all_noout[,i], fill = Dx))+
    scale_fill_manual(values = cbPalette)+
    geom_boxplot(outlier.color = "black")+
    scale_x_discrete(breaks = levels(all_noout$Site)[c(T,F,F)])+
    xlab("Site")+
    ylab(paste0("Units for ", current))+
    lims(y = range_axis)
  
  # histogram of raw values before

  hist1 <- ggplot(all3, aes(x=all3[,i])) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666")+
    xlab(paste0(current))+
    lims(x = range_axis)+
    ylab("Raw Value Density")
  
  # histogram of raw values now
  hist2 <- ggplot(all_noout, aes(x = all_noout[,i]))+
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666")+
    xlab(paste0(current, " no outlier sites"))+
    lims(x = range_axis)
  
  # histogram of the means before (in #C1D0D0)
  hist3 <- ggplot(site_mean, aes(x = site_mean[,current]))+
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.6, fill="#C1D0D0")+
    xlab(paste0(current))+
    lims(x = range_axis)+
    ylab("Mean Density")
    
  # histogram of the means after
  hist4 <- ggplot(site_mean_noout, aes(x = site_mean_noout[,current]))+
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.6, fill="#C1D0D0")+
    xlab(paste0(current, " no outlier sites"))+
    lims(x = range_axis)
  
  filename <- paste0(wrkdir_iqr, current, '_outliers.png')
  filename2 <- paste0(wrkdir_iqr, current, '_NOoutliers.png')
  
  ggsave(filename = filename, plot = plot, width = 16, height = 9)
  ggsave(filename = filename2, plot = plot2, width = 16, height = 9)
  
  filename3 <- paste0(wrkdir_iqr, current, '_histograms.png')
  histograms <- ggarrange(hist3, hist4, hist1, hist2, nrow = 2, ncol = 2)
  ggsave(filename = filename3, plot = histograms, width = 16, height = 16)
  
} # finish the else statement
} # finish the for loop

# save out dataframe that has all of the site exclusion information
write.csv(site_out_df, paste0(wrkdir_iqr, iqr_range, '_SiteOutlier_List.csv'))

site_out_concat <- site_out_df %>% 
  group_by(ROI) %>%
  #summarise(OutlierSites = toString(Site), OutliersDx = toString(Dx)) %>% 
  summarise(OutlierSites = toString(Site)) %>%
  ungroup()

write.csv(site_out_concat, paste0(wrkdir_iqr, iqr_range, '_SiteToCompare.csv'))

####### VARIABILITY WITHIN SITE ######
novar <- as.data.frame(matrix(NA, ncol = 1))

dir.create(paste0(wrkdir, 'Variance_Within_Site/')) # create new directory for this analysis
  dir.create(paste0(wrkdir, "Variance_Within_Site/1.5IQR"))
  dir.create(paste0(wrkdir, "Variance_Within_Site/2.5IQR"))

if(iqr_range == 1.5){
  var_iqr <- paste0(wrkdir, "Variance_Within_Site/1.5IQR/")
}else if(iqr_range == 2.5){
  var_iqr <- paste0(wrkdir, "Variance_Within_Site/2.5IQR/")
}

filename_novar <- paste0(var_iqr, iqr_range, '_NoVarianceList.txt') # create a text file that will contain sites with ZERO variance (which is abnormal)
file.create(filename_novar) # create the no variance file

# function to create spaces in the text file 
tab_newline <- function(){
  write(paste0(" "), file = filename_novar, append = TRUE)
}

#site_out_df_sd <- data.frame(matrix(NA, ncol = 3)) # dataframe to save sites that have wild variability
#colnames(site_out_df_sd) <- c("Site", "Dx", "ROI") # change column names

site_out_df_sd <- data.frame(matrix(NA, ncol = 2))
colnames(site_out_df_sd) <- c("Site", "ROI")

for(stdev in start:end){
  
  tab_newline()
  
  current <- colnames(all3[stdev]) # current ROI

  # weighted_sd <- all3 %>% group_by(Site, Dx) %>% 
    #summarise_at(current, funs(weighted.sd(., w = Weights, na.rm = TRUE))) # take the weighed SD
  # colnames(weighted_sd) <- c("Site", "Dx", "weight_site_var")
  
  weighted_sd <- as.data.frame(all3 %>% group_by(Site) %>% summarise_at(current, sd, na.rm = TRUE))
  colnames(weighted_sd) <- c("Site", current)
  
  # define outliers outside IQR
  sd_outliers <- unique(boxplot(weighted_sd[,current], plot = FALSE, range = iqr_range)$out)

if(length(sd_outliers) == 0){ # if there are no outliers
  hist5 <- ggplot(weighted_sd, aes(x=weighted_sd[,current])) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.6, fill="#C1D0D0")+
    xlab(paste0(current))
  
  filename_hist5 <- paste0(var_iqr, current, '_variance.png')
  ggsave(filename = filename_hist5, plot = hist5, width = 9, height = 9)
    
  
  }else{
 
  no_outlier_sd <- weighted_sd

  for(m in 1:length(sd_outliers)){ # indices of the sd outliers
    outlier_idx <- which(weighted_sd[,current] == sd_outliers[m])
    #site_out_df_sd <- site_out_df_sd %>% add_row(Site = weighted_sd$Site[outlier_idx], Dx = weights_grouped$Dx[outlier_idx], ROI = current)
    site_out_df_sd <- site_out_df_sd %>% add_row(Site = weighted_sd$Site[outlier_idx], ROI = current)
    no_outlier_sd[outlier_idx,current] <- NA
    
  }

  
  if(any(sd_outliers == 0)){ # if any of the variances were 0 write out to a text file
    index <- which(weighted_sd[,current] == 0)
    sentence <- paste0("For the ROI ", current, ", site ", weighted_sd$Site[index], " with dx ", weighted_sd$Dx[index], " had variance 0 and should be checked.")
    
    novar <- write(paste0(sentence), file = filename_novar, append = TRUE)
    
  } # end of if statement for variances of 0

site_out_ROI <- site_out_df_sd %>% filter(ROI == current)
    
# we can either graph the distribution of variances OR of all of the values.... I am going to go with distribution of variances for this one. Previously we were graphing the overall distribution of the values

noout_all3 <- all3

    # take out sites that were marked as outliers in the entire dataset to plot raw values
    for(all in 1:dim(site_out_ROI)[1]){
      
    #tidx <- which(noout_all3$Site == site_out_ROI$Site[all] & noout_all3$Dx == site_out_ROI$Dx[all])
    tidx <- which(noout_all3$Site == site_out_ROI$Site[all])
    noout_all3[tidx,i] <- NA
    
    }

    
range_x <- range(weighted_sd[,current], na.rm = TRUE) 
range_x2 <- range(all3[,stdev], na.rm = TRUE)

   hist_sd_uncorrect <-  ggplot(weighted_sd, aes(x=weighted_sd[,current])) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.6, fill="#C1D0D0")+
      xlab(paste0(current))+
      ylab("Variance Density")+
      lims(x = range_x)
    # histogram of weighted SD uncorrected
    
    hist_sd_correct <-  ggplot(no_outlier_sd, aes(x=no_outlier_sd[,current])) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.6, fill="#C1D0D0")+
      xlab(paste0(current, " No Outliers"))+
      lims(x = range_x)
    # histogram of the variances between sites
    
    hist_raw1 <-  ggplot(all3, aes(x=all3[,stdev])) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666")+
      xlab(paste0(current))+
      lims(x = range_x2)+
      ylab("Raw Values Density")

   hist_raw2 <- ggplot(noout_all3, aes(x=noout_all3[,stdev]))+
     geom_histogram(aes(y=..density..), color = "black", fill = "white")+
     geom_density(alpha=.2, fill="#FF6666")+
     xlab(paste0(current, " No Outliers"))+
     lims(x = range_x2)
   
histplots <- ggarrange(hist_sd_uncorrect, hist_sd_correct, hist_raw1, hist_raw2,  nrow = 2, ncol = 2)
    filename_hist <- paste0(var_iqr, current, "_variance.png")
    ggsave(filename = filename_hist, plot = histplots, width = 16, height = 16)
    
  } # end the else statement
} # end the for loop

write.csv(site_out_df_sd, paste0(var_iqr, iqr_range, "_VarianceOutliers.csv"))

site_out_concat_sd <- site_out_df_sd %>% 
  group_by(ROI) %>%
  #summarise(OutlierSites = toString(Site), OutliersDx = toString(Dx)) %>% 
  summarise(OutlierSites = toString(Site)) %>%
  ungroup()

write.csv(site_out_concat_sd, paste(var_iqr, iqr_range, "_VarianceOutliers_ToCompare.csv"))


##### NOTES #####
# run both 1.5 and 2.5 to see how the numbers being taken out compare: will take a while, so start early.
# the weighted SD and weighted mean analyses are SEPARATE, not building upon one another

###### PROBLEMS ######


##### HOW TO IMPROVE THE CODE #####
# combine everything into one loop to make code a bit more efficient

