# Cross Dx Outlier Script

###### SETTING UP ######
library(ggplot2)
library(dplyr)
library(psych)
library(tidyverse)
library(ggpubr)

options(pillar.sigfig=20)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")
dir.create("~/Desktop/Site_Outlier_Detection/1.5IQR_Sites")
dir.create("~/Desktop/Site_Outlier_Detection/2.5IQR_Sites")

all3 <- read.csv("~/Desktop/ALL2_13July2020.csv")
all3$Site_new <- as.factor(all3$Site_new)
all3$Dx <- as.factor(all3$Dx)
all3$Site <- as.factor(all3$Site)

####### FUNCTIONS ######

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

####### OUTLIER ######

# first, get rid of individual outliers within the region of interest & plot
start <- which(colnames(all3)=="LLatVent")
end <- which(colnames(all3)=="FullSurfArea")

# As a note: this is going by 1.5*IQR WITHIN sites, through EACH of the regions of interest.
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

site_out_index <- NA
count <- 1
site_out_df <- as.data.frame(matrix(NA, ncol = 3))
colnames(site_out_df) <- c("Site", "Dx", "ROI")


for(i in start:end){
  
  current <- colnames(all3[i])
  
  site_mean <- all3 %>% select(2, i, 238, 239) %>% na.omit() %>%
                    group_by(Site, Dx) %>% summarise_at(current, mean)

  
  site_outliers <- boxplot(site_mean[,3], plot = FALSE, range = 1.5)$out
  
  filename1 <- paste0(current, "_SiteOutliers.csv")
  
  
  if(length(site_outliers) == 0){
   # noout <- paste0("There were no site outliers detected in the ", current)
  #  write.csv(noout, paste0("~/Desktop/Site_Outlier_Detection/1.5IQR_Sites/", filename1), row.names = FALSE, col.names = FALSE)
    
    histzero <- ggplot(all3, aes(x=all3[,i])) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666")+
      xlab(paste0(current))
    
    filenamezero <- paste0('~/Desktop/Site_Outlier_Detection/1.5IQR_Sites/', colnames(all3)[i], 'hist_nooutliers.png')
    ggsave(filename = filenamezero, plot = histzero, width = 16, height = 9)
    
    plotzero <- ggplot(all3, aes(x=Site, y=all3[,i], fill=Dx)) + 
      scale_fill_manual(values=cbPalette)+
      geom_boxplot(outlier.colour="black") +
      scale_x_discrete(breaks = levels(all3$Site)[c(T,F,F)]) +
      xlab("Site")+
      ylab(paste0("Units for ", current))

      filenamezero2 <- paste0('~/Desktop/Site_Outlier_Detection/1.5IQR_Sites/', colnames(all3)[i], '_nooutliers.png')
      ggsave(filename = filenamezero2, plot = plotzero, width = 16, height = 9)
    
      }else{

  
  site_out_index <- NA 
  
  for(j in (1:length(site_outliers))){
    site_out_index[j] <- which(site_mean[,3] == site_outliers[[j]])
    
  }
  
  # append to the current dataframe
  site_out_df <- site_out_df %>% add_row(Site = site_mean$Site[site_out_index], Dx = site_mean$Dx[site_out_index], ROI = current)
  

  all_noout <- all3
  
  # take out the outliers to plot
  for(k in 1:dim(site_out_df)[1]){
    takeout <- which(all3$Site == site_out_df[k,1] & all3$Dx == site_out_df[k,2])
    all_noout[takeout,] <- NA
  }


  # print what it looked like before
  plot <- ggplot(all3, aes(x=Site, y=all3[,i], fill=Dx)) + 
    scale_fill_manual(values=cbPalette)+
    geom_boxplot(outlier.colour="black") +
    scale_x_discrete(breaks = levels(all3$Site)[c(T,F,F)]) +
    xlab("Site")+
    ylab(paste0("Units for ", current))
  
  # print what it looks like now
  plot2 <- ggplot(all_noout, aes(x = Site, y = all_noout[,i], fill = Dx))+
    scale_fill_manual(values = cbPalette)+
    geom_boxplot(outlier.color = "black")+
    scale_x_discrete(breaks = levels(all_noout$Site)[c(T,F,F)])+
    xlab("Site")+
    ylab(paste0("Units for ", current))
  
  hist1 <- ggplot(all3, aes(x=all3[,i])) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666")+
    xlab(paste0(current))
  
  hist2 <- ggplot(all_noout, aes(x = all_noout[,i]))+
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666")+
    xlab(paste0(current, " no outlier sites"))
  
  filename = paste0('~/Desktop/Site_Outlier_Detection/1.5IQR_Sites/', colnames(all3)[i], '_outmarked.png')
  filename2 = paste0('~/Desktop/Site_Outlier_Detection/1.5IQR_Sites/', colnames(all3)[i], '.png')
  ggsave(filename = filename, plot = plot, width = 16, height = 9)
  ggsave(filename = filename2, plot = plot2, width = 16, height = 9)
  
  filename3 <- paste0('~/Desktop/Site_Outlier_Detection/1.5IQR_Sites/', colnames(all3)[i], '_histograms.png')
  histograms <- ggarrange(hist1, hist2, nrow = 1, ncol = 2)
  ggsave(filename = filename3, plot = histograms, width = 16, height = 9)
} # finish the else statement
} # finish the for loop

# save out tibble that has all of the site exclusion information
write.csv(site_out_df, "~/Desktop/Site_Outlier_Detection/1.5_MeanOutsideIQR_OutlierList.csv")


# here we define outliers as 1.5 outside of the IQR for the means of all of the sites

##### NOTES #####
# check the math on all of these concepts....
# do we need to weight based on sample size? (FOR MEAN)
# keep this as an option for a weighted mean if we want it. May want to replace the raw mean measure above with this to account for differences in the size of the groups
#weighted_mean <- all3 %>% group_by(Site, Dx) %>% summarise(site_mean = weighted.mean(FullSurfArea, Weights, na.rm = TRUE))
# Maybe we should go back out to 2.5 IQR just to make sure that we aren't being too stringent.

##### HOW TO IMPROVE #####
# combine the two plots into one PDF if that is feasible
# standardize the x-axis of the plots to make them look a bit better
# flip 1.5_MeanOutsideIQR_OutlierList.csv so that it prints out all the sites that are getting taken out to match with the visual QC

####### VARIABILITY WITHIN SITE ######
weights_grouped <- all3 %>% count(Site, Dx)
all3$Weights <- NA

novar <- as.data.frame(matrix(NA, ncol = 1))
filename_novar <- paste0("~/Desktop/Site_Outlier_Detection/NoVarianceList.txt")
file.create(filename_novar)

tab_newline <- function(){
  write(paste0(" "), file = filename_novar, append = TRUE)
}


for(i in 1:length(weights_grouped$n)){
  index <- which((all3$Site == weights_grouped$Site[i]) & (all3$Dx == weights_grouped$Dx[i]))
  all3[index,which(colnames(all3) == "Weights")] <- weights_grouped$n[i]
}

site_out_df_sd <- data.frame(matrix(NA, ncol = 3))
colnames(site_out_df_sd) <- c("Site", "Dx", "ROI")

for(stdev in start:end){
  
  tab_newline()
  
  current <- colnames(all3[stdev])

  weighted_sd <- all3 %>% group_by(Site, Dx) %>% summarise_at(current, funs(weighted.sd(., w = Weights, na.rm = TRUE)))
  colnames(weighted_sd) <- c("Site", "Dx", "weight_site_var")
  
  # define outliers outside 1.5 IQR
  sd_outliers <- unique(boxplot(weighted_sd$weight_site_var, plot = FALSE, range = 2.5)$out)

  outlier_idx <- NA
  
  if(length(sd_outliers) > 0)
  for(m in 1:length(sd_outliers)){
    outlier_idx[m] <- which(weighted_sd$weight_site_var == sd_outliers[m])
    
  }

  
  if(any(sd_outliers == 0)){
    index <- which(weighted_sd$weight_site_var == 0)
    sentence <- paste0("For the ROI ", current, ", site ", weighted_sd$Site[index], " with dx ", weighted_sd$Dx[index], " had variance 0 and should be checked.")
    
    novar <- write(paste0(sentence), file = filename_novar, append = TRUE)
    
  } # end of if statement for variances of 0
  
  if(length(outlier_idx) == 0){
    # I don't think this needs to be written out to a file as it will just be missing from the spreadsheet
    sentence <- paste0("The ", current, " had no sites/dx that had variances that were outside 2.5 IQR.")
    # write.csv(sentence, paste0("~/Desktop/Site_Outlier_Detection/", current, "sd_out.png"))
    
  }else{
    #for(out in (1:length(sd_outliers))){
    #  outlier_idx[out] <- which(weighted_sd$weight_site_var == sd_outliers)
    #}
    # DO WE NEED TO HAVE THIS FOR LOOP HERE?
    
    for(idx in 1:length(outlier_idx)){
      site_out_df_sd <- site_out_df_sd %>% add_row(Site = weighted_sd$Site[outlier_idx], Dx = weighted_sd$Dx[outlier_idx], ROI = current)
      
    }
    
    site_out_df_sd <- site_out_df_sd %>% add_row(Site = weights_grouped$Site[outlier_idx], Dx = weights_grouped$Dx[outlier_idx], ROI = current)
    
    no_outlier_sd <- weighted_sd
    # we can either graph the distribution of variances OR of all of the values.... I am going to go with distribution of   variances for this one. Previously we were graphing the overall distribution of the values
    
    for(k in 1:dim(site_out_df_sd)[1]){
      takeout_sd <- which(weighted_sd$Site == site_out_df_sd[k,1] & weighted_sd$Dx == site_out_df_sd[k,2])
      no_outlier_sd[takeout_sd,] <- NA
    }
    
    hist_sd_correct <-  ggplot(no_outlier_sd, aes(x=weight_site_var)) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.6, fill="#C1D0D0")+
      xlab(paste0(current))
    # histogram of the variances between sites
    
    hist_sd_uncorrect <-  ggplot(weighted_sd, aes(x=weight_site_var)) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      geom_density(alpha=.6, fill="#C1D0D0")+
      xlab(paste0(current))
    # histogram of weighted SD uncorrected
    
    histplots <- ggarrange(hist_sd_uncorrect, hist_sd_correct, nrow = 1, ncol = 2)
    filename_hist <- paste0("~/Desktop/Site_Outlier_Detection/", current, "_variance.png")
    ggsave(filename = filename_hist, plot = histplots, width = 16, height = 9)
  } # end the else statement
  
  
} # end the for loop

write.csv(site_out_df_sd, "~/Desktop/Site_Outlier_Detection/SD_2.5OutlierSites.csv")

