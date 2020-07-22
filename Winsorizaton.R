# Winsorization of Outliers for ENIGMA Cross Disorder Dataset

##### SETTING UP #####
library(ggplot2)
library(dplyr)
library(DescTools)
library(psych)
library(ggpubr)
library(tidyr)
library(reshape2)

wrkdir <- "~/Desktop/XDisorder_Outlier/"

options(pillar.sigfig=20) # sets dplyr to 20 significant figures to increase accuracy
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7")

current_iqr <- 1.5

all3 <- read.csv(paste0(wrkdir, "ALL2_13July2020.csv")) # read in the latest data sheet

all3$Site_new <- as.factor(all3$Site_new) # change everything to factors for plotting later
all3$Dx <- as.factor(all3$Dx)
all3$Site <- as.factor(all3$Site)

all3_sorted <- all3 %>% arrange(Site) # sort by site to match up with the way the loops are ran below

start <- which(colnames(all3)=="LLatVent")
end <- which(colnames(all3)=="FullSurfArea")

wrkdir_plots <- "~/Desktop/WinsorizedPlots/"


###### DEMONSTRATION OF WINSORIZATION #######
# Example ROI: Full Surface Area, across all sites
# the following is to show how different thresholds affect the distribution of values differently
fsa <- all3 %>% select(FullSurfArea)

fsa_raw <- ggplot(fsa, aes(x=fsa[,1])) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 50)+
  geom_density(alpha=.5, fill="#D1DDE6")+
  xlab(paste0("Full Surface Area"))+
  lims(x = range(fsa, na.rm = TRUE))+
  ylab("Raw Value Density")+
  ggtitle("Raw Values")

ggsave(filename = paste0(wrkdir_plots, "fsa_raw.png"), plot = fsa_raw, width = 9, height = 9)

fsa_winsor_10 <- as.data.frame(winsor(fsa, trim = 0.05, na.rm = TRUE))

fsa_10 <- ggplot(fsa_winsor_10, aes(x=fsa_winsor_10[,1])) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 50)+
  geom_density(alpha=.5, fill="#838B96")+
  xlab(paste0("Full Surface Area"))+
  lims(x = range(fsa, na.rm = TRUE))+
  ylab("Raw Value Density")+
  ggtitle("Winsorized 10%")

ggsave(filename = paste0(wrkdir_plots, "fsa_10.png"), plot = fsa_10, width = 9, height = 9)

fsa_winsor_5 <- as.data.frame(winsor(fsa, trim = 0.025, na.rm = TRUE))

fsa_5 <- ggplot(fsa_winsor_5, aes(x=fsa_winsor_5[,1])) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 50)+
  geom_density(alpha=.2, fill="#79D7F7")+
  xlab(paste0("Full Surface Area"))+
  lims(x = range(fsa, na.rm = TRUE))+
  ylab("Raw Value Density")+
  ggtitle("Winsorized 5%")

ggsave(filename = paste0(wrkdir_plots, "fsa_5.png"), plot = fsa_5, width = 9, height = 9)

fsa_winsor_1 <- as.data.frame(winsor(fsa, trim = 0.005, na.rm = TRUE))

fsa_1 <- ggplot(fsa_winsor_1, aes(x=fsa_winsor_1[,1])) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 50)+
  geom_density(alpha=.7, fill="#F5F5F5")+
  xlab(paste0("Full Surface Area"))+
  lims(x = range(fsa, na.rm = TRUE))+
  ylab("Raw Value Density")+
  ggtitle("Winsorized 1%")

ggsave(filename = paste0(wrkdir_plots, "fsa_1.png"), plot = fsa_1, width = 9, height = 9)

###### PERCENT OUTLIERS: WITHIN SITE ######
# Decide what percentage of data points are outside set IQR. 
# Used as a data driven method for winsorization

perc_out_df <- as.data.frame(matrix(NA, ncol = 2)) # data frame that holds the percent of points that are outliers within site for each ROI
colnames(perc_out_df) <- c("Percent", "ROI")


for(i in start:end){
  roi <- colnames(all3)[i]
  
  for(j in levels(all3$Site)){
    
    site_specific <- all3 %>% filter(Site == as.numeric(j)) %>% select(roi)
    outliers <- boxplot(site_specific[,1], plot = FALSE, range = current_iqr)$out
    perc_outlier <- (length(outliers) / length(site_specific[,1]))*100 
    
    # now append onto the dataframe w/ a label for ROI
    perc_out_df <- perc_out_df %>% add_row(Percent = perc_outlier, ROI = roi)
    
    
  } # end site loop
} # end ROI group

perc_out_df

range(perc_out_df$Percent, na.rm = TRUE)
mean(perc_out_df$Percent, na.rm = TRUE)
median(perc_out_df$Percent, na.rm = TRUE)

perc_out_df$Index <- 1:length(perc_out_df$Percent)

percent_within_index <- ggplot(perc_out_df, aes(x = Index, y = Percent, label = ROI, col = as.factor(ROI)))+
  geom_point()+
  geom_line()+
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(data = subset(perc_out_df, Percent > 25))+
  geom_hline(yintercept = mean(perc_out_df$Percent, na.rm = TRUE))+
  geom_hline(yintercept = median(perc_out_df$Percent, na.rm = TRUE), linetype = "dashed")+
  scale_color_grey()

ggsave(filename = paste0(wrkdir_plots, "percent_within_index.png"), plot = percent_across_index, height = 9, width = 16)

####### WINSORIZE WITHIN SITE #########
winsorized_values <- as.data.frame(matrix(NA, ncol = 3))
colnames(winsorized_values) <- c("Winsorized", "Site", "ROI")

for(i in start:end){ # start the ROI loop
  roi <- colnames(all3)[i]
  
  for(j in levels(all3$Site)){ # start the site loop 
    
    site <- as.numeric(j)
    
    site_current <- all3_sorted %>% filter(Site == site) %>% select(roi)
    winsorized <- cbind.data.frame(winsor(site_current[,1], trim = 0.01, na.rm = TRUE), site, roi)
    colnames(winsorized) <- c("Winsorized", "Site", "ROI")
    
    winsorized_values <- winsorized_values %>% add_row(winsorized)
    
    # save the winsorized values to a new dataframe (add columns with a new ROI and add rows with a new site each time)
    
  }
}

winsorized_values <- winsorized_values %>% group_by(ROI) %>% mutate(grouped_id = row_number())

all_win_within <- as.data.frame(winsorized_values %>% spread(ROI, Winsorized) %>% select(-grouped_id))

all_win_within <- all_win_within[-dim(all_win_within)[1],]
all_win_within$Dx <- as.factor(all3_sorted$Dx)

####### EXAMPLE BOX PLOT #########
# For example: Full Surface Area
fsa_before <- all3 %>% select(Site, Dx, FullSurfArea)

before_fsa <- ggplot(fsa_before, aes(x = Site, y = FullSurfArea, fill = Dx))+
  scale_fill_manual(values = cbPalette)+
  geom_boxplot(outlier.colour = "black")+
  scale_x_discrete(breaks = levels(fsa_before$Site)[c(T,F,F)])+
  xlab("Site")+
  ylim(c(50000, 250000))+
  ylab("Units for Full Surface Area")

ggsave(filename = paste0(wrkdir_plots, "before_fsa.png"), plot = before_fsa, height = 9, width = 16)

fsa_after <- all_win_within %>% select(Site, Dx, FullSurfArea)

after_fsa <- ggplot(fsa_before, aes(x = Site, y = FullSurfArea, fill = Dx))+
  scale_fill_manual(values = cbPalette)+
  geom_boxplot(outlier.colour = "black")+
  scale_x_discrete(breaks = levels(fsa_before$Site)[c(T,F,F)])+
  xlab("Site")+
  ylim(c(50000, 250000))+
  ylab("Units for Full Surface Area")

ggsave(filename = paste0(wrkdir_plots, "after_fsa.png"), plot = after_fsa, height = 9, width = 16)

####### DEMONSTRATIONS FOR WINSORIZATION WITHIN #######
fsa_site_113 <- all3 %>% filter(Site == 113) %>% select(FullSurfArea)
x_range <- range(fsa_site_113$FullSurfArea, na.rm = TRUE)

threeraw <- ggplot(fsa_site_113, aes(x=fsa_site_113[,1])) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 50)+
  geom_density(alpha=.4, fill="#A8B6C8")+
  xlab(paste0("Full Surface Area"))+
  ylab("Raw Value Density")+
  xlim(x_range)+
  ggtitle("Raw Values, Site 113")

ggsave(filename = paste0(wrkdir_plots, "threeraw.png"), plot = threeraw, height = 9, width = 9)

fsa_wins_113 <- as.data.frame(winsor(fsa_site_113, trim = ((mean(perc_out_df$Percent, na.rm = TRUE)/2)/100), na.rm = TRUE))

threewin <- ggplot(fsa_wins_113, aes(x=fsa_wins_113[,1])) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 50)+
  geom_density(alpha=.6, fill="#2A5078")+
  xlab(paste0("Full Surface Area"))+
  ylab("Wisorized Density")+
  xlim(x_range)+
  ggtitle("Winsorized at 1.3%, Site 113")

combined_three <- ggarrange(threeraw, threewin, nrow = 1, ncol = 2)

ggsave(filename = paste0(wrkdir_plots, "threecombo.png"), plot = combined_three, height = 12, width = 12)

######### PERCENT OUTLIERS: ACROSS SITES #########
perc_out_df_across <- as.data.frame(matrix(NA, ncol = 2)) # data frame that holds the percent of points that are outliers within site for each ROI
colnames(perc_out_df_across) <- c("Percent", "ROI")

# need to test if this works
for(i in start:end){
  
  roi <- colnames(all3)[i]
  roi_specific <- all3 %>% select(roi)
  outliers <- boxplot(roi_specific[,1], plot = FALSE, range = current_iqr)$out
  perc_outlier <- (length(outliers) / length(roi_specific[,1]))*100 
  
  # now append onto the dataframe w/ a label for ROI
  perc_out_df_across <- perc_out_df_across %>% add_row(Percent = perc_outlier, ROI = roi)
  
  
} # end ROI group

head(perc_out_df_across)

range(perc_out_df_across$Percent, na.rm = TRUE)
mean(perc_out_df_across$Percent, na.rm = TRUE)
median(perc_out_df_across$Percent, na.rm = TRUE)

perc_out_df_across$Index <- 1:length(perc_out_df_across$Percent)

perc_out_df_across_plot <- ggplot(perc_out_df_across, aes(x = Index, y = Percent, label = ROI))+
  geom_point()+
  geom_line()+
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(data = subset(perc_out_df_across, Percent > 3))+
  geom_hline(yintercept = mean(perc_out_df_across$Percent, na.rm = TRUE))+
  geom_hline(yintercept = median(perc_out_df_across$Percent, na.rm = TRUE), linetype = "dashed")

ggsave(filename = paste0(wrkdir_plots, "perc_out_df_across_plot.png"), plot = perc_out_df_across_plot, height = 12, width = 16)

####### WINSORIZE ACROSS SITE #########
winsorized_values_x <- as.data.frame(matrix(NA, ncol = 2))
colnames(winsorized_values_x) <- c("Winsorized", "ROI")

for(i in start:end){ # start the ROI loop
  roi <- colnames(all3)[i]
  
  site_current <- all3_sorted %>% select(roi)
  winsorized <- cbind.data.frame(winsor(site_current[,1], trim = 0.01, na.rm = TRUE), roi)
  colnames(winsorized) <- c("Winsorized", "ROI")
  
  winsorized_values_x <- winsorized_values_x %>% add_row(winsorized)
  
}

winsorized_values_x <- winsorized_values_x %>% group_by(ROI) %>% mutate(grouped_id = row_number())

all_win_across <- as.data.frame(winsorized_values_x %>% spread(ROI, Winsorized) %>% select(-grouped_id))

all_win_across$Site <- as.factor(all3_sorted$Site)
all_win_across$Dx <- as.factor(all3_sorted$Dx)

####### EXAMPLE BOX PLOT: ACROSS #########
# For example: Full Surface Area
fsa_before # rendered previously

fsa_after <- all_win_across %>% select(Site, Dx, FullSurfArea)

after_fsa <- ggplot(fsa_after, aes(x = Site, y = FullSurfArea, fill = Dx))+
  scale_fill_manual(values = cbPalette)+
  geom_boxplot(outlier.colour = "black")+
  scale_x_discrete(breaks = levels(fsa_after$Site)[c(T,F,F)])+
  xlab("Site")+
  ylim(c(50000, 250000))+
  ylab("Units for Full Surface Area")

ggsave(filename = paste0(wrkdir_plots, "after_fsa2.png"), plot = after_fsa, height = 9, width = 16)

######## DEMONSTRATIONS FOR WINSORIZE ACROSS SITE ########
fsa_raw # raw plot: created above

fsa_winsor_across <- as.data.frame(winsor(fsa, trim = ((mean(perc_out_df_across$Percent, na.rm = TRUE)/2)/100), na.rm = TRUE))
  
fsa_1 <- ggplot(fsa_winsor_across, aes(x=fsa_winsor_across[,1])) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 50)+
    geom_density(alpha=.7, fill="#2A5078")+
    xlab(paste0("Full Surface Area"))+
    lims(x = range(fsa, na.rm = TRUE))+
    ylab("Winsorized Value Density")+
    ggtitle("Winsorized Values")

across <- ggarrange(fsa_raw, fsa_1, nrow = 1, ncol = 2)

ggsave(filename = paste0(wrkdir_plots, "across.png"), plot = across, width = 12, height = 9)

# pick a site with REALLY low outliers to show what would happen
# Lpal
lpal <- all3 %>% select(Lpal)
x_range <- range(lpal, na.rm = TRUE)

lpal_raw <- ggplot(lpal, aes(x=lpal[,1])) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 50)+
  geom_density(alpha=.7, fill="#D1DDE6")+
  xlab("Left Pallidum")+
  ylab("Raw Value Denisty")+
  xlim(x_range)+
  ggtitle("Raw Values")

lpal_win <- as.data.frame(winsor(lpal, trim = ((mean(perc_out_df_across$Percent, na.rm = TRUE)/2)/100), na.rm = TRUE))

lpal_adjust <- ggplot(lpal_win, aes(x=lpal_win[,1])) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 50)+
  geom_density(alpha=.7, fill="#2A5078")+
  xlab("Left Pallidum")+
  ylab("Winsor Value Denisty")+
  xlim(x_range)+
  ggtitle("Winsorized Values")

lpal_combo <- ggarrange(lpal_raw, lpal_adjust, nrow = 1, ncol = 2)

ggsave(filename = paste0(wrkdir_plots, "lpal_combo.png"), plot = lpal_combo, height = 9, width = 12)


######## FLAG OUTLIER MEAN SITES ########
# all_win_within
# all_win_across

site_out_df_raw <- as.data.frame(matrix(NA, ncol = 2))
colnames(site_out_df_raw) <- c("Site", "ROI")

# original flagging
for(i in start:end){
  
current <- colnames(all3_sorted[i])
  
site_mean <- as.data.frame(all3_sorted %>% select(2,i,238,239) %>% na.omit() %>% group_by(Site) %>% summarise_at(current, mean))
  
site_outliers <- boxplot(site_mean[,current], plot = FALSE, range = current_iqr)$out

if(length(site_outliers) == 0){
  
}else{
  
  site_out_index <- NA
  
  for(j in 1:length(site_outliers)){ # find index site mean data frame has outliers
    site_out_index[j] <- which(site_mean[,current] == site_outliers[[j]])
  }
  
}# end if/else loop for no site outliers

site_out_df_raw <- site_out_df_raw %>% add_row(Site = site_out_index, ROI = current)

}# end loop for ROI

# within flagging

site_out_df_within <- as.data.frame(matrix(NA, ncol = 2))
colnames(site_out_df_within) <- c("Site", "ROI")

start_within <- which(colnames(all_win_within)=="FullSurfArea")
end_within <- which(colnames(all_win_within)=="RThickness")

for(i in start_within:end_within){
  
  current <- colnames(all_win_within[i])
  
  site_mean <- as.data.frame(all_win_within %>% select("Dx",current,"Site") %>% na.omit() %>% group_by(Site) %>% summarise_at(current, mean))

  site_outliers <- boxplot(site_mean[,current], plot = FALSE, range = current_iqr)$out
  
  if(length(site_outliers) == 0){
    
  }else{
    
    site_out_index <- NA
    
    for(j in 1:length(site_outliers)){ # find index site mean data frame has outliers
      site_out_index[j] <- which(site_mean[,current] == site_outliers[[j]])
    }
    
  }# end if/else loop for no site outliers
  
  site_out_df_within <- site_out_df_within %>% add_row(Site = site_out_index, ROI = current)
  
}# end loop for ROI

# across flagging

site_out_df_across <- as.data.frame(matrix(NA, ncol = 2))
colnames(site_out_df_across) <- c("Site", "ROI")

start_across <- which(colnames(all_win_across)=="FullSurfArea")
end_across <- which(colnames(all_win_within)=="RThickness")

for(i in start_across:end_across){
  
  current <- colnames(all_win_across[i])
  
  site_mean <-  site_mean <- as.data.frame(all_win_across %>% select("Dx",current,"Site") %>% na.omit() %>% group_by(Site) %>% summarise_at(current, mean))
  
  site_outliers <- boxplot(site_mean[,current], plot = FALSE, range = current_iqr)$out
  
  if(length(site_outliers) == 0){
    
  }else{
    
    site_out_index <- NA
    
    for(j in 1:length(site_outliers)){ # find index site mean data frame has outliers
      site_out_index[j] <- which(site_mean[,current] == site_outliers[[j]])
    }
    
    site_out_df_across <- site_out_df_across %>% add_row(Site = site_out_index, ROI = current)
    
    
  }# end if/else loop for no site outliers

}# end loop for ROI

######## FLAG OUTLIER SD SITES ########
# raw flagging
site_out_df_sd <- data.frame(matrix(NA, ncol = 2))
colnames(site_out_df_sd) <- c("Site", "ROI")

for(stdev in start:end){
  
current <- colnames(all3_sorted[stdev])

sd_dist <- as.data.frame(all3_sorted %>% group_by(Site) %>% summarise_at(current, sd, na.rm = TRUE))
colnames(sd_dist) <- c("Site", current)

sd_outliers <- boxplot(sd_dist[,current], plot = FALSE, range = current_iqr)$out

if(length(sd_outliers) == 0){
  
} else {
  
  outlier_idx <- NA
  
  for(m in 1:length(sd_outliers)){
    outlier_idx[m] <- which(sd_dist[,current] == sd_outliers[m])
  }# end for loop for finding outlier idx

site_out_df_sd <- site_out_df_sd %>% add_row(Site = outlier_idx, ROI = current)
  
  }# end if/else statement
} # end ROI loop

# flag after within winsor
site_out_df_sd_within <- data.frame(matrix(NA, ncol = 2))
colnames(site_out_df_sd_within) <- c("Site", "ROI")

for(stdev in start_within:end_within){
  
  current <- colnames(all_win_within[stdev])
  
  sd_dist <- as.data.frame(all_win_within %>% group_by(Site) %>% summarise_at(current, sd, na.rm = TRUE))
  colnames(sd_dist) <- c("Site", current)
  
  sd_outliers <- boxplot(sd_dist[,current], plot = FALSE, range = current_iqr)$out
  
  if(length(sd_outliers) == 0){
    
  } else {
    
    outlier_idx <- NA
    
    for(m in 1:length(sd_outliers)){
      outlier_idx[m] <- which(sd_dist[,current] == sd_outliers[m])
    }# end for loop for finding outlier idx
    
    site_out_df_sd_within <- site_out_df_sd_within %>% add_row(Site = outlier_idx, ROI = current)
    
  }# end if/else statement
} # end ROI loop



# flag after across winsor
site_out_df_sd_across <- data.frame(matrix(NA, ncol = 2))
colnames(site_out_df_sd_across) <- c("Site", "ROI")

for(stdev in start_within:end_within){
  
  current <- colnames(all_win_within[stdev])
  
  sd_dist <- as.data.frame(all_win_within %>% group_by(Site) %>% summarise_at(current, sd, na.rm = TRUE))
  colnames(sd_dist) <- c("Site", current)
  
  sd_outliers <- boxplot(sd_dist[,current], plot = FALSE, range = current_iqr)$out
  
  if(length(sd_outliers) == 0){
    
  } else {
    
    outlier_idx <- NA
    
    for(m in 1:length(sd_outliers)){
      outlier_idx[m] <- which(sd_dist[,current] == sd_outliers[m])
    }# end for loop for finding outlier idx
    
    site_out_df_sd_across <- site_out_df_sd_across %>% add_row(Site = outlier_idx, ROI = current)
    
  }# end if/else statement
} # end ROI loop

###### MERGE AND PRINT OUTLIER DATA FRAME  ######

one <- site_out_df_raw %>% group_by(ROI) %>% summarise(OutlierSites_MeanRaw = toString(Site)) %>% ungroup()
two <- site_out_df_within %>% group_by(ROI) %>% summarise(OutlierSites_MeanWithin = toString(Site)) %>% ungroup()
three <- site_out_df_across %>% group_by(ROI) %>% summarise(OutlierSites_MeanAcross = toString(Site)) %>% ungroup()
four <- site_out_df_sd %>% group_by(ROI) %>% summarise(OutlierSites_SDRaw = toString(Site)) %>% ungroup()
five <- site_out_df_sd_within %>% group_by(ROI) %>% summarise(OutlierSites_SDWithin = toString(Site)) %>% ungroup()
six <- site_out_df_sd_across %>% group_by(ROI) %>% summarise(OutlierSites_SDAcross = toString(Site)) %>% ungroup()

outlier_sites <- one %>% merge(two, by = "ROI", all = TRUE) %>% merge(three, by = "ROI", all = TRUE) %>% merge(four, by = "ROI", all = TRUE) %>% merge(five, by = "ROI", all = TRUE) %>% merge(six, by = "ROI", all = TRUE)

write.csv(outlier_sites, paste0(wrkdir, "Table_OutlierSites.csv"), row.names = F)

####### NOTES ######










