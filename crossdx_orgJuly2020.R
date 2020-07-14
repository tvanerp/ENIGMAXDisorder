setwd("/Users/Elena/Documents/Work/ENIGMA/crossdx")

library(plyr)
library(psych)
library(reshape)

##MDD###

MDD_freesurfer <- read.csv("Subcort_Cort_11June2020.csv")
MDD_cov <- read.csv("ENIGMA_MDD_Covariates_08July2020.csv")

MDD_cov2 <- MDD_cov[,c(1:4, 7, 9)] #select only relevant variables to crossdx

#!!!check site column if new version of db!!!

MDD_freesurfer2 <- MDD_freesurfer[,c(1:236, 301)] #select only freesurfer variables + site

MDD = merge(MDD_freesurfer2, MDD_cov2, by="SubjID") #merge freesurfer and cov

#keep only data from sites that joined the crossdx project

MDD_xd <- MDD[grep("Episca|Munster|LOND|Gron|Stanf|SF|Edi|SaoPaulo|Melb|Novo|FOR2107|Barc|AMC|MPIP|CLING|CODE|SHIP|Houst|SYD|Magd|Dub|QTIM|McMas|Oxford|Singapore|AFFDIS|NESDA|ETPB|Minnesota|Calgary|Moraldilemma|BiDirect|Socat|FIDMAG|CARDIFF|SanRaffaele|Oxfordyoung|TAD|TIGER|Hiroshima", MDD$SubjID),]

#recode AD variable

MDD_xd$AD[MDD_xd$AD=="2"] <- "1"

#rename IDs to match BP
MDD_xd$SubjID <- gsub("SaoPaulo_", "", MDD_xd$SubjID)
MDD_xd$SubjID <- gsub("Houstyoung_", "Houst_", MDD_xd$SubjID)

#Add working group
MDD_xd['WG']='1'

###BD###

BD <- read.csv("BD_MasterData_04_23_2020.csv")

BD2 <- rename(BD, c("AgeofOnset" = "AO")) #rename variable for consistency
BD2 <- rename(BD, c("AntiDep" = "AD")) #rename variable for consistency
BD2 <- BD2[,-c(173:180, 182:250)] #remove means - will recalculate them later
BD2$Dx[BD2$Dx=="1"] <- "2" #recode Dx

#add constant to Site

BD2$Site <- BD2$Site + 100

#calculate means. Not super efficient way but it does not rely on particular columns numbers

BD2$Mvent <- rowMeans(BD2[,c("LLatVent","RLatVent")]);  #calculate mean Ventricle
BD2$Mthal <- rowMeans(BD2[,c("Lthal","Rthal")]);  #calculate mean Thalamus
BD2$Mcaud <- rowMeans(BD2[,c("Lcaud","Rcaud")]);  #calculate mean Caudate
BD2$Mput <- rowMeans(BD2[,c("Lput","Rput")]);   #calculate mean Putamen
BD2$Mpal <- rowMeans(BD2[,c("Lpal","Rpal")]);   #calculate mean Pallidum
BD2$Mhippo <- rowMeans(BD2[,c("Lhippo","Rhippo")]);   #calculate mean Hippocampus
BD2$Mamyg <- rowMeans(BD2[,c("Lamyg","Ramyg")]);   #calculate mean Amygdala
BD2$Maccumb <- rowMeans(BD2[,c("Laccumb","Raccumb")]);   #calculate mean Accumbens

BD2$M_caudalanteriorcingulate_thickavg <- rowMeans(BD2[,c("L_caudalanteriorcingulate_thickavg" , "R_caudalanteriorcingulate_thickavg")]);
BD2$M_caudalmiddlefrontal_thickavg <- rowMeans(BD2[,c("L_caudalmiddlefrontal_thickavg" , "R_caudalmiddlefrontal_thickavg")]);
BD2$M_cuneus_thickavg <- rowMeans(BD2[,c("L_cuneus_thickavg" , "R_cuneus_thickavg")]);
BD2$M_entorhinal_thickavg <- rowMeans(BD2[,c("L_entorhinal_thickavg" , "R_entorhinal_thickavg")]);
BD2$M_fusiform_thickavg <- rowMeans(BD2[,c("L_fusiform_thickavg" , "R_fusiform_thickavg")]);
BD2$M_inferiorparietal_thickavg <- rowMeans(BD2[,c("L_inferiorparietal_thickavg" , "R_inferiorparietal_thickavg")]);
BD2$M_inferiortemporal_thickavg <- rowMeans(BD2[,c("L_inferiortemporal_thickavg" , "R_inferiortemporal_thickavg")]);
BD2$M_isthmuscingulate_thickavg <- rowMeans(BD2[,c("L_isthmuscingulate_thickavg" , "R_isthmuscingulate_thickavg")]);
BD2$M_lateraloccipital_thickavg <- rowMeans(BD2[,c("L_lateraloccipital_thickavg" , "R_lateraloccipital_thickavg")]);
BD2$M_lateralorbitofrontal_thickavg <- rowMeans(BD2[,c("L_lateralorbitofrontal_thickavg" , "R_lateralorbitofrontal_thickavg")]);
BD2$M_lingual_thickavg <- rowMeans(BD2[,c("L_lingual_thickavg" , "R_lingual_thickavg")]);
BD2$M_medialorbitofrontal_thickavg <- rowMeans(BD2[,c("L_medialorbitofrontal_thickavg" , "R_medialorbitofrontal_thickavg")]);
BD2$M_middletemporal_thickavg <- rowMeans(BD2[,c("L_middletemporal_thickavg" , "R_middletemporal_thickavg")]);
BD2$M_parahippocampal_thickavg <- rowMeans(BD2[,c("L_parahippocampal_thickavg" , "R_parahippocampal_thickavg")]);
BD2$M_paracentral_thickavg 	<- rowMeans(BD2[,c("L_paracentral_thickavg" , "R_paracentral_thickavg")]);
BD2$M_parsopercularis_thickavg <- rowMeans(BD2[,c("L_parsopercularis_thickavg" , "R_parsopercularis_thickavg")]);
BD2$M_parsorbitalis_thickavg <- rowMeans(BD2[,c("L_parsorbitalis_thickavg" , "R_parsorbitalis_thickavg")]);
BD2$M_parstriangularis_thickavg <- rowMeans(BD2[,c("L_parstriangularis_thickavg" , "R_parstriangularis_thickavg")]);
BD2$M_pericalcarine_thickavg <- rowMeans(BD2[,c("L_pericalcarine_thickavg" , "R_pericalcarine_thickavg")]);
BD2$M_postcentral_thickavg <- rowMeans(BD2[,c("L_postcentral_thickavg" , "R_postcentral_thickavg")]);
BD2$M_posteriorcingulate_thickavg <- rowMeans(BD2[,c("L_posteriorcingulate_thickavg" , "R_posteriorcingulate_thickavg")]);
BD2$M_precentral_thickavg <- rowMeans(BD2[,c("L_precentral_thickavg" , "R_precentral_thickavg")]);
BD2$M_precuneus_thickavg <- rowMeans(BD2[,c("L_precuneus_thickavg" , "R_precuneus_thickavg")]);
BD2$M_rostralanteriorcingulate_thickavg <- rowMeans(BD2[,c("L_rostralanteriorcingulate_thickavg" , "R_rostralanteriorcingulate_thickavg")]);
BD2$M_rostralmiddlefrontal_thickavg <- rowMeans(BD2[,c("L_rostralmiddlefrontal_thickavg" , "R_rostralmiddlefrontal_thickavg")]);
BD2$M_superiorfrontal_thickavg <- rowMeans(BD2[,c("L_superiorfrontal_thickavg" , "R_superiorfrontal_thickavg")]);
BD2$M_superiorparietal_thickavg <- rowMeans(BD2[,c("L_superiorparietal_thickavg" , "R_superiorparietal_thickavg")]);
BD2$M_superiortemporal_thickavg <- rowMeans(BD2[,c("L_superiortemporal_thickavg" , "R_superiortemporal_thickavg")]);
BD2$M_supramarginal_thickavg <- rowMeans(BD2[,c("L_supramarginal_thickavg" , "R_supramarginal_thickavg")]);
BD2$M_frontalpole_thickavg <- rowMeans(BD2[,c("L_frontalpole_thickavg" , "R_frontalpole_thickavg")]);
BD2$M_temporalpole_thickavg <- rowMeans(BD2[,c("L_temporalpole_thickavg" , "R_temporalpole_thickavg")]);
BD2$M_transversetemporal_thickavg <- rowMeans(BD2[,c("L_transversetemporal_thickavg" , "R_transversetemporal_thickavg")]);
BD2$M_insula_thickavg <- rowMeans(BD2[,c("L_insula_thickavg" , "R_insula_thickavg")]);
BD2$M_transversetemporal_thickavg <- rowMeans(BD2[,c("L_transversetemporal_thickavg" , "R_transversetemporal_thickavg")]);
BD2$M_bankssts_thickavg <- rowMeans(BD2[,c("L_bankssts_thickavg" , "R_bankssts_thickavg")]);
BD2$MThickness <- rowMeans(BD2[,c("LThickness" , "RThickness")]);

BD2$M_caudalanteriorcingulate_surfavg <- rowMeans(BD2[,c("L_caudalanteriorcingulate_surfavg" , "R_caudalanteriorcingulate_surfavg")]);
BD2$M_caudalmiddlefrontal_surfavg <- rowMeans(BD2[,c("L_caudalmiddlefrontal_surfavg" , "R_caudalmiddlefrontal_surfavg")]);
BD2$M_cuneus_surfavg <- rowMeans(BD2[,c("L_cuneus_surfavg" , "R_cuneus_surfavg")]);
BD2$M_entorhinal_surfavg <- rowMeans(BD2[,c("L_entorhinal_surfavg" , "R_entorhinal_surfavg")]);
BD2$M_fusiform_surfavg <- rowMeans(BD2[,c("L_fusiform_surfavg" , "R_fusiform_surfavg")]);
BD2$M_inferiorparietal_surfavg <- rowMeans(BD2[,c("L_inferiorparietal_surfavg" , "R_inferiorparietal_surfavg")]);
BD2$M_inferiortemporal_surfavg <- rowMeans(BD2[,c("L_inferiortemporal_surfavg" , "R_inferiortemporal_surfavg")]);
BD2$M_isthmuscingulate_surfavg <- rowMeans(BD2[,c("L_isthmuscingulate_surfavg" , "R_isthmuscingulate_surfavg")]);
BD2$M_lateraloccipital_surfavg <- rowMeans(BD2[,c("L_lateraloccipital_surfavg" , "R_lateraloccipital_surfavg")]);
BD2$M_lateralorbitofrontal_surfavg <- rowMeans(BD2[,c("L_lateralorbitofrontal_surfavg" , "R_lateralorbitofrontal_surfavg")]);
BD2$M_lingual_surfavg <- rowMeans(BD2[,c("L_lingual_surfavg" , "R_lingual_surfavg")]);
BD2$M_medialorbitofrontal_surfavg <- rowMeans(BD2[,c("L_medialorbitofrontal_surfavg" , "R_medialorbitofrontal_surfavg")]);
BD2$M_middletemporal_surfavg <- rowMeans(BD2[,c("L_middletemporal_surfavg" , "R_middletemporal_surfavg")]);
BD2$M_parahippocampal_surfavg <- rowMeans(BD2[,c("L_parahippocampal_surfavg" , "R_parahippocampal_surfavg")]);
BD2$M_paracentral_surfavg 	<- rowMeans(BD2[,c("L_paracentral_surfavg" , "R_paracentral_surfavg")]);
BD2$M_parsopercularis_surfavg <- rowMeans(BD2[,c("L_parsopercularis_surfavg" , "R_parsopercularis_surfavg")]);
BD2$M_parsorbitalis_surfavg <- rowMeans(BD2[,c("L_parsorbitalis_surfavg" , "R_parsorbitalis_surfavg")]);
BD2$M_parstriangularis_surfavg <- rowMeans(BD2[,c("L_parstriangularis_surfavg" , "R_parstriangularis_surfavg")]);
BD2$M_pericalcarine_surfavg <- rowMeans(BD2[,c("L_pericalcarine_surfavg" , "R_pericalcarine_surfavg")]);
BD2$M_postcentral_surfavg <- rowMeans(BD2[,c("L_postcentral_surfavg" , "R_postcentral_surfavg")]);
BD2$M_posteriorcingulate_surfavg <- rowMeans(BD2[,c("L_posteriorcingulate_surfavg" , "R_posteriorcingulate_surfavg")]);
BD2$M_precentral_surfavg <- rowMeans(BD2[,c("L_precentral_surfavg" , "R_precentral_surfavg")]);
BD2$M_precuneus_surfavg <- rowMeans(BD2[,c("L_precuneus_surfavg" , "R_precuneus_surfavg")]);
BD2$M_rostralanteriorcingulate_surfavg <- rowMeans(BD2[,c("L_rostralanteriorcingulate_surfavg" , "R_rostralanteriorcingulate_surfavg")]);
BD2$M_rostralmiddlefrontal_surfavg <- rowMeans(BD2[,c("L_rostralmiddlefrontal_surfavg" , "R_rostralmiddlefrontal_surfavg")]);
BD2$M_superiorfrontal_surfavg <- rowMeans(BD2[,c("L_superiorfrontal_surfavg" , "R_superiorfrontal_surfavg")]);
BD2$M_superiorparietal_surfavg <- rowMeans(BD2[,c("L_superiorparietal_surfavg" , "R_superiorparietal_surfavg")]);
BD2$M_superiortemporal_surfavg <- rowMeans(BD2[,c("L_superiortemporal_surfavg" , "R_superiortemporal_surfavg")]);
BD2$M_supramarginal_surfavg <- rowMeans(BD2[,c("L_supramarginal_surfavg" , "R_supramarginal_surfavg")]);
BD2$M_frontalpole_surfavg <- rowMeans(BD2[,c("L_frontalpole_surfavg" , "R_frontalpole_surfavg")]);
BD2$M_temporalpole_surfavg <- rowMeans(BD2[,c("L_temporalpole_surfavg" , "R_temporalpole_surfavg")]);
BD2$M_transversetemporal_surfavg <- rowMeans(BD2[,c("L_transversetemporal_surfavg" , "R_transversetemporal_surfavg")]);
BD2$M_insula_surfavg <- rowMeans(BD2[,c("L_insula_surfavg" , "R_insula_surfavg")]);
BD2$M_transversetemporal_surfavg <- rowMeans(BD2[,c("L_transversetemporal_surfavg" , "R_transversetemporal_surfavg")]);
BD2$M_bankssts_surfavg <- rowMeans(BD2[,c("L_bankssts_surfavg" , "R_bankssts_surfavg")]);

#rename inconsistent cols

BD2 <- rename(BD2, c("L_entorhinal_thickavg" = "L_entorhil_thickavg", "R_entorhinal_thickavg" = "R_entorhil_thickavg", 
                     "L_entorhinal_surfavg" = "L_entorhil_surfavg", "R_entorhinal_surfavg" = "R_entorhil_surfavg",
                     "L_supramarginal_thickavg" = "L_supramargil_thickavg", "R_supramarginal_thickavg" = "R_supramargil_thickavg",
                     "L_supramarginal_surfavg" = "L_supramargil_surfavg", "R_supramarginal_surfavg" = "R_supramargil_surfavg",
                     "M_entorhinal_thickavg" = "M_entorhil_thickavg", "M_entorhinal_surfavg" = "M_entorhil_surfavg",
                     "M_supramarginal_thickavg" = "M_supramargil_thickavg", "M_supramarginal_surfavg" = "M_supramargil_surfavg"))

#Add working group
BD2['WG']='2'

#rename IDs to make it easier for later (add prefix)
BD2$SubjID <- gsub("^(\\d+-1_for2107)$", "FOR2107_\\1", BD2$SubjID)
BD2$SubjID <- gsub("_for2107", "", BD2$SubjID)
BD2$SubjID <- gsub("^(HHR\\d+.*)", "Halifax_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(2HHR\\d+.*)", "Halifax_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(DM\\d+.*)", "Halifax_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(70\\d+.*)", "Halifax_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(S\\d+\\w+)$", "Yale_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(0\\d+_\\d+_\\w+)$", "Barcelona_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(1\\d+_\\d+_\\w+)$", "Barcelona_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(2\\d+_\\d+_\\w+)$", "Barcelona_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(con_.*)", "moodinflame_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(pat_.*)", "moodinflame_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(1\\d+\\.*)", "Penn_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(8\\d+\\.*)", "Penn_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(9\\d+\\.*)", "Penn_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(C00\\d+.*)", "Brazil_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(P00\\d+.*)", "Brazil_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(T000\\d+.*)", "Brazil_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(A\\d+_.*)", "Munster_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(A\\d+-.*)", "Munster_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(C\\d+_.*)", "Munster_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(C\\d+-.*)", "Munster_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(N\\d+_.*)", "Munster_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(N\\d+.*)", "Munster_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(s\\d+.*)", "struct_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(2\\d\\d\\d\\d)$", "Houst_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(3\\d\\d\\d\\d)$", "Houst_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(ENBD\\d+.*)", "Houst_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(ENBDO\\d+.*)", "Houst_\\1", BD2$SubjID)
BD2$SubjID <- gsub("^(\\w\\w\\d\\d\\d\\d\\d\\d)$", "Paris_\\1", BD2$SubjID)


###SCZ###

SCZ <- read.csv("ENIGMA_SZ_16_ES_CrossDisorder_20200709.csv")

SCZ$Dx[SCZ$Dx=="1"] <- "3" #recode Dx

SCZ$Dx[SCZ$Dx=="2"] <- "3" #recode Dx
SCZ$Dx[SCZ$Dx=="5"] <- "3" #recode Dx
SCZ$Dx[SCZ$Dx=="6"] <- "3" #recode Dx

#remove subj with 0

#SCZ <- SCZ[!SCZ$SubjID == "UNIMAAS_SMF055", ]
#SCZ <- SCZ[!SCZ$SubjID == "UNIMAAS_SMF075", ]

#replace 0 with NA in freesurfer measures

SCZ[, 18:174][SCZ[, 18:174] == 0] <- NA

#calculate means

SCZ$Mvent <- rowMeans(SCZ[,c("LLatVent","RLatVent")]);  #calculate mean Ventricle
SCZ$Mthal <- rowMeans(SCZ[,c("Lthal","Rthal")]);  #calculate mean Thalamus
SCZ$Mcaud <- rowMeans(SCZ[,c("Lcaud","Rcaud")]);  #calculate mean Caudate
SCZ$Mput <- rowMeans(SCZ[,c("Lput","Rput")]);   #calculate mean Putamen
SCZ$Mpal <- rowMeans(SCZ[,c("Lpal","Rpal")]);   #calculate mean Pallidum
SCZ$Mhippo <- rowMeans(SCZ[,c("Lhippo","Rhippo")]);   #calculate mean Hippocampus
SCZ$Mamyg <- rowMeans(SCZ[,c("Lamyg","Ramyg")]);   #calculate mean Amygdala
SCZ$Maccumb <- rowMeans(SCZ[,c("Laccumb","Raccumb")]);   #calculate mean Accumbens

SCZ$M_caudalanteriorcingulate_thickavg <- rowMeans(SCZ[,c("L_caudalanteriorcingulate_thickavg" , "R_caudalanteriorcingulate_thickavg")]);
SCZ$M_caudalmiddlefrontal_thickavg <- rowMeans(SCZ[,c("L_caudalmiddlefrontal_thickavg" , "R_caudalmiddlefrontal_thickavg")]);
SCZ$M_cuneus_thickavg <- rowMeans(SCZ[,c("L_cuneus_thickavg" , "R_cuneus_thickavg")]);
SCZ$M_entorhinal_thickavg <- rowMeans(SCZ[,c("L_entorhinal_thickavg" , "R_entorhinal_thickavg")]);
SCZ$M_fusiform_thickavg <- rowMeans(SCZ[,c("L_fusiform_thickavg" , "R_fusiform_thickavg")]);
SCZ$M_inferiorparietal_thickavg <- rowMeans(SCZ[,c("L_inferiorparietal_thickavg" , "R_inferiorparietal_thickavg")]);
SCZ$M_inferiortemporal_thickavg <- rowMeans(SCZ[,c("L_inferiortemporal_thickavg" , "R_inferiortemporal_thickavg")]);
SCZ$M_isthmuscingulate_thickavg <- rowMeans(SCZ[,c("L_isthmuscingulate_thickavg" , "R_isthmuscingulate_thickavg")]);
SCZ$M_lateraloccipital_thickavg <- rowMeans(SCZ[,c("L_lateraloccipital_thickavg" , "R_lateraloccipital_thickavg")]);
SCZ$M_lateralorbitofrontal_thickavg <- rowMeans(SCZ[,c("L_lateralorbitofrontal_thickavg" , "R_lateralorbitofrontal_thickavg")]);
SCZ$M_lingual_thickavg <- rowMeans(SCZ[,c("L_lingual_thickavg" , "R_lingual_thickavg")]);
SCZ$M_medialorbitofrontal_thickavg <- rowMeans(SCZ[,c("L_medialorbitofrontal_thickavg" , "R_medialorbitofrontal_thickavg")]);
SCZ$M_middletemporal_thickavg <- rowMeans(SCZ[,c("L_middletemporal_thickavg" , "R_middletemporal_thickavg")]);
SCZ$M_parahippocampal_thickavg <- rowMeans(SCZ[,c("L_parahippocampal_thickavg" , "R_parahippocampal_thickavg")]);
SCZ$M_paracentral_thickavg 	<- rowMeans(SCZ[,c("L_paracentral_thickavg" , "R_paracentral_thickavg")]);
SCZ$M_parsopercularis_thickavg <- rowMeans(SCZ[,c("L_parsopercularis_thickavg" , "R_parsopercularis_thickavg")]);
SCZ$M_parsorbitalis_thickavg <- rowMeans(SCZ[,c("L_parsorbitalis_thickavg" , "R_parsorbitalis_thickavg")]);
SCZ$M_parstriangularis_thickavg <- rowMeans(SCZ[,c("L_parstriangularis_thickavg" , "R_parstriangularis_thickavg")]);
SCZ$M_pericalcarine_thickavg <- rowMeans(SCZ[,c("L_pericalcarine_thickavg" , "R_pericalcarine_thickavg")]);
SCZ$M_postcentral_thickavg <- rowMeans(SCZ[,c("L_postcentral_thickavg" , "R_postcentral_thickavg")]);
SCZ$M_posteriorcingulate_thickavg <- rowMeans(SCZ[,c("L_posteriorcingulate_thickavg" , "R_posteriorcingulate_thickavg")]);
SCZ$M_precentral_thickavg <- rowMeans(SCZ[,c("L_precentral_thickavg" , "R_precentral_thickavg")]);
SCZ$M_precuneus_thickavg <- rowMeans(SCZ[,c("L_precuneus_thickavg" , "R_precuneus_thickavg")]);
SCZ$M_rostralanteriorcingulate_thickavg <- rowMeans(SCZ[,c("L_rostralanteriorcingulate_thickavg" , "R_rostralanteriorcingulate_thickavg")]);
SCZ$M_rostralmiddlefrontal_thickavg <- rowMeans(SCZ[,c("L_rostralmiddlefrontal_thickavg" , "R_rostralmiddlefrontal_thickavg")]);
SCZ$M_superiorfrontal_thickavg <- rowMeans(SCZ[,c("L_superiorfrontal_thickavg" , "R_superiorfrontal_thickavg")]);
SCZ$M_superiorparietal_thickavg <- rowMeans(SCZ[,c("L_superiorparietal_thickavg" , "R_superiorparietal_thickavg")]);
SCZ$M_superiortemporal_thickavg <- rowMeans(SCZ[,c("L_superiortemporal_thickavg" , "R_superiortemporal_thickavg")]);
SCZ$M_supramarginal_thickavg <- rowMeans(SCZ[,c("L_supramarginal_thickavg" , "R_supramarginal_thickavg")]);
SCZ$M_frontalpole_thickavg <- rowMeans(SCZ[,c("L_frontalpole_thickavg" , "R_frontalpole_thickavg")]);
SCZ$M_temporalpole_thickavg <- rowMeans(SCZ[,c("L_temporalpole_thickavg" , "R_temporalpole_thickavg")]);
SCZ$M_transversetemporal_thickavg <- rowMeans(SCZ[,c("L_transversetemporal_thickavg" , "R_transversetemporal_thickavg")]);
SCZ$M_insula_thickavg <- rowMeans(SCZ[,c("L_insula_thickavg" , "R_insula_thickavg")]);
SCZ$M_transversetemporal_thickavg <- rowMeans(SCZ[,c("L_transversetemporal_thickavg" , "R_transversetemporal_thickavg")]);
SCZ$M_bankssts_thickavg <- rowMeans(SCZ[,c("L_bankssts_thickavg" , "R_bankssts_thickavg")]);
SCZ$MThickness <- rowMeans(SCZ[,c("LThickness" , "RThickness")]);

SCZ$M_caudalanteriorcingulate_surfavg <- rowMeans(SCZ[,c("L_caudalanteriorcingulate_surfavg" , "R_caudalanteriorcingulate_surfavg")]);
SCZ$M_caudalmiddlefrontal_surfavg <- rowMeans(SCZ[,c("L_caudalmiddlefrontal_surfavg" , "R_caudalmiddlefrontal_surfavg")]);
SCZ$M_cuneus_surfavg <- rowMeans(SCZ[,c("L_cuneus_surfavg" , "R_cuneus_surfavg")]);
SCZ$M_entorhinal_surfavg <- rowMeans(SCZ[,c("L_entorhinal_surfavg" , "R_entorhinal_surfavg")]);
SCZ$M_fusiform_surfavg <- rowMeans(SCZ[,c("L_fusiform_surfavg" , "R_fusiform_surfavg")]);
SCZ$M_inferiorparietal_surfavg <- rowMeans(SCZ[,c("L_inferiorparietal_surfavg" , "R_inferiorparietal_surfavg")]);
SCZ$M_inferiortemporal_surfavg <- rowMeans(SCZ[,c("L_inferiortemporal_surfavg" , "R_inferiortemporal_surfavg")]);
SCZ$M_isthmuscingulate_surfavg <- rowMeans(SCZ[,c("L_isthmuscingulate_surfavg" , "R_isthmuscingulate_surfavg")]);
SCZ$M_lateraloccipital_surfavg <- rowMeans(SCZ[,c("L_lateraloccipital_surfavg" , "R_lateraloccipital_surfavg")]);
SCZ$M_lateralorbitofrontal_surfavg <- rowMeans(SCZ[,c("L_lateralorbitofrontal_surfavg" , "R_lateralorbitofrontal_surfavg")]);
SCZ$M_lingual_surfavg <- rowMeans(SCZ[,c("L_lingual_surfavg" , "R_lingual_surfavg")]);
SCZ$M_medialorbitofrontal_surfavg <- rowMeans(SCZ[,c("L_medialorbitofrontal_surfavg" , "R_medialorbitofrontal_surfavg")]);
SCZ$M_middletemporal_surfavg <- rowMeans(SCZ[,c("L_middletemporal_surfavg" , "R_middletemporal_surfavg")]);
SCZ$M_parahippocampal_surfavg <- rowMeans(SCZ[,c("L_parahippocampal_surfavg" , "R_parahippocampal_surfavg")]);
SCZ$M_paracentral_surfavg 	<- rowMeans(SCZ[,c("L_paracentral_surfavg" , "R_paracentral_surfavg")]);
SCZ$M_parsopercularis_surfavg <- rowMeans(SCZ[,c("L_parsopercularis_surfavg" , "R_parsopercularis_surfavg")]);
SCZ$M_parsorbitalis_surfavg <- rowMeans(SCZ[,c("L_parsorbitalis_surfavg" , "R_parsorbitalis_surfavg")]);
SCZ$M_parstriangularis_surfavg <- rowMeans(SCZ[,c("L_parstriangularis_surfavg" , "R_parstriangularis_surfavg")]);
SCZ$M_pericalcarine_surfavg <- rowMeans(SCZ[,c("L_pericalcarine_surfavg" , "R_pericalcarine_surfavg")]);
SCZ$M_postcentral_surfavg <- rowMeans(SCZ[,c("L_postcentral_surfavg" , "R_postcentral_surfavg")]);
SCZ$M_posteriorcingulate_surfavg <- rowMeans(SCZ[,c("L_posteriorcingulate_surfavg" , "R_posteriorcingulate_surfavg")]);
SCZ$M_precentral_surfavg <- rowMeans(SCZ[,c("L_precentral_surfavg" , "R_precentral_surfavg")]);
SCZ$M_precuneus_surfavg <- rowMeans(SCZ[,c("L_precuneus_surfavg" , "R_precuneus_surfavg")]);
SCZ$M_rostralanteriorcingulate_surfavg <- rowMeans(SCZ[,c("L_rostralanteriorcingulate_surfavg" , "R_rostralanteriorcingulate_surfavg")]);
SCZ$M_rostralmiddlefrontal_surfavg <- rowMeans(SCZ[,c("L_rostralmiddlefrontal_surfavg" , "R_rostralmiddlefrontal_surfavg")]);
SCZ$M_superiorfrontal_surfavg <- rowMeans(SCZ[,c("L_superiorfrontal_surfavg" , "R_superiorfrontal_surfavg")]);
SCZ$M_superiorparietal_surfavg <- rowMeans(SCZ[,c("L_superiorparietal_surfavg" , "R_superiorparietal_surfavg")]);
SCZ$M_superiortemporal_surfavg <- rowMeans(SCZ[,c("L_superiortemporal_surfavg" , "R_superiortemporal_surfavg")]);
SCZ$M_supramarginal_surfavg <- rowMeans(SCZ[,c("L_supramarginal_surfavg" , "R_supramarginal_surfavg")]);
SCZ$M_frontalpole_surfavg <- rowMeans(SCZ[,c("L_frontalpole_surfavg" , "R_frontalpole_surfavg")]);
SCZ$M_temporalpole_surfavg <- rowMeans(SCZ[,c("L_temporalpole_surfavg" , "R_temporalpole_surfavg")]);
SCZ$M_transversetemporal_surfavg <- rowMeans(SCZ[,c("L_transversetemporal_surfavg" , "R_transversetemporal_surfavg")]);
SCZ$M_insula_surfavg <- rowMeans(SCZ[,c("L_insula_surfavg" , "R_insula_surfavg")]);
SCZ$M_transversetemporal_surfavg <- rowMeans(SCZ[,c("L_transversetemporal_surfavg" , "R_transversetemporal_surfavg")]);
SCZ$M_bankssts_surfavg <- rowMeans(SCZ[,c("L_bankssts_surfavg" , "R_bankssts_surfavg")]);

#sum area SurfArea
SCZ$FullSurfArea <- SCZ$LSurfArea + SCZ$RSurfArea

#rename columns

SCZ2 <- rename(SCZ, c("L_entorhinal_thickavg" = "L_entorhil_thickavg", "R_entorhinal_thickavg" = "R_entorhil_thickavg", 
                      "L_entorhinal_surfavg" = "L_entorhil_surfavg", "R_entorhinal_surfavg" = "R_entorhil_surfavg",
                      "L_supramarginal_thickavg" = "L_supramargil_thickavg", "R_supramarginal_thickavg" = "R_supramargil_thickavg",
                      "L_supramarginal_surfavg" = "L_supramargil_surfavg", "R_supramarginal_surfavg" = "R_supramargil_surfavg",
                      "M_entorhinal_thickavg" = "M_entorhil_thickavg", "M_entorhinal_surfavg" = "M_entorhil_surfavg",
                      "M_supramarginal_thickavg" = "M_supramargil_thickavg", "M_supramarginal_surfavg" = "M_supramargil_surfavg"))

##recode site with strings

levels(SCZ2$Site)[match("MGH",levels(SCZ2$Site))] <- 11
levels(SCZ2$Site)[match("UMN",levels(SCZ2$Site))] <- 12
levels(SCZ2$Site)[match("UNM",levels(SCZ2$Site))] <- 13

write.csv(SCZ2, "scz_temp.csv")
SCZ3 <- read.csv("scz_temp.csv")

##recode site with numbers so they don't overlap but keep an order

SCZ3$Site <- ifelse(grepl("GROUP", SCZ3$AMC),
                    gsub(0, 6, SCZ3$Site),
                    SCZ3$Site)
SCZ3$Site <- ifelse(grepl("GROUP", SCZ3$AMC),
                    gsub(1, 7, SCZ3$Site),
                    SCZ3$Site)
SCZ3$Site <- ifelse(grepl("Huilong", SCZ3$AMC),
                    gsub(1, 8, SCZ3$Site),
                    SCZ3$Site)
SCZ3$Site <- ifelse(grepl("Huilong", SCZ3$AMC),
                    gsub(3, 9, SCZ3$Site),
                    SCZ3$Site)

#recode 2 NA in MCIC dataset
#!!!check the cell coordinates if new version of db!!!

SCZ3[3025, 176][is.na(SCZ3[3025, 176])] <- 11
SCZ3[3069, 176][is.na(SCZ3[3069, 176])] <- 12

#add missing site in order
SCZ3$Site = ifelse(grepl("(CASSI)",SCZ3$SubjID),14,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(COBRE)",SCZ3$SubjID),15,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(EONK)",SCZ3$SubjID),16,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(ESO)",SCZ3$SubjID),17,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(Frankfurt)",SCZ3$SubjID),18,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(GAP)",SCZ3$SubjID),19,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(GIPSI)",SCZ3$SubjID),20,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(Singapore)",SCZ3$SubjID),21,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(CIAM)",SCZ3$SubjID),22,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(IGP)",SCZ3$SubjID),23,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(Madrid)",SCZ3$SubjID),24,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(MPRC)",SCZ3$SubjID),25,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(Oxford_)",SCZ3$SubjID),26,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(PAFIP1.5T)",SCZ3$SubjID),27,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(PAFIP3T)",SCZ3$SubjID),28,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(RSCZ)",SCZ3$SubjID),29,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(RomeSL)",SCZ3$SubjID),30,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(SCORE)",SCZ3$SubjID),31,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(UCISZ)",SCZ3$SubjID),32,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(UNIBA)",SCZ3$SubjID),33,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(UNIMAAS)",SCZ3$SubjID),34,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(UPENN)",SCZ3$SubjID),35,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(Zurich)",SCZ3$SubjID),36,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(FOR2107Muenster)",SCZ3$SubjID),37,SCZ3$Site)
SCZ3$Site = ifelse(grepl("(FOR2107Marburg)",SCZ3$SubjID),38,SCZ3$Site)
SCZ3$Site[is.na(SCZ3$Site)] <- 10  #NA left Huilong

#add constant to Site
SCZ3$Site <- as.numeric(SCZ3$Site) + 200

#rename IDs to match BP and MDD
SCZ3$SubjID <- gsub("FOR2107Muenster", "FOR2107", SCZ3$SubjID)
SCZ3$SubjID <- gsub("FOR2107Marburg", "FOR2107", SCZ3$SubjID)
SCZ3$SubjID <- gsub("_for2107", "", SCZ3$SubjID)
SCZ3$SubjID <- gsub("CIAM_CIAM", "CIAM", SCZ3$SubjID)
SCZ3$SubjID <- gsub("(^Singapore_.*\\d+)_1", "\\1", SCZ3$SubjID)
SCZ3$SubjID <- ifelse(grepl("0", SCZ3$Dx),
                     gsub("Singapore_Parietal", "U0", SCZ3$SubjID),
                     SCZ3$SubjID)
#Add working group
SCZ3['WG']='3'

#merge BD, MDD and SCZ databases

ALL <- rbind.fill(MDD_xd, BD2, SCZ3)

#Check for duplicated SubjIDs that may cause issues with merging data sets.
if(anyDuplicated(ALL[,c("SubjID")]) != 0) { stop(paste0('You have duplicate SubjIDs in your file.\nMake sure there are no repeat SubjIDs.\n')) }

#identify which IDs occurs more than 1

n_occur <- data.frame(table(ALL$SubjID))
n_occur_more1 <- n_occur[n_occur$Freq > 1,]
n_occur_more1_IDs <- rename(n_occur_more1, c("Var1" = "SubjID"))
duplicated_IDs <- merge(ALL, n_occur_more1_IDs, by="SubjID") 

library(dplyr)

#if any, remove duplicated SubjIDs rows
ALL2 <- distinct(ALL, SubjID, .keep_all= TRUE, na.rm = TRUE)

#add a new Site variable for each cohort (regardless of scan)
ALL2$Site_new <- 999

#MDD
ALL2$Site_new = ifelse(grepl("(Episca)",ALL2$SubjID),1,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(AFFDIS)",ALL2$SubjID),2,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(AMC)",ALL2$SubjID),3,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Barc_)",ALL2$SubjID),4,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(BiDirect)",ALL2$SubjID),5,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Calgary)",ALL2$SubjID),6,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(CARDIFF)",ALL2$SubjID),7,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(CLING)",ALL2$SubjID),8,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(CODE)",ALL2$SubjID),9,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Dub)",ALL2$SubjID),10,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Edin)",ALL2$SubjID),11,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(ETPB)",ALL2$SubjID),12,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(FIDMAG)",ALL2$SubjID),13,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(FOR2107)",ALL2$SubjID),14,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Gron)",ALL2$SubjID),15,ALL2$Site_new) 
ALL2$Site_new = ifelse(grepl("(Houst)",ALL2$SubjID),16,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Hiroshima)",ALL2$SubjID),17,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(LOND)",ALL2$SubjID),18,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Magd)",ALL2$SubjID),19,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(McMas)",ALL2$SubjID),20,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(MD_SFB)",ALL2$SubjID),21,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Melb)",ALL2$SubjID),22,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Minnesota)",ALL2$SubjID),23,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Moraldilemma)",ALL2$SubjID),24,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(MPIP)",ALL2$SubjID),25,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Munster)",ALL2$SubjID),26,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(NESDA)",ALL2$SubjID),27,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Novo)",ALL2$SubjID),28,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Oxford)",ALL2$SubjID),29,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Oxfordyoung)",ALL2$SubjID),30,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(QTIM)",ALL2$SubjID),31,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(SanRaffaele)",ALL2$SubjID),32,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Well)",ALL2$SubjID),33,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(SF_)",ALL2$SubjID),34,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(SHIP)",ALL2$SubjID),35,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(SHIPtrend)",ALL2$SubjID),36,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Singapore)",ALL2$SubjID),37,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Socat)",ALL2$SubjID),38,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(StanfFAA)",ALL2$SubjID),39,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(StanfT1wAggr)",ALL2$SubjID),40,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(SYD)",ALL2$SubjID),41,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(TAD)",ALL2$SubjID),42,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(TIGER)",ALL2$SubjID),43,ALL2$Site_new)

#BP
ALL2$Site_new = ifelse(grepl("(CIAM)",ALL2$SubjID),44,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Barcelona)",ALL2$SubjID),45,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(BAR_)",ALL2$SubjID),46,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(BD_)",ALL2$SubjID),47,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(BF)",ALL2$SubjID),48,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(BDII)",ALL2$SubjID),49,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Brazil)",ALL2$SubjID),50,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(ESNA)",ALL2$SubjID),51,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Halifax)",ALL2$SubjID),52,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(moodinflame)",ALL2$SubjID),53,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Paris)",ALL2$SubjID),54,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Penn)",ALL2$SubjID),55,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(SBP)",ALL2$SubjID),56,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(struct)",ALL2$SubjID),57,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Yale)",ALL2$SubjID),58,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(TAB)",ALL2$SubjID),59,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("^(U00)",ALL2$SubjID),37,ALL2$Site_new)

#SCZ
ALL2$Site_new = ifelse(grepl("(ASRB)",ALL2$SubjID),60,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(CASSI)",ALL2$SubjID),61,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(COBRE)",ALL2$SubjID),62,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(EONK)",ALL2$SubjID),63,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(ESO)",ALL2$SubjID),64,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Frankfurt)",ALL2$SubjID),65,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(GAP)",ALL2$SubjID),66,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(GIPSI)",ALL2$SubjID),67,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(GROUP)",ALL2$SubjID),68,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Huilong)",ALL2$SubjID),69,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(IGP)",ALL2$SubjID),70,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Madrid)",ALL2$SubjID),71,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(MPRC)",ALL2$SubjID),72,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(MCIC)",ALL2$SubjID),73,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Oxford_subj)",ALL2$SubjID),74,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Oxford_con)",ALL2$SubjID),75,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(PAFI)",ALL2$SubjID),76,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(RSCZ)",ALL2$SubjID),77,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(RomeSL)",ALL2$SubjID),78,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(SCORE)",ALL2$SubjID),79,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(UCISZ)",ALL2$SubjID),80,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(UNIBA)",ALL2$SubjID),81,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(UNIMAAS)",ALL2$SubjID),82,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(UPENN)",ALL2$SubjID),83,ALL2$Site_new)
ALL2$Site_new = ifelse(grepl("(Zurich)",ALL2$SubjID),84,ALL2$Site_new)

#remove subjs with no dx

ALL2 <- ALL2[!(is.na(ALL2$Dx) | ALL2$Dx==""), ]

##remove duplicated freesurfer subjs

ALL2 <- ALL2[!ALL2$SubjID == "U00304", ]
ALL2 <- ALL2[!ALL2$SubjID =="Singapore_SBIC0094", ] 
ALL2 <- ALL2[!ALL2$SubjID =="Singapore_Parietal0055", ]
ALL2 <- ALL2[!ALL2$SubjID =="Singapore_Parietal0058", ]
ALL2 <- ALL2[!ALL2$SubjID == "U00203", ]  
ALL2 <- ALL2[!ALL2$SubjID == "Singapore_SBIC0061", ]
ALL2 <- ALL2[!ALL2$SubjID == "U00279", ]
ALL2 <- ALL2[!ALL2$SubjID == "Singapore_SBIC0066", ]
ALL2 <- ALL2[!ALL2$SubjID == "U00232", ]
ALL2 <- ALL2[!ALL2$SubjID == "Singapore_SBIC0017", ]
ALL2 <- ALL2[!ALL2$SubjID == "U00278", ]
ALL2 <- ALL2[!ALL2$SubjID == "Singapore_SBIC0065", ]
ALL2 <- ALL2[!ALL2$SubjID == "U00228", ]
ALL2 <- ALL2[!ALL2$SubjID == "Singapore_SBIC0013", ]
ALL2 <- ALL2[!ALL2$SubjID =="U00220", ] 
ALL2 <- ALL2[!ALL2$SubjID =="Singapore_SBIC0005", ]
ALL2 <- ALL2[!ALL2$SubjID =="U00243", ]
ALL2 <- ALL2[!ALL2$SubjID == "Singapore_SBIC0029", ]
ALL2 <- ALL2[!ALL2$SubjID =="U00234", ]
ALL2 <- ALL2[!ALL2$SubjID =="Singapore_SBIC0020", ]
ALL2 <- ALL2[!ALL2$SubjID =="con_11", ]
ALL2 <- ALL2[!ALL2$SubjID =="pat_13", ]
ALL2 <- ALL2[!ALL2$SubjID =="GAP_93", ]
ALL2 <- ALL2[!ALL2$SubjID =="GAP_94", ]
ALL2 <- ALL2[!ALL2$SubjID =="CIAM_CON_39", ]
ALL2 <- ALL2[!ALL2$SubjID =="CIAM_SCHIZ_01", ]
ALL2 <- ALL2[!ALL2$SubjID =="struct628", ]
ALL2 <- ALL2[!ALL2$SubjID =="struct710", ]
ALL2 <- ALL2[!ALL2$SubjID =="PAFIP1.5T_80277_001", ]
ALL2 <- ALL2[!ALL2$SubjID =="PAFIP1.5T_80330_001", ]
ALL2 <- ALL2[!ALL2$SubjID =="UPENN_4911", ]
ALL2 <- ALL2[!ALL2$SubjID =="UPENN_4956", ]
ALL2 <- ALL2[!ALL2$SubjID =="MPRC_Z90562", ]
ALL2 <- ALL2[!ALL2$SubjID =="MPRC_Z90582",]
ALL2 <- ALL2[!ALL2$SubjID =="AFFDIS_K004", ]
ALL2 <- ALL2[!ALL2$SubjID =="AFFDIS_U003", ]
ALL2 <- ALL2[!ALL2$SubjID == "MPIP_WG_413", ]
ALL2 <- ALL2[!ALL2$SubjID =="MPIP_WG_444", ]

############

write.csv(ALL2, "ALL2.csv")

####QC

cnames=colnames(ALL2);

#Check to make sure none of the age values are negative
agezero=which(ALL2$Age==0);
if(length(agezero)>0){
  stop("Some of your age values are zero This does not makes sense.");
}

#Check to make sure none of the sex values are negative
sexzero=which(ALL2$Sex==0);
if(length(sexzero)>0){
  stop("Some of your Sex values are zero This does not makes sense");
}

#Check to make sure none of the sex values are negative
sex2=which(ALL2$Sex>2);
if(length(sex2)>0){
  stop("Some of your Sex values are larger than 2. This does not makes sense");
}

#Check to make sure none of the Dx values are negative
Dx2=which(ALL2$Dx>3);
if(length(Dx2)>0){
  stop("Some of your Dx values are larger than 3. This does not makes sense");
}

ALL2_MDD = ALL2[which(ALL$Dx==1),]

AOzero=which(ALL2$AO==0);
if(length(AOzero)>0){
  stop("Some of your age of onset values are zero This does not makes sense.");
}

AOyoung=which(ALL2$AO<4);
if(length(AOyoung)>0){
  stop("Some of your age of onset values are younger than 4 This does not makes sense.");
}

#Loop through the different structures
for(x in 2:13){
  
  #Check to make sure none of the values are negative
  negs=which(ALL2[,x]<0);
  if(length(negs)>0){
    stop("Some of your volume values are negative. This does not makes sense. Open your CorticalMeasuresENIGMA_ThickAvg.csv file \
         in Excel and set negative volume values and poorly segmented values to NA in the file.\n\n");
  }
} 
#Loop through the different structures
for(x in 3:237){
  
  badsegs=0;
  ind=which(ALL[,x]=="x")
  ind2=which(ALL[,x]=="X")
  index=c(ind,ind2);
  if(length(index) > 0){
    stop("There were values marked with x or X detected in your CorticalMeasuresENIGMA_ThickAvg.csv file. All poorly segmented \
         or otherwise missing values should be marked with the letters NA in the CorticalMeasuresENIGMA_ThickAvg.csv file.\n\n");
  }
  
  #Check to make sure there are not any missing values
  miss=which(ALL[,x]=="");
  if(length(miss)>0){
    stop("There were missing values detected in your CorticalMeasuresENIGMA_ThickAvg.csv file. Open your CorticalMeasuresENIGMA_ThickAvg.csv file \
         in Excel and locate any blank cells. Missing values should be marked with NA. \n\n");
  }
  
  #Check to make sure none of the values are negative
  negs=which(ALL[,x]<0);
  if(length(negs)>0){
    stop("Some of your volume values are negative. This does not makes sense. Open your CorticalMeasuresENIGMA_ThickAvg.csv file \
         in Excel and set negative volume values and poorly segmented values to NA in the file.\n\n");
  }
  
  #Find which subjects are marked as NA for a given structure in the loop
  nas=which(is.na(ALL[,x]));
  if(length(nas)>0){
    interm=ALL[-nas,x];
    badsegs=badsegs+length(nas)
    cat(paste("You marked ", as.character(badsegs), " subjects as poorly segmented in the ", cnames[x], "\n", sep=''));
  } else {
    interm=ALL[,x];
    cat(paste("None of the subjects in the ", cnames[x], " were marked as poorly segmented\n", sep=''));
  }
}

###plots###

library(ggplot2)

ALL2 <- read.csv("ALL2.csv")

ALL2$Site <- as.factor(ALL2$Site) #factorize site column for plotting
ALL2$WG <- as.factor(ALL2$WG) #factorize WG column for plotting
ALL2$Dx <- as.factor(ALL2$Dx) #factorize Dx column for plotting

# colour blind palette

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

###boxplots###

##replace outliers with NA for plotting purposes

start <- which(colnames(ALL2)=="LLatVent")
end <- which(colnames(ALL2)=="FullSurfArea")
#create vector with these row numbers
regions <- start:end
#divide by site
df_list <- split(ALL2, as.factor(ALL2$Site))
#loop through regions 
for (j in df_list){
  outliers <-list()
  #run for loop for each region, get a vector of row numbers that are outliers
  for (i in regions){
    outliers[[(i-(start-1))]] <- boxplot(j[,i], plot=FALSE, range=2.5)$out
    #for each outlier, get the row number in the variable excl and then set that row to NA
    for (h in 1:length(outliers[[(i-(start-1))]])){
      excl <- which(ALL2[,regions[(i-(start-1))]]==outliers[[(i-(start-1))]][h])
      ALL2[excl,regions[(i-(start-1))]] <- NA
    }
  }
}

##boxplots

start <- which(colnames(ALL2)=="LLatVent")
end <- which(colnames(ALL2)=="FullSurfArea")
for (i in start:end){
  plot <- ggplot(ALL2, aes(x=Site, y=ALL2[,i], fill=Dx)) + 
    scale_fill_manual(values=cbPalette)+
    geom_boxplot(outlier.colour="black") +
    scale_x_discrete(breaks = levels(ALL2$Site)[c(T,F,F)]) +
    xlab("Site")
  filename = paste0('/Users/Elena/Documents/Work/ENIGMA/crossdx/boxplot/', colnames(ALL2)[i], '.png')
  ggsave(filename = filename, plot = plot, width = 16, height = 9)
} 

###outliers detection###

## print list of outliers ACROSS sites
#this will print a csv file with the outliers (IDs) detected for each region (each region is a column)

ALL2 <- read.csv("ALL2.csv")

#get rows numbers for areas that we want outliers for (assuming the columns are all next to each other)
start <- which(colnames(ALL2)=="LLatVent")
end <- which(colnames(ALL2)=="FullSurfArea")
#create vector with these row numbers
regions <- start:end
outliers_subjID_list <- list()
count <- 0
#loop through regions and save subjID in ALL2frame outliers_subjID
for (i in regions){
  count <- count + 1
  OutVals <- boxplot(ALL2[,i], plot=FALSE, range=1.5)$out
  outliers_subjID_list[[count]] <- ALL2$SubjID[which(ALL2[,i] %in% OutVals)]
}
n.obs <- sapply(outliers_subjID_list, length)
seq.max <- seq_len(max(n.obs))
outliers_subjID <- as.data.frame(sapply(outliers_subjID_list, "[", i = seq.max))
colnames(outliers_subjID) <- colnames(ALL2)[regions]

write.csv(outliers_subjID, "list_outliers1.5IQC.csv")

## print list of outliers WITHIN sites 
#this will print csv file for each site, with the outliers (IDs) for each region

#start
ALL2 <- read.csv("ALL2.csv")
#get rows
start <- which(colnames(ALL2)=="LLatVent")
end <- which(colnames(ALL2)=="FullSurfArea")
#create vector with these row numbers
regions <- start:end
  #divide by site
  df_list <- split(ALL2, as.factor(ALL2$Site))
  #loop through regions and save subjID in ALL2frame outliers_subjID
  for (j in df_list){
    outliers_subjID_list <- list()
    count <- 0
    for (i in regions){
    count <- count + 1
    OutVals <- boxplot(j[,i], plot=FALSE, range=2.5)$out
    outliers_subjID_list[[count]] <- j$SubjID[which(j[,i] %in% OutVals)]
  }
  n.obs <- sapply(outliers_subjID_list, length)
  seq.max <- seq_len(max(n.obs))
  outliers_subjID <- as.data.frame(sapply(outliers_subjID_list, "[", i = seq.max))
  colnames(outliers_subjID) <- colnames(j)[regions]
    write.csv(outliers_subjID, paste0(unique(j$Site), ".csv"))
  }
  
#merge all csv from above
  
  filenames <- list.files(path="/Users/Elena/Documents/Work/ENIGMA/crossdx/sites_outliers/",pattern=".csv")
  fullpath=file.path("/Users/Elena/Documents/Work/ENIGMA/crossdx/sites_outliers/",filenames)
  dataset <- do.call("rbind",lapply(fullpath,FUN=function(files){ read.csv(files)}))
  
#counting how many outlier per region (only actual values)
  
 dataset %>%
   select_if(function(x) any(!is.na(x))) %>%
   summarise_all(funs(sum(!is.na(.)))) -> dataset
 
 write.csv(dataset, "sum_outliers_1.5.csv")
 