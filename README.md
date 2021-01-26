# Heavitree Quartzite detrital zircon data

![alt text][Heavitree_photo]

[Heavitree_photo]: https://github.com/cverdel/Heavitree_Quartzite_DZ_data/blob/main/IMG_1680.jpg?raw=true

The Neoproterozoic (Tonian) Heavitree Quartzite forms the lower part of the Neoproterozoic stratigraphic succession of the Amadeus Basin. Detrital zircon data have been collected from the Heavitree Quartzite as part of [several studies that go back almost 30 years.](https://github.com/cverdel/Heavitree_Quartzite_DZ_data/blob/main/Heavitree_DZ_references) I've compiled and rearranged the data into eight data tables in this repository. The following R code illustrates some ways of organising and plotting the data. 

Note: at several intermediate stages this script writes data tables to a temporary directory. To see the location of this directory, use tempdir()

```
require(tidyverse)
require(data.table)
require(IsoplotR)
require(provenance)

#Evaluates data from the following 8 datasets
#1  Blades_etal_2021      #2 Camacho_etal_2015    #3 Haines_etal_2016
#4  Hollis_etal_2013      #5 Kositcin_etal_2014   #6 Maidment_etal_2007
#7  Normington_etal_2016  #8 Zhao_etal_1992           

#The first part downloads the data
#Download data
data_url<-"https://github.com/cverdel/Heavitree_Quartzite_DZ_data/raw/main/Blades_etal_2021.csv"
temp<-tempfile()
download.file(data_url, temp, mode="wb")
#Read data
Blades_etal_2021<-read.csv(temp)

#Download data
data_url<-"https://github.com/cverdel/Heavitree_Quartzite_DZ_data/raw/main/Camacho_etal_2015_Dean.csv"
temp<-tempfile()
download.file(data_url, temp, mode="wb")
#Read data
Camacho_etal_2015<-read.csv(temp)

#Download data
data_url<-"https://github.com/cverdel/Heavitree_Quartzite_DZ_data/raw/main/Haines_etal_2016_Dean.csv"
temp<-tempfile()
download.file(data_url, temp, mode="wb")
#Read data
Haines_etal_2016<-read.csv(temp)

#Download data
data_url<-"https://github.com/cverdel/Heavitree_Quartzite_DZ_data/raw/main/Hollis_etal_2013_Heavitree.csv"
temp<-tempfile()
download.file(data_url, temp, mode="wb")
#Read data
Hollis_etal_2013<-read.csv(temp)

#Download data
data_url<-"https://github.com/cverdel/Heavitree_Quartzite_DZ_data/raw/main/Kositcin_etal_2014.csv"
temp<-tempfile()
download.file(data_url, temp, mode="wb")
#Read data
Kositcin_etal_2014<-read.csv(temp)

#Download data
data_url<-"https://github.com/cverdel/Heavitree_Quartzite_DZ_data/raw/main/Maidment_etal_2007_Heavitree.csv"
temp<-tempfile()
download.file(data_url, temp, mode="wb")
#Read data
Maidment_etal_2007<-read.csv(temp)

#Download data
data_url<-"https://github.com/cverdel/Heavitree_Quartzite_DZ_data/raw/main/Normington_etal_2016_Heavitree.csv"
temp<-tempfile()
download.file(data_url, temp, mode="wb")
#Read data
Normington_etal_2016<-read.csv(temp)

#Download data
data_url<-"https://github.com/cverdel/Heavitree_Quartzite_DZ_data/raw/main/Zhao_etal_1992_Heavitree.csv"
temp<-tempfile()
download.file(data_url, temp, mode="wb")
#Read data
Zhao_etal_1992<-read.csv(temp)

#Ages are calculated for each dataset using the Isoplot R package. Ages are calculated for each datset individually because the errors are reported slightly differently.
#1
subset <- Blades_etal_2021[c("U238Pb206", "U238Pb206_error_2s_abs", "Pb207Pb206", "Pb207Pb206_error_2s_abs")] #Subsets data for age calculation
data<-read.data(subset, method="U-Pb", format=2, ierr=2) #Reads subset in IsoplotR format
tUPb <- age(data) #Calculates ages and errors from isotope ratios. Errors are always 1 sigma, regardless of input format.
combined<-cbind(tUPb,Blades_etal_2021) #Combines ages back with original data
Blades_etal_2021<-combined %>% 
  rename(t.conc_error = "s[t.conc]",
         t.76_error = "s[t.76]",
         t.68_error = "s[t.68]")

#2
subset <- Camacho_etal_2015[c("U238Pb206", "U238Pb206_error_1s_rel", "Pb207Pb206", "Pb207Pb206_error_1s_rel")] #Subsets data for age calculation
data<-read.data(subset, method="U-Pb", format=2, ierr=3) #Reads subset in IsoplotR format
tUPb <- age(data) #Calculates ages and errors from isotope ratios. Errors are always 1 sigma, regardless of input format.
combined<-cbind(tUPb,Camacho_etal_2015) #Combines ages back with original data
Camacho_etal_2015<-combined %>% 
  rename(t.conc_error = "s[t.conc]",
         t.76_error = "s[t.76]",
         t.68_error = "s[t.68]")

#3
subset <- Haines_etal_2016[c("U238Pb206", "U238Pb206_error_1s_abs", "Pb207Pb206", "Pb207Pb206_error_1s_abs")] #Subsets data for age calculation
data<-read.data(subset, method="U-Pb", format=2, ierr=1) #Reads subset in IsoplotR format
tUPb <- age(data) #Calculates ages and errors from isotope ratios. Errors are always 1 sigma, regardless of input format.
combined<-cbind(tUPb,Haines_etal_2016) #Combines ages back with original data
Haines_etal_2016<-combined %>% 
  rename(t.conc_error = "s[t.conc]",
         t.76_error = "s[t.76]",
         t.68_error = "s[t.68]")

#4
subset <- Hollis_etal_2013[c("U238Pb206", "U238Pb206_error_1s_abs", "Pb207Pb206", "Pb207Pb206_error_1s_abs")] #Subsets data for age calculation
data<-read.data(subset, method="U-Pb", format=2, ierr=1) #Reads subset in IsoplotR format
tUPb <- age(data) #Calculates ages and errors from isotope ratios. Errors are always 1 sigma, regardless of input format.
combined<-cbind(tUPb,Hollis_etal_2013) #Combines ages back with original data
Hollis_etal_2013<-combined %>% 
  rename(t.conc_error = "s[t.conc]",
         t.76_error = "s[t.76]",
         t.68_error = "s[t.68]")

#5
subset <- Kositcin_etal_2014[c("U238Pb206", "U238Pb206_error_1s_rel", "Pb207Pb206", "Pb207Pb206_error_1s_rel")] #Subsets data for age calculation
data<-read.data(subset, method="U-Pb", format=2, ierr=3) #Reads subset in IsoplotR format
tUPb <- age(data) #Calculates ages and errors from isotope ratios. Errors are always 1 sigma, regardless of input format.
combined<-cbind(tUPb,Kositcin_etal_2014) #Combines ages back with original data
Kositcin_etal_2014<-combined %>% 
  rename(t.conc_error = "s[t.conc]",
         t.76_error = "s[t.76]",
         t.68_error = "s[t.68]")

#6
subset <- Maidment_etal_2007[c("U238Pb206", "U238Pb206_error_1s_abs", "Pb207Pb206", "Pb207Pb206_error_1s_abs")] #Subsets data for age calculation
data<-read.data(subset, method="U-Pb", format=2, ierr=1) #Reads subset in IsoplotR format
tUPb <- age(data) #Calculates ages and errors from isotope ratios. Errors are always 1 sigma, regardless of input format.
combined<-cbind(tUPb,Maidment_etal_2007) #Combines ages back with original data
Maidment_etal_2007<-combined %>% 
  rename(t.conc_error = "s[t.conc]",
         t.76_error = "s[t.76]",
         t.68_error = "s[t.68]")

#7
subset <- Normington_etal_2016[c("U238Pb206", "U238Pb206_error_1s_rel", "Pb207Pb206", "Pb207Pb206_error_1s_rel")] #Subsets data for age calculation
data<-read.data(subset, method="U-Pb", format=2, ierr=3) #Reads subset in IsoplotR format
tUPb <- age(data) #Calculates ages and errors from isotope ratios. Errors are always 1 sigma, regardless of input format.
combined<-cbind(tUPb,Normington_etal_2016) #Combines ages back with original data
Normington_etal_2016<-combined %>% 
  rename(t.conc_error = "s[t.conc]",
         t.76_error = "s[t.76]",
         t.68_error = "s[t.68]")

#8
subset <- Zhao_etal_1992[c("U238Pb206", "U238Pb206_error_1s_rel", "Pb207Pb206", "Pb207Pb206_error_1s_rel")] #Subsets data for age calculation
data<-read.data(subset, method="U-Pb", format=2, ierr=3) #Reads subset in IsoplotR format
tUPb <- age(data) #Calculates ages and errors from isotope ratios. Errors are always 1 sigma, regardless of input format.
combined<-cbind(tUPb,Zhao_etal_1992) #Combines ages back with original data
Zhao_etal_1992<-combined %>% 
  rename(t.conc_error = "s[t.conc]",
         t.76_error = "s[t.76]",
         t.68_error = "s[t.68]")

#This part merges the 8 datasets into 1
merged<-merge(Blades_etal_2021, Camacho_etal_2015, all = TRUE)
merged<-merge(merged, Haines_etal_2016, all = TRUE)
merged<-merge(merged, Hollis_etal_2013, all = TRUE)
merged<-merge(merged, Kositcin_etal_2014, all = TRUE)
merged<-merge(merged, Maidment_etal_2007, all = TRUE)
merged<-merge(merged, Normington_etal_2016, all = TRUE)
merged<-merge(merged, Zhao_etal_1992, all = TRUE)

write.csv(merged, paste0(tempdir(), "/", "merged_data.csv"), row.names=T) #Writes the merged dataset to a temp directory
tempdir() #Location of the temp directory

#Make a plot showing the sources of data
merged$Reference <- fct_infreq(merged$Reference) #Sorts Reference by frequency
p<-ggplot(merged, aes(Reference, ..count..)) +
  geom_bar(colour="black", fill="lightblue")+
  ylab("Number of zircons")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1))+
  theme(axis.title.x = element_blank())
p
```
![alt text][Heavitree_plot1]

[Heavitree_plot1]: https://github.com/cverdel/Heavitree_Quartzite_DZ_data/blob/main/Heavitree_Rplot1.jpg?raw=true


#Calculates discordance
merged$ratio86<-as.numeric(merged$t.68/merged$t.76) #238/206 age divided by #207/206 age

#Discordance filter
Concordance_co=0.10 #Chooses concordance cutoff (10%, 15%, etc.)
upper_limit=1+Concordance_co
lower_limit=1-Concordance_co

#Creates a new column called "Concordance" with y for concordant data, n for discordant data
merged$Concordance<-ifelse(merged$ratio86 <=upper_limit & merged$ratio86 >=lower_limit, "y", "n")

#Plots a simple (i.e., no error ellipses) Tera-Wasserburg diagram
p<-ggplot(data=merged, aes(x=U238Pb206,y=Pb207Pb206, colour=Concordance))+
  geom_point()+
  theme_bw()
p

#Creates 2 new datatables for concordant and discordant results
merged_concordant<-subset(merged, Concordance=="y")
merged_discordant<-subset(merged, Concordance=="n")
write.csv(merged_concordant, paste0(tempdir(), "/", "concordant_analyses.csv"), row.names=T) #Writes the concordant data table to a temp directory

#Changes name of merged_concordant to df
df<-merged_concordant

#Adds a column with 1:250k mapsheet based on sample location
df$map_sheet<-ifelse(df$Longitude >127.5 & df$Longitude <129   & df$Latitude < -22 & df$Latitude > -23, 'WEBB', 
              ifelse(df$Longitude >127.5 & df$Longitude <129   & df$Latitude < -23 & df$Latitude > -24, 'MACDONALD', 
              ifelse(df$Longitude >127.5 & df$Longitude <129   & df$Latitude < -24 & df$Latitude > -25, 'RAWLINSON', 
              ifelse(df$Longitude >129   & df$Longitude <130.5 & df$Latitude < -23 & df$Latitude > -24, 'MOUNT RENNIE', 
              ifelse(df$Longitude >129   & df$Longitude <130.5 & df$Latitude < -24 & df$Latitude > -25, 'BLOODS RANGE', 
              ifelse(df$Longitude >129   & df$Longitude <130.5 & df$Latitude < -25 & df$Latitude > -26, 'PETERMANN RANGES', 
              ifelse(df$Longitude >130.5 & df$Longitude <132   & df$Latitude < -23 & df$Latitude > -24, 'MOUNT LIEBIG', 
              ifelse(df$Longitude >130.5 & df$Longitude <132   & df$Latitude < -24 & df$Latitude > -25, 'LAKE AMADEUS', 
              ifelse(df$Longitude >130.5 & df$Longitude <132   & df$Latitude < -25 & df$Latitude > -26, 'AYERS ROCK', 
              ifelse(df$Longitude >132   & df$Longitude <133.5 & df$Latitude < -23 & df$Latitude > -24, 'HERMANNSBURG', 
              ifelse(df$Longitude >132   & df$Longitude <133.5 & df$Latitude < -24 & df$Latitude > -25, 'HENBURY', 
              ifelse(df$Longitude >132   & df$Longitude <133.5 & df$Latitude < -25 & df$Latitude > -26, 'KULGERA', 
              ifelse(df$Longitude >133.5 & df$Longitude <135   & df$Latitude < -23 & df$Latitude > -24, 'ALICE SPRINGS', 
              ifelse(df$Longitude >133.5 & df$Longitude <135   & df$Latitude < -24 & df$Latitude > -25, 'RODINGA', 
              ifelse(df$Longitude >133.5 & df$Longitude <135   & df$Latitude < -25 & df$Latitude > -26, 'FINKE', 
              ifelse(df$Longitude >135   & df$Longitude <136.5 & df$Latitude < -23 & df$Latitude > -24, 'ILLOGWA CREEK', 
              ifelse(df$Longitude >135   & df$Longitude <136.5 & df$Latitude < -24 & df$Latitude > -25, 'HALE RIVER', 
              'NA')))))))))))))))))

#Make a plot showing which mapsheets the samples are from
df$map_sheet <- fct_infreq(df$map_sheet) #Sorts Region by frequency
p<-ggplot(df, aes(map_sheet, ..count..)) +
  geom_bar(colour="black", fill="lightblue")+
  ylab("Number of concordant zircons")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10, angle=45, hjust=1))+
  theme(axis.title.x = element_blank())
p

#Make density plots for each sample
p<-ggplot(data=merged_concordant, aes(x=t.conc))+
  geom_density(colour="black", alpha=0.4, fill="#3d5fff")+
  geom_rug(outside=FALSE)+
  scale_x_continuous(breaks = seq(1000,3000, by=2000))+
  theme_bw()+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  ylab("Density")+
  xlab("Äge (Ma)")+
  facet_wrap(~Sample_number,  scales="free_y", ncol=7)
dev.new()
p

#Calculate MDS coordinates
#Rearrange table first
df2 <- df[c("Sample_number", "t.conc")] #Subsets data 
df2<-df2%>% dplyr::mutate(ID=row_number()) #Re-shapes data
df3<-spread(df2, Sample_number, t.conc, fill = NA, convert = FALSE)
df3 <- df3[ -c(1) ]
df3 <- data.table(df3)[, lapply(.SD, function(x) x[order(is.na(x))])]
df4<-df3[!df3[, Reduce(`&`, lapply(.SD, is.na))]] #df4 is in MDS format with sample numbers at the top

write.csv(df4, paste0(tempdir(), "/", "concordant_MDS.csv"), row.names=T) #Writes to a temp directory

#Read data
MDS_data<-read.distributional(paste0(tempdir(), "/", "concordant_MDS.csv")) 

#Could skip this next section if you didn't want to add synthetic spectra
#Function needed to combine synthetic spectra
cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

#Create synthetic spectra
s600<-as.data.frame(runif(100, min=550, max=650)) #100 analyses between 550 and 650 Ma
colnames(s600) <- c("600")
s1050<-as.data.frame(runif(100, min=1000, max=1100))
colnames(s1050) <- c("1050")
s1200<-as.data.frame(runif(100, min=1150, max=1250))
colnames(s1200) <- c("1200")
s1600<-as.data.frame(runif(100, min=1550, max=1650))
colnames(s1600) <- c("1600")
s1800<-as.data.frame(runif(100, min=1750, max=1850))
colnames(s1800) <- c("1800")

alldata<-as.data.frame(cbind.fill(df4, s600, s1050, s1200, s1600, s1800)) #Combines MDS coordinates of samples and synthetics
write.csv(alldata, paste0(tempdir(), "/", "alldata_MDS.csv"), row.names=F) 

#Read data
MDS_data<-read.distributional(paste0(tempdir(), "/", "alldata_MDS.csv")) 

#Calculate MDS
MDS_Heavitree<-MDS(MDS_data)

###Put MDS coordinates into a dataframe
all_coord<-as.data.frame(MDS_Heavitree$points) 
all_coord<-setDT(all_coord, keep.rownames = "Sample_number")[] #3 columns: 1st is sample number, then the two MDS coordinates, V1 and V2

row.names.remove <- as.character(c("600", "1050", "1200", "1600", "1800"))

Heavitree_coord_wo_synthetics<-all_coord[!(all_coord$Sample_number %in% row.names.remove), ] #Dataset without synthetics
synthetics_coord<-all_coord[(all_coord$Sample_number  %in% row.names.remove), ] #Dataset with only synthetics

#Adds MDS coordinates to the data table of concordant analyses that was created before (df)
df<-merge(df, Heavitree_coord_wo_synthetics, all = TRUE) 

sample_parameters<-df %>% distinct(Sample_number, .keep_all = TRUE) #Condenses dataframe down to 1 row per sample
sample_parameters<-subset(sample_parameters, select=c(Sample_number, Formation, map_sheet, Reference, Latitude, Longitude, 
                                                      V1, V2))#Selects important parameters for sample dataframe

#Creates a sample summary table in the temp directory
write.csv(sample_parameters, paste0(tempdir(), "/", "Heavitree_sample_summary.csv"), row.names=F) 
tempdir()


#MDS plotting by mapsheet
max_x<-max(abs(all_coord$V1)) #MDS V1 (horizontal) coordinate with the greatest absolute value
max_y<-max(abs(all_coord$V2)) ##MDS V2 (vertical) coordinate with the greatest absolute value
buffer=0.1
factor_x<-0.75/(max_x+buffer)
factor_y<-0.5/(max_y+buffer)

sample_parameters$rescale_V1<-sample_parameters$V1*factor_x
sample_parameters$rescale_V2<-sample_parameters$V2*factor_y

#Adds columns with latitude and longitude of 1:250k mapsheet centres
sample_parameters$mapsheet_c_lat<-as.numeric(ifelse(sample_parameters$map_sheet=='HERMANNSBURG', -23.5,
                     ifelse(sample_parameters$map_sheet=='ALICE SPRINGS', -23.5,
                            ifelse(sample_parameters$map_sheet=='ILLOGWA CREEK', -23.5, 
                                   ifelse(sample_parameters$map_sheet=='RAWLINSON', -24.5, 
                                          ifelse(sample_parameters$map_sheet=='WEBB', -22.5 , 
                                                 ifelse(sample_parameters$map_sheet=='BLOODS RANGE', -24.5 ,
                                                 'NA')))))))

sample_parameters$mapsheet_c_long<-as.numeric(ifelse(sample_parameters$map_sheet=='HERMANNSBURG', 132.75,
                                         ifelse(sample_parameters$map_sheet=='ALICE SPRINGS', 134.25,
                                                ifelse(sample_parameters$map_sheet=='ILLOGWA CREEK', 135.75, 
                                                       ifelse(sample_parameters$map_sheet=='RAWLINSON', 128.25, 
                                                              ifelse(sample_parameters$map_sheet=='WEBB', 128.25 , 
                                                                     ifelse(sample_parameters$map_sheet=='BLOODS RANGE', 129.75 ,
                                                                            'NA')))))))

#Adds re-scaled MDS values to map centre coordinates
sample_parameters$y<-sample_parameters$mapsheet_c_lat+sample_parameters$rescale_V2
sample_parameters$x<-sample_parameters$mapsheet_c_long+sample_parameters$rescale_V1

write.csv(sample_parameters, paste0(tempdir(), "/", "Heavitree_sample_summary_w_geographic_MDS.csv"), row.names=F) 

ggplot(sample_parameters, aes(x=x, y=y))+
  geom_point()

#Does the same for synthetic spectra
synthetics_coord$rescale_V1<-synthetics_coord$V1*factor_x
synthetics_coord$rescale_V2<-synthetics_coord$V2*factor_y

#HERMANNSBURG synthetics
synthetics_coord$y<-synthetics_coord$rescale_V2 + (-23.5)
synthetics_coord$x<-synthetics_coord$rescale_V1 + 132.75
Hermannsburg_synthetics<-synthetics_coord

#ALICE SPRINGS synthetics
synthetics_coord$y<-synthetics_coord$rescale_V2 + (-23.5)
synthetics_coord$x<-synthetics_coord$rescale_V1 + 134.25
Alice_Springs_synthetics<-synthetics_coord

#ILLOGWA CREEK synthetics
synthetics_coord$y<-synthetics_coord$rescale_V2 + (-23.5)
synthetics_coord$x<-synthetics_coord$rescale_V1 + 135.75
Illogwa_Creek_synthetics<-synthetics_coord

#RAWLINSON synthetics
synthetics_coord$y<-synthetics_coord$rescale_V2 + (-24.5)
synthetics_coord$x<-synthetics_coord$rescale_V1 + 128.25
Rawlinson_synthetics<-synthetics_coord

#WEBB synthetics
synthetics_coord$y<-synthetics_coord$rescale_V2 + (-22.5)
synthetics_coord$x<-synthetics_coord$rescale_V1 + 128.25
Webb_synthetics<-synthetics_coord

#BLOODS RANGE synthetics
synthetics_coord$y<-synthetics_coord$rescale_V2 + (-24.5)
synthetics_coord$x<-synthetics_coord$rescale_V1 + 129.75
Bloods_Range_synthetics<-synthetics_coord

#Merge
synthetics_rescaled<-merge(Hermannsburg_synthetics, Alice_Springs_synthetics, all = TRUE)
synthetics_rescaled<-merge(synthetics_rescaled, Illogwa_Creek_synthetics, all = TRUE)
synthetics_rescaled<-merge(synthetics_rescaled, Rawlinson_synthetics, all = TRUE)
synthetics_rescaled<-merge(synthetics_rescaled, Webb_synthetics, all = TRUE)
synthetics_rescaled<-merge(synthetics_rescaled, Bloods_Range_synthetics, all = TRUE)

write.csv(synthetics_rescaled, paste0(tempdir(), "/", "Synthetic_spectra_w_geographic_MDS.csv"), row.names=F) 
tempdir()
```
