#Full Rethomics tutorial @ https://rethomics.github.io/

install.packages("ggpubr", "damr", "ggetho")

DATA_DIR<- "/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data" 			#create working directory --> set to path with activity monitor files in it
list.files(DATA_DIR, pattern= "*.txt|*.csv")	#check files in working directory --> make sure they included monitor files
setwd(DATA_DIR) 								#set working directory

library(ggplot2)
library(damr)									#use to read in DAM data
library(ggetho)								#use to plot and check data
library(ggpubr)								#use to mulitple plot ggplots

metadata <- fread("Dmel+Dsec_all_1212LD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])						#load data into behavior table
summary(dt)

#library(zeitgebr)
per_dt <- periodogram(activity, dt, FUN = chi_sq_periodogram, resample_rate = 1/mins(5))
#per_dt
#ggperio(per_dt, aes(period, power, colour=strain)) + 
#  stat_pop_etho(colour = c(rgb(66, 199, 244, 100, maxColorValue = 255))) +
#  theme_classic()

per_dt <- find_peaks(per_dt)

#ggperio(per_dt) + geom_line(aes(group = id, colour=strain))

ggperio(per_dt) + 
  geom_line(aes(group = id, colour = strain)) +
  geom_peak(col = "black") +
  geom_line(aes(y = signif_threshold)) +
  facet_wrap(~ region_id, ncol = 8)

#change strain name to check for dead flies, change status of dead flies in metadata file and re-run analysis

ggetho(dt[xmv(strain) == c("BT410-3")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
       stat_bar_tile_etho() +
       facet_wrap(~ region_id) +
       stat_ld_annotations(ld_colours =  c("lightyellow", "black"), outline  = NA, l_duration = hours(12), height = 1, alpha = 0.05)

#plot average activity by strain- tiles:      
       
ggetho(dt, aes(x = t, y = strain, z=activity)) +
       stat_tile_etho() +
       stat_ld_annotations(ld_colours =  c("lightyellow", "black"), outline  = NA, l_duration = hours(12))

#double actograms of individuals binned by 30 minutes
       
ggetho(dt[xmv(strain) == c("Dsim.196")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
       stat_bar_tile_etho() +
       #facet_wrap(~ region_id) +
       stat_ld_annotations(ld_colours =  c("lightyellow", "black"), l_duration = hours(12))

#double actograms of average activity binned by 30 minutes
       
ggetho(dt[xmv(strain) == c("Dsec.07")], aes(x = t, z=activity), summary_FUN = mean, summary_time_window = mins(30) , time_wrap = hours(24), time_offset = hours(6)) +
       stat_bar_tile_etho(method = mean) +
       theme_classic() 
       #stat_ld_annotations(ld_colours =  c("lightyellow", "black"), outline = NA, height = 1, alpha = 0.05, l_duration = hours(12))

#average activity by strain - line standardized by max activity of Dmel.CS              

DmelOR <- ggetho(dt[xmv(strain) == c("Dmel.OR")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-DmelOR$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmelOR<-max(colMeans(df_x))
DmelOR<- ggetho(dt[xmv(strain) == c("Dmel.OR")], aes(x = t, y = activity/maxDmelOR, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

DmelCS <- ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-DmelCS$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmelCS<-max(colMeans(df_x))
DmelCS<- ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, y = activity/maxDmelCS, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dsec07 <- ggetho(dt[xmv(strain) == c("Dsec.07")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dsec07$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec07<-max(colMeans(df_x))
Dsec07 <- ggetho(dt[xmv(strain) == c("Dsec.07")], aes(x = t, y = activity/maxDsec07, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dsec28 <- ggetho(dt[xmv(strain) == c("Dsec.28")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Dsec28$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec28<-max(colMeans(df_x))
Dsec28 <- ggetho(dt[xmv(strain) == c("Dsec.28")], aes(x = t, y = activity/maxDsec28, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.452, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

ggarrange(DmelCS+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), DmelOR+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dsec07+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dsec28+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          heights = c(3, 3, 3, 3), ncol = 1, nrow = 4)

dt[, MA := ifelse(t %% hours(24) < hours(21), "R", "MA")]
dt[, S := ifelse(t %% hours(24) < hours(18), "R", "S")]


summary_dt <- 
  rejoin(dt[,
               .(
                 # this is where the computation happens
                 activity_MA = mean(activity[MA == "MA"]),
                 activity_S = mean(activity[S == "S"]) 
               ),
               ,by=id])
summary_dt

p<- ggplot(summary_dt, aes(x=strain, y=activity_MA/activity_S, fill=strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "") + 
  theme_classic() + 
  scale_x_discrete(name = "", labels = c("","","","")) + 
  theme(legend.position = "none")
p
p$data$strain
p$data$activity_MA
x<-cbind(p$data$strain, p$data$activity_MA)
write.csv(x, "~/Desktop/x.csv")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
           ncol = 1, nrow = 1)

pairwise.wilcox.test(summary_dt[, activity_MA/max(activity_MA)], 
                     summary_dt[, strain], p.adjust = "bonf")

x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data/Extended photoperiod data/Automated_peakdetection/Dmel_Dsec_MPeak_1212LD.csv")

p<- ggplot(x, aes(x=Order, y=M1212, fill=Order), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.75),alpha=0.3) +
  scale_y_continuous(name= "Morning peak timing (hours from 'lights-on')") + 
  theme_classic() + 
  scale_x_discrete(name = "", labels = c(expression(italic("Dmel.CS")), expression(italic("Dmel.OR")), expression(italic("Dsec.07")), expression(italic("Dsec.28")))) + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")
p

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(x$M1212, 
                     x$Strain, p.adjust = "bonf")

##########################################################LD2DD analysis############################################################

metadata <- fread("metadataLD2DD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])						#load data into behavior table
summary(dt)

Dsec07<-ggetho(dt[xmv(strain) == c("Dsec.07")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho(colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("grey", "black"), outline  = NA, l_duration = hours(12), height = .02, alpha = 1) +
  theme_classic()

Dsec28<-ggetho(dt[xmv(strain) == c("Dsec.28")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho(colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("grey", "black"), outline  = NA, l_duration = hours(12), height = .02, alpha = 1) +
  theme_classic()

DmelCS<-ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho(colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("grey", "black"), outline  = NA, l_duration = hours(12), height = .02, alpha = 1) +
  theme_classic()

DmelOR<-ggetho(dt[xmv(strain) == c("Dmel.OR")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho(colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("grey", "black"), outline  = NA, l_duration = hours(12), height = .02, alpha = 1) +
  theme_classic()

ggarrange(DmelCS+ rremove("x.text")+rremove("xlab")+rremove("ylab")+ rremove("y.text"), DmelOR+ rremove("x.text")+rremove("xlab")+rremove("ylab")+ rremove("y.text"), Dsec07+ rremove("x.text")+rremove("xlab")+rremove("ylab")+ rremove("y.text"), Dsec28+ rremove("x.text")+rremove("xlab")+rremove("ylab")+ rremove("y.text"), 
          heights = c(3, 3, 3, 3), ncol = 1, nrow = 4)

##########################################################DD PERIODOGRAM ANALYSIS############################################################

metadata <- fread("metadataDD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[species == "Dmel" & status == "OK" | species == "Dsec" & status == "OK"])						#load data into behavior table
summary(dt)

library(zeitgebr)
per_dt <- periodogram(activity, dt, FUN = chi_sq_periodogram, resample_rate = 1/mins(5))
per_dt
ggperio(per_dt, aes(period, power, colour=strain)) + 
        stat_pop_etho(colour = c(rgb(66, 199, 244, 100, maxColorValue = 255))) +
        theme_classic()
         
per_dt <- find_peaks(per_dt)

ggperio(per_dt) + geom_line(aes(group = id, colour=strain))

ggperio(per_dt) + 
  geom_line(aes(group = id, colour = strain)) +
  geom_peak(col = "black") +
  geom_line(aes(y = signif_threshold)) +
  facet_wrap(~ region_id, ncol = 8)

summary_dt <- rejoin(per_dt[peak==1])
summary_dt

x<-cbind(summary_dt$strain, summary_dt$region_id,  summary_dt$period)
x
write.csv(x, "~/Desktop/periods.csv")
x<-read.csv("~/Desktop/periods.csv")
x

pairwise.wilcox.test(x$period, x$strain, p.adjust.method = "bonferroni")

#          Dmel.CS  Dmel.OR  Dsec.07
# Dmel.OR  1.5e-09  -        -      
# Dsec.07  9.8e-10  0.0113   -      
# Dsec.28  4.7e-06  0.5395   0.0041 

strain_means <- 
    summary_dt[,
             .(
               # this is where the computation happens
               average_peak = mean(period)/60/60
               ),
             by=strain]

strain_means  
#DD       
#1:  Dsec.28     23.56667
#2:  Dsec.07     23.16296
#3:  Dmel.CS     24.25517
#4:  Dmel.OR     23.45312
#5: Dsim.196     24.58636
#6:  Dsim.04     23.36538

################################################################EXTENDED PHOTOPERIOD ANALYSIS##################################################################

metadata <- fread("metadata_20-4LD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])						#load data into behavior table
summary(dt)

#library(zeitgebr)
per_dt <- periodogram(activity, dt, FUN = chi_sq_periodogram, resample_rate = 1/mins(5))
#per_dt
#ggperio(per_dt, aes(period, power, colour=strain)) + 
#  stat_pop_etho(colour = c(rgb(66, 199, 244, 100, maxColorValue = 255))) +
#  theme_classic()

per_dt <- find_peaks(per_dt)

#ggperio(per_dt) + geom_line(aes(group = id, colour=strain))

ggperio(per_dt) + 
  geom_line(aes(group = id, colour = strain)) +
  geom_peak(col = "black") +
  geom_line(aes(y = signif_threshold)) +
  facet_wrap(~ region_id, ncol = 8)

ggetho(dt[xmv(strain) == c("Dsec.28")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho() +
  facet_wrap(~ region_id) +
  stat_ld_annotations(ld_colours =  c("lightyellow", "black"), outline  = NA, l_duration = hours(12), height = 1, alpha = 0.05)

DmelOR <- ggetho(dt[xmv(strain) == c("Dmel.OR")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-DmelOR$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmelOR<-max(colMeans(df_x))
DmelOR<- ggetho(dt[xmv(strain) == c("Dmel.OR")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(20)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,1,2,3,4), limits = c(0,4)) 
#scale_y_continuous(name= "") 

DmelCS <- ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-DmelCS$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmelCS<-max(colMeans(df_x))
DmelCS<- ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(20)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,1,2,3,4), limits = c(0,4)) 
#scale_y_continuous(name= "") 

Dsec07 <- ggetho(dt[xmv(strain) == c("Dsec.07")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dsec07$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec07<-max(colMeans(df_x))
Dsec07 <- ggetho(dt[xmv(strain) == c("Dsec.07")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(20)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,1,2,3,4), limits = c(0,4)) 
#scale_y_continuous(name= "") 

Dsec28 <- ggetho(dt[xmv(strain) == c("Dsec.28")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Dsec28$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec28<-max(colMeans(df_x))
Dsec28 <- ggetho(dt[xmv(strain) == c("Dsec.28")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.452, outline  = NA, l_duration = hours(20)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,1,2,3,4), limits = c(0,4)) 
  #scale_y_continuous(name= "") 

ggarrange(DmelCS+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), DmelOR+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dsec07+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dsec28+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          heights = c(3, 3, 3, 3), ncol = 1, nrow = 4)

x<-read.csv("/Users/michaelshahandeh/Desktop/168Epeak_rep2.csv", header = TRUE)
p<- ggplot(x, aes(x=Order, y=E168, fill=Order), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "Evening peak time 14:10h LD(hours)", limits = c(10,21), breaks = c(10, 12, 14, 16, 18, 20)) + 
  theme_classic() + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(x$Epeak, x$Strain, p.adjust.method = "bonferroni")
wilcox.test(x$Epeak, x$Species)

x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data/Extended photoperiod data/Automated_peakdetection/14-1Epeak/DeltaME.csv", header = TRUE)
p<- ggplot(x, aes(x=Order, y=Epeak, fill =Order), inherit.aes = TRUE) +
  geom_boxplot( fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "Evening peak time 16:8h LD(hours)", limits = c(10,21), breaks = c(10, 12, 14, 16, 18, 20)) + 
  theme_classic() + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(x$Epeak, x$Strain, p.adjust.method = "bonferroni")
wilcox.test(x$Epeak, x$Species)

x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data/Extended photoperiod data/Automated_peakdetection/Epeak_1410.csv", header = TRUE)
p<- ggplot(x, aes(x=Order, y=Epeak, fill=Order), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "Evening peak time 16:8h LD(hours)", limits = c(10,21), breaks = c(10, 12, 14, 16, 18, 20)) + 
  theme_classic() + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)
pairwise.wilcox.test(x$Epeak, x$Strain, p.adjust.method = "bonferroni")


x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data/Extended photoperiod data/Automated_peakdetection/Epeak_168.csv", header = TRUE)
p<- ggplot(x, aes(x=Order, y=Epeak, fill=Order), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "Evening peak time 16:8h LD(hours)", limits = c(10,21), breaks = c(10, 12, 14, 16, 18, 20)) + 
  theme_classic() + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(x$Epeak, x$Strain, p.adjust.method = "bonferroni")
wilcox.test(x$Epeak, x$Species)

x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data/Extended photoperiod data/Automated_peakdetection/Epeak_186.csv", header = TRUE)
p<- ggplot(x, aes(x=Order, y=Epeak, fill=Order), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "Evening peak time 18:6h LD(hours)", limits = c(10,21), breaks = c(10, 12, 14, 16, 18, 20)) + 
  theme_classic() + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(x$Epeak, x$Strain, p.adjust.method = "bonferroni")
wilcox.test(x$Epeak, x$Species)

x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data/Extended photoperiod data/Automated_peakdetection/Epeak_204.csv", header = TRUE)
p<- ggplot(x, aes(x=Order, y=Epeak, fill=Order), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "Evening peak time 20:4h LD(hours)", limits = c(10,21), breaks = c(10, 12, 14, 16, 18, 20)) + 
  theme_classic() + 
  geom_hline(yintercept = 0) +
  theme(legend.position = "none")

pairwise.wilcox.test(x$Epeak, x$Strain, p.adjust.method = "bonferroni")
wilcox.test(x$Epeak, x$Species)

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

################################################################tropical Dmel/Dsim/Dmau 12h LD ANALYSIS##################################################################

metadata <- fread("metadata-12LDmau.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[replicate == "5"])

ggetho(dt[xmv(strain) == c("Dmau.91")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho(colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), outline  = NA, l_duration = hours(12), height = .02, alpha = 1) +
  theme_classic() +
  facet_wrap(~region_id)

metadata <- fread("metadata-12LDmau.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])


LZV76 <- ggetho(dt[xmv(strain) == c("Dmel.LZV_L76")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-LZV76$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxLZV76<-max(colMeans(df_x))
LZV76<- ggetho(dt[xmv(strain) == c("Dmel.LZV_L76")], aes(x = t, y = activity/maxLZV76, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

LZV72 <- ggetho(dt[xmv(strain) == c("Dmel.LZV_L72")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-LZV72$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxLZV72<-max(colMeans(df_x))
LZV72<- ggetho(dt[xmv(strain) == c("Dmel.LZV_L72")], aes(x = t, y = activity/maxLZV72, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

MD221 <- ggetho(dt[xmv(strain) == c("Dsim.MD221")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-MD221$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxMD221<-max(colMeans(df_x))
MD221<- ggetho(dt[xmv(strain) == c("Dsim.MD221")], aes(x = t, y = activity/maxMD221, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

MD242 <- ggetho(dt[xmv(strain) == c("Dsim.MD242")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-MD242$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxMD242<-max(colMeans(df_x))
MD242<- ggetho(dt[xmv(strain) == c("Dsim.MD242")], aes(x = t, y = activity/maxMD242, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dmau90 <- ggetho(dt[xmv(strain) == c("Dmau.90")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dmau90$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmau90<-max(colMeans(df_x))
Dmau90<- ggetho(dt[xmv(strain) == c("Dmau.90")], aes(x = t, y = activity/maxDmau90, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dmau91 <- ggetho(dt[xmv(strain) == c("Dmau.91")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dmau91$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmau91<-max(colMeans(df_x))
Dmau91<- ggetho(dt[xmv(strain) == c("Dmau.91")], aes(x = t, y = activity/maxDmau91, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

ggarrange(LZV72+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), LZV76+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),MD221+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),MD242+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),Dmau90+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dmau91+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),
          heights = c(3), ncol = 1, nrow = 6)

dt[, MA := ifelse(t %% hours(24) < hours(21), "R", "MA")]
dt[, S := ifelse(t %% hours(24) < hours(18), "R", "S")]

summary_dt <- 
  rejoin(dt[,
            .(
              # this is where the computation happens
              activity_MA = mean(activity[MA == "MA"]),
              activity_S = mean(activity[S == "S"]) 
            ),
            ,by=id])
summary_dt
order <-c("Dmel.LZV_L72", "Dmel.LZV_L76", "Dsim.MD221", "Dsim.MD242", "Dmau.90", "Dmau.91")
p<- ggplot(summary_dt, aes(x=factor(strain, level = order), y=activity_MA/max(activity_MA), fill=strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "", breaks = c(0, 0.5, 1), labels = c("", "", "")) + 
  theme_classic() + 
  scale_x_discrete(name = "", labels = c("","")) + 
  geom_hline(yintercept = 0.016508465, lwd = 2 , col = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(summary_dt[, activity_MA/max(activity_MA)], 
                     summary_dt[, strain], p.adjust = "bonf")


x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data/Extended photoperiod data/DeltaME_Dmau_afromelsim.csv", header = TRUE)
p<- ggplot(x, aes(x=Order, y=E12, fill=Order), inherit.aes = TRUE) +
  geom_boxplot(fill= rgb(66, 199, 244, 100, maxColorValue = 255), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 1),alpha=0.3) +
  scale_y_continuous(name= "Evening peak time 16:8h LD(hours)", limits = c(11,17), breaks = c(10, 12, 14, 16, 18, 20)) + 
  theme_classic() + 
  geom_hline(yintercept = 12, lwd = 2 , col = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

################################################################tropical Dmel/Dsim/Dmau 16:8h LD ANALYSIS##################################################################

metadata <- fread("metadata-168LDmau.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[replicate == "5" ])

ggetho(dt[xmv(strain) == c("Dmau.91")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho(colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), outline  = NA, l_duration = hours(16), height = .02, alpha = 1) +
  theme_classic() +
  facet_wrap(~region_id)

metadata <- fread("metadata-168LDmau.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])

LZV76 <- ggetho(dt[xmv(strain) == c("Dmel.LZV_L76")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-LZV76$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxLZV76<-max(colMeans(df_x))
LZV76<- ggetho(dt[xmv(strain) == c("Dmel.LZV_L76")], aes(x = t, y = activity/maxLZV76, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

LZV72 <- ggetho(dt[xmv(strain) == c("Dmel.LZV_L72")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-LZV72$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxLZV72<-max(colMeans(df_x))
LZV72<- ggetho(dt[xmv(strain) == c("Dmel.LZV_L72")], aes(x = t, y = activity/maxLZV72, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

A10 <- ggetho(dt[xmv(strain) == c("Dmel.LZV_A10")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-A10$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxA10<-max(colMeans(df_x))
A10<- ggetho(dt[xmv(strain) == c("Dmel.LZV_A10")], aes(x = t, y = activity/maxA10, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 


MD221 <- ggetho(dt[xmv(strain) == c("Dsim.MD221")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-MD221$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxMD221<-max(colMeans(df_x))
MD221<- ggetho(dt[xmv(strain) == c("Dsim.MD221")], aes(x = t, y = activity/maxMD221, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

MD242 <- ggetho(dt[xmv(strain) == c("Dsim.MD242")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-MD242$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxMD242<-max(colMeans(df_x))
MD242<- ggetho(dt[xmv(strain) == c("Dsim.MD242")], aes(x = t, y = activity/maxMD242, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dmau90 <- ggetho(dt[xmv(strain) == c("Dmau.90")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dmau90$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmau90<-max(colMeans(df_x))
Dmau90<- ggetho(dt[xmv(strain) == c("Dmau.90")], aes(x = t, y = activity/maxDmau90, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dmau91 <- ggetho(dt[xmv(strain) == c("Dmau.91")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dmau91$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmau91<-max(colMeans(df_x))
Dmau91<- ggetho(dt[xmv(strain) == c("Dmau.91")], aes(x = t, y = activity/maxDmau91, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

ggarrange(LZV72+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), LZV76+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),MD221+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),MD242+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),Dmau90+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dmau91+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),
          heights = c(3), ncol = 1, nrow = 6)

x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data/Extended photoperiod data/DeltaME_Dmau_afromelsim.csv", header = TRUE)
p<- ggplot(x, aes(x=Order, y=E168, fill=Order), inherit.aes = TRUE) +
  geom_boxplot(fill= rgb(66, 199, 244, 100, maxColorValue = 255), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 1),alpha=0.3) +
  scale_y_continuous(name= "Evening peak time 16:8h LD(hours)", limits = c(11,17), breaks = c(10, 12, 14, 16, 18, 20)) + 
  theme_classic() + 
  geom_hline(yintercept = 13, lwd = 2 , col = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(x$E168, x$strain, p.adjust.method = "none")

########################################Longevity Cumulative survival plots##############################################

x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Longevity_assay/rep1/Longevity_168hLD.csv")

plot(x$Day, x$Dmel_CSP12, type = "l", col = rgb(66, 199, 244, 100, maxColorValue = 255), lwd = 2,  xlab = "Day", ylab = "Cumulative survival probability", axes = FALSE, ylim = c(0,1), xlim = c(1,55))
points(x$Day, x$Dmel_CSP12, pch = 16, col = rgb(66, 199, 244, 255, maxColorValue = 255), cex = 1)
lines(x$Day, x$Dmel_CSP168, lty = 1, col = rgb(66, 199, 244, 100, maxColorValue = 255), lwd = 2)
points(x$Day, x$Dmel_CSP168, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 15, cex = 1)
lines(x$Day, x$Dsec_CSP12, lty = 1, col = rgb(232, 164, 35, 100, maxColorValue = 255), lwd = 2)
points(x$Day, x$Dsec_CSP12, pch = 16, col = rgb(232, 164, 35, 255, maxColorValue = 255), cex = 1)
lines(x$Day, x$Dsec_CSP168, lty = 1, col = rgb(232, 164, 35, 100, maxColorValue = 255), lwd = 2)
points(x$Day, x$Dsec_CSP168, col = rgb(232, 164, 35, 255, maxColorValue = 255), pch = 15, cex = 1)
axis(1, lty=1, at = c(0,10,20,30,40,50), labels = FALSE)
axis(2)


########################################Longevity logrank analaysis##############################################

library(survival)
x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Longevity_assay/rep1/pooled_cps.csv")

## compare Dsecpooled 12hLD to 168h LD   p= 3e-04 
Dsecpooled<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Longevity_assay/rep1/pooled_cps.csv")
Dsecpooledsurv<-Surv(Dsecpooled$Time_Dsec, Dsecpooled$Status_Dsec)
Dsecpooledfit<-survfit(Dsecpooledsurv ~ Dsecpooled$Treatment_Dsec)
plot(Dsecpooledfit, col = c("blue", "green"))
survdiff(Dsecpooledsurv ~ Dsecpooled$Treatment_Dsec)

## compare Dmelpooled 12hLD to 168h LD    p= 1e-05 
Dmelpooled<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Longevity_assay/rep1/pooled_cps.csv")
Dmelpooledsurv<-Surv(Dmelpooled$Time_Dsec, Dmelpooled$Status_Dmel)
Dmelpooledfit<-survfit(Dmelpooledsurv ~ Dmelpooled$Treatment_Dmel)
plot(Dmelpooledfit, col = c("blue", "green"))
survdiff(Dmelpooledsurv ~ Dmelpooled$Treatment_Dmel)

## compare DmelCS 12h LD to 168h LD rep1: p= 3e-05 rep2: p= 9e-05 
rep1pooled<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Longevity_assay/rep2/pooled_cps_rep2.csv")
rep1pooledsurv<-Surv(rep1pooled$Time_DmelCScomp, rep1pooled$Status_DmelCScomp)
rep1pooledfit<-survfit(rep1pooledsurv ~ rep1pooled$Treatment_DmelCScomp)
plot(rep1pooledfit, col = c("blue", "green"))
survdiff(rep1pooledsurv ~ rep1pooled$Treatment_DmelCScomp)

## compare DmelOR 12h LD to 168h LD rep1: p= 8e-04 rep2: p= 6e-05 
rep1pooled<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Longevity_assay/rep2/pooled_cps_rep2.csv")
rep1pooledsurv<-Surv(rep1pooled$Time_DmelORcomp, rep1pooled$Status_DmelORcomp)
rep1pooledfit<-survfit(rep1pooledsurv ~ rep1pooled$Treatment_DmelORcomp)
plot(rep1pooledfit, col = c("blue", "green"))
survdiff(rep1pooledsurv ~ rep1pooled$Treatment_DmelORcomp)

## compare Dsec07 12h LD to 168h LD rep1: p= 0.05  rep2: p= 0.008 
rep1pooled<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Longevity_assay/rep2/pooled_cps_rep2.csv")
rep1pooledsurv<-Surv(rep1pooled$Time_Dsec07comp, rep1pooled$Status_Dsec07comp)
rep1pooledfit<-survfit(rep1pooledsurv ~ rep1pooled$Treatment_Dsec07comp)
plot(rep1pooledfit, col = c("blue", "green"))
survdiff(rep1pooledsurv ~ rep1pooled$Treatment_Dsec07comp)

## compare Dsec28 12h LD to 168h LD rep1: p= 0.003 rep2: p= 0.02 
rep1pooled<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Longevity_assay/rep2/pooled_cps_rep2.csv")
rep1pooledsurv<-Surv(rep1pooled$Time_Dsec28comp, rep1pooled$Status_Dsec28comp)
rep1pooledfit<-survfit(rep1pooledsurv ~ rep1pooled$Treatment_Dsec28comp)
plot(rep1pooledfit, col = c("blue", "green"))
survdiff(rep1pooledsurv ~ rep1pooled$Treatment_Dsec28comp)

################################################################KO HYBRID ANALYSIS##################################################################

metadata <- fread("metadata_KODFs+PDF-1212LD.csv")				#read in metadata table --> change to .csv for your experiments
metadata <- fread("metadata_KODFs+PDF-168LD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])
dt <- load_dam(metadata[status == "OK" & strain == "Dsec28" |
                          status == "OK" & strain == "Dsec07" | 
                          status == "OK" & strain == "w1118" |
                          status == "OK" & strain == "Dsec28xw1118" |
                          status == "OK" & strain == "Dsec07xw1118"|
                          status == "OK" & strain == "Dsec07xpdf01(CS)" | 
                          status == "OK" & strain == "Dsec28xpdf01(CS)" | 
                          status == "OK" & strain == "w1118xPdf01(CS)" ], FUN = sleepr::sleep_dam_annotation)
summary(dt)

dt <- load_dam(metadata[status == "OK" & strain == "Dsec28" |
                             status == "OK" & strain == "Dsec07" | 
                             status == "OK" & strain == "w1118" |
                             status == "OK" & strain == "Dsec28xw1118" |
                             status == "OK" & strain == "Dsec07xw1118" |
                             status == "OK" & strain == "Dsec28xcry02" |
                             status == "OK" & strain == "Dsec07xcry02" |
                             status == "OK" & strain == "Dsec07xcyc01" |
                             status == "OK" & strain == "Dsec28xcyc01"|
                             status == "OK" & strain == "Dsec07xHr38^56" |                           
                             status == "OK" & strain == "Dsec28xHr38^56" |
                             status == "OK" & strain == "Dsec07xdfexel6011" |
                             status == "OK" & strain == "Dsec28xdfexel6011"] , FUN = sleepr::sleep_dam_annotation)   

#change strain name to check for dead flies, change status of dead flies in metadata file and re-run analysis

ggetho(dt[xmv(strain) == c("Dsec07")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho() +
  facet_wrap(~ region_id) +
  stat_ld_annotations(ld_colours =  c("lightyellow", "black"), outline  = NA, l_duration = hours(16), height = 1, alpha = 0.05)

#average activity by strain - line standardized by max activity         

Dsec07 <- ggetho(dt[xmv(strain) == "Dsec07"], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dsec07$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec07<-max(colMeans(df_x))
Dsec07 <- ggetho(dt[xmv(strain) == c("Dsec07")], aes(x = t, y = activity/maxDsec07, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dsec07xw1118 <- ggetho(dt[xmv(strain) == c("Dsec07xw1118")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dsec07xw1118$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec07xw1118<-max(colMeans(df_x))
Dsec07xw1118 <- ggetho(dt[xmv(strain) == c("Dsec07xw1118")], aes(x = t, y = activity/maxDsec07xw1118, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(102, 150, 0, 255, maxColorValue = 255), fill = rgb(102, 150, 0, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(102, 150, 0, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

Dsec28 <- ggetho(dt[xmv(strain) == c("Dsec28")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Dsec28$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec28<-max(colMeans(df_x))
Dsec28 <- ggetho(dt[xmv(strain) == c("Dsec28")], aes(x = t, y = activity/maxDsec28, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.452, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dsec28xw1118 <- ggetho(dt[xmv(strain) == c("Dsec28xw1118")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Dsec28xw1118$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec28xw1118<-max(colMeans(df_x))
Dsec28xw1118 <- ggetho(dt[xmv(strain) == c("Dsec28xw1118")], aes(x = t, y = activity/2.62, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(102, 150, 0, 255,  maxColorValue = 255), fill = rgb(102, 150, 0, 100,maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(102, 150, 0, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

w1118 <- ggetho(dt[xmv(strain) == c("w1118")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-w1118$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxw1118<-max(colMeans(df_x))
w1118 <- ggetho(dt[xmv(strain) == c("w1118")], aes(x = t, y = activity/maxw1118, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

ggarrange(Dsec07+ rremove("x.text")+rremove("xlab"), Dsec28+ rremove("x.text")+rremove("xlab"), Dsec07xw1118+ rremove("x.text")+rremove("xlab"),  Dsec28xw1118+ rremove("x.text")+rremove("xlab"), eya2+ rremove("x.text")+rremove("xlab"), eya228+ rremove("x.text")+rremove("xlab"), Clkout+ rremove("x.text")+rremove("xlab"), w1118+ rremove("x.text")+rremove("xlab"), 
          heights = c(3, 3, 3, 3, 4.1), ncol = 2, nrow = 4)

#Dsec07 crosses

Clkout <- ggetho(dt[xmv(strain) == c("Dsec07xClk[out]")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Clkout$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxClkout<-max(colMeans(df_x))
Clkout <- ggetho(dt[xmv(strain) == c("Dsec07xClk[out]")], aes(x = t, y = activity/maxClkout, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

CCHa <- ggetho(dt[xmv(strain) == c("Dsec07xCCHa-")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-CCHa$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxCCHa<-max(colMeans(df_x))
CCHa <- ggetho(dt[xmv(strain) == c("Dsec07xCCHa-")], aes(x = t, y = activity/maxCCHa, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Clk <- ggetho(dt[xmv(strain) == c("Dsec07xClk^JRK")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Clk$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxClk<-max(colMeans(df_x))
Clk <- ggetho(dt[xmv(strain) == c("Dsec07xClk^JRK")], aes(x = t, y = activity/maxClk, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

cwo <- ggetho(dt[xmv(strain) == c("Dsec07xcwo-")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-cwo$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxcwo<-max(colMeans(df_x))
cwo <- ggetho(dt[xmv(strain) == c("Dsec07xcwo-")], aes(x = t, y = activity/maxcwo, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

cyc <- ggetho(dt[xmv(strain) == c("Dsec07xcyc01")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-cyc$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxcyc<-max(colMeans(df_x))
cyc <- ggetho(dt[xmv(strain) == c("Dsec07xcyc01")], aes(x = t, y = activity/maxcyc, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

cry <- ggetho(dt[xmv(strain) == c("Dsec07xcry02")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-cry$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxcry<-max(colMeans(df_x))
cry <- ggetho(dt[xmv(strain) == c("Dsec07xcry02")], aes(x = t, y = activity/maxcry, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

DfED7853 <- ggetho(dt[xmv(strain) == c("Dsec07xDfED7853")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-DfED7853$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDfED7853<-max(colMeans(df_x))
DfED7853 <- ggetho(dt[xmv(strain) == c("Dsec07xDfED7853")], aes(x = t, y = activity/maxDfED7853, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

Hr38 <- ggetho(dt[xmv(strain) == c("Dsec07xHr38^56")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Hr38$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxHr38<-max(colMeans(df_x))
Hr38 <- ggetho(dt[xmv(strain) == c("Dsec07xHr38^56")], aes(x = t, y = activity/maxHr38, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

scro <- ggetho(dt[xmv(strain) == c("Dsec07xscro-")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-scro$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxscro<-max(colMeans(df_x))
scro <- ggetho(dt[xmv(strain) == c("Dsec07xscro-")], aes(x = t, y = activity/maxscro, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

Rh70 <- ggetho(dt[xmv(strain) == c("Dsec07xRh70")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Rh70$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxRh70<-max(colMeans(df_x))
Rh70 <- ggetho(dt[xmv(strain) == c("Dsec07xRh70")], aes(x = t, y = activity/maxRh70, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

Fer2 <- ggetho(dt[xmv(strain) == c("Dsec07xFer2-")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Fer2$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxFer2<-max(colMeans(df_x))
Fer2 <- ggetho(dt[xmv(strain) == c("Dsec07xFer2-")], aes(x = t, y = activity/maxFer2, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

pdp1 <- ggetho(dt[xmv(strain) == c("Dsec07xPDP1E")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-pdp1$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxpdp1<-max(colMeans(df_x))
pdp1 <- ggetho(dt[xmv(strain) == c("Dsec07xPDP1E")], aes(x = t, y = activity/maxpdp1, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

vri <- ggetho(dt[xmv(strain) == c("Dsec07xdfexel6011")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-vri$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxvri<-max(colMeans(df_x))
vri <- ggetho(dt[xmv(strain) == c("Dsec07xdfexel6011")], aes(x = t, y = activity/maxvri, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

ITP <- ggetho(dt[xmv(strain) == c("Dsec07xITP-")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-ITP$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxITP<-max(colMeans(df_x))
ITP <- ggetho(dt[xmv(strain) == c("Dsec07xITP-")], aes(x = t, y = activity/maxITP, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

Dsec07xpdf01CS <- ggetho(dt[xmv(strain) == c("Dsec07xpdf01(CS)")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dsec07xpdf01CS$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec07xpdf01CS<-max(colMeans(df_x))
Dsec07xpdf01CS <- ggetho(dt[xmv(strain) == c("Dsec07xpdf01(CS)")], aes(x = t, y = activity/maxDsec07xpdf01CS, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

#Morning peak 12hLD
ggarrange(w1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Dsec07xw1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),   CCHa+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),  Clkout+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),  cwo+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cyc+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cry + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Fer2+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Hr38+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          ITP+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),  DfED7853+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Dsec07xpdf01CS+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Rh70+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), scro+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),   vri+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          heights = c(3), nrow = 2, ncol = 8)

metadata <- fread("metadata_KODFs+PDF-1212LD.csv")				#read in metadata table --> change to .csv for your experiments
metadata <- fread("metadata_KODFs+PDF-168LD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK" & strain == "Dsec07xw1118" | 
                            status == "OK" & strain == "Dsec07xcry02" |
                            status == "OK" & strain == "Dsec07xscro-" |
                            status == "OK" & strain == "Dsec07xcwo-" |     
                            status == "OK" & strain == "Dsec07xcyc01" |
                            status == "OK" & strain == "Dsec07xClk(out)" |
                            status == "OK" & strain == "Dsec07xHr38^56" |                           
                            status == "OK" & strain == "Dsec07xFer2-" |                           
                            status == "OK" & strain == "Dsec07xRh70" |  
                            status == "OK" & strain == "Dsec07xdfexel6011" |
                            status == "OK" & strain == "Dsec07xDfED7853" |
                            status == "OK" & strain == "Dsec07xCCHa-" |
                            status == "OK" & strain == "Dsec07xITP-" |
                            status == "OK" & strain == "Dsec07xpdf01(CS)" |
                            status == "OK" & strain == "w1118"] , FUN = sleepr::sleep_dam_annotation)

dt[, MA := ifelse(t %% hours(24) < hours(21), "R", "MA")]
dt[, S := ifelse(t %% hours(24) < hours(18), "R", "S")]

summary_dt <- 
  rejoin(dt[,
            .(
              # this is where the computation happens
              activity_MA = mean(activity[MA == "MA"]),
              activity_S = mean(activity[S == "S"]) 
            ),
            ,by=id])
summary_dt$strain

order <- c("w1118", "Dsec07xw1118", "Dsec07xcyc01", "Dsec07xcry02",  "Dsec07xHr38^56", "Dsec07xdfexel6011", "Dsec07", "Dsec28xw1118", "Dsec28xcyc01", "Dsec28xcry02",  "Dsec28xHr38^56", "Dsec28xdfexel6011", "Dsec28")
order <- c("w1118", "Dsec07xw1118", "Dsec07xCCHa-", "Dsec07xClK(out)", "Dsec07xcwo-", "Dsec07xcyc01", "Dsec07xcry02", "Dsec07xFer2-", "Dsec07xHr38^56", "Dsec07xITP-", "Dsec07xDfED7853",   "Dsec07xpdf01(CS)", "Dsec07xRh70", "Dsec07xscro-",  "Dsec07xdfexel6011")

p<- ggplot(summary_dt, aes(x= factor(strain, level = order), y=activity_MA/max(activity_MA), fill=strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(102, 150, 0, 100, maxColorValue = 255),  rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(232, 164, 35, 255, maxColorValue = 255), rgb(102, 150, 0, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(232, 164, 35, 255, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1),alpha=0.3) +
  scale_y_continuous(name= "Increase in pre-dawn activity (%)", labels = scales::percent) + 
  theme_classic() + 
  scale_x_discrete(name = "") + 
 # geom_hline(yintercept = 0.016508465, lwd = 2 , col = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(summary_dt[, activity_MA/max(activity_MA)], 
                     summary_dt[, strain], p.adjust = "none")

#Evening peak plasticity 168hLD
ggarrange(Dsec07+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), CCHa + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Clkout + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cwo + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cyc+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cry + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), DfED7853+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          Hr38+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), scro+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Rh70+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Fer2+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), pdp1+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),  vri+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Dsec07xw1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), w1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          heights = c(3), ncol = 2, nrow = 7)

x<-read.csv("/Users/michaelshahandeh/Desktop/peakdetection/KO_DF_hybrids/DeltaME_Dsec07.csv")
order <- c("w1118", "Dsec07xw1118", "Dsec07xCCHa-", "Dsec07xClK(out)", "Dsec07xcwo-", "Dsec07xcyc01", "Dsec07xcry02", "Dsec07xFer2-", "Dsec07xHr38^56", "Dsec07xITP-", "Dsec07xDfED7853",   "Dsec07xpdf01(CS)", "Dsec07xRh70", "Dsec07xscro-",  "Dsec07xdfexel6011")
p<- ggplot(x, aes(x= factor(Strain, level = order), y=E168, fill=Strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(102, 150, 0, 100,  maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(232, 164, 35, 255, maxColorValue = 255), rgb(102, 150, 0, 100,  maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(232, 164, 35, 255, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 1),alpha=0.3) +
  scale_y_continuous(name= "", limits = c(10,16)) + 
  theme_classic() + 
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(x$E168, x$Strain, p.adjust.method = "none")

x<-read.csv("/Users/michaelshahandeh/Desktop/peakdetection/KO_DF_hybrids/DeltaME_Dsec07.csv")
pairwise.wilcox.test(x$E168, x$Strain, p.adjust = "none")

#w1118 crosses
dt <- load_dam(metadata[status == "OK"])

PDP1 <- ggetho(dt[xmv(strain) == c("w1118xPDP1E")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-PDP1$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxPDP1<-max(colMeans(df_x))
PDP1 <- ggetho(dt[xmv(strain) == c("w1118xPDP1E")], aes(x = t, y = activity/maxPDP1, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 100, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

cwo <- ggetho(dt[xmv(strain) == c("cwo-xw1118")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-cwo$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxcwo<-max(colMeans(df_x))
cwo <- ggetho(dt[xmv(strain) == c("cwo-xw1118")], aes(x = t, y = activity/maxcwo, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 100, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

ED7853 <- ggetho(dt[xmv(strain) == c("DfED7853xw1118")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-ED7853$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxED7853<-max(colMeans(df_x))
ED7853 <- ggetho(dt[xmv(strain) == c("DfED7853xw1118")], aes(x = t, y = activity/maxED7853, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 100, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

Hr38 <- ggetho(dt[xmv(strain) == c("Hr38^56xw1118")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Hr38$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxHr38<-max(colMeans(df_x))
Hr38 <- ggetho(dt[xmv(strain) == c("Hr38^56xw1118")], aes(x = t, y = activity/maxHr38, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 100, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

Clk <- ggetho(dt[xmv(strain) == c("Clk^JRKxw1118")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Clk$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxClk<-max(colMeans(df_x))
Clk <- ggetho(dt[xmv(strain) == c("Clk^JRKxw1118")], aes(x = t, y = activity/maxClk, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

vri <- ggetho(dt[xmv(strain) == c("DfExel6011xw1118")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-vri$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxvri<-max(colMeans(df_x))
vri <- ggetho(dt[xmv(strain) == c("DfExel6011xw1118")], aes(x = t, y = activity/maxvri, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 100, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

cyc <- ggetho(dt[xmv(strain) == c("cyc01xw1118")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-cyc$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxcyc<-max(colMeans(df_x))
cyc <- ggetho(dt[xmv(strain) == c("cyc01xw1118")], aes(x = t, y = activity/maxcyc, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 100, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

cry <- ggetho(dt[xmv(strain) == c("cry02xw1118")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-cry$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxcry<-max(colMeans(df_x))
cry <- ggetho(dt[xmv(strain) == c("cry02xw1118")], aes(x = t, y = activity/maxcry, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 100, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

pdf <- ggetho(dt[xmv(strain) == c("w1118xPdf01(CS)")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-pdf$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxpdf<-max(colMeans(df_x))
pdf <- ggetho(dt[xmv(strain) == c("w1118xPdf01(CS)")], aes(x = t, y = activity/maxpdf, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

ggarrange(w1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), pdf+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), 
          heights = c(3), ncol = 1, nrow = 2)

dt <- load_dam(metadata[  status == "OK" & strain == "Hr38^56xw1118" |
                            status == "OK" & strain == "Clk^JRKxw1118" | 
                            status == "OK" & strain == "dfexel6011xw1118" |
                            status == "OK" & strain == "cyc01xw1118" |
                            status == "OK" & strain == "cry02xw1118" |
                            status == "OK" & strain == "w1118"], FUN = sleepr::sleep_dam_annotation)

dt[, MA := ifelse(t %% hours(24) < hours(21), "R", "MA")]
dt[, S := ifelse(t %% hours(24) < hours(18), "R", "S")]

summary_dt <- 
  rejoin(dt[,
            .(
              # this is where the computation happens
              activity_MA = mean(activity[MA == "MA"]),
              activity_S = mean(activity[S == "S"]) 
            ),
            ,by=id])
summary_dt

order <- c("w1118", "Clk^JRKxw1118" , "cry02xw1118", "cyc01xw1118", "Hr38^56xw1118", "dfexel6011xw1118")
p<- ggplot(summary_dt, aes(x= factor(strain, level = order), y= factor(activity_MA, level = order), fill=strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1),alpha=0.3) +
  theme_classic() + 
  scale_x_discrete(name = "") + 
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(summary_dt[, activity_MA/max(activity_MA)], 
                     summary_dt[, strain], p.adjust = "none")

#168hLD

ggarrange(w1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Clk+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          heights = c(3), ncol = 1, nrow = 2)

x<-read.csv("/Users/michaelshahandeh/Desktop/peakdetection/KO_DF_hybrids/DeltaME_w1118.csv", header = TRUE)
order <- c("w1118", "Clk^JRKxw1118", "w1118xPdf01(CS)")
p<- ggplot(x, aes(x= factor(Strain, level = order), y=E168, fill=Strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 1),alpha=0.3) +
  scale_y_continuous(name= "E peak time (hours)", limits = c(10,16)) + 
  theme_classic() + 
  scale_x_discrete(name = "") + 
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(x$E168, x$Strain, p.adjust.method = "none")

#Dsec28 crosses

dt <- load_dam(metadata[status == "OK"])

ITP <- ggetho(dt[xmv(strain) == c("Dsec28xITP-")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-ITP$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxITP<-max(colMeans(df_x))
ITP <- ggetho(dt[xmv(strain) == c("Dsec28xITP-")], aes(x = t, y = activity/maxITP, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Clk <- ggetho(dt[xmv(strain) == c("Dsec28xClk^JRK")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Clk$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxClk<-max(colMeans(df_x))
Clk <- ggetho(dt[xmv(strain) == c("Dsec28xClk^JRK")], aes(x = t, y = activity/maxClk, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

cwo <- ggetho(dt[xmv(strain) == c("Dsec28xcwo-")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-cwo$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxcwo<-max(colMeans(df_x))
cwo <- ggetho(dt[xmv(strain) == c("Dsec28xcwo-")], aes(x = t, y = activity/maxcwo, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

cyc28 <- ggetho(dt[xmv(strain) == c("Dsec28xcyc01")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-cyc28$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxcyc28<-max(colMeans(df_x))
cyc28 <- ggetho(dt[xmv(strain) == c("Dsec28xcyc01")], aes(x = t, y = activity/maxcyc28, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

cry28 <- ggetho(dt[xmv(strain) == c("Dsec28xcry02")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-cry28$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxcry28<-max(colMeans(df_x))
cry28 <- ggetho(dt[xmv(strain) == c("Dsec28xcry02")], aes(x = t, y = activity/maxcry28, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

DfED7853 <- ggetho(dt[xmv(strain) == c("Dsec28xDfED7853")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-DfED7853$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDfED7853<-max(colMeans(df_x))
DfED7853 <- ggetho(dt[xmv(strain) == c("Dsec28xDfED7853")], aes(x = t, y = activity/maxDfED7853, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

Hr3828 <- ggetho(dt[xmv(strain) == c("Dsec28xHr38^56")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Hr3828$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxHr3828<-max(colMeans(df_x))
Hr3828 <- ggetho(dt[xmv(strain) == c("Dsec28xHr38^56")], aes(x = t, y = activity/maxHr3828, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

scro <- ggetho(dt[xmv(strain) == c("Dsec28xscro-")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-scro$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxscro<-max(colMeans(df_x))
scro <- ggetho(dt[xmv(strain) == c("Dsec28xscro-")], aes(x = t, y = activity/maxscro, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

Fer2 <- ggetho(dt[xmv(strain) == c("Dsec28xFer2-")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Fer2$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxFer2<-max(colMeans(df_x))
Fer2 <- ggetho(dt[xmv(strain) == c("Dsec28xFer2-")], aes(x = t, y = activity/maxFer2, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

pdp1 <- ggetho(dt[xmv(strain) == c("Dsec28xPDP1E")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-pdp1$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxpdp1<-max(colMeans(df_x))
pdp1 <- ggetho(dt[xmv(strain) == c("Dsec28xPDP1E")], aes(x = t, y = activity/maxpdp1, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

vri28 <- ggetho(dt[xmv(strain) == c("Dsec28xdfexel6011")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-vri28$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxvri28<-max(colMeans(df_x))
vri28 <- ggetho(dt[xmv(strain) == c("Dsec28xdfexel6011")], aes(x = t, y = activity/maxvri28, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 225, 198, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

CCHa <- ggetho(dt[xmv(strain) == c("Dsec28xCCHa-")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-CCHa$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxCCHa<-max(colMeans(df_x))
CCHa <- ggetho(dt[xmv(strain) == c("Dsec28xCCHa-")], aes(x = t, y = activity/maxCCHa, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 100, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none")  +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1))

Dsec28xpdf01CS <- ggetho(dt[xmv(strain) == c("Dsec28xpdf01(CS)")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dsec28xpdf01CS$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec28xpdf01CS<-max(colMeans(df_x))
Dsec28xpdf01CS <- ggetho(dt[xmv(strain) == c("Dsec28xpdf01(CS)")], aes(x = t, y = activity/maxDsec28xpdf01CS, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 225, 198, 255, maxColorValue = 255), fill = rgb(191, 225, 198, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

#12hLD
ggarrange(w1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Dsec28xw1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), CCHa+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Clk+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cwo + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cyc+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cry + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Fer2+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Hr38+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          ITP+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),  DfED7853+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Dsec07xpdf01CS+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), pdp1+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), scro+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),   vri+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          heights = c(3), nrow = 2, ncol = 8)

#12hLD both 07 and 28
ggarrange(w1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Dsec07xw1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cyc+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cry + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Hr38+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          vri+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),  Dsec07+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Dsec28xw1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Dsec28xw1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cyc28+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), cry28 + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Hr3828+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          vri28+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),  Dsec28+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          heights = c(3), nrow = 2, ncol = 7)


dt <- load_dam(metadata[ 
                            status == "OK" & strain == "Dsec07xw1118" | 
                            status == "OK" & strain == "Dsec07xcry02" |
                            status == "OK" & strain == "Dsec07" |
                            status == "OK" & strain == "Dsec07xcyc01" |
                            status == "OK" & strain == "Dsec07xHr38^56" |                           
                            status == "OK" & strain == "Dsec07xdfexel6011" |
                            status == "OK" & strain == "Dsec28xw1118" | 
                            status == "OK" & strain == "Dsec28xcry02" |
                            status == "OK" & strain == "Dsec28" |
                            status == "OK" & strain == "Dsec28xcyc01" |
                            status == "OK" & strain == "Dsec28xHr38^56" | 
                            status == "OK" & strain == "Dsec28xdfexel6011" |
                            status == "OK" & strain == "w1118"] , FUN = sleepr::sleep_dam_annotation)

dt[, MA := ifelse(t %% hours(24) < hours(21), "R", "MA")]
dt[, S := ifelse(t %% hours(24) < hours(18), "R", "S")]

summary_dt <- 
  rejoin(dt[,
            .(
              # this is where the computation happens
              activity_MA = mean(activity[MA == "MA"]),
              activity_S = mean(activity[S == "S"]) 
            ),
            ,by=id])
summary_dt[,strain]

order <- c("w1118", "Dsec07xw1118", "Dsec07xcyc01", "Dsec07xcry02", "Dsec07xHr38^56", "Dsec07xdfexel6011", "Dsec07", "Dsec28xw1118",  "Dsec28xcyc01", "Dsec28xcry02", "Dsec28xHr38^56", "Dsec28xdfexel6011", "Dsec28")
p<- ggplot(summary_dt, aes(x= factor(strain, level = order), y=activity_MA/max(activity_MA), fill=strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(102, 150, 0, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255),  rgb(102, 150, 0, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255),rgb(232, 164, 35, 100, maxColorValue = 255) ), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1),alpha=0.3) +
  scale_y_continuous(name= "Increase in pre-dawn activity (%)", labels = scales::percent) + 
  theme_classic() + 
  geom_hline(yintercept = 0.016508465, lwd = 2 , col = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  scale_x_discrete(name = "") + 
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)


dt <- load_dam(metadata[ 
  status == "OK" & strain == "Dsec28xw1118" | 
    status == "OK" & strain == "Dsec28xcry02" |
    status == "OK" & strain == "Dsec28xscro-" |
    status == "OK" & strain == "Dsec28xcwo-" |     
    status == "OK" & strain == "Dsec28xcyc01" |
    status == "OK" & strain == "Dsec28xClk^JRK" |
    status == "OK" & strain == "Dsec28xHr38^56" |                           
    status == "OK" & strain == "Dsec28xFer2-" |                           
    status == "OK" & strain == "Dsec28xPDP1E" |
    status == "OK" & strain == "Dsec28xdfexel6011" |
    status == "OK" & strain == "Dsec28xDfED7853" |
    status == "OK" & strain == "Dsec28xCCHa-" |
    status == "OK" & strain == "Dsec28xITP-" |
    status == "OK" & strain == "Dsec28xpdf01(CS)" |
    status == "OK" & strain == "w1118"] , FUN = sleepr::sleep_dam_annotation)

dt[, MA := ifelse(t %% hours(24) < hours(21), "R", "MA")]
dt[, S := ifelse(t %% hours(24) < hours(18), "R", "S")]

summary_dt <- 
  rejoin(dt[,
            .(
              # this is where the computation happens
              activity_MA = mean(activity[MA == "MA"]),
              activity_S = mean(activity[S == "S"]) 
            ),
            ,by=id])
summary_dt

order <- c("w1118", "Dsec28xw1118", "Dsec28xCCHa-", "Dsec28xClk^JRK", "Dsec28xcwo-", "Dsec28xcyc01", "Dsec28xcry02", "Dsec28xFer2-", "Dsec28xHr38^56", "Dsec28xITP-", "Dsec28xDfED7853",   "Dsec28xpdf01(CS)", "Dsec28xPDP1E", "Dsec28xscro-",  "Dsec28xdfexel6011")
p<- ggplot(summary_dt, aes(x= factor(strain, level = order), y=activity_MA/max(activity_MA), fill=strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(102, 150, 0, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.1),alpha=0.3) +
  scale_y_continuous(name= "Increase in pre-dawn activity (%)", labels = scales::percent) + 
  theme_classic() + 
  scale_x_discrete(name = "") + 
  geom_hline(yintercept = 0.016508465, lwd = 2 , col = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

#168hLD
ggarrange(w1118+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Dsec28xw1118 + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Clk + rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), 
          Hr38+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), scro+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"), Fer2+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),  vri+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),  Dsec28+ rremove("x.text")+rremove("xlab") + rremove("y.text")+rremove("ylab"),
          heights = c(3), ncol = 2, nrow = 6)

x<-read.csv("/Users/michaelshahandeh/Desktop/peakdetection/KO_DF_hybrids/DeltaME_Dsec28.csv", header = TRUE)
order <- c("w1118", "Dsec28xw1118", "Dsec28xCCHa-", "Dsec28xClk^JRK", "Dsec28xcwo-", "Dsec28xcyc01", "Dsec28xcry02", "Dsec28xFer2-", "Dsec28xHr38^56", "Dsec28xITP-", "Dsec28xDfED7853", "Dsec28xpdf01(CS)", "Dsec28xPDP1E", "Dsec28xRh70",  "Dsec28xscro-",  "Dsec28xdfexel6011")
p<- ggplot(x, aes(x= factor(Strain, level = order), y=E168, fill=Strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(102, 150, 0, 100,  maxColorValue = 255), rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255),rgb(191, 225, 198, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 1),alpha=0.3) +
  scale_y_continuous(name= "Increase in pre-dawn activity (%)", labels = scales::percent) + 
  theme_classic() + 
  geom_hline(yintercept = 13, lwd = 2 , col = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  scale_x_discrete(name = "") + 
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

x<-read.csv("/Users/michaelshahandeh/Desktop/peakdetection/KO_DF_hybrids/DeltaME_Dsec28.csv")
pairwise.wilcox.test(x$E168, x$Strain, p.adjust = "none")

###################################################EnricoCS comparison#############################################                        

metadata <- fread("metadata-EnricoCS_16-8.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])						#load data into behavior table
summary(dt)

DmelCS <- ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-DmelCS$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmelCS<-max(colMeans(df_x))
DmelCS<- ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, y = activity/maxDmelCS, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

EnricoCS <- ggetho(dt[xmv(strain) == c("Enrico.CS")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-EnricoCS$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxEnricoCS<-max(colMeans(df_x))
EnricoCS<- ggetho(dt[xmv(strain) == c("Enrico.CS")], aes(x = t, y = activity/maxEnricoCS, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

ggarrange(DmelCS+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), EnricoCS+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),
          heights = c(3, 3), ncol = 1, nrow = 2)

dt[, MA := ifelse(t %% hours(24) < hours(21), "R", "MA")]
dt[, S := ifelse(t %% hours(24) < hours(18), "R", "S")]

summary_dt <- 
  rejoin(dt[,
            .(
              # this is where the computation happens
              activity_MA = mean(activity[MA == "MA"]),
              activity_S = mean(activity[S == "S"]) 
            ),
            ,by=id])
summary_dt

p<- ggplot(summary_dt, aes(x=strain, y=activity_MA/max(activity_MA), fill=strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "") + 
  theme_classic() + 
  scale_x_discrete(name = "", labels = c("","","","")) + 
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(summary_dt[, activity_MA/max(activity_MA)], 
                     summary_dt[, strain], p.adjust = "bonf")

metadata <- fread("metadata-EnricoCS_16-8.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])						#load data into behavior table
summary(dt)

EnricoCS <- ggetho(dt[xmv(strain) == c("Enrico.CS")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-EnricoCS$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxEnricoCS<-max(colMeans(df_x))
EnricoCS<- ggetho(dt[xmv(strain) == c("Enrico.CS")], aes(x = t, y = activity/maxEnricoCS, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

metadata <- fread("metadata_16-8LD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])						#load data into behavior table
summary(dt)

DmelCS <- ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-DmelCS$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmelCS<-max(colMeans(df_x))
DmelCS<- ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, y = activity/maxDmelCS, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

ggarrange(DmelCS+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), EnricoCS+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),
          heights = c(3, 3), ncol = 1, nrow = 2)

x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data/Enrico.CS/DeltaME_EnricoCS.csv")

p<- ggplot(x, aes(x=Strain, y=E1212, fill=Strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "Evening peak time 16:8h LD(hours)", limits = c(10,18), breaks = c(10, 12, 14, 16, 18)) + 
  theme_classic() + 
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

wilcox.test(x$E168~ x$Strain)

##########################################Pdf RNAi ################################################

metadata <- fread("metadata_PDFRNAi_12hLD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[replicate == "9"])						#load data into behavior table
summary(dt)

ggetho(dt[xmv(strain) == c("BT411-3")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho() +
  facet_wrap(~ region_id) +
  stat_ld_annotations(ld_colours =  c("lightyellow", "black"), outline  = NA, l_duration = hours(16), height = 1, alpha = 0.05)


metadata <- fread("metadata_PDFRNAi_168hLD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])						#load data into behavior table
summary(dt)

BT4103 <- ggetho(dt[xmv(strain) == c("BT410-3")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT4103$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT4103<-max(colMeans(df_x))
BT4103 <- ggetho(dt[xmv(strain) == c("BT410-3")], aes(x = t, y = activity/maxBT4103, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

BT4113 <- ggetho(dt[xmv(strain) == c("BT411-3")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT4113$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT4113<-max(colMeans(df_x))
BT4113 <- ggetho(dt[xmv(strain) == c("BT411-3")], aes(x = t, y = activity/maxBT4113, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

BT4103xRNAi <- ggetho(dt[xmv(strain) == c("BT410-3xRNAi")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT4103xRNAi$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT4103xRNAi<-max(colMeans(df_x))
BT4103xRNAi <- ggetho(dt[xmv(strain) == c("BT410-3xRNAi")], aes(x = t, y = activity/maxBT4103xRNAi, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 191, 191, 255, maxColorValue = 255), fill = rgb(0, 0, 0, 0, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 191, 191, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

BT4113xRNAi <- ggetho(dt[xmv(strain) == c("BT411-3xRNAi")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT4113xRNAi$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT4113xRNAi<-max(colMeans(df_x))
BT4113xRNAi <- ggetho(dt[xmv(strain) == c("BT411-3xRNAi")], aes(x = t, y = activity/maxBT4113xRNAi, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 191, 191, 255, maxColorValue = 255), fill = rgb(191, 191, 191, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 191, 191, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

w1118xRNAi <- ggetho(dt[xmv(strain) == c("w1118xRNAi")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-w1118xRNAi$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxw1118xRNAi<-max(colMeans(df_x))
w1118xRNAi <- ggetho(dt[xmv(strain) == c("w1118xRNAi")], aes(x = t, y = activity/maxw1118xRNAi, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

ggarrange(BT4103 + rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), BT4103xRNAi+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), BT4113 + rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), BT4113xRNAi + rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), w1118xRNAi + rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          heights = 3, ncol = 1, nrow = 5)

x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM_data/PDF_RNAi/DeltaME_RNAi.csv", header = TRUE)

p<- ggplot(x, aes(x=Strain, y=E168, fill=Strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(0, 0, 0, 0, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(191, 191, 191, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = .5),alpha=0.3) +
  scale_y_continuous(name= "", limits = c(10,16)) + 
  theme_classic() + 
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(x$E1212, 
                     x$Strain, p.adjust = "none")

metadata <- fread("metadata_PDFRNAi_12hLD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[replicate != "9" & status == "OK"])						#load data into behavior table
dt <- load_dam(metadata[strain == "BT410-3" & status == "OK" | strain == "BT411-3" & status == "OK" | strain == "BT410-3xRNAi" & status == "OK" | strain == "BT411-3xRNAi" & status == "OK" | strain == "w1118xRNAi" & status == "OK"])						#load data into behavior table
summary(dt)

BT4103 <- ggetho(dt[xmv(strain) == c("BT410-3")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT4103$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT4103<-max(colMeans(df_x))
BT4103 <- ggetho(dt[xmv(strain) == c("BT410-3")], aes(x = t, y = activity/maxBT4103, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

BT4113 <- ggetho(dt[xmv(strain) == c("BT411-3")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT4113$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT4113<-max(colMeans(df_x))
BT4113 <- ggetho(dt[xmv(strain) == c("BT411-3")], aes(x = t, y = activity/maxBT4113, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

BT4103xRNAi <- ggetho(dt[xmv(strain) == c("BT410-3xRNAi")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT4103xRNAi$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT4103xRNAi<-max(colMeans(df_x))
BT4103xRNAi <- ggetho(dt[xmv(strain) == c("BT410-3xRNAi")], aes(x = t, y = activity/maxBT4103xRNAi, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 191, 191, 255, maxColorValue = 255), fill = rgb(0, 0, 0, 0, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 191, 191, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

BT4113xRNAi <- ggetho(dt[xmv(strain) == c("BT411-3xRNAi")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT4113xRNAi$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT4113xRNAi<-max(colMeans(df_x))
BT4113xRNAi <- ggetho(dt[xmv(strain) == c("BT411-3xRNAi")], aes(x = t, y = activity/maxBT4113xRNAi, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(191, 191, 191, 255, maxColorValue = 255), fill = rgb(191, 191, 191, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(191, 191, 191, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

w1118xRNAi <- ggetho(dt[xmv(strain) == c("w1118xRNAi")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-w1118xRNAi$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxw1118xRNAi<-max(colMeans(df_x))
w1118xRNAi <- ggetho(dt[xmv(strain) == c("w1118xRNAi")], aes(x = t, y = activity/maxw1118xRNAi, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.025, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

ggarrange(BT4103 + rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), BT4103xRNAi+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), BT4113 + rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), BT4113xRNAi + rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), w1118xRNAi + rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          heights = 3, ncol = 1, nrow = 5)

dt[, MA := ifelse(t %% hours(24) < hours(21), "R", "MA")]
dt[, S := ifelse(t %% hours(24) < hours(18), "R", "S")]

summary_dt <- 
  rejoin(dt[,
            .(
              # this is where the computation happens
              activity_MA = mean(activity[MA == "MA"]),
              activity_S = mean(activity[S == "S"]) 
            ),
            ,by=id])
summary_dt

p<- ggplot(summary_dt, aes(x=strain, y=activity_MA/max(activity_MA), fill=strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(0, 0, 0, 0, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(191, 191, 191, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.3),alpha=0.3) +
  scale_y_continuous(name= "") + 
  theme_classic() + 
  scale_x_discrete(name = "", labels = c("", "", "", "", "", "", "", "")) + 
  #geom_hline(yintercept = 0) +
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)


pairwise.wilcox.test(summary_dt[, activity_MA/max(activity_MA)], 
                     summary_dt[, strain], p.adjust = "none")

##########################################Pdf rescue data##########################################

metadata <- fread("metadata_PDFswap-12hLD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[replicate == "pdf01-2" & status == "OK"])						#load data into behavior table
summary(dt)

#change strain name to check for dead flies, change status of dead flies in metadata file and re-run analysis

ggetho(dt[xmv(strain) == c("pdf01(CS)")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho() +
  facet_wrap(~ region_id) +
  stat_ld_annotations(ld_colours =  c("lightyellow", "black"), outline  = NA, l_duration = hours(12), height = 1, alpha = 0.05)

metadata <- fread("metadata_PDFswap-168hLD.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])						#load data into behavior table
summary(dt)


BT421.2pdf01 <- ggetho(dt[xmv(strain) == c("BT421.2-pdf01")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT421.2pdf01$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT421.2pdf01<-max(colMeans(df_x))
BT421.2pdf01<- ggetho(dt[xmv(strain) == c("BT421.2-pdf01")], aes(x = t, y = activity/maxBT421.2pdf01, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

BT421.3pdf01 <- ggetho(dt[xmv(strain) == c("BT421.3-pdf01")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT421.3pdf01$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT421.3pdf01<-max(colMeans(df_x))
BT421.3pdf01<- ggetho(dt[xmv(strain) == c("BT421.3-pdf01")], aes(x = t, y = activity/maxBT421.3pdf01, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

BT420.3pdf01 <- ggetho(dt[xmv(strain) == c("BT420.3-pdf01")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT420.3pdf01$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT420.3pdf01<-max(colMeans(df_x))
BT420.3pdf01<- ggetho(dt[xmv(strain) == c("BT420.3-pdf01")], aes(x = t, y = activity/maxBT420.3pdf01, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

BT420.5pdf01 <- ggetho(dt[xmv(strain) == c("BT420.5-pdf01")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-BT420.5pdf01$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxBT420.5pdf01<-max(colMeans(df_x))
BT420.5pdf01<- ggetho(dt[xmv(strain) == c("BT420.5-pdf01")], aes(x = t, y = activity/maxBT420.5pdf01, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

DmelrescueattP2 <- ggetho(dt[xmv(species) == c("DmelrescueattP2")], aes(x = t, y = activity, colour = species), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-DmelrescueattP2$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmelrescueattP2<-max(colMeans(df_x))
DmelrescueattP2<- ggetho(dt[xmv(species) == c("DmelrescueattP2")], aes(x = t, y = activity/maxDmelrescueattP2, colour = species), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

DsecrescueattP2 <- ggetho(dt[xmv(species) == c("DsecrescueattP2")], aes(x = t, y = activity, colour = species), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-DsecrescueattP2$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsecrescueattP2<-max(colMeans(df_x))
DsecrescueattP2<- ggetho(dt[xmv(species) == c("DsecrescueattP2")], aes(x = t, y = activity/maxDsecrescueattP2, colour = species), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

pdf01 <- ggetho(dt[xmv(strain) == c("pdf01(CS)")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-pdf01$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxpdf01<-max(colMeans(df_x))
pdf01<- ggetho(dt[xmv(strain) == c("pdf01(CS)")], aes(x = t, y = activity/maxpdf01, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(190, 190, 190, 255, maxColorValue = 255), fill = rgb(190, 190, 190, 100, maxColorValue = 255)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(190, 190, 190, 255, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(16)) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

ggarrange(BT420.3pdf01+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), BT421.2pdf01+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), BT420.5pdf01+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          BT421.3pdf01+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), DmelrescueattP2+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), DsecrescueattP2+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),
          ncol = 2, nrow = 3)

ggarrange(pdf01+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), DmelrescueattP2+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), DsecrescueattP2+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),
          ncol = 1, nrow = 3)

dt <- load_dam(metadata[ 
  status == "OK" & strain == "BT420.3-pdf01" | 
    status == "OK" & strain == "BT420.5-pdf01" |
    status == "OK" & strain == "BT421.2-pdf01" |
    status == "OK" & strain == "BT421.3-pdf01" |     
    status == "OK" & strain == "pdf01(CS)" ], FUN = sleepr::sleep_dam_annotation)

dt[, MA := ifelse(t %% hours(24) < hours(21), "R", "MA")]
dt[, S := ifelse(t %% hours(24) < hours(18), "R", "S")]

summary_dt <- 
  rejoin(dt[,
            .(
              # this is where the computation happens
              activity_MA = mean(activity[MA == "MA"]),
              activity_S = mean(activity[S == "S"]) 
            ),
            ,by=id])
summary_dt

p<- ggplot(summary_dt, aes(x=species, y=activity_MA/max(activity_MA), fill=species), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "") + 
  theme_classic() + 
  scale_x_discrete(name = "") + 
  theme(legend.position = "none")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(summary_dt[, activity_MA/max(activity_MA)], 
                     summary_dt[, strain], p.adjust = "bonf")

p$data$strain
p$data$activity_MA
x<-cbind(p$data$strain, p$data$activity_MA)
write.csv(x, "~/Desktop/rescue_MA.csv")

x<-read.csv("~/Desktop/rescue_MA.csv")
p<- ggplot(x, aes(x=Order, y=MA/max(MA), fill=Order), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(190, 190, 190, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.6),alpha=0.3) +
  scale_y_continuous(name= "") + 
  theme_classic() + 
  scale_x_discrete(name = "") + 
  theme(legend.position = "none")
ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(x$MA/max(x$MA), 
                     x$Strain, p.adjust = "bonf")

x<-read.csv("~/Desktop/Epeak_attP2pdf01_rescue.csv")

p<- ggplot(x, aes(x=Order, y=E1212, fill=Order), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(190, 190, 190, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.8),alpha=0.3) +
  scale_y_continuous(name= "", limits = c(10,16)) + 
  theme_classic() + 
  scale_x_discrete(name = "") + 
  theme(legend.position = "none")
ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)
p
pairwise.wilcox.test(x$E1212, x$genotype, p.adjust.method = "bonferroni" )
pairwise.wilcox.test(x$E168, x$genotype, p.adjust.method = "bonferroni" )

pairwise.wilcox.test(x$E1212, x$rescue, p.adjust.method = "bonferroni" )
pairwise.wilcox.test(x$E168, x$rescue, p.adjust.method = "bonferroni" )

##########################################M peak 12h LD Pdf immuno ##########################################

#Points plus lines connecting means
mel_means<-c(4.30406296, 4.66706875, 3.818580975, 4.74091806)
sec07_means<-c(2.03241638, 3.36690474, 4.9458239, 5.7924907)

x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/Immunostains/PDF-morning peak time points/fixed at time points/PDF sLNv immunostains/Threshold/s-LNV_terminals_PDF_CS.csv")
y<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/Immunostains/PDF-morning peak time points/fixed at time points/PDF sLNv immunostains/Threshold/s-LNV_terminals_PDF_07.csv")
xaxis<-c(1,2,3,4)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$AU~x$Order ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 8), xlim = c(0.8, 4.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$AU~y$Order ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 8), xlim = c(0.5, 4.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec07_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4), labels = c("6", "8", "10", "12"),  col.axis = "white", col.lab = "white")

#cell bodies
x<-read.csv("/Users/michaelshahandeh/Documents/Circadian_rhythm/Immunostains/PDF-morning peak time points/fixed at time points/Threshold/s-LNV_cellbodies_re_PDF.csv")
xaxis<-c(1,2,3,4)
par(mar = c(1,1,1,1), bty = "l")
mel_means<-c(19.912809, 18.50518, 18.86567025, 16.96022)
sec07_means<-c(20.0844715, 16.5138985, 17.406542, 9.029429)
stripchart(x$AU1/5~x$Order ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 5), xlim = c(0.8, 4.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(x$AU3/5~x$Order ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 5), xlim = c(0.5, 4.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means/5, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec07_means/5, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4), labels = c("6", "8", "10", "12"),  col.axis = "white", col.lab = "white")

##########################################E peak 12h LD Pdf immuno ##########################################

mel_means<-c(0.3212139, 1.0624952, 1.834007, 1.8531314, 2.04104875)
sec07_means<-c(1.807869475, 1.98941995, 1.5335451, 1.259125, 1.9402565)

x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Pdf_immunos/Dmel_Dsec_E_peak_12hLD/CS_Epeak_12hLD.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Pdf_immunos/Dmel_Dsec_E_peak_12hLD/07_Epeak_12hLD.csv")
xaxis<-c(1,2,3,4,5)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$AU~x$Order ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 3), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$AU~y$Order ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 3), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec07_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4, 5), labels = c("6", "8", "10", "12", "14"),  col.axis = "white", col.lab = "white")
#legend(1, 4, col = c(rgb(66, 199, 244, 255, maxColorValue = 255)), pch = c(16), c(expression(italic(""))), bty = "n", pt.cex = 2, xpd = TRUE )
#legend(1, 3.4, col = c(rgb(232, 164, 35, 255, maxColorValue = 255)), pch = c(16), c(expression(italic(""))), bty = "n", pt.cex = 2, xpd = TRUE )

##########################################E peak 16:8h LD Pdf immuno ##########################################

mel_means<-c(7.744539, 19.3977485, 4.380659, 5.122052, 1.831151)
sec07_means<-c(22.1973185, 18.2687765, 15.063515, 11.8607155, 2.219081)

x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Pdf_immunos/Dmel_Dsec_E_peak_168hLD/CS_Epeak_168hLD.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/Pdf_immunos/Dmel_Dsec_E_peak_168hLD/07_Epeak_168hLD.csv")
xaxis<-c(1,2,3,4,5)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$AU~x$Order ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 30), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$AU~y$Order ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 30), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec07_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4, 5), labels = c("6", "8", "10", "12", "14"),  col.axis = "white", col.lab = "white")
#legend(1, 4, col = c(rgb(66, 199, 244, 255, maxColorValue = 255)), pch = c(16), c(expression(italic(""))), bty = "n", pt.cex = 2, xpd = TRUE )
#legend(1, 3.4, col = c(rgb(232, 164, 35, 255, maxColorValue = 255)), pch = c(16), c(expression(italic(""))), bty = "n", pt.cex = 2, xpd = TRUE )

mel_means<-c(1.467532, 19.2912874, 3.1951675, 2.2569248, 6.4599372)
sec07_means<-c(50.7359442, 64.9477357, 42.0004419, 72.3015658, 0.9263199)

x<-read.csv("/Users/michaelshahandeh/Desktop/Pdf_immunos/Dmel_Dsec_E_peak_168hLD/CS_Epeak_168hLD_dendrites.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/Pdf_immunos/Dmel_Dsec_E_peak_168hLD/07_Epeak_168hLD_dendrites.csv")
xaxis<-c(1,2,3,4,5)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$AU~x$Order ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 125), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$AU~y$Order ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 125), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec07_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4, 5), labels = c("6", "8", "10", "12", "14"),  col.axis = "white", col.lab = "white")
#legend(1, 4, col = c(rgb(66, 199, 244, 255, maxColorValue = 255)), pch = c(16), c(expression(italic(""))), bty = "n", pt.cex = 2, xpd = TRUE )
#legend(1, 3.4, col = c(rgb(232, 164, 35, 255, maxColorValue = 255)), pch = c(16), c(expression(italic(""))), bty = "n", pt.cex = 2, xpd = TRUE )

##########################################M peak 12h LD Pdf smFISH ##########################################

x<-read.csv("/Users/michaelshahandeh/Desktop/smFISH/smRNA_FISH_sLNVs_ZT20-ZT4/Mpeak_RNAspots_DmelCS.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/smFISH/smRNA_FISH_sLNVs_ZT20-ZT4/Mpeak_RNAspots_Dsec07.csv")
mel_means<-c(14.1, 21, 22.2, 23)
sec07_means<-c(2.8, 8.8, 7.625, 17.3)

xaxis<-c(1,2,3,4)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$Spots~x$Order ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 30), xlim = c(0.8,4.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$Spots~y$Order ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 30), xlim = c(0.8,4.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec07_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4), labels = c("20", "22", "0", "2"),  col.axis = "white", col.lab = "white")


##########################################E peak 12h LD Pdf smFISH ##########################################

mel_means<-c(40, 37.5, 44, 43, 15)
sec07_means<-c(51, 53, 45, 52, 59.5)

#Points plus lines connecting means
x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/smFISH/63xsmFISH_images/smFISH_RS_12hLD_Epeak_DmelCS.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/smFISH/63xsmFISH_images/smFISH_RS_12hLD_Epeak_Dsec07.csv")
xaxis<-c(1,2,3,4,5)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$LRavg~x$Order ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 80), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$LRavg~y$Order ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 80), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec07_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4, 5), labels = c("6", "8", "10", "12", "14"),  col.axis = "white", col.lab = "white")
#legend(1, 118, col = c(rgb(66, 199, 244, 255, maxColorValue = 255)), pch = c(16), c(expression(italic(""))), bty = "n", pt.cex = 2, xpd = TRUE )
#legend(1, 110, col = c(rgb(232, 164, 35, 255, maxColorValue = 255)), pch = c(16), c(expression(italic(""))), bty = "n", pt.cex = 2, xpd = TRUE )

##########################################E peak 168h LD Pdf smFISH ##########################################

mel_means<-c(39, 35, 43.5, 47.5, 54.5)
sec07_means<-c(45, 41, 47.5, 43, 44)

x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/smFISH/63xsmFISH_images/smFISH_RS_168hLD_Epeak_DmelCS.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/smFISH/63xsmFISH_images/smFISH_RS_168hLD_Epeak_Dsec07.csv")
xaxis<-c(1,2,3,4,5)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$AvgLR~x$Order ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 70), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$AvgLR~y$Order ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 70), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec07_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4, 5), labels = c("6", "8", "10", "12", "14"),  col.axis = "white", col.lab = "white")
#legend(1, 78, col = c(rgb(66, 199, 244, 255, maxColorValue = 255)), pch = c(16), c(expression(italic(""))), bty = "n", pt.cex = 2, xpd = TRUE )
#legend(1, 70, col = c(rgb(232, 164, 35, 255, maxColorValue = 255)), pch = c(16), c(expression(italic(""))), bty = "n", pt.cex = 2, xpd = TRUE )


##########################################M peak 12h LD transcriptional reporters##########################################

x<-read.csv("/Users/michaelshahandeh/Desktop/Pdf transc. reporter immunos/BT412:413_Mpeak 12hLD/re_projections.csv")
mel_means<-c(36.1274218, 43.825737, 34.478029, 28.3213493)
sec_means<-c(3.0530163, 18.2006696, 10.6542492, 9.2224295)
xaxis<-c(1,2,3,4)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$AU1~x$ZT ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 60), xlim = c(0.8, 4.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(x$AU2~x$ZT ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 60), xlim = c(0.8, 4.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255), lty = 2)
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4), labels = c("-4", "-2", "0", "2"),  col.axis = "white", col.lab = "white")

##########################################E peak 12h LD transcriptional reporters##########################################
mel_means<-c(2.823852113, 0.8920811, 2.75466508, 3.07777854, 1.77719945)
sec_means<-c(1.37580293, 0.90809911, 1.02829739, 1.36068639, 1.18932182)
Order<-c(-4, -2, 0, 2, 4)
melSE<-c(0.246180392, 0.064727613, 0.20608038, 0.094658973, 0.18589523)
secSE<-c(0.055591853, 0.072917192, 0.112160109, 0.142319233, 0.095572217)

#Points plus lines connecting means
x<-read.csv("/Users/michaelshahandeh/Desktop/Pdf transc. reporter immunos/BT412:413 E peak 12hLD/BT412_Epeak_12hLD.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/Pdf transc. reporter immunos/BT412:413 E peak 12hLD/BT413_Epeak_12hLD.csv")
xaxis<-c(1,2,3,4,5)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$AU~x$Order ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 4), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$AU~y$Order ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 4), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255), lty = 2)
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4, 5), labels = c("6", "8", "10", "12", "14"),  col.axis = "white", col.lab = "white")

##########################################E peak 16:8h LD transcriptional reporters##########################################

mel_means<-c(4.17871374, 3.65880354, 4.00978433, 1.80095522, 3.14484452)
sec_means<-c(2.80788305, 2.3021594, 1.97167121, 1.55092419, 1.82343425)
Order<-c(-4, -2, 0, 2, 4)
melSE<-c(0.11576698, 0.240396145, 0.205912497, 0.359412079, 0.281851879)
secSE<-c(0.142710656, 0.080008729, 0.163605159, 0.128367273, 0.117209398)


#Points plus lines connecting means
x<-read.csv("/Users/michaelshahandeh/Desktop/Pdf transc. reporter immunos/BT412:413 E peak 168hLD/BT412_Epeak_168hLD.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/Pdf transc. reporter immunos/BT412:413 E peak 168hLD/BT413_Epeak_168hLD.csv")
xaxis<-c(1,2,3,4,5)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$AU~x$Order ,pch = 16, cex = 3, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 5.2), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$AU~y$Order ,pch = 16, cex = 3, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 5.2), xlim = c(0.8, 5.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255), lty = 2)
axis(2, las = 2,  col.axis = "white", col.lab = "white")
axis(1, at = c(1, 2, 3, 4, 5), labels = c("6", "8", "10", "12", "14"),  col.axis = "white", col.lab = "white")

##############################################Fecundity Assay######################################################
library(boot)
cop2hr<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/fecundity assay/2hr_bootsraps.csv")
cop3d<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/fecundity assay/3d_bootsraps.csv")
mean.fun <- function(dat, idx) mean(dat[idx], na.rm = TRUE)

boot.out.cop<- boot(data =cop2hr$Dsec16, statistic = mean.fun, R = 1000)
boot_cop<-boot.out.cop$t
plot(boot.out.cop)
boot.ci(boot.out.cop)

#2 hr barplots strains
par(bty = "l", mar = c(0.5,0.75,0.5, 0.5))
xaxis<-c(1,2,4,5,6,7,8)
yaxis<-c(96.15,100,92,96,51.61,19.23,54.84,29.03,50,34.62,51.61,16.67)
barplot(yaxis, space = c(0,0,1,0,1,0,1,0,1,0,1,0), col = c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), ylim = c(0,101), axes = FALSE, ann = FALSE)
axis(2, lty = 1, labels = FALSE)
segments(0.5,88.46,0.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(1.5,100,1.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(3.5,80.77,3.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(4.5,88.00,4.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(6.5,32.26,6.5,70.97, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(7.5,3.85,7.5,34.62, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(9.5,38.71,9.5,70.97, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(10.5,12.90,10.5,45.16, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(12.5,29.17,12.5,69.57, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(13.5,15.38,13.5,53.85, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(15.5,32.34,15.5,67.74, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(16.5,3.70,16.5,32.00, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)

#3 d barplots strains
par(bty = "l", mar = c(0.5,0.75,0.5, 0.5))
xaxis<-c(1,2,4,5,6,7,8)
yaxis<-c(96.00, 88.46, 100.00, 96.15, 92.31, 45.83, 72.00, 56.00, 96.00, 84.62, 87.50, 95.45)
barplot(yaxis, space = c(0,0,1,0,1,0,1,0,1,0,1,0), col = c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), ylim = c(0,101), axes = FALSE, ann = FALSE)
axis(2, lty = 1, labels = FALSE)
segments(0.5,88.46,0.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(1.5,76.92,1.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(3.5,100,3.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(4.5,88.46,4.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(6.5,80.77,6.5,100.00, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(7.5,27.27,7.5,65.21, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(9.5,88.46,9.5,53.85, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(10.5,34.81,10.5,76, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(12.5,87.50,12.5,100, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(13.5,69.23,13.5,96.15, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(15.5,72.73,15.5,100, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(16.5,85.71,16.5,100, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)

#2 hr barplots species
par(bty = "l", mar = c(0.5,0.75,0.5, 0.5))
xaxis<-c(1,2,3,4)
yaxis<-c(94.12,98.04,52.14,25.23)
barplot(yaxis, space = c(0,0,1,0), col = c( rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), ylim = c(0,101), axes = FALSE, ann = FALSE)
segments(0.5,86.36,0.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(1.5,93.34,1.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(3.5,42.74,3.5,60.68, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(4.5,21.84,4.5,40.62, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
axis(2, lty = 1, labels = FALSE)

#3 d barplots species
par(bty = "l", mar = c(0.5,0.75,0.5, 0.5))
xaxis<-c(1,2,3,4)
yaxis<-c(98.04,92.31,87.00,70.10)
barplot(yaxis, space = c(0,0,1,0), col = c( rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), ylim = c(0,101), axes = FALSE, ann = FALSE)
segments(0.5,93.34,0.5,100, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(1.5,84.21,1.5,98.30, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
segments(3.5,80.00,3.5,93, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
segments(4.5,60.21,4.5,78.99, col = rgb(232, 164, 35, 255, maxColorValue = 255), lwd = 2)
axis(2, lty = 1, labels = FALSE)

##########################################Scholl analysis########################################################

x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/ZT2-ZT14_Scholl/Scholl.csv")

pairwise.wilcox.test(x$Avg, x$Order, p.adjust.method = "none")

#                      pairwise wilcox	  corrected	
#DmelCS	2	  vs 	DmelCS	14		0.00578	  0.01734	*
#Dsec07	2	  vs 	Dsec07	14		0.02206	  0.04412	*
#DmelCS	2	  vs 	Dsec07	2		  0.07332	  0.07332	n.s.
#DmelCS	14	vs 	Dsec07	14		0.00047	  0.00188	**

x<-read.csv("/Users/michaelshahandeh/Desktop/ZT2-ZT14_Scholl/Scholl_CS.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/ZT2-ZT14_Scholl/Scholl_07.csv")
sec07_means<-c(17.92857143, 15.88888889)
mel_means<-c(14.71428571, 10.27777778)

xaxis<-c(1,2)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$Avg~x$Order ,pch = 16, cex = 4, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 28), xlim = c(0.8,2.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$Avg~y$Order ,pch = 16, cex = 4, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 28), xlim = c(0.8,2.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec07_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")

x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/ZT2-ZT14_Scholl/Scholl OR&28/Scholl.csv")

pairwise.wilcox.test(x$Avg, x$Order, p.adjust.method = "none")


x<-read.csv("/Users/michaelshahandeh/Desktop/Scholl OR&28/Scholl_OR.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/Scholl OR&28/Scholl_28.csv")
sec28_means<-c(15.3, 14.8)
mel_means<-c(12.6, 9.1)

xaxis<-c(1,2)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$Avg~x$Order ,pch = 16, cex = 4, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 28), xlim = c(0.8,2.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$Avg~y$Order ,pch = 16, cex = 4, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 28), xlim = c(0.8,2.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec28_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")

x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/ZT2:ZT14/168hLD_ZT2:14_PDFimmunos/168hLD_Scholl/Scholl_CS.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/ZT2:ZT14/168hLD_ZT2:14_PDFimmunos/168hLD_Scholl/Scholl_07.csv")
sec07_means<-c(17.5, 13.9)
mel_means<-c(11.9, 11.9)
xaxis<-c(1,2)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$Avg~x$Order ,pch = 16, cex = 4, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 28), xlim = c(0.8,2.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$Avg~y$Order ,pch = 16, cex = 4, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 28), xlim = c(0.8,2.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec07_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")

pairwise.wilcox.test(y$Avg, y$Order, p.adjust.method = "none")

#pairwise                    wilcov	corrected
#DmelCS	2	vs 	DmelCS	14		1	      1
#Dsec07	2	vs 	Dsec07	14		0.289	  0.867
#DmelCS	2	vs 	Dsec07	2	  	0.012	  0.048*
#DmelCS	14	vs 	Dsec07	14	0.527	  1

x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/ZT2:ZT14/OR:28_168hLD_ZT2:14_PDFimmunos/Scholl_OR:28_168hLD/Scholl_OR.csv")
y<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/ZT2:ZT14/OR:28_168hLD_ZT2:14_PDFimmunos/Scholl_OR:28_168hLD/Scholl_28.csv")
sec_means<-c(13,11)
mel_means<-c(9.5, 8)

xaxis<-c(1,2)
par(mar = c(1,1,1,1), bty = "l")
stripchart(x$Avg~x$Order ,pch = 16, cex = 4, col = rgb(66, 199, 244, 100, maxColorValue = 255), ylim = c(0, 28), xlim = c(0.8,2.2), axes = FALSE, ann = FALSE, vertical = TRUE, method = "jitter")
stripchart(y$Avg~y$Order ,pch = 16, cex = 4, col = rgb(232, 164, 35, 100, maxColorValue = 255), ylim = c(0, 28), xlim = c(0.8,2.2), axes = FALSE, ann = FALSE, vertical = TRUE, add = TRUE, method = "jitter")
lines(xaxis, mel_means, lwd =2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
lines(xaxis, sec_means, lwd =2, col = rgb(232, 164, 35, 255, maxColorValue = 255))
axis(2, las = 2,  col.axis = "white", col.lab = "white")

x<-read.csv("/Users/michaelshahandeh/Desktop/circ plasticity pub/ZT2:ZT14/OR:28_168hLD_ZT2:14_PDFimmunos/Scholl_OR:28_168hLD/Scholl_ZT2_ZT14_168hLD_OR&28.csv", header = TRUE)

######################################replicate/females for revisions################################################
metadata <- fread("metadata_12-12LD_replicate.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[status == "OK"])						#load data into behavior table
summary(dt)

#change strain name to check for dead flies, change status of dead flies in metadata file and re-run analysis

ggetho(dt[xmv(strain) == c("Dsec.28")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
  stat_bar_tile_etho() +
  facet_wrap(~ region_id) +
  stat_ld_annotations(ld_colours =  c("lightyellow", "black"), outline  = NA, l_duration = hours(12), height = 1, alpha = 0.05)

metadata <- fread("metadata_16-8LD_replicate.csv")				#read in metadata table --> change to .csv for your experiments
metadata										#check metadata
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
metadata
dt <- load_dam(metadata[strain == "Dmel.CS" & status == "OK" |
                          strain == "Dmel.OR" & status == "OK" |
                          strain == "Dsec.07" & status == "OK" |
                          strain == "Dsec.28" & status == "OK"])	
summary(dt)

DmelOR <- ggetho(dt[xmv(strain) == c("Dmel.OR")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-DmelOR$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmelOR<-max(colMeans(df_x))
DmelOR<- ggetho(dt[xmv(strain) == c("Dmel.OR")], aes(x = t, y = activity/maxDmelOR, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

DmelCS <- ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-DmelCS$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmelCS<-max(colMeans(df_x))
DmelCS<- ggetho(dt[xmv(strain) == c("Dmel.CS")], aes(x = t, y = activity/maxDmelCS, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dsec07 <- ggetho(dt[xmv(strain) == c("Dsec.07")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dsec07$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec07<-max(colMeans(df_x))
Dsec07 <- ggetho(dt[xmv(strain) == c("Dsec.07")], aes(x = t, y = activity/maxDsec07, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  #facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dsec28 <- ggetho(dt[xmv(strain) == c("Dsec.28")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Dsec28$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec28<-max(colMeans(df_x))
Dsec28 <- ggetho(dt[xmv(strain) == c("Dsec.28")], aes(x = t, y = activity/maxDsec28, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.452, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dsec13 <- ggetho(dt[xmv(strain) == c("Dsec.13")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Dsec13$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec13<-max(colMeans(df_x))
Dsec13 <- ggetho(dt[xmv(strain) == c("Dsec.13")], aes(x = t, y = activity/maxDsec13, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.452, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dsec32 <- ggetho(dt[xmv(strain) == c("Dsec.32")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6))
x<-Dsec32$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDsec32<-max(colMeans(df_x))
Dsec32 <- ggetho(dt[xmv(strain) == c("Dsec.32")], aes(x = t, y = activity/maxDsec32, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(232, 164, 35, 255, maxColorValue = 255), fill = rgb(232, 164, 35, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.452, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(232, 164, 35, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  #scale_y_continuous(name= "", breaks = c(0,1,2,3,4)) 
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

LZV76 <- ggetho(dt[xmv(strain) == c("Dmel.LZV_L76")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-LZV76$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxLZV76<-max(colMeans(df_x))
LZV76<- ggetho(dt[xmv(strain) == c("Dmel.LZV_L76")], aes(x = t, y = activity/maxLZV76, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

LZV72 <- ggetho(dt[xmv(strain) == c("Dmel.LZV_L72")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-LZV72$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxLZV72<-max(colMeans(df_x))
LZV72<- ggetho(dt[xmv(strain) == c("Dmel.LZV_L72")], aes(x = t, y = activity/maxLZV72, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(66, 199, 244, 255, maxColorValue = 255), fill = rgb(66, 199, 244, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(66, 199, 244, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

MD221 <- ggetho(dt[xmv(strain) == c("Dsim.MD221")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-MD221$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxMD221<-max(colMeans(df_x))
MD221<- ggetho(dt[xmv(strain) == c("Dsim.MD221")], aes(x = t, y = activity/maxMD221, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(169, 196, 255, 255, maxColorValue = 255), fill = rgb(169, 196, 255, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(169, 196, 255, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

MD242 <- ggetho(dt[xmv(strain) == c("Dsim.MD242")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-MD242$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxMD242<-max(colMeans(df_x))
MD242<- ggetho(dt[xmv(strain) == c("Dsim.MD242")], aes(x = t, y = activity/maxMD242, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(169, 196, 255, 255, maxColorValue = 255), fill = rgb(169, 196, 255, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(169, 196, 255, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dmau90 <- ggetho(dt[xmv(strain) == c("Dmau.90")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dmau90$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmau90<-max(colMeans(df_x))
Dmau90<- ggetho(dt[xmv(strain) == c("Dmau.90")], aes(x = t, y = activity/maxDmau90, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(197, 168, 197, 255, maxColorValue = 255), fill = rgb(197, 168, 197, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(197, 168, 197, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

Dmau91 <- ggetho(dt[xmv(strain) == c("Dmau.91")], aes(x = t, y = activity, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) 
x<-Dmau91$data$activity
df_x <- as.data.frame(matrix(x, ncol = 48,  byrow = TRUE))
maxDmau91<-max(colMeans(df_x))
Dmau91<- ggetho(dt[xmv(strain) == c("Dmau.91")], aes(x = t, y = activity/maxDmau91, colour = strain), summary_time_window = mins(30), time_wrap = hours(24), time_offset = hours(6)) +
  stat_pop_etho(geom = "bar", colour = rgb(197, 168, 197, 255, maxColorValue = 255), fill = rgb(197, 168, 197, 100, maxColorValue = 255)) +
  stat_ld_annotations(ld_colours =  c("yellow", "black"), height = 0.04, alpha = 0.45, outline  = NA, l_duration = hours(12)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", col = rgb(197, 168, 197, 255, maxColorValue = 255)) +
  # facet_grid(strain ~ .) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(name= "", breaks = c(0,.25,0.5,0.75,1)) 

ggarrange(DmelCS+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), DmelOR+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), LZV72+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), LZV76+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),MD221+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),MD242+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),Dmau90+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dmau91+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dsec07+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dsec13+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dsec28+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dsec32+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"),
          heights = c(3, 3, 3, 3, 3), ncol = 6, nrow = 2)

ggarrange(DmelCS+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), DmelOR+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dsec07+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), Dsec28+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          heights = c(3, 3, 3, 3), ncol = 1, nrow = 4)

dt[, MA := ifelse(t %% hours(24) < hours(21), "R", "MA")]
dt[, S := ifelse(t %% hours(24) < hours(18), "R", "S")]

summary_dt <- 
  rejoin(dt[,
            .(
              # this is where the computation happens
              activity_MA = mean(activity[MA == "MA"]),
              activity_S = mean(activity[S == "S"]) 
            ),
            ,by=id])
summary_dt

order <-c("Dmel.CS", "Dmel.OR", "Dmel.LZV_L72", "Dmel.LZV_L76", "Dsim.MD221", "Dsim.MD242", "Dmau.90", "Dmau.91", "Dsec.07", "Dsec.13", "Dsec.28", "Dsec.32")
order <-c("Dmel.CS", "Dmel.OR",  "Dsec.07", "Dsec.28")

p<- ggplot(summary_dt, aes(x=factor(strain, level = order), y=activity_MA/max(activity_MA), fill=strain), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(169, 196, 255, 100, maxColorValue = 255), rgb(169, 196, 255, 100, maxColorValue = 255), rgb(197, 168, 197, 100, maxColorValue = 255), rgb(197, 168, 197, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.2),alpha=0.3) +
  scale_y_continuous(name= "") + 
  theme_classic() + 
  scale_x_discrete(name = "") + 
  theme(legend.position = "none")
p$data$strain
p$data$activity_MA
x<-cbind(p$data$strain, p$data$activity_MA)
write.csv(x, "~/Desktop/Females_MA.csv")

ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

pairwise.wilcox.test(summary_dt[, activity_MA/max(activity_MA)], 
                     summary_dt[, strain], p.adjust = "bonf")
x<-read.csv("/Users/michaelshahandeh/Dropbox/CircadEvo/Nature REVISION/DATA/Data and statistics/EDFig4/Epeak_replicate_168hLD.csv", header = TRUE)
p<- ggplot(x, aes(x=Order, y=E168, fill =Order), inherit.aes = TRUE) +
  geom_boxplot(fill= c(rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(66, 199, 244, 100, maxColorValue = 255), rgb(169, 196, 255, 100, maxColorValue = 255), rgb(169, 196, 255, 100, maxColorValue = 255), rgb(197, 168, 197, 100, maxColorValue = 255), rgb(197, 168, 197, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255), rgb(232, 164, 35, 100, maxColorValue = 255)), color="black", outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 1.2),alpha=0.3) +
  scale_y_continuous(name= "Evening peak time 16:8h LD(hours)", limits = c(11,16), breaks = c(10, 12, 14, 16, 18, 20)) + 
  theme_classic() + 
  theme(legend.position = "none")
p
ggarrange(p+ rremove("x.text")+rremove("xlab")+ rremove("y.text")+rremove("ylab"), 
          ncol = 1, nrow = 1)

PT<-pairwise.wilcox.test(x$E1212, x$strain, p.adjust.method = "bonferroni")
pairwise.wilcox.test(x$E1212, x$strain)

###############################Pdf enhancer correlation##########################################

latenhancer <- c( 27.6648,	32.1656,	33.8361,	35.7596,	41.2033,	45.2538, 39.3999,	47.5162,	41.8719, 9.145,	7.3697,	1.9403, 9.082)

MAFenhancer<-c(0.256789404,	0.242007872,	0.255764755,	0.218248066,	0.258537463,	0.22414767, 0.244833719,	0.301817974,	0.334667539, 0.134920635,	0.042857143,	0.078817734, 0.095238095)

cor.test(lat, MAFCDS, method = "spearman")

cor.test(lat, MAFenhancer, method = "spearman")

cor.test(lat, MAFgene, method = "spearman")

par(pty = "s")
#plot(lat, MAFCDS, type = "p", pch = 19, col = "darkolivegreen4", axes = FALSE, ann = FALSE, ylim = c(0,.6), xlim = c(0,50))
plot(lat, MAFenhancer, pch = 19, col = "goldenrod", axes = FALSE, ann = FALSE, ylim = c(0,.6), xlim = c(0,50))
points(latitude_ctrl, SNV_avg, pch = 19, col = "black" )
#points(lat, MAFgene, pch = 19, col = "black")
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAFenhancer ~ lat), lty = 2, col = "goldenrod")
#abline(lm(MAFCDS ~ lat), lty = 2, col = "darkolivegreen4")
abline(lm(SNV_avg ~ latitude_ctrl), lty = 2, col = "black")
#abline(lm(MAFgene ~ lat), lty = 2, col = "black")
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")
#text(12, .27, "rho = 0.659 p = 0.017", col = "darkolivegreen4")
text(25, .58, "enhancer: rho = 0.769 p = 0.007", col = "goldenrod")
#text(40, .17, "rho = 0.742 p = 0.005", col = "black")
text(24.5, .53, "paralytic: rho = 0.775 p = 0.003", col = "black")


#enhancer
p22282195<- c(	0.183673469387755,	0.222222222222222,	0.140845070422535,	0.175,	0.188118811881188,	0.0697674418604651,	0.0728476821192053,	NA,	NA,	0,	0.1,	0.0344827586206897,	0.166666666666667)																
p22283060<- c(	0.131578947368421,	0.230769230769231,	0.247787610619469,	0.125,	0.220183486238532,	0.128787878787879,	0.150537634408602,	0.221311475409836,	0.254950495049505,	0.25,	0,	0.0344827586206897,	0)																
p22283066<- c(	0.11965811965812,	0.222972972972973,	0.241379310344828,	0.12280701754386,	0.207207207207207,	0.136363636363636,	0.194444444444444,	0.239915074309979,	0.254408060453401,	0.25,	0,	0.0344827586206897,	0)																
p22283068<- c(	0.686440677966102,	0.428571428571429,	0.398305084745763,	0.421052631578947,	0.459459459459459,	0.424242424242424,	0.338624338624339,	0.326226012793177,	0.385,	0,	0.2,	0.137931034482759,	0.166666666666667)																
p22283104<- c(	0.271317829457364,	0.153409090909091,	0.226277372262774,	0.115942028985507,	0.184210526315789,	0.211678832116788,	NA,	0.335341365461847,	0.421717171717172,	0,	0,	0.0689655172413793,	0)																
p22283133<- c(	0.336363636363636,	0.145569620253165,	0.297297297297297,	0.130434782608696,	0.205128205128205,	0.242424242424242,	0.303191489361702,	0.393684210526316,	0.433701657458564,	0.444444444444444,	0,	0.172413793103448,	0.166666666666667)																
p22284354<- c(	0.0684931506849315,	0.290540540540541,	0.238461538461538,	0.4375,	0.345454545454545,	0.355769230769231,	0.409356725146199,	0.294429708222812,	0.258227848101266,	0,	0,	0.0689655172413793,	0.166666666666667)																
cor.test(lat, p22282195, method = "spearman")
cor.test(lat, p22283060, method = "spearman")
cor.test(lat, p22283066, method = "spearman")
cor.test(lat, p22283068, method = "spearman")
cor.test(lat, p22283104, method = "spearman") #p=0.03 BUT 1241 bp upstream of the start codon
cor.test(lat, p22283133, method = "spearman")
cor.test(lat, p22284354, method = "spearman") #p=0.03 BUT 2491 bp upstream of the start codon

enhancer_avg_minus_sigs<-c(	0.29154297,	0.250021095,	0.265122875,	0.194858886,	0.256019434,	0.200317125,	0.211929118,	0.295284193,	0.332015053,	0.188888889,	0.06,	0.082758621)
cor.test(lat, enhancer_avg_minus_sigs, method = "spearman")


###############################Czn enhancer correlatio##########################################
lat <- c( 27.6648,	32.1656,	33.8361,	35.7596,	41.2033,	45.2538,	39.3999,	47.5162,	41.8719,	7.3697,	9.145,	9.082,	1.9403)

MAF<-c(0.53535939,	0.498899313,	0.557696028,	0.601349379,	0.546203077,	0.559577669, 0.62075191,	0.596554604,	0.538049876, 0.777777778,	0.511111111,	0.613793103,	0.466666667)
MAF2.4CRZ<-c(0.369407898,	0.373694439,	0.403906527,	0.369110475,	0.354029857,	0.389418964,	0.457374551,	0.370405069,	0.355834006,	0.418253968,	0.416666667,	0.404761905,	0.468050366)
y<-cor.test(lat, MAF2.4CRZ, method = "spearman")
y

#	Spearman's rank correlation rho
#data:  lat and MAF2.4CRZ
#S = 584, p-value = 0.03212
#alternative hypothesis: true rho is not equal to 0
#sample estimates:  rho = -0.6043956 


par(mar=c(1,1,1,1))
plot(lat, MAF2.4CRZ, type = "p", pch = 19, col = "darkgrey", axes = FALSE, ann = FALSE, ylim = c(0,0.65), xlim = c(0,50))
points(lat, MAFenhancer, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 19)
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAF2.4CRZ ~ lat), lty = 1, col = "darkgrey")
abline(lm(MAFenhancer ~ lat), lty = 1, col = rgb(66, 199, 244, 255, maxColorValue = 255))
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")

text(27, .01, "Czn rho = 0.1758242 p = 0.566", col = "darkgrey")
text(27, .09, "Pdf enhancer: rho = 0.769 p = 0.007", col = rgb(66, 199, 244, 255, maxColorValue = 255))

###############################Ast-A enhancer correlatio##########################################
lat <- c( 27.6648,	32.1656,	33.8361,	35.7596,	41.2033,	45.2538,	39.3999,	47.5162,	41.8719,	7.3697,	9.145,	9.082,	1.9403)

MAF2.1ASTA<-c(0.445399757,	0.398801497,	0.531814438,	0.360317299,	0.467556959,	0.409468892,	0.434026126,	0.421810552,	0.513204648,	0.598290598,	0.575579976,	0.570524724,	0.602564103)
MAF2.4ASTA<-c(0.438717764,	0.389686136,	0.506814706,	0.33878203,	0.44058058,	0.386617702,	0.434026126,	0.399734448,	0.482456106,	0.582086168,	0.626984127,	0.613095238,	0.596439625)
y<-cor.test(lat, MAF2.4ASTA, method = "spearman")
y
#Spearmans rank correlation rho
#data:  lat and MAF2.4
#S = 604, p-value = 0.01713
#alternative hypothesis: true rho is not equal to 0
#sample estimates: rho = -0.6593407

par(pty = "m")
plot(lat, MAF2.4ASTA, type = "p", pch = 19, col = "darkgrey", axes = FALSE, ann = FALSE, ylim = c(0,0.65), xlim = c(0,50))
points(lat, MAFenhancer, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 19)
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAF2.4ASTA ~ lat), lty = 1, col = "darkgrey")
abline(lm(MAFenhancer ~ lat), lty = 1, col = rgb(66, 199, 244, 255, maxColorValue = 255))
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")

text(27, .05, "Ast-A enhancer: rho = -0.6263736 p = 0.02527", col = "darkgrey")
text(27, .09, "Pdf enhancer: rho = 0.769 p = 0.007", col = rgb(66, 199, 244, 255, maxColorValue = 255))


###############################CCHa1 enhancer correlatio##########################################
lat <- c( 27.6648,	32.1656,	33.8361,	35.7596,	41.2033,	45.2538,	39.3999,	47.5162,	41.8719,	7.3697,	9.145,	9.082,	1.9403)

MAF2.4CCHa1<-c(0.47748498,	0.483965078,	0.473332165,	0.457529853,	0.429485539,	0.396924429,	0.357137124,	0.320628381,	0.35853746,	0.406060606,	0.402777778,	0.481481481,	0.456510922)

y<-cor.test(lat, MAF2.4CCHa1, method = "spearman")
y

#	Spearman's rank correlation rho
#data:  lat and MAF2.4CCHa1
#S = 572, p-value = 0.04489
#alternative hypothesis: true rho is not equal to 0
#sample estimates:  rho = -0.5714286

par(pty = "m")
plot(lat, MAF2.4CCHa1, type = "p", pch = 19, col = "darkgrey", axes = FALSE, ann = FALSE, ylim = c(0,0.65), xlim = c(0,50))
points(lat, MAFenhancer, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 19)
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAF2.4CCHa1 ~ lat), lty = 1, col = "darkgrey")
abline(lm(MAFenhancer ~ lat), lty = 1, col = rgb(66, 199, 244, 255, maxColorValue = 255))
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")

text(27, .05, "CCHa1 enhancer: rho = -0.8736264 p = 6.817e-05", col = "darkgrey")
text(27, .09, "Pdf enhancer: rho = 0.769 p = 0.007", col = rgb(66, 199, 244, 255, maxColorValue = 255))


###############################sNPF enhancer correlatio##########################################
lat <- c( 27.6648,	32.1656,	33.8361,	35.7596,	41.2033,	45.2538,	39.3999,	47.5162,	41.8719,	7.3697,	9.145,	9.082,	1.9403)

MAF2.4sNPF<-c(0.437976708,	0.487079003,	0.526268542,	0.462699832,	0.488284337,	0.541842386,	0.390479213,	0.405403543,	0.4058489,	0.472727273,	0.402173913,	0.449275362,	0.412538544)

y<-cor.test(lat, MAF2.4sNPF, method = "spearman")
y

#Spearman's rank correlation rho

#data:  lat and MAF2.4sNPF
#S = 342, p-value = 0.849
#alternative hypothesis: true rho is not equal to 0
#sample estimates: rho = 0.06043956 

par(pty = "m")
plot(lat, MAF2.4sNPF, type = "p", pch = 19, col = "darkgrey", axes = FALSE, ann = FALSE, ylim = c(0,0.65), xlim = c(0,50))
points(lat, MAFenhancer, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 19)
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAF2.4sNPF ~ lat), lty = 1, col = "darkgrey")
abline(lm(MAFenhancer ~ latenhancer), lty = 1, col = rgb(66, 199, 244, 255, maxColorValue = 255))
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")

text(25,0.8, "sNPF enhancer: rho = 0.1153846 p = 0.7096", col = "darkgrey")
text(25,.7, "Pdf enhancer: rho = 0.769 p = 0.007", col = rgb(66, 199, 244, 255, maxColorValue = 255))


###############################ITP enhancer correlatio##########################################
lat <- c( 27.6648,	32.1656,	33.8361,	35.7596,	41.2033,	45.2538,	39.3999,	47.5162,	41.8719,	7.3697,	9.145,	9.082,	1.9403)

MAF2.4ITP<-c(0.402533093,	0.432175864,	0.394969701,	0.326079929,	0.428163805,	0.400353036,	0.483908704,	0.446294796,	0.393528845,	0.466666667,	0.518518519,	0.416666667,	0.500410509)

y<-cor.test(lat, MAF2.4ITP, method = "spearman")
y

#Spearman's rank correlation rho

#data:  lat and MAF2.4ITP
#S = 518, p-value = 0.1516
#alternative hypothesis: true rho is not equal to 0
#sample estimates: rho = -0.4230769 

par(pty = "m")
plot(lat, MAF2.4ITP, type = "p", pch = 19, col = "darkgrey", axes = FALSE, ann = FALSE, ylim = c(0,0.6), xlim = c(0,50))
points(lat, MAFenhancer, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 19)
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAF2.4ITP ~ lat), lty = 1, col = "darkgrey")
abline(lm(MAFenhancer ~ lat), lty = 1, col = rgb(66, 199, 244, 255, maxColorValue = 255))
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")

text(25,0.8, "sNPF enhancer: rho = 0.1153846 p = 0.7096", col = "darkgrey")
text(25,.7, "Pdf enhancer: rho = 0.769 p = 0.007", col = rgb(66, 199, 244, 255, maxColorValue = 255))

###############################ASTC enhancer correlatio##########################################
lat <- c( 27.6648,	32.1656,	33.8361,	35.7596,	41.2033,	45.2538,	39.3999,	47.5162,	41.8719,	7.3697,	9.145,	9.082,	1.9403)

MAF2.4ASTC<-c(0.467944306,	0.486835945,	0.476908575,	0.442300533,	0.45055653,	0.438268643,	0.461430948,	0.480939534,	0.495954514,	0.35380117,	0.393274854,	0.443859649,	0.416743556)

y<-cor.test(lat, MAF2.4ASTC, method = "spearman")
y

#Spearman's rank correlation rho

#data:  lat and MAF2.4ASTC
#S = 170, p-value = 0.06416
#alternative hypothesis: true rho is not equal to 0
#sample estimates:  rho = 0.532967  

par(mar =c(1,1,1,1))
plot(lat, MAF2.4ASTC, type = "p", pch = 19, col = "darkgrey", axes = FALSE, ann = FALSE, ylim = c(0,0.6), xlim = c(0,50))
points(lat, MAFenhancer, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 19)
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAF2.4ASTC ~ lat), lty = 1, col = "darkgrey")
abline(lm(MAFenhancer ~ lat), lty = , col = rgb(66, 199, 244, 255, maxColorValue = 255))
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")

###############################Tk enhancer correlatio##########################################
lat <- c( 27.6648,	32.1656,	33.8361,	35.7596,	41.2033,	45.2538,	39.3999,	47.5162,	41.8719,	7.3697,	9.145,	9.082,	1.9403)

MAF2.4Tk<-c(0.288796733,	0.293201238,	0.378867567,	0.377096564,	0.362776831,	0.285282142,	0.397222409,	0.39746843,	0.210582011,	0.149659864,	0.266666667,	0.281787389, 0.379820867)

y<-cor.test(lat, MAF2.4Tk, method = "spearman")
y

#	Spearman's rank correlation rho
#data:  lat and MAF2.4Tk
#S = 258, p-value = 0.3341
#alternative hypothesis: true rho is not equal to 0
#sample estimates: rho = 0.2912088 


par(pty = "m")
plot(lat, MAF2.4Tk, type = "p", pch = 19, col = "darkgrey", axes = FALSE, ann = FALSE, ylim = c(0,0.8), xlim = c(0,50))
points(lat, MAFenhancer, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 19)
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAF2.4Tk ~ lat), lty = 2, col = "darkgrey")
abline(lm(MAFenhancer ~ lat), lty = 2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")

###############################NPF enhancer correlatio##########################################
lat <- c( 27.6648,	32.1656,	33.8361,	35.7596,	41.2033,	45.2538,	39.3999,	47.5162,	41.8719,	7.3697,	9.145,	9.082,	1.9403)

MAF2.4NPF<-c(0.48056673,	0.466027721,	0.305519041,	0.322729734,	0.360190864,	0.441851547,	0.44386108,	0.381013542,	0.313971352,	0.52,	0.955555556,	0.866666667,	0.664778325)

y<-cor.test(lat, MAF2.4NPF, method = "spearman")
y

#	Spearman's rank correlation rho
#data:  lat and MAF2.4NPF
#S = 628, p-value = 0.006892
#alternative hypothesis: true rho is not equal to 0
#sample estimates: rho = -0.7252747  

par(pty = "m")
plot(lat, MAF2.4NPF, type = "p", pch = 19, col = "darkgrey", axes = FALSE, ann = FALSE, ylim = c(0,1), xlim = c(0,50))
points(lat, MAFenhancer, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 19)
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAF2.4NPF ~ lat), lty = 2, col = "darkgrey")
abline(lm(MAFenhancer ~ lat), lty = 2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")

###############################DILP2 enhancer correlatio##########################################
lat <- c( 27.6648,	32.1656,	33.8361,	35.7596,	41.2033,	45.2538,	39.3999,	47.5162,	41.8719,	7.3697,	9.145,	9.082,	1.9403)

MAF2.4DILP2<-c(0.384230855,	0.258853131,	0.362681936,	0.277794008,	0.209944637,	0.323422175,	0.16709751,	0.22073603,	0.193178313,	0.785119048,	0.6,	0.746666667,	0.538679424)

y<-cor.test(lat, MAF2.4DILP2, method = "spearman")
y

#	Spearman's rank correlation rho
#data:  lat and MAF2.4DILP2
#S = 654, p-value = 0.001844
#alternative hypothesis: true rho is not equal to 0
#sample estimates: rho = -0.7967033 

par(pty = "m")
plot(lat, MAF2.4DILP2, type = "p", pch = 19, col = "darkgrey", axes = FALSE, ann = FALSE, ylim = c(0,1), xlim = c(0,50))
points(lat, MAFenhancer, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 19)
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAF2.4DILP2 ~ lat), lty = 2, col = "darkgrey")
abline(lm(MAFenhancer ~ lat), lty = 2, col = rgb(66, 199, 244, 255, maxColorValue = 255))
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")


############all together############################all together############################all together################

par(mar = c(1,1,1,1))
plot(latenhancer, MAFenhancer, col = rgb(66, 199, 244, 255, maxColorValue = 255), pch = 19, type ="p", axes = FALSE, ann = FALSE, ylim = c(0,0.6), xlim = c(0,50))
axis(1, lty = 1)
axis(2, lty = 1)
abline(lm(MAF2.4ITP ~ lat), lty = 1, col = "darkgrey", lwd = 2)
abline(lm(MAF2.4sNPF ~ lat), lty = 1, col = "darkgrey", lwd = 2)
abline(lm(MAF2.4CCHa1 ~ lat), lty = 1, col = "darkgrey", lwd = 2)
abline(lm(MAF2.1ASTA ~ lat), lty = 1, col = "darkgrey", lwd = 2)
abline(lm(MAF2.4CRZ ~ lat), lty = 1, col = "darkgrey", lwd = 2)
#abline(lm(MAF2.4NPF ~ lat), lty = 1, col = "darkgrey", lwd = 2)
#abline(lm(MAF2.4DILP2 ~ lat), lty = 1, col = "darkgrey", lwd = 2)
abline(lm(MAF2.4ASTC ~ lat), lty = 1, col = "darkgrey", lwd = 2)
#abline(lm(MAF2.4Tk ~ lat), lty = 1, col = "darkgrey", lwd = 2)
abline(lm(MAFenhancer ~ latenhancer), lty = 1, col = rgb(66, 199, 244, 255, maxColorValue = 255), lwd = 2)
title(ylab = "Mean alt allele frequency", xlab = "Distance from equator")



