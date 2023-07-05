library(damr)
library(ggetho)
library(peakPick)
library(zoo)

#construct metadata file for last 3 days of long photophase

#set directory to where your monitor files and metadata files are
DATA_DIR<- "/Users/michaelshahandeh/Documents/Circadian_rhythm/DAM data" 			#create working directory
list.files(DATA_DIR, pattern= "*.txt|*.csv")	#check files in working directory
setwd(DATA_DIR) 

metadata <- fread("metadata_KODFs-1212LD.csv")				#read in metadata table with updated statuses
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
dt <- load_dam(metadata[replicate == "20" & status == "OK"])						#load data into behavior table excluding "DEAD" flies
summary(dt)

#Plot double actograms for all individuals of a single strain to identify dead flies; repeat for all strains by changing strain name in quotes

ggetho(dt[xmv(strain) == c("Dsec28xHr38^56")], aes(x = t, z=activity), summary_time_window = mins(30), multiplot = 2) +
       stat_bar_tile_etho() +
       facet_wrap(~ region_id)
              
#change status of dead flies in metadata file ("OK" --> "DEAD") and re-run analysis

metadata <- fread("metadata_KODFs-1212LD.csv")				#read in metadata table with updated statuses
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)	#link metadata to monitorfiles
dt <- load_dam(metadata[status == "OK" & replicate == "9"])						#load data into behavior table excluding "DEAD" flies
# | strain == "Dsec.28noni"]) add to end of metadata to sort by XX or YY
summary(dt)

p<- ggetho(dt[xmv(strain) == "Dsec.XF1" & xmv(region_id) == MPS], aes(x = t, y = activity), summary_time_window = mins(10), time_wrap = hours(24), time_offset = hours(6)) +
       stat_pop_etho() +
       stat_ld_annotations(ld_colours =  c("lightyellow", "black"), height = 1, alpha = 0.3, outline  = NA, l_duration = hours(12))
       
x<-p$data$activity
last6<-x[109:144]
first18<-x[1:108]
ts<-append(last6, first18) # rearrange timeseries so that morning and evening peaks are more centered
smooth<-rollmean(rollmean(ts, 3), 3) #traingular (double) rolling average witha 30 minute sliding window

morning<-mean(which(grepl(max(smooth[18:42]), smooth[18:42]))) +17   #takes max value position in morning from 18 to 1
evening<-mean(which(grepl(max(smooth[102:130]), smooth[102:132]))) + 101	#takes max value position in evening from 11 to 16

tiff("~/Desktop/peakdetection/Dsec.XF1_MPS.tiff", width = 5, height = 2.5, units = "in", res = 600)
par(bty = "l", mar = c(4,4,1,1))
plot(ts, type = "l", axes = FALSE, ann = FALSE, lwd = 2)
axis(1, at = c(0, 18, 36, 54, 72, 90, 108, 126, 144), labels = c(-6, -3, 0, 3, 6, 9, 12, 15, 18))
axis(2)
title(xlab = "ZT (hours)", ylab = "Activity")
lines(smooth, col = "blue", lwd = 2)
#points((1:144)[peakspicked], smooth[peakspicked], col="red", pch = 16)
points(morning, max(smooth[morning]), col = "red", pch = 16)
points(evening, max(smooth[evening]), col = "red", pch = 16)
box()
dev.off()

peakloc<-c(morning/6-6, evening/6-6)
peakloc <-append("Dsec.XF1_MPS", peakloc)
write(peakloc, file = "~/Desktop/peakdetection/peak_locations_Dsec.XF1.txt", ncolumns = 3, append = TRUE, sep = ",")

