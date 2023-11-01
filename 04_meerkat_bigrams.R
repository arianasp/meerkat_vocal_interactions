library(gtools)
library(stringr)
library(stringi) 
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(chron)
library(lubridate) 
library(gridExtra)

#load data
load("04_meerkat_bigrams.RData")

#calculate inter call intercals
ICI <- quantile(call_data$lag,  probs = seq(.1, .9, by = .1), na.rm = T)[9]

#recode some inconsistencies in call type names 
call_data$stn_call_type <- str_replace(call_data$stn_call_type, "lc", "ld")
call_data$type_group <- str_replace(call_data$type_group, "lc", "ld")

#check sample sizes and call type names
table(call_data$type_group)

## the columns containing call type info:
## hyb - indicator of fused or sequential calls
## stn_call_type - details the call type elements in each call
## type_group - only main calls, hybrids were collapsed into the parent type acording to overal call frequency 
options(digits.secs = 3)

#reformat times and dates
call_data$year <- str_sub(call_data$date, 1, 4)
call_data <- call_data[which(call_data$isCall == 1),]

call_data$t0GPS_UTC <- as.POSIXct(call_data$t0GPS_UTC, tz = "UTC")
call_data$tMidGPS_UTC <- as.POSIXct(call_data$tMidGPS_UTC, tz = "UTC")
call_data$tendGPS_UTC <- as.POSIXct(call_data$tendGPS_UTC, tz = "UTC")

 
#reformatting sequential calls into two row format
for (row in which(call_data$hyb == "sq")) {
  print(row)
  ext_row <- call_data[row, ]
  call_data$type_group[row] <- stri_extract_first(call_data[row, "stn_call_type"], regex="\\w+")
  call_data$tendGPS_UTC[row] <- as.POSIXct(call_data$tMidGPS_UTC[row], tz = "UTC")-0.05
  ext_row$type_group <- stri_extract_last(call_data[row, "stn_call_type"], regex="\\w+")
  ext_row$t0GPS_UTC <- as.POSIXct(ext_row$tMidGPS_UTC,  tz = "UTC")+0.05
  call_data <- rbind(call_data,  ext_row)
  }


all_pairs <- data.frame()
for (year in unique(call_data$year)) {
  year_select <- call_data[which(call_data$year == year), ]
  
  for (ind in unique(year_select$ind)) {
    ind_select <- year_select[which(year_select$ind == ind) ,] 
    
    for (day in unique(ind_select$date)) {
      day_select <- ind_select[which(ind_select$date == day) ,]
      day_select <-  day_select[order(day_select$t0GPS_UTC),]
      day_select$brake <- as.POSIXct(day_select$t0GPS_UTC , tz = "UTC") - 
        lag(as.POSIXct(day_select$tendGPS_UTC,  tz = "UTC"))
      pairs_per_day <-
        cbind(
          paste(day_select$type_group[-length(day_select$type_group)], "_",   day_select$type_group[-1], sep = ""),
          paste(day_select$hyb[-length(day_select$hyb)], "_",   day_select$hyb[-1], sep = ""),
          day_select$brake[-1],
          day_select$ind[1],
          day_select$wavFileName[1])
      
      pairs_per_day <- as.data.frame(pairs_per_day)
      pairs_per_day$V3 <- as.numeric(pairs_per_day$V3)
      pairs_per_day <- pairs_per_day[which(pairs_per_day$V3 < ICI), ] #here we filter by minimum lag
      
      all_pairs <- rbind(all_pairs, pairs_per_day)
    }
  }
}


colnames(all_pairs) <- c("pair", "hyb", "brake", "ID", "file")




all_pairs <- all_pairs[-which(grepl("unk", all_pairs$pair)), ]
all_pairs <- all_pairs[-which(grepl("NA", all_pairs$pair)), ]

pair_counts <- as.data.frame(table(all_pairs$pair))

frequent_call_pairs <- pair_counts[which(pair_counts$Freq > 50),]

frequent_pairs <- all_pairs[which(all_pairs$pair %in% frequent_call_pairs$Var1) ,]

pairs_only <- separate(frequent_pairs, pair, c("call_1", "call_2"), sep = "_") [ , 1:2]

input.matrix <- pairs_only

#input.matrix <- read.table(file.choose(), header=TRUE, sep="\t", quote="", comment.char="")
names(input.matrix) <- c("W_C", "Coll_Word")
input.matrix <- table(input.matrix$Coll_Word, input.matrix$W_C)

# computation
pearson.residuals <- as.data.frame.matrix(chisq.test(input.matrix, correct=FALSE)$residuals)
output.table <- data.frame(COLLOCATE=rownames(pearson.residuals), pearson.residuals,
                           SUMABSDEV=apply(pearson.residuals, 1, \(af) sum(abs(af))),
                           LARGESTPREF=colnames(pearson.residuals)[apply(pearson.residuals, 1, \(af) which.max(af))])


longData<-melt(output.table[ , -(ncol(output.table) - 1)])
longData<-longData[longData$value!=0,]


longData$r_val <- round(longData$value, 2)



self <- ggplot(longData, aes(x = variable, y = COLLOCATE)) + 
  geom_raster(aes(fill= value)) + 
  scale_fill_gradient2(low="blue", high="red", mid = "white") +
  labs(x="Response", y="Initiation", title="self-reply") +
  geom_text(  aes(x=variable, y=COLLOCATE, label = r_val), color="black", size=5) + 
  theme_bw() + theme(axis.text.x=element_text(size=12, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=12),
                     plot.title=element_text(size=14))


all_pairs <- data.frame()
for (year in unique(call_data$year)) {
  year_select <- call_data[which(call_data$year == year), ]
   
    for (day in unique(year_select$date)) {
      day_select <- year_select[which(year_select$date == day) ,]
      day_select <-  day_select[order(day_select$t0GPS_UTC),]
      day_select$brake <- as.POSIXct(day_select$t0GPS_UTC , tz = "UTC") - 
        lag(as.POSIXct(day_select$tendGPS_UTC,  tz = "UTC"))
      day_select$resp <- lag(day_select$ind)
      pairs_per_day <-
        cbind(
          paste(day_select$type_group[-length(day_select$type_group)], "_",   day_select$type_group[-1], sep = "") ,
          day_select$brake[-1],
          lag(day_select$ind), 
          day_select$ind) 
          
      
          pairs_per_day <- as.data.frame(pairs_per_day)
          pairs_per_day$V2 <- as.numeric(pairs_per_day$V2)
          pairs_per_day <- pairs_per_day[which(pairs_per_day$V2 < ICI), ] #here we filter by minimum lag
          pairs_per_day <- pairs_per_day[which(pairs_per_day$V3 != pairs_per_day$V4), ]
          
          all_pairs <- rbind(all_pairs, pairs_per_day)
    }
  }



colnames(all_pairs) <- c("pair","brake", "ID1", "ID2")




all_pairs <- all_pairs[-which(grepl("unk", all_pairs$pair)), ]
all_pairs <- all_pairs[-which(grepl("NA", all_pairs$pair)), ]

pair_counts <- as.data.frame(table(all_pairs$pair))

frequent_call_pairs <- pair_counts[which(pair_counts$Freq > 50),]

frequent_pairs <- all_pairs[which(all_pairs$pair %in% frequent_call_pairs$Var1) ,]

pairs_only <- separate(frequent_pairs, pair, c("call_1", "call_2"), sep = "_") [ , 1:2]

input.matrix <- pairs_only


names(input.matrix) <- c("W_C", "Coll_Word")
input.matrix <- table(input.matrix$Coll_Word, input.matrix$W_C)

# computation
pearson.residuals <- as.data.frame.matrix(chisq.test(input.matrix, correct=FALSE)$residuals)
output.table <- data.frame(COLLOCATE=rownames(pearson.residuals), pearson.residuals,
                           SUMABSDEV=apply(pearson.residuals, 1, \(af) sum(abs(af))),
                           LARGESTPREF=colnames(pearson.residuals)[apply(pearson.residuals, 1, \(af) which.max(af))])

 



longData<-melt(output.table[ , -(ncol(output.table) - 1)])
longData<-longData[longData$value!=0,]


longData$r_val <- round(longData$value, 2)



non_self <- ggplot(longData, aes(x = variable, y = COLLOCATE)) + 
       geom_raster(aes(fill= value)) + 
       scale_fill_gradient2(low="blue", high="red", mid = "white", limits = c(-70,180), name = "") +
       labs(x="Response", y="", title="non-self-reply") +
       geom_text(  aes(x=variable, y=COLLOCATE, label = r_val), color="black", size=5) + 
       theme_bw() + theme(axis.text.x=element_text(size=12, angle=0, vjust=0.3),
                          axis.text.y=element_text(size=12),
                          plot.title=element_text(size=14))
 
grid.arrange(self, non_self, nrow = 1)

