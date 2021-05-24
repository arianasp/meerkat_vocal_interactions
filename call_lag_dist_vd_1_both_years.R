library(lubridate)
library(dplyr)
library(splitstackshape) 
library(varhandle)
library(ggplot2)
library(gridExtra)
library(markovchain) 
library(asnipe)
library(igraph)
library(CINNA)
library(reshape2)
library(gtools)

library(stringr)

# Analyses of meerkat call sequences

# parameters

#uncomment to use on Vlad's machine
codedir <- "C:/Users/vdemartsev/Dropbox/meerkats_shared/code"
datadir <- "C:/Users/vdemartsev/Dropbox/meerkats_shared/data"

#uncomment to use on Ari's machine
#codedir <- '~/Dropbox/meerkats/meerkats_shared/code'
#datadir <- '~/Dropbox/meerkats/meerkats_shared/data'

#select field season
years <- c(2019)

#setup of analysis variables
dist_thresh <- 5 #distance to nearest neighbor
response_time <- 1:6 #response time cutoff
sequence_lenght <- 3 #length of call sequence
seq_start = 10 #silence time before initiating a sequence

#set kernell variables
max.lag <- 10
hop.time <- 0.01
bandwidth <- 0.2
jitter <- F
jitter.time <- 10


# set working directory
setwd(codedir)
# source library
source("call_interactions_functions.R")
both_years <- data.frame()
 for (year in years){
   # load call data
   #calls.all <- read.delim(paste(datadir, "/HM" , year ,"_ALL_CALLS_SYNCHED.csv",sep=""))
   calls.all <- read.csv("V:/meerkat/working/processed/acoustic/resolve_conflicts/2021-05-17/2021-05-17_HM2019_ALL_CALLS_SYNCHED_conflicts_resolved.csv", sep = "\t")
   # recode call type to lowercase for consistency between years
   calls.all$callType <- tolower(calls.all$callType)
   
   # load movement data
   load("V:/meerkat/working/processed/movement/HM2019_COORDINATES_all_sessions.RData")
   
   # load individual info
   indInfo <- read.delim(paste(datadir, "/HM_INDIVIDUAL_INFO_", year,".txt", sep = ""))
   individuals <- unique(calls.all$ind)
#caller_list <- unique(calls.all$ind)
date_list <- as.character(unique(calls.all$date))

all_calls_seq <- data.frame() 


for (date in date_list) {
  
  #get current data
  calls <- calls.all[which(calls.all$date == date),]
  #select only the middle hour
  if (year == 2019)
  {calls <- calls[which(substr(calls$t0GPS_UTC,12,13)==13) , ]}
  
  # change times to numeric
  options(digits.secs=3)
  calls$t0.numeric <- as.numeric(as.POSIXlt(calls$t0GPS_UTC,tz="UTC", format = "%Y-%m-%d %H:%M:%OS"))  #NOTE: Looking at these times, they seem rounded to the second - this would be a big problem! LOOK INTO THIS -- FIX
  calls$tf.numeric <- as.numeric(as.POSIXlt(calls$tEndGPS_UTC,tz="UTC", format  = "%Y-%m-%d %H:%M:%OS")) 
  
  # remove some miss labled  non calls
  calls[ which(calls$entryName == "x"),"isCall"] <- 0
  calls[ which(calls$entryName == "eating"),"isCall"] <- 0
  calls[ which(calls$entryName == "chew"),"isCall"] <- 0
  calls[ which(calls$entryName == "digging"),"isCall"] <- 0
  calls[ which(calls$entryName == "?"),"isCall"] <- 0
  calls[ which(calls$entryName == "na"),"isCall"] <- 0
  calls[ which(calls$entryName == "//"),"isCall"] <- 0
  calls[ which(calls$callType == "skip"),"isCall"] <- 0
  calls[ which(calls$callType == "syn"),"isCall"] <- 0
  calls[ which(calls$callType == "sync"),"isCall"] <- 0
  
  calls[ which(grepl("Mar", calls$entryName)), "isCall"] <- 0
  calls[ which(grepl("bark", calls$entryName)), "isCall"] <- 0
  calls[ which(grepl("bird", calls$entryName)), "isCall"] <- 0
  calls[ which(grepl("#", calls$entryName)), "isCall"] <- 0
  calls[ which(grepl("be", calls$entryName)), "isCall"] <- 0
  
  calls[ which(grepl("1:", calls$entryName)), "isCall"] <- 0
  
  calls$entryName <- droplevels(calls$entryName)
  
  # remove nonfocal calls and non calls
  calls <- calls[which(calls$pred_focalType == "F" &  calls$isCall == 1), ]
  calls <- calls[order(calls$t0.numeric ),]
  
  #get caller distance
  calls$c_dist <- NA
  for(i in 2:nrow(calls)){
    dist <- sqrt((calls$x_emitted[i]-calls$x_emitted[i-1])^2 + (calls$y_emitted[i]-calls$y_emitted[i-1])^2 ) 
    if (!is.na(dist)) {calls$c_dist[i] <- dist}
  }
  
  #duplicate all calls data frame
  all_calls <- calls
  all_calls$lag <- NA
  for(i in 2:nrow(all_calls)){
    all_calls$lag[i] <- abs(all_calls$t0.numeric[i-1] - all_calls$t0.numeric[i])}
  
  
  all_calls$seq <- NA
  seq_start_points <- which(all_calls$lag > seq_start)+1
  
  tseq <- 1
  for (r in seq_start_points){
    l <- r+1
    if (l >= nrow(all_calls)-1) {next}
    all_calls$seq[l-1] = tseq
    while(all_calls$lag[l-1] < 5) {
      all_calls$seq[l] <- tseq
      l <- l+1}
    
    tseq <- tseq + 1
  }
  
  
  #####recode some calls due to between year mismatches######
  
  all_calls$entryName<-  droplevels(all_calls$entryName)
  all_calls <- all_calls[c(names(all_calls)[-1], "entryName")]
  
  
  all_calls[which(all_calls$callType == "s"), "callType"] <- "sn"  
  all_calls[which(all_calls$callType == "social"), "callType"] <- "soc"  
  all_calls[which(all_calls$callType == "c"), "callType"] <- "cc"  
  all_calls[which(all_calls$callType == "aggress"), "callType"] <- "agg"  
  all_calls[which(grepl("chat",all_calls$callType)), "callType"] <- "agg" 
  all_calls[which(grepl("ala",all_calls$callType)), "callType"] <- "al" 
  
 all_calls$ini <- NA
 all_calls$self <- NA
 for (i in 2:nrow(all_calls))
 {all_calls$ini[i] <- as.character(all_calls$callType[i-1])
 all_calls$self[i] <- as.character(all_calls$ind [i-1])}
  
  print(paste("end", date))
  
  all_calls_seq <- rbind(all_calls_seq, all_calls)
  
}


both_years <- smartbind(both_years, all_calls_seq)
 }

both_years$callType <- both_years$entryName ##### at this point we have all focal enteries together 


### recode some of the calls due to between years inconsistency ###
### this is the part that does the actual name recoding ###
#If we encounter more naming issues which are not resolved, the folowing lines could be modified

#here set the vectors for unifying modifier and call names
hybrid <- c("hybrid|hyb")
sequence <- c("seq|sq")
move <- c("move|mov")
agression <- c("aggress|agress|chat|growl|grunt|bark")
alarm <- c("alarm|alrm|alert|ala")
lost <- c("lost|loc|lc")

both_years$callType <- str_replace_all(both_years$callType, hybrid, "hyb")
both_years$callType <- str_replace_all(both_years$callType, sequence, "sq")
both_years$callType <- str_replace_all(both_years$callType, move, "mo")
both_years$callType <- str_replace_all(both_years$callType, c("lead"), "ld")
both_years$callType <- str_replace_all(both_years$callType, c("social"), "soc")
both_years$callType <- str_replace_all(both_years$callType, c("ukn"), "unk")
both_years$callType <- str_replace_all(both_years$callType, agression, "agg")
both_years$callType <- str_replace_all(both_years$callType, alarm, "al")
both_years$callType <- str_replace_all(both_years$callType, lost, "lc")
both_years[which(both_years$callType == "s"), "callType"] <- "sn"  
both_years[which(both_years$callType == "c"), "callType"] <- "cc"  


both_years$callType <- str_replace(both_years$callType, fixed("s+"), "sn ") #renaming s as the first element
both_years$callType <- str_replace(both_years$callType, "\\+s$", " sn ")   #renaming s as the second element

### get hyb calls ###

both_years$hyb <- NA
both_years[grepl("fu|hyb", both_years$callType), "hyb"] <- "fu"
both_years[grepl("sq", both_years$callType), "hyb"] <- "sq"



### get call type elements ####
call_types <- c("agg|al|cc|ld|mo|sn|soc|lc|unk") #call types that we need
both_years <- cbind(both_years , str_extract_all(both_years$callType, call_types, simplify = TRUE)) #getting the call type elements 
both_years$`1` <- as.character(both_years$`1`) #removing factors
both_years$`2` <- as.character(both_years$`2`) #removing factors

#collecting the call type elements in an alphabetic order
both_years$final <- ifelse(both_years$`1` < both_years$`2`, paste(both_years$`1`, both_years$`2`), paste(both_years$`2`, both_years$`1`)) 
# keeping the original order for sequential calls
both_years[which(both_years$hyb == "sq") , "final"] <- paste(both_years[which(both_years$hyb == "sq") , "1"], both_years[which(both_years$hyb == "sq") , "2"])


#looking at the frequencies of the main call types as a self check
call_types <- c("cc", "soc", "al", "agg", "sn", "ld", "mo", "lc", "unk")


freq <- data.frame()
for (i in 1:length (call_types))
{
  freq[i,1] <- sum(str_count(rbind(both_years$`1`, both_years$`2`) , pattern = call_types[i]), na.rm = T)
  freq[i,2] <- call_types[i]
}
freq<- freq[order(freq$V1, decreasing = F),] ### if we decide to recode by the rarest category decreasing is to be set to TRUE


#get sample size plot for sanity check
par(mfrow=c(1,2))
bp <- barplot(freq$V1, names.arg = freq$V2, main = "Sample sizes per call type") 
text(bp, 0, freq$V1 ,cex=1,pos=3) 
#recode the call-types into main call categories. The order of recoding is frequency based 
# hybrid call_types are collapsed to the more frequent type. 

both_years$type_group <- NA  

for ( i in 1:nrow(freq))
{ type <- freq$V2[i]

both_years[which(grepl(type, both_years$final)), "type_group"] <- type 
}


#sample sizes per call category

freq <- data.frame()
for (i in 1:length (call_types))
{
  freq[i,1] <- sum(str_count(both_years$type_group , pattern = call_types[i]), na.rm = T)
  freq[i,2] <- call_types[i]
}
freq<- freq[order(freq$V1, decreasing = F),] 

bp <- barplot(freq$V1, names.arg = freq$V2, main = "Sample sizes per call category") #get sample size plot for sanity check
text(bp, 0, freq$V1 ,cex=1,pos=3) 


#cleaning
both_years <- subset(both_years, select=-c(`1`, `2`, `3`))
names(both_years)[names(both_years) == "final"] <- "stn_call_type"



all_calls_seq <- both_years
all_calls_seq$pair <- NA
all_calls_seq$pair <- paste(lag(all_calls_seq$type_group, n=1), all_calls_seq$type_group) #end data arranging bit



####making call pair transition tables for core call types####

##check for reply lag times
#remove fast self replies
replieT <-  all_calls_seq[-(which( all_calls_seq$lag < 0.1 &  all_calls_seq$ind ==all_calls_seq$self)) , ]  #get rid of fast self replies
quantile(replieT$lag, na.rm = T, probs = c(0.25 ,0.5, 0.75, 0.90))    # get response lag IQR and 90th percentiles

#duplicate data for easier handling
all_pairs <- all_calls_seq


all_pairs[which(all_pairs$lag >= 10) , "trigger"] <- "trigger" #mark trigger calls after 10 sec of group silence

responce_idx <- which(all_pairs$trigger == "trigger")+1 #getting the index of potential responce calls

for (i in responce_idx){
  if (all_pairs [ i , "lag"] <= 5) {all_pairs[i , "trigger"] <- "responce"}  #getting responce calls
}


#all_pairs <- all_pairs [which(all_pairs$c_dist < 20), ] #filter by distance here
all_pairs <- all_pairs[-(which(all_pairs$lag < 0.1 &  all_pairs$ind == all_pairs$self)) , ] #remove fast self replies

nrow( all_pairs) #get sample sizes
table( all_pairs$type_group)

all_pairs <- all_pairs[which(all_pairs$lag <= 5), ] #remove slow replies

all_pairs <- all_pairs[which(!is.na(all_pairs$ini) & !is.na(all_pairs$type_group)), ] #remove any accidental NAs
#all_pairs <- all_pairs[which(all_pairs$trigger == "responce") ,]

allcalls <- c(as.character(all_pairs$ini ), as.character(all_pairs$type_group ))
all_pairs$lag <- as.numeric(as.character(all_pairs$lag))

all_pairs$null <- sample(all_pairs$type_group, size = nrow(all_pairs), replace = F)
all_pairs$rand_lag <- sample(all_pairs$lag, size = nrow(all_pairs), replace = F)

#only focal-focal pairs
all_focal_pairs <- all_pairs[which(all_pairs$ind == all_pairs$self), ]
#all call pairs
all_foc_NF_pairs <-  all_pairs
#only focal-non_focal pairs
all_NF_pairs <- all_pairs[which(all_pairs$ind != all_pairs$self), ]
#all call pairs no same pairs
inter_type <- all_pairs[which(all_pairs$type_group != all_pairs$ini) ,]


call_pairs_list <- list(all_NF_pairs, all_focal_pairs, all_foc_NF_pairs , inter_type)
names(call_pairs_list) <- c("caller exchange", "self reply", "all replies", "inter_type")

#make heat maps with real probability values and color code corresponding to null subtraction
p <- list()
for (x in 1:length(call_pairs_list))
{
  #make empty matrix with call types
  gbiMat <- matrix(0,nrow=8,ncol=8)
  row.names(gbiMat) <- c("sn",  "soc",  "cc",  "mo", "ld", "agg", "al", "unk")
  colnames(gbiMat) <- c("sn",  "soc",  "cc",  "mo", "ld", "agg", "al", "unk")
  
  #make empty matrix with call types for null
  gbiMat_null <- gbiMat
  
  all_pairs <- call_pairs_list[[x]]
  
  #transform the raw assosiation into  matrix
  for (j in 1:nrow(all_pairs))
  {
    gbiMat[which(row.names(gbiMat) == all_pairs$ini[j]), which(colnames(gbiMat)== all_pairs$type_group[j])] <- gbiMat[which(row.names(gbiMat) == all_pairs$ini[j]), which(colnames(gbiMat)== all_pairs$type_group[j])]+1
  }
  
  prob_matrix <- gbiMat/rowSums(gbiMat)
  if (x == 4) {prob_matrix <- sna::diag.remove(prob_matrix, remove.val=NA)} #for the intercall only remove diagonal
  melted_prob_matrix <- melt(prob_matrix) #melt table into long format
  
  
  colnames(melted_prob_matrix) <- c("initiation", "reply", "value")
  
  
  #transform the null call association into  matrix
  for (j in 1:nrow(all_pairs))
  {
    gbiMat_null[which(row.names(gbiMat_null) == all_pairs$ini[j]), which(colnames(gbiMat_null)== all_pairs$null[j])] <- gbiMat_null[which(row.names(gbiMat_null) == all_pairs$ini[j]), which(colnames(gbiMat_null)== all_pairs$null[j])]+1
  }
  
  
  null_prob_matrix <- gbiMat_null/rowSums(gbiMat_null)
  melted_net_null <- melt(null_prob_matrix) #melt table into long format
  colnames(melted_net_null) <- c("initiation", "reply", "value")
  
  
  #substracting data from null
  delta_probs <- prob_matrix - null_prob_matrix 
  if (x == 4) {delta_probs <- sna::diag.remove(delta_probs, remove.val=NA)} #for the intercall only remove diagonal
  
  melted_delta_probs <- melt(delta_probs) #melt table into long format
  colnames(melted_delta_probs) <- c("initiation", "reply", "value")
  
  p[[x]] <- ggplot(data = melted_delta_probs, aes(y=initiation, x=reply, fill=value)) + 
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", name="Null substracion" ) +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed() + ggtitle(names(call_pairs_list[x])) + geom_text(data = melted_prob_matrix ,aes(label = round(value, 2))) #get the values from the real prob matrix
  
  
}

grid.arrange(p[[1]], p[[2]],p[[3]], p[[4]] ,  nrow = 2, top = paste("Delta call_type transition matrix (20m max range) "))



### calcuate self reply time vs caller exchange reply time



###calculate proportions of self replies vs caller exchange per call type

tmp_pairs <- all_calls_seq[which(!is.na(all_calls_seq$ini) & !is.na(all_calls_seq$type_group)) ,]
tmp_pairs$caller_match <-ifelse( tmp_pairs$ind == tmp_pairs$self, T, F) 


call_pairs <- c("cc cc", "sn sn", "agg agg", "al al", "soc soc", "mo mo")
tmp_pairs <- tmp_pairs[which(tmp_pairs$pair %in% call_pairs), ] #select only same call type pairs
tmp_pairs <- tmp_pairs[-(which(tmp_pairs$lag < 0.1 & tmp_pairs$caller_match ==T)) , ] #remove quick self replies
tmp_pairs <- tmp_pairs[which(tmp_pairs$lag <= 5), ] #remove slow replies
tmp_pairs <- tmp_pairs [which(tmp_pairs$c_dist < 20), ] #filter by distance here
tmp_pairs$rand_lag <- sample(tmp_pairs$lag, nrow(tmp_pairs), replace = F)

#plot proportions of self and non-self reply




##### look at self reply probability by the number of neighbors ######
#collect only calls with neighbor data

pos_neigh <- tmp_pairs[which(!is.na(tmp_pairs$nbr_neigh_below_10m)) , ]

pos_neigh$self_prob <- 1/(pos_neigh$nbr_neigh_below_10m+1)
pos_neigh$neigh_prob <- pos_neigh$nbr_neigh_below_10m/(pos_neigh$nbr_neigh_below_10m+1)

#get self reply random probs for all call types
self_means <- data.frame()
for (type in call_pairs){
  self_means[1 , type] <- mean(pos_neigh[which(pos_neigh$pair == type) , "self_prob"])}

#get self reply rates for each call type
self_rates <- table(pos_neigh$caller_match, pos_neigh$pair)

self_rates <- self_rates[2 , ]/(self_rates[1 , ] + self_rates[2 , ])

self_rates <- c(self_rates[3], self_rates[5], self_rates[1], self_rates[2], self_rates[6], self_rates[4])


self_means <- rbind(self_means, self_rates)
self_means$condition <- c("rand_rate", "rate") #get values for self replies per call type (rate vs random probability)
self_means$condition <- as.factor(self_means$condition)
data_long <- reshape2::melt(self_means, id.vars= "condition") 

data_long$variable <- ordered(data_long$variable, levels = c("agg agg", "al al", "cc cc", "mo mo", "sn sn", "soc soc")) #get the factor into a proper order

# Grouped bar plot
p2 <- ggplot(data_long, aes(fill=condition, y=value, x=variable)) + 
  geom_bar(position="dodge", stat="identity") + xlab("call_pair")+ ylab("Self reply") + ggtitle("Mean self reply") + 
  labs(fill = "Self reply") + scale_fill_discrete(labels=c( "By chance", "Rate"))

# general proportion plot
p1 <- ggplot(tmp_pairs, aes(pair , fill=caller_match)) +
  geom_bar(position="fill", width=0.7) + xlab("call_pair")+ ylab("Proportion") + ggtitle("Self reply / caller exchange proportion") + labs(fill = "Self reply")  
#Plot all together
grid.arrange(p1, p2, nrow = 2)


### box plots #### reconsider this , perhaps not needed ####
par(mfrow=c(1,1))
b <- boxplot(lag ~ pair, data=tmp_pairs[which(tmp_pairs$caller_match == F) , ], plot=0) #get names and sample sizes
boxplot( lag ~ pair , data = tmp_pairs[which(tmp_pairs$caller_match == F) , ],main = "Response time lag for caller exchange dyads", ylab = "Time lag (sec)", xlab = "Call type pair",  names=paste(b$names, "(n=", b$n, ")"))  #plot

#####boxplot(tmp_pairs[which(tmp_pairs$caller_match == F),"lag"] , tmp_pairs[which(tmp_pairs$caller_match == T),"lag"],  names = c("non-self", "self"), main = "all_calls")
#### par(mfrow=c(2,3))
####
####  for (n in call_pairs)
#### {
#### #hist(tmp_pairs[which(tmp_pairs$caller_match == F & tmp_pairs$pair == n),"lag"], xlim=c(0,30), col="red", breaks = 100,freq=FALSE, main = paste(n,"calls trainsition response time"), xlab = "Time lag (sec)")
#### #hist(tmp_pairs[which(tmp_pairs$caller_match == T & tmp_pairs$pair == n),"lag"], add=T, col=rgb(0, 1, 0, 0.5),breaks = 100, freq=FALSE )
#### 
#### boxplot(tmp_pairs[which(tmp_pairs$caller_match == F & tmp_pairs$pair == n),"lag"] , tmp_pairs[which(tmp_pairs$caller_match == T & tmp_pairs$pair == n),"lag"],   names = c("non-self", "self"), main = n)
#### 
#### #hist(tmp_pairs[which(tmp_pairs$caller_match == F & tmp_pairs$pair == n),"rand_lag"], xlim=c(0,30), col="red", breaks = 100,freq=FALSE, main = paste(n, "calls trainsition null response time"), xlab = "Time lag (sec)")
#### #hist(tmp_pairs[which(tmp_pairs$caller_match == T & tmp_pairs$pair == n),"rand_lag"], add=T, col=rgb(0, 1, 0, 0.5),breaks = 100, freq=FALSE )
#### }
#### 



#####################################################################

test <- all_calls_seq[which(!is.na(all_calls_seq$c_dist)), ] #only get pairs with distance
test <- test[which(test$self != test$ind), ] #only get caller exchange
test <- test[which(test$c_dist < 20), ] #limit distance to 20m

test <- test[which(test$lag < 5), ] #limit time lag to 5sec


call_pairs <- c("cc cc", "sn sn", "agg agg", "al al", "soc soc", "mo mo")
test <- test[which(test$pair %in% call_pairs), ] #select only same call type pairs
table(test$pair)

#run plots and models for call types with reasonable sample sizes only

par(mfrow=c(1,1))

call_pairs <- c("cc cc"  , "sn sn" ,  "agg agg",   "soc soc")
colors <- c("#000000", "#E69F00", "#56B4E9", "#009E73")

plot(1, type="n", xlab=bquote(bold("Distance(m)")), ylab=bquote(bold("Response time (sec)")), xlim=c(0, 20), ylim=c(0, 2.5), main = "Response time by distance", bty = "l")
legend("topleft", 
       legend = call_pairs, 
       col = adjustcolor(colors, alpha.f=0.4), 
       pch = 20, 
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.05, 0.01))

for (call in call_pairs)
{
  clr <- colors[which(call_pairs == call)]
  test_plot <- test[which(test$pair == call), ]
  
  library(lme4)
  
  model<- lm(lag~c_dist, data = test_plot)
  
  summary(model)
  
  
  
  pred.data=data.frame(c_dist= seq(from=min(test_plot$c_dist),
                                   to=max(test_plot$c_dist),
                                   length.out=100))
  ci.plot=predict.lm(object=model, newdata=pred.data,
                     interval="confidence") 
  
  #plot(1, type="n", xlab=bquote(bold("Distance(m)")), ylab=bquote(bold("Response time (sec)")), xlim=c(min(test_plot$c_dist), max(test_plot$c_dist)), ylim=c(0, 3.5), main = call)
  
  polygon(
    x=c(pred.data$c_dist, rev(pred.data$c_dist)),
    y=c(ci.plot[, "lwr"], rev(ci.plot[, "upr"])),
    border=NA, col = adjustcolor(clr, alpha.f=0.2))
  lines(x=pred.data$c_dist, y=ci.plot[, "fit"], lwd=2, lty=2,col=clr )
  
}
par(mfrow=c(2,2))
 hist(test[which(test$pair == "cc cc") ,"c_dist"], col = adjustcolor("#000000", alpha.f=0.4), main = "Close Call response distance", xlab = "meters")
 hist(test[which(test$pair == "sn sn") ,"c_dist"], col = adjustcolor("#E69F00", alpha.f=0.4), main = "Short note response distance", xlab = "meters")
 hist(test[which(test$pair == "agg agg") ,"c_dist"],col = adjustcolor("#56B4E9", alpha.f=0.4), main = "Aression response distance", xlab = "meters")
 hist(test[which(test$pair == "soc soc") ,"c_dist"],col = adjustcolor("#009E73", alpha.f=0.4), main = "Social response distance", xlab = "meters")

###    #### this is an old version with null model generation #####
### ##loop over all call-type pairs and plot with  null models##
### 
### call_types <- c("cc", "sn", "agg", "al", "soc", "mo")
### g <- list()
### for (type in call_types)
### {
### same_type <- test[which(test$pair == paste(type, type)),]
### #same_type <- same_type[which(same_type$lag < 60), ]
### null <- data.frame(matrix(NA, ncol = 1000, nrow = nrow(same_type)))
### for (i in 1:1000){
###   null[i] <- sample(same_type$c_dist,nrow(same_type), replace = F)}
### null$lag <- same_type$lag
### #null_long <- tidyr::gather(null, condition, measurement, X1:X1000, factor_key=TRUE)
### 
### g[[type]] <- ggplot()   
### 
### for (i in 1:1000)
### { g[[type]] <- g[[type]] + stat_smooth(geom='line', data=null, aes_string(x = paste("X",i, sep = ""), y = "lag"), method = "lm", se = F,  col = "red", alpha = 0.01, size = 1) }
### 
### g[[type]]<- g[[type]] + geom_smooth(data=same_type, aes(x = c_dist, y = lag),  method = "lm", col = "blue") + 
###   #geom_point(data=same_type, aes(x = c_dist, y = lag)) + 
###   xlab("Distance (meters)") + ylab("Response time (sec)")  +
###   ggtitle(paste(type, "-", type, "reponse time vs distance")) 
### 
### }
### do.call(grid.arrange,g)





### ##### looking at acoustic structure of CC ini vs CC responce ######
### 
### CCCC <- all_calls_seq
### CCCC$trigger <- NA
### CCCC$responce <- NA
### CCCC[which(CCCC$type_group == "cc" & CCCC$lag >10), "trigger"] <- "trigger"  #find CC trigger calls produced after 10sec silence
### CCCC[which(CCCC$trigger == "trigger")+1 ,  "responce"] <- "responce"
### CCCC[which(CCCC$responce == "responce" & CCCC$type_group != "cc"), "responce"] <- NA
### CCCC[which(CCCC$responce == "responce" & CCCC$ind == CCCC$self), "responce"] <- "self"
### 
### 
### 
### CCCC[which(CCCC$fundamental.meanentire. < 0), "fundamental.meanentire."] <- NA
### 
### 
###
### metrix <- c("duration" ,"rms" , "fundamental.meanentire.", "entropy.meanentire.",  "hnr.meanentire.", "T_Jitter_rel.nf","T_Jitter_rel.f" ,  "Shimmer_rel.nf" , "Shimmer_rel.f")
###
###
###
###  CCCC_i <- CCCC
###  CCCC_i <- CCCC_i[which(CCCC_i$trigger == "trigger" | CCCC_i$responce == "responce" ), ]
###  CCCC_i <- CCCC_i[-which(CCCC_i$trigger == "trigger" & CCCC_i$responce == "responce"), ]
###  CCCC_i <- CCCC_i[which(CCCC_i$c_dist < 50), ]
###  
###  dist <- list()
###  met <- list()
###  for (m in metrix) {
###    dat <- scale(CCCC_i[, m])
###    
###    met[[m]] <- ggplot(data= CCCC_i[which(CCCC_i$trigger == "trigger"),]) + geom_boxplot(aes_string(x = "trigger", y = dat[which(CCCC_i$trigger == "trigger")])) + geom_boxplot( data= CCCC_i[which(CCCC_i$responce == "responce"),], aes_string(x = "responce", y = dat[which(CCCC_i$responce == "responce")])) +
###      ylab(m) + xlab("")
###   dist[[m]] <- ggplot(data= CCCC_i[which(CCCC_i$responce == "responce"),]) + geom_smooth(aes_string(x = "c_dist", y = dat[which(CCCC_i$responce == "responce")]))  +
###     geom_point(aes_string(x = "c_dist", y = dat[which(CCCC_i$responce == "responce")])) +  ylab(m) + xlab("distance")
###    
###  }
###  grid.arrange(grobs = met,  nrow = 3, top = paste("trigger n=",length(which(CCCC_i$trigger == "trigger")), "responce n=", length(which(CCCC_i$responce == "responce"))))
### grid.arrange(grobs = dist,  nrow = 3, top = paste("responce n=", length(which(CCCC_i$responce == "responce"))))
###      
###  
###      ###### looking at self reply timing #####
###     CCCC <- all_calls_seq
###
###alldata_self_replies <- data.frame(stringsAsFactors = F)
### for (date in dates) {
###   print(date)
###   curr <- CCCC[which(CCCC$date == date), ]
###   CC_trigger <- which(curr$type_group == "cc" & curr$lag > 10)
###    if (length (which(CC_trigger == nrow(curr))) != 0){CC_trigger <- CC_trigger[-which(CC_trigger == nrow(curr))]}
###   
###   self_replies <- data.frame(stringsAsFactors = F)
###   for (t in CC_trigger)
###   { 
###     r=t
###     id <- curr[r, "ind"]
###     while (curr[r+1, "ind"] != id | curr[r+1, "type_group"] != "cc") 
###     {r<- r+1
###     if (r == nrow(curr)) break }
###     next_self <- r+1
###     ifelse(next_self-t == 1, reply <- "self", reply <- "non_self")
###     
###    if (reply == "non_self") {
###    first<-  curr[(min(which(curr[(t+1):(next_self-1), "type_group"] == "cc")) + t), "t0.numeric"] - curr[t, "tf.numeric"]
###    reply_lag <- curr[next_self, "t0.numeric"] -  curr[(min(which(curr[(t+1):(next_self-1), "type_group"] == "cc")) + t), "tf.numeric"]
###    last_lag <-  curr[r+1, "t0.numeric"] -  curr[(max(which(curr[(t+1):(next_self-1), "type_group"] == "cc")) +t) , "tf.numeric"]
###    reply_count <- length(which(curr[(t+1):(next_self-1), "type_group"] == "cc"))
###     
###    }else{
###      
###     first <- "none" 
###     last <- "none"
###     reply_lag <- "none" 
###     last_lag <- "none"
###     reply_count <- "none" }
###     loc <- data.frame(date , as.character(id), curr$type_group[t] ,curr$t0.numeric[next_self] - curr$tf.numeric[t], next_self, t, first, reply_lag, last_lag, reply, reply_count, stringsAsFactors = F)
###     self_replies <- rbind(self_replies, loc)
###     }
###  alldata_self_replies <- rbind(alldata_self_replies, self_replies)
###  
### }
###
### colnames (alldata_self_replies) <- c("date", "ID", "Call_type", "self_lag" , "next_self_row", "trigger_row", "neighbour_first_reply_lag", "self_lag_after_first_reply", "self_lag_after_last_reply","reply_type", "#_of_non_self_reply_calls")
###
### ggplot(data = alldata_self_replies, aes(x = reply_type, y = self_lag)) + geom_boxplot()+ ylim(0, 500)
### 

###I have lost all track of this bit , it needs to be redesigned completely ##

 ###### looking at probability with  kernell function ########
 dates <- as.character(unique(all_calls_seq$date)) #get relevant dates
 individuals <- as.character(unique(all_calls_seq$ind))
 #empty lists to fill with data
 all_ind_self_kernels<-list()
 all_ind_non_self_kernels <- list()
 
 #choose individual
 for (main in individuals) {
   #empty lists to fill with data
   all_self_kernels<- list()
   all_non_self_kernels <- list()
   
   #choose date  
   for (date in dates) {
     
     print(date)
     print(main)
     
     
     CCCC <- all_calls_seq 
     #get calls from relevant date
     curr <- CCCC[which(CCCC$date == date), ]
     
     #identify self replies and non_self replies
     curr$self_trigger <- NA
     curr$main <- 0
     curr[which(curr$ind == main), "main"] = 1
     
     ###getting trigger tims for self and no self replies
     ##for (row in 2:(nrow(curr)-1)) {
     ##  
     ##   if (curr$type_group[row] == "cc" & curr$lag[row] > 5 & curr$ind[row] == main & curr$ind[row+1] == main & curr$type_group[row+1] == "cc" ) {curr$self_trigger[row] = 1} 
     ##   if (curr$type_group[row] == "cc" & curr$lag[row] > 5 & curr$ind[row] == main & curr$ind[row+1] != main & curr$type_group[row+1] == "cc" & curr$lag[row+1] < 5) {curr$self_trigger[row] = 0} 
     ##}
     
     #getting trigger tims for self and no self replies
     for (row in 2:(nrow(curr)-1)) {
       
       if (curr$type_group[row] == "cc" & curr$lag[row] > 5 & curr$ind[row] == main) {curr$self_trigger[row] = 1} 
       if (curr$type_group[row] == "cc" & curr$lag[row] > 5 & curr$ind[row] == main & curr$ind[row+1] != main & curr$type_group[row+1] == "cc" & curr$lag[row+1] < 5) {curr$self_trigger[row] = 0} 
     }
     
     
     
     
     
     #getting main call times and trigger times
     self_reply_trigger <-  curr[which(curr$self_trigger == 1), "tf.numeric"]
     non_self_reply_trigger <-  curr[which(curr$self_trigger == 0), "tf.numeric"]
     
     #get call times and remove trigger calls from it
     call_times <-  curr[which(curr$ind == main & curr$type_group == "cc"), "t0.numeric"] 
     
     if (length(self_reply_trigger) >0)
     {call_times <- call_times[-which(call_times %in% curr[which(curr$self_trigger == 1), "t0.numeric"])]} 
     if (length(non_self_reply_trigger) >0)
     {call_times <- call_times[- which(call_times %in% curr[which(curr$self_trigger == 0), "t0.numeric"])]}
     
     # get latest start time and earliest end time of all individuals in table
     start.times <- tapply(curr$t0.numeric, curr$ind, min, na.rm = T)
     latest.start <- max(start.times, na.rm = T)
     
     # get start and end times
     end.times <- tapply(curr$t0.numeric, curr$ind, max, na.rm = T)
     earliest.finish <- min(end.times, na.rm = T)
     
     #run the kernell function while skipping non existing conditions
     
     if (length(self_reply_trigger >0))
     { kernel_self <- kernel.call.rates(call.times = call_times, 
                                        trigger.times = self_reply_trigger, start = latest.start, 
                                        end = earliest.finish, max.lag = max.lag, hop.time = hop.time, kernel.size = bandwidth, 
                                        jitter = jitter, jitter.time = jitter.time)
     }else{
       kernel_self <- NULL }
     
     
     #run the kernell function for non-self replies
     if (length(non_self_reply_trigger >0))
     {kernel_non_self <- kernel.call.rates(call.times = call_times, 
                                           trigger.times = non_self_reply_trigger, start = latest.start, 
                                           end = earliest.finish, max.lag = max.lag, hop.time = hop.time, kernel.size = bandwidth, 
                                           jitter = jitter, jitter.time = jitter.time)
     
     }else{
       kernel_non_self <- NULL } 
     
     #collect all outputs
     all_self_kernels[[date]] <- kernel_self
     all_non_self_kernels[[date]] <- kernel_non_self
   }
   
   all_ind_self_kernels[[main]] <-  all_self_kernels
   all_ind_non_self_kernels[[main]] <- all_non_self_kernels
 }
 
 # a data frame to collect all calls from all days and individuals
 all_self_call_rates <- (data.frame(matrix(NA, nrow = 2001, ncol = 0)))
 all_nonself_call_rates <- (data.frame(matrix(NA, nrow = 2001, ncol = 0)))
 
 #summarise call probabilities
 
 for (i in 1:length(all_ind_self_kernels)) {
   test <- sapply(all_ind_self_kernels[[i]], colMeans)
   
   group_mean <- rowMeans(test)
   
   
   all_self_call_rates <- cbind(all_self_call_rates, group_mean)
 }
 
 
 for (i in 1:length(all_ind_non_self_kernels)) {
   nonself_test <- sapply(all_ind_non_self_kernels[[i]], colMeans)
   
   
   
   group_nonself_mean <- rowMeans(nonself_test)
   
   
   all_nonself_call_rates <- cbind( all_nonself_call_rates,group_nonself_mean)
 }
 
 #plot
 par(mfrow=c(1,1))
 
 group_mean <- apply(all_self_call_rates, 1, mean)
 group_nonself_mean <- apply(all_nonself_call_rates,1,mean)
 plot(seq(-max.lag, max.lag, hop.time), group_mean, type = "l",ylim = c(0,0.2), xlim = c(-1, max.lag), main = paste("HM_",year))
 lines(x =seq(-max.lag, max.lag, hop.time),  y = group_nonself_mean, col = "blue")
 abline(v = 0, col = "red")
 
 
 #### windowed call rate of non focal responses bs no non focal responses #######
 
 plots<- list()
 dist <- c() 
 par(mfrow=c(2,2))
 date_list <- as.character(unique(all_calls_seq$date)) #get relevant dates
 distances <- seq(5, 20, by = 5) #the distance cutoff is controled here
 for (dist_thresh in distances){
   #empty lists to fill with data
   all_ind_foc_reply<-list()
   all_ind_non_foc_reply <- list()
   
   #some settings
   window_bins <- seq(from = 0, to = 10, by = 0.5) #setting the bins for counting calls
   bin_duration <- window_bins[2]-window_bins[1] #getting bin lag
   bin_count <- length(window_bins) #counting bins
   
   
   #choose individual
   for (main in individuals) {
     #empty lists to fill with data
     all_foc_reply<- list()
     all_non_foc_reply <- list()
     
     #choose date  
     for (date in date_list) {
       
       print(date)
       print(main)
       
       #get calls from relevant date
       curr <- CCCC[which(CCCC$date == date), ]
       
       
       #finding trigger CCs by the focal given after 10 sec of silence and having position data
       test_trigger_x <-curr[which(curr$ind == main & !is.na(curr$x_emitted)), "t0.numeric"] #find times of CC calls of main ind with position data
       
      trigger_times <- c()
        for (i in test_trigger_x) {
          
         calls_b_4 <- curr[which(curr$tf.numeric < i & curr$tf.numeric > i-10 ),] #get all calls in 10 sec range
         
         if(nrow(calls_b_4) == 0) {trigger_times = c(trigger_times , i)} #if no calls add time as trigger
        else{  #find dustances to calls
           distance_2_trigger <- c()
           for(r in 1:nrow(calls_b_4)){
             distance_2_trigger[r] <- sqrt((curr[which(curr$t0.numeric ==i),"x_emitted"] -calls_b_4$x_emitted[r])^2 + (curr[which(curr$t0.numeric ==i),"y_emitted"]-calls_b_4$y_emitted[r])^2 ) 
           }
          
          
         if (min(distance_2_trigger, na.rm = T) > 30) { trigger_times <- c(trigger_times , i)} #if all calls further than 30m dd to trigger
         }
          
        }
      
       #test_trigger <- curr[which(curr$type_group == "cc" & curr$lag > seq_start & curr$ind == main & !is.na(curr$x_emitted)), "tf.numeric"]
       test_trigger <- curr[which(curr$t0.numeric %in% trigger_times) , "tf.numeric"]
       
       #if no suitable trigger calls go to next day
       if (length(test_trigger) == 0) {next}
       
       
       #making empty matrices to collect focal call events under two conditions
       non_foc <- matrix(NA, length(test_trigger), bin_count) #non focal is calling
       foc <- matrix(NA, length(test_trigger), bin_count)     #non focal is not calling
       
       for (call in 1:length(test_trigger)) {
         for (bin in 1:length(window_bins)) {
           
           t0 <- test_trigger[call] + window_bins[bin] #setting start time of a bin relative to the spesific call
           
           tf <- t0+bin_duration                       #setting end time of the bin
           no_dist <- F
           
           #selecting the relevant calls within the bin
           wind <- curr[which(curr$t0.numeric > t0 & curr$t0.numeric < tf & curr$type_group == "cc"), ]
           
           #if non focal calls are given starting from the trigger call and up to the end of current bin  
           if (any(curr[which(curr$t0.numeric >= test_trigger[call] & curr$t0.numeric < tf & curr$type_group == "cc" ), "ind"] != main)){
             trigger_type <- "non_focal"
           } 
           else {
             trigger_type <- "focal"
           }
           
           if (trigger_type == "non_focal") {
             
             ##calculate distance between the focal and non_focal 
             fc <- curr[which(curr$tf.numeric == test_trigger[call]) , ]
             nfc<- curr[which(curr$t0.numeric >= test_trigger[call] & curr$t0.numeric < tf & curr$type_group == "cc" & curr$ind != main ), ] 
             
             #if no distance could be measured
             if (all(is.na(nfc$x_emitted))) {
               no_dist <- T
             }
             
             #calculate distanses between focal and non focal   
             dist <- c()  
             for (i in 1:nrow(nfc)) {
               dist[i] <-  sqrt((fc$x_emitted[1]-nfc$x_emitted[i])^2 + (fc$y_emitted[1]-nfc$y_emitted[i])^2 ) 
             } 
           }
           
           #if there is a valid distance measurment mark relevant cell in non_focal matrix if not change the trigger type to "focal"
           ifelse (no_dist != T & length(which(dist < dist_thresh)) > 0 & trigger_type == "non_focal", non_foc[call, bin] <- length(which(wind$ind == main)), trigger_type <-  "focal") #mark the relevant bin in the non_focal matrix
           
           #if trigger type is focal and there are is valid distace calculation mark as "focal"
           if (trigger_type == "focal" & no_dist != T)  {
             foc[call,bin] <- length(which(wind$ind == main)) 
           }  #mark the relevant bin in the focal matrix
         }
       }
       #collect all outputs into lists
       all_foc_reply[[date]] <- foc
       all_non_foc_reply[[date]] <- non_foc
     }
     
     all_ind_foc_reply[[main]] <-  all_foc_reply
     all_ind_non_foc_reply[[main]] <-  all_non_foc_reply
   }
   
   # a data frame to collect all calls from all days and individuals
   all_foc_rates <- (data.frame(matrix(NA, nrow = bin_count, ncol = 0)))
   all_nonfoc_rates <- (data.frame(matrix(NA, nrow = bin_count, ncol = 0)))
   
   #summarise call probabilities ### not sure if I am making a mess here????? - Old version
   
   # for (i in 1:length(all_ind_foc_reply)) {
   #   test <- sapply(all_ind_foc_reply[[i]],   colMeans, na.rm = T)
   #   if (length(test) == 0) next
   #   group_mean <- rowMeans(test, na.rm = T)
   #    
   #   
   #   all_foc_rates <- cbind(all_foc_rates, group_mean)
   # }
   # 
   # 
   # for (i in 1:length(all_ind_non_foc_reply)) {
   #   nonfoc_test <- sapply(all_ind_non_foc_reply[[i]], colMeans, na.rm = T)
   #   if (length(nonfoc_test) == 0) next
   #   
   #   
   #   group_nonfoc_mean <- rowMeans(nonfoc_test, na.rm = T)
   #   
   #   
   #   all_nonfoc_rates <- cbind(all_nonfoc_rates, group_nonfoc_mean)
   # }
   
   #get a big matrix of all trigger calls, also make a data frame that matches its rows with metadata (individual + date)
   nonfoc_all <- foc_all <- matrix(NA, nrow = 0, ncol = bin_count)
   metadata_nonfoc <- metadata_foc <- data.frame()
   for(i in 1:length(all_ind_non_foc_reply)){
     for(j in 1:length(all_ind_non_foc_reply[[i]])){
       if(length(all_ind_non_foc_reply[[i]])>0){
         
         curr_mat <- all_ind_non_foc_reply[[i]][[j]]
         nonfoc_all <- rbind(nonfoc_all, curr_mat)
         
         #get metadata for those rows
         metadat_rows <- data.frame(ind = rep(individuals[i], nrow(curr_mat)), date = rep(date_list[j], nrow(curr_mat)))
         metadata_nonfoc <- rbind(metadata_nonfoc, metadat_rows)
       }
     }
   }
   
   for(i in 1:length(all_ind_foc_reply)){
     for(j in 1:length(all_ind_foc_reply[[i]])){
       if(length(all_ind_foc_reply[[i]])>0){
         
         curr_mat <- all_ind_foc_reply[[i]][[j]]
         foc_all <- rbind(foc_all, curr_mat)
         
         #get metadata for those rows
         metadat_rows <- data.frame(ind = rep(individuals[i], nrow(curr_mat)), date = rep(date_list[j], nrow(curr_mat)))
         metadata_foc <- rbind(metadata_foc, metadat_rows)
       }
     }
   }
   
   #now nonfoc_all and foc_all are matrices where each row is a trigger call and each column is a bin
   #the rows in metadata_foc and metadata_nonfoc match with the rows o nonfoc_all and foc_all
   
   #par(mfrow=c(1,1))
   
   foc_mean <- colMeans(foc_all, na.rm = T)
   nonfoc_mean <- colMeans(nonfoc_all, na.rm = T)
   plot(x = window_bins, y = foc_mean, type = "l",ylim = c(0,0.1), xlim = c(0, max(window_bins)), main = dist_thresh)
   lines(x = window_bins,  y = nonfoc_mean, col = "blue")
   
###   
###  ##### plotting the output with CIs in ggplot ####
###  #all_foc_rates <- cbind(all_foc_rates, window_bins) #add bins as a grouping factors
###  #all_foc_rates[ ,"window_bins"] <- factor(all_foc_rates[ ,"window_bins"])
###  #colnames(all_foc_rates) <- c(1:(ncol(all_foc_rates)-1), "bin")
###  #melted_all_focal_rates <- melt(all_foc_rates, id.vars="bin") #transpose to a long format
###  #
###  #all_nonfoc_rates <- cbind(all_nonfoc_rates, window_bins)
###  #all_nonfoc_rates[ ,"window_bins"] <- factor(all_nonfoc_rates[ ,"window_bins"])
###  #colnames(all_nonfoc_rates) <- c(1:(ncol(all_nonfoc_rates)-1), "bin")
###  #melted_all_nonfoc_rates <- melt(all_nonfoc_rates, id.vars="bin")
###  #
###  ##plot
###  #plots[[which(distances == dist_thresh)]] <- ggplot(data = melted_all_focal_rates, aes(x = as.numeric(as.character(bin)), y = as.numeric(value)) ) + geom_smooth(method = "loess", se = F)  + geom_jitter(col = "blue") +
###  #  geom_smooth(data = melted_all_nonfoc_rates, aes(x = as.numeric(as.character(bin)), y = as.numeric(value)), method = "loess", col = "red", se = F) + geom_jitter(data = melted_all_nonfoc_rates, aes(x = as.numeric(as.character(bin)), y = as.numeric(value)), col = "red") +
###  #  xlab("sec") + ylab("probability") + theme_bw()+ ggtitle(dist_thresh)
###  #
}
### 
### 
### 
### #grid.arrange(grobs = plots,  nrow = 3, top = paste("distance cutoff (meters)"))
### 
###  