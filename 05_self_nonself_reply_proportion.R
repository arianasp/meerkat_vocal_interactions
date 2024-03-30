library(tidyverse)

# Analyses of Meerkat call pair dynamics
# finding the proportion of self vs non self reply

load("05_self_nonself_reply_proportions_anon.RData")

#find intercall interval destribution
ICI <- quantile(both_years$lag,  probs = seq(.1, .9, by = .1), na.rm = T)[9]

#get the pairs for all calls
all_calls_seq <- both_years
all_calls_seq$pair <- NA
all_calls_seq$pair <-
  paste(lag(all_calls_seq$type_group, n = 1), all_calls_seq$type_group)
for (i in 2:nrow(all_calls_seq)) {
  all_calls_seq[i, "ini"] <- all_calls_seq[i - 1, "type_group"]
}
length(unique(all_calls_seq$date))
table(all_calls_seq$ind, all_calls_seq$date)

all_calls_seq$self <- lag(all_calls_seq$ind)
#remove first call of each day
all_calls_seq[which(is.na(all_calls_seq$lag)), "ini"] <- NA


###calculate proportions of self replies vs caller exchange per call type

tmp_pairs <-
  all_calls_seq[which(!is.na(all_calls_seq$ini) &
                        !is.na(all_calls_seq$type_group)) , ]
tmp_pairs$caller_match <-
  ifelse(tmp_pairs$ind == tmp_pairs$self, T, F)

tmp_pairs <-
  tmp_pairs[-(which(tmp_pairs$lag < 0.1 &
                      tmp_pairs$caller_match == T)) ,] #remove quick self replies
#remove slow replies
tmp_pairs <-
  tmp_pairs[which(tmp_pairs$lag <= ICI)  ,]


#get sample sizes
sample_sizes <- data.frame(table(tmp_pairs$pair))

#get call pairs of reasonable sample size
#call_pairs <-
#sample_sizes[which(sample_sizes$Freq > pair_sample_size), 1]


call_pairs <-
  c("agg agg", "al al", "cc cc", "sn sn", "soc soc")
tmp_pairs <-
  tmp_pairs[which(tmp_pairs$pair %in% call_pairs),] #select  call type pairs of interest


# general proportion plot self vs non self reply
p1 <- ggplot(tmp_pairs, aes(pair , fill = caller_match)) +
  geom_bar(position = "fill", width = 0.7) + xlab("call type pair") + ylab("Proportion") +
  geom_text(
    aes(label=signif(..count.. / tapply(..count.., ..x.., sum)[as.character(..x..)], digits=2)),
    stat="count", size=7,
    position=position_fill(vjust=0.5))+
  labs(fill = "self-reply") + scale_fill_manual(labels=c('non-self-reply', 'self-reply'), values=c("grey70", "grey40"), name="") +
  theme_minimal(base_size = 25) + theme( legend.position='top')  + scale_y_continuous(breaks=c(0,0.5, 1))
 

p1
