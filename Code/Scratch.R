# Scratch code 
summary_stats_WT <- as.data.frame(summary_stats_WT)
parameters <- rownames(summary_stats_WT[2:5,])
ggplot(summary_stats_WT, aes(x = rownames(summary_stats_WT[2:5], y = summary_stats_WT[2:5,4], fill=parameters)))