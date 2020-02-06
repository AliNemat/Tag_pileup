setwd("/Users/labadmin/Ali_Files/GalaxyProjects/bioinformatics/TagPileUp_development/R_scripts")
library(gridExtra)
library(ggplot2)

infile_1 <- paste("output_Ali_tagPileup_Chexmix_001_6th.txt",sep="")
infile_2 <- paste("output_Ali_tagPileup_Genetrack_005_6th.txt",sep="  ") 
df_1 <- read.table(infile_1, sep = "",header = FALSE)
df_2 <- read.table(infile_2, sep = "",header = FALSE)


bps_loc_1<-df_1$V1
tags_p_1<-df_1$V2
tags_n_1<-df_1$V3

bps_loc_2<- df_2$V1
tags_p_2<- df_2$V2
tags_n_2<- df_2$V3


bps_dist_1=abs (bps_loc_1)
bps_dist_2=abs (bps_loc_2)

avg_tags_dist_p_1=sum (bps_dist_1*tags_p_1)/sum (tags_p_1)
avg_tags_dist_n_1=sum (bps_dist_1*tags_n_1)/sum (tags_n_1)

avg_tags_dist_p_2=sum (bps_dist_2*tags_p_2)/sum (tags_p_2)
avg_tags_dist_n_2=sum (bps_dist_2*tags_n_2)/sum (tags_n_2)

sd_tags_dist_p_1= sqrt ( sum(as.numeric ( (bps_dist_1 - avg_tags_dist_p_1) * (bps_dist_1 - avg_tags_dist_p_1)*tags_p_1))/(sum (tags_p_1)-1) )
sd_tags_dist_n_1= sqrt ( sum(as.numeric ( (bps_dist_1 - avg_tags_dist_n_1) * (bps_dist_1 - avg_tags_dist_n_1)*tags_n_1))/(sum (tags_n_1)-1) )

sd_tags_dist_p_2= sqrt ( sum(as.numeric ( (bps_dist_2 - avg_tags_dist_p_2) * (bps_dist_2 - avg_tags_dist_p_2)*tags_p_2))/(sum (tags_p_2)-1) )
sd_tags_dist_n_2= sqrt ( sum(as.numeric ( (bps_dist_2 - avg_tags_dist_n_2) * (bps_dist_2 - avg_tags_dist_n_2)*tags_n_2))/(sum (tags_n_2)-1) )

tags_total_dist_1 = rep(0, max (bps_dist_1) +1 ) # + 1 is to include 0 distance.
tags_total_dist_2 = rep(0, max (bps_dist_2) +1 )

for (i in 1:length(tags_p_1)) {
  index=bps_dist_1[i]+1
  tags_total_dist_1[index] =  tags_total_dist_1[index] + tags_p_1[i] + tags_n_1[i]
}

for (i in 1:length(tags_p_2)) {
  index = bps_dist_2[i] + 1
  tags_total_dist_2[index] = tags_total_dist_2[index] + tags_p_2[i] + tags_n_2[i]
}

tags_prob_dist_1= tags_total_dist_1/(sum (tags_p_1) + sum (tags_n_1))
tags_prob_dist_2= tags_total_dist_2/(sum (tags_p_2) + sum (tags_n_2))


composites=cbind (bps_loc_1, tags_p_1, tags_n_1, bps_loc_2, tags_p_2, tags_n_2)
stats=cbind (0:max(bps_dist_1), tags_prob_dist_1, cumsum(tags_prob_dist_1), 0:max(bps_dist_2), tags_prob_dist_2, cumsum(tags_prob_dist_2))

composites_data=as.data.frame(composites)
stats_data=as.data.frame(stats)
pdf("ChexmixVSGenetrack_002.pdf",width=12, height=4, useDingbats=FALSE, bg="white")
p_1 = ggplot() + 
  geom_line(data = composites_data, aes(x = bps_loc_1, y = tags_p_1), color = "blue") +
  geom_line(data = composites_data, aes(x = bps_loc_1, y = -tags_n_1), color = "red") +
  geom_line(data = composites_data, aes(x = bps_loc_2, y = tags_p_2), color = "green") +
  geom_line(data = composites_data, aes(x = bps_loc_2, y = -tags_n_2), color = "black") +
  xlab('Distance from center of motif') +
  ylab('Number of tags')


p_2 = ggplot() +
  geom_line(data = stats_data, aes(x = V1, y = tags_prob_dist_1), color = "blue") +
  geom_line(data = stats_data, aes(x = V4, y = tags_prob_dist_2), color = "red") +
  xlab('Number of tags') +
  ylab('Probability') 
  

p_3 = ggplot() + 
  geom_line(data = stats_data, aes(x = V1, y = V3), color = "blue") +
  geom_line(data = stats_data, aes(x = V4, y = V6), color = "red") +
  xlab('Number of tags') +
  ylab('Cumulative probability') 
  
grid.arrange(p_1,p_2,p_3, ncol=3,top="Main Title")
dev.off()

