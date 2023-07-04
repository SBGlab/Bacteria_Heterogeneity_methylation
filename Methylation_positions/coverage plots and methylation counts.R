library(ggplot2)
library(gridExtra)
library(grid)
setwd("E:/Universidad/MBB/TFM/bisulfite-seq")

# Read the data from the file
data1 <- read.table(file="SRR15242564_deduplicated_filtered_coverage", header = FALSE)
data2<- read.table(file="SRR15242565_deduplicated_filtered_coverage", header = FALSE)
data3 <- read.table(file="SRR15242568_deduplicated_filtered_coverage", header = FALSE)
data4 <- read.table(file="SRR15242569_deduplicated_filtered_coverage", header = FALSE)

# Rename the columns
colnames(data1) <- c("Chromosome", "Position", "Coverage")
interval_size <- 10000


# Aggregate the coverage by interval
aggregated_data1 <- data1 %>%
  mutate(Position = (Position - 1) %/% interval_size * interval_size) %>%
  group_by(Chromosome, Position) %>%
  summarize(Coverage = mean(Coverage))

# Plot the aggregated coverage
p1=ggplot(aggregated_data1, aes(x = Position, y = Coverage)) +
  geom_line() +
  geom_area(fill = "skyblue", color = "steelblue", alpha = 0.5)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.title.position = "plot") +
  ggtitle("SRR15242564")+
  xlab(NULL) +
  ylab(NULL)+
  scale_x_continuous(labels = scales::comma_format(scale = 1e-6, big.mark = ",", accuracy = 1))+
  coord_cartesian(ylim = c(min(aggregated_data1$Coverage), NA))


colnames(data2) <- c("Chromosome", "Position", "Coverage")
aggregated_data2 <- data2 %>%
  mutate(Position = (Position - 1) %/% interval_size * interval_size) %>%
  group_by(Chromosome, Position) %>%
  summarize(Coverage = mean(Coverage))
p2=ggplot(aggregated_data2, aes(x = Position, y = Coverage)) +
  geom_line() +
  geom_area(fill = "skyblue", color = "steelblue", alpha = 0.5)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.title.position = "plot") +
  ggtitle("SRR15242565")+
  xlab(NULL) +
  ylab(NULL)+
  scale_x_continuous(labels = scales::comma_format(scale = 1e-6, big.mark = ",", accuracy = 1))


colnames(data3) <- c("Chromosome", "Position", "Coverage")
aggregated_data3 <- data3 %>%
  mutate(Position = (Position - 1) %/% interval_size * interval_size) %>%
  group_by(Chromosome, Position) %>%
  summarize(Coverage = mean(Coverage))

# Plot the aggregated coverage
p3=ggplot(aggregated_data3, aes(x = Position, y = Coverage)) +
  geom_line() +
  geom_area(fill = "skyblue", color = "steelblue", alpha = 0.5)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.title.position = "plot") +
  ggtitle("SRR15242568")+
  xlab(NULL) +
  ylab(NULL)+
  scale_x_continuous(labels = scales::comma_format(scale = 1e-6, big.mark = ",", accuracy = 1))+
  coord_cartesian(ylim = c(min(aggregated_data3$Coverage), NA))




colnames(data4) <- c("Chromosome", "Position", "Coverage")
aggregated_data4 <- data4 %>%
  mutate(Position = (Position - 1) %/% interval_size * interval_size) %>%
  group_by(Chromosome, Position) %>%
  summarize(Coverage = mean(Coverage))

# Plot the aggregated coverage
p4=ggplot(aggregated_data4, aes(x = Position, y = Coverage)) +
  geom_line() +
  geom_area(fill = "skyblue", color = "steelblue", alpha = 0.5)+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5),
        plot.title.position = "plot") +
  ggtitle("SRR15242569")+
  xlab(NULL) +
  ylab(NULL)+
  scale_x_continuous(labels = scales::comma_format(scale = 1e-6, big.mark = ",", accuracy = 1))


margins <- unit(1, "cm")  # Specify the margins in centimeters

arranged_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, top = "Aggregated Coverage along Chromosome (Interval Size: 1000)",bottom="Chromosome position (million bp)", left="Average coverage")


methylated_positions=data.frame(full=c(14659,953,2614,6850,11430,664,16492,21416,46866,15075,16709),partial=c(630,25,80,294,491,11,741,1045,494,1739,924))
rownames(methylated_positions)=c("0-0","1-1","2-2","3-3","4-4","5-5","6-6","7-7"," SRR1536433","BSeq_1","BSeq_2")
par(mar = c(5, 8, 4, 2))
par(oma=c(3,3,3,3))
par(cex.axis = 0.8)
barplot(methylated_positions$full,names.arg=rownames(methylated_positions),ylim=c(0,max(methylated_positions$full)+10000),ylab="",xlab="",main="Methylation postions
found per dataset",col=c("#2ca25f","#fc9272","#fc9272",
"#fc9272","#2ca25f","#fc9272","#2ca25f","#2ca25f","#2ca25f","#2ca25f","#2ca25f"),cex.axis=1.2,las=2,cex.names=1.2)
text(x = par("usr")[1] - 6, y = 30000, labels = "Number of methylations", pos = 1, cex = 1.2, xpd = TRUE, srt = 90)
abline(h = 10000, col = "red", lwd = 2)
barplot(methylated_positions$partial,names.arg=rownames(methylated_positions),ylim=c(0,max(methylated_positions$partial)+500),ylab="",xlab="",main="Partial methylation 
postions per dataset",col=c("#2ca25f","#fc9272","#fc9272","#fc9272","#2ca25f","#fc9272","#2ca25f","#2ca25f","#2ca25f","#2ca25f","#2ca25f"),las=2,cex.axis=1.2,las=2,cex.names=1.2)
text(x = par("usr")[1] - 5, y = 1200, labels = "Number of methylations", pos = 1, cex = 1.2, xpd = TRUE, srt = 90)
abline(h = 450, col = "red", lwd = 2)

         


