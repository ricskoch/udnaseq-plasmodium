## Import and clean, calculate new data
library(tidyverse)

g4_distance_merged_0p5 <- read.csv("nooverlap_distance_merged_10bp.bed", sep = "\t",
                                   header = FALSE)
colnames(g4_distance_merged_0p5) <- c("chr_U", "start_U", "end_U", "log2U", 
                                      "chr_G4", "start_G4", "end_G4", "G4H", "distance")
IP.vs.INPUT.merged.log2.bin10bp <-  read.csv("IP.vs.INPUT.merged.log2.bin10bp.bedGraph", sep = "\t",
                                             header = FALSE)
colnames(IP.vs.INPUT.merged.log2.bin10bp) <- c("chr", "start", "end", "log2U")

merged_above_0p5_bin10bp <- subset(IP.vs.INPUT.merged.log2.bin10bp,
                                   IP.vs.INPUT.merged.log2.bin10bp[,4] > 0.5)

write.table(x = merged_above_0p5_bin10bp, file = "../../../Merged/bin10bp/merged_above_0p5_bin10bp.bed", 
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Calculate G4_size
g4_distance_merged_0p5 <- mutate(g4_distance_merged_0p5, G4_size = end_G4 - start_G4)

# Categorization based on G4Hunter threshold score for later visualization
g4_distance_merged_0p5 <- mutate(g4_distance_merged_0p5, G4H_abs = abs(G4H))

# Filtering of the data based on categories in the Ana article
g4_distance_merged_0p5 <- mutate(g4_distance_merged_0p5, 
                                 filter = ifelse(G4H_abs >= 1.75, "1.75",
                                                 ifelse(G4H_abs >= 1.5, "1.5", 
                                                        ifelse(G4H_abs >= 1.2, "1.2", "Unfiltered"))))

g4_distance_merged_0p5 %>% 
  ggplot()+
  geom_boxplot(mapping = aes(y = G4_size, color = filter), show.legend = FALSE)+
  scale_y_log10()+
  theme_bw()+
  facet_wrap(~filter, nrow = 1)

# Distribution of log2(U-enrichment) in the blacklisted sample
u_enrichment_distribution <- IP.vs.INPUT.merged.log2.bin10bp %>% 
  filter(chr != "Pf3D7_MIT_v3" & chr != "Pf3D7_API_v3") %>% 
  ggplot(mapping = aes(x = log2U))+
  geom_histogram(mapping = aes(fill = chr, color = "black"), 
                 binwidth = 0.05, show.legend = FALSE)+
  scale_fill_discrete()+
  scale_color_manual(values = c("black"))+
  geom_vline(xintercept = 0, color = "red", linetype = "dotted", linewidth = 0.5)+
  coord_cartesian(xlim = c(-1, 1))+
  theme_bw(base_size = 10)+
  xlab("log2(U-enrichment)")+
  ylab("Count")+
  facet_wrap(~chr, nrow = 2)+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"))
u_enrichment_distribution

# Descriptive analysis of G4 statistics
summary(nooverlap_distance_merged_10bp$G4_size)

g4_distance_merged_0p5 %>% 
  ggplot()+
  geom_boxplot(mapping = aes(y = G4_size, color = filter), show.legend = FALSE)+
  scale_y_log10()+
  theme_bw()+
  facet_wrap(~filter, nrow = 1)

ggplot(g4_distance_merged_0p5, mapping = aes(x = G4_size)) +
  stat_ecdf(geom = "step", color = "blue") +
  stat_density(mapping = aes(y = ..scaled..), geom = "line", color = "red")+
  theme_bw() +
  xlab("Predicted G4 size range (bp)") +
  ylab("Cumulative Probability")

# Filtering for a size as well, because some large outliers, based on summary 
# it should be less than 200 bp, which contains around 80% of the data

g4_distance_merged_0p5 %>% 
  filter(G4H_abs > 1.2) %>% 
  filter(G4_size <= 200) %>%
  summarize(mean_distance = round(mean(distance), digits = 2), 
            median_distance = median(distance),
            range_distance = range(distance))
g4_distance_merged_0p5 %>% 
  filter(G4H_abs > 1.2) %>% 
  filter(G4_size <= 200) %>%
summary(distance)

g4_distance_uenrichment <- g4_distance_merged_0p5 %>% 
  filter(G4H_abs > 1.2) %>% 
  filter(G4_size <= 200) %>% 
  ggplot(aes(x = distance))+
  geom_histogram(fill = "grey50", binwidth = 30, 
                 color = "white", linewidth = 0.1)+
  coord_cartesian(xlim = c(-2000, 2000))+
  theme_bw()+
  xlab("Distance to closest G4-coordinates (bp)")+
  ylab("Count")+
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
g4_distance_uenrichment

g4_distance_merged_0p5 %>% 
  distinct(G4_size, .keep_all = TRUE) %>% 
  ggplot()+
  geom_point(mapping = aes(x = G4_size, y = G4H_abs))+
  coord_cartesian(ylim = c(0.5,3))+
  xlab("G4Hunter prediction score")+
  ylab("G4 size (bp)")

#G4 prediction score distribution 
g4_prediction_hist <- ggplot(nooverlap_distance_merged_10bp)+
  geom_density(mapping = aes(x = V8, 
               fill = ifelse(V8 > 0, "green", "red"), color = "grey50"), show.legend = FALSE)+
  xlab("G4 prediction score")+
  ylab("Density")+
  scale_fill_manual(values = c("green2", "red2"))+
  scale_color_manual(values = c("grey40"))+
  theme_bw()+
  scale_x_continuous(breaks = seq(-3, 3, by = 1),
                     minor_breaks = seq(-3, 3, by = 0.5))
g4_prediction_hist

ggsave("u_enrichment_distribution.eps", u_enrichment_distribution, device = 'eps', dpi = 300, width = 12, height = 4)


