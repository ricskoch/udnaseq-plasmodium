library(ggplot2)
library(dplyr)

GCcontent <- ggplot(data = merged_blacklisted_GCcontent, mapping = aes(x = V5, y = V4))+
  geom_jitter(mapping = aes(color = case_when(
                              V5 > 0.4 & V4 > 0 ~ "green2",
                              V5 > 0.4 & V4 < 0 ~ "red2",
                              V4 > 0 ~ "green4",
                              V4 < 0 ~ "red4",
                              TRUE ~ "black"
  )), show.legend = FALSE)+
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x)+
  geom_text(aes(x = max(V5), y = max(V4), label = paste("Correlation =", round(cor(V5, V4), 2))),
            hjust = 1, vjust = 1) +
  scale_color_identity()+
  theme_minimal()+
  coord_cartesian(ylim=c(-0.2,1))+
  labs(x = "GC%", y = "log2(U-enrichment)")
GCcontent
ggsave("GCcontent.svg")

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow2", "#FF7F00", "red", "#7F0000"))


GCcontent <- ggplot(data = merged_blacklisted_GCcontent, mapping = aes(x = V5 * 100, y = ifelse(V4 > -0.2, V4, NA), color = V4))+
  geom_jitter(show.legend = FALSE)+
  scale_color_gradientn(colours = jet.colors(9),
                        limits = c(-0.8,0.5), na.value = '#7F0000')+
  coord_cartesian(ylim=c(-0.1,1), xlim = c(0,80))+
  labs(x = "GC%", y = "log2(U-enrichment)")+
  theme(text = element_text(family = 'arial'),
        panel.border = element_rect(color = 'black', fill = NA),
        panel.grid = element_line(color = 'grey90'),
        panel.background = element_rect(fill = 'white'),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  geom_hline(yintercept = 0, color = 'grey50')+
  geom_vline(xintercept = 0, color = 'grey50')
GCcontent

ggsave("GCcontent_jet.png", width = 10, height = 8, dpi = 300, device = "png")

# Correlation analysis using imported table columns

cor(merged_blacklisted_GCcontent$V4, merged_blacklisted_GCcontent$V5, method = 'pearson')
cor(merged_blacklisted_GCcontent$V4, merged_blacklisted_GCcontent$V5, method = 'spearman')

# Plotting the dataset

plot(merged_blacklisted_GCcontent$V5, merged_blacklisted_GCcontent$V4, 
     xlim = c(0, 0.8), 
     ylim = c(-2, 2), 
     xlab = 'GC%', 
     ylab = 'log2 U enrichment', 
     main = 'Correlation between GC% and U-enrichment')

# Draw correlation line

abline(lm(merged_blacklisted_GCcontent$V4 ~ merged_blacklisted_GCcontent$V5), 
       col = 'red', 
       lwd = 3)
