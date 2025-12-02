library(ggplot2)
library(dplyr)
library("ggsignif")

gene_categories <- c("VAR", "RIF", "stevor")

datasets <- c("dataset_1", "dataset_2", "dataset_3", "dataset_4")

dataset_4 <- transform(dataset_4,
                      Category = ifelse(name %in% gene_categories, name, "Others"))


combined <- transform(combined,
                      SurfaceProtein = ifelse(is.element(set = surfaceproteins$V4, el = combined$id), "Surface", "Not surface"))
#Individual curves

combined %>%
  filter(dataset == 4) %>%
    ggplot( mapping = aes(x = category, y = mean)) +
    geom_violin() +
    labs(x = "Gene Category", y = "Mean log2(U-enrichment)")+
    coord_cartesian(ylim = c(-1.1,1.1))+
    stat_signif(comparisons = list(c("RIF", "Others"), c("stevor", "Others"), c("VAR", "Others")),
                test = "t.test", step_increase = 0.05)+
    theme_minimal()

# Genewise comparison within each dataset

ggplot(data = combined, mapping = aes(x = category, y = mean)) +
  geom_violin(aes(colour = category)) +
  stat_signif(comparisons = list(c("RIF", "Others"), c("stevor", "Others"), c("VAR", "Others")),
              test = "t.test", step_increase = 0.05, map_signif_level = TRUE)+
  labs(y = "Mean log2(U-enrichment)")+
  coord_cartesian(ylim = c(-0.6,1.1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 1))+
  theme_minimal()+
  scale_color_manual(values = c("red4","blue4","cyan4", "green4"))+
  facet_wrap(vars(dataset))

ggsave("UDNAenrichment_violinplot_sign.eps", dpi = 300, width = 10, height = 15)

# Comparison of the same group of genes across the datasets

ggplot(data = combined, mapping = aes(x = dataset, y = mean)) +
  geom_violin(aes(colour = dataset)) +
  stat_signif(comparisons = list(c(1,3),c(2,3),c(3,4)),
              test = "t.test", step_increase = 0.12, map_signif_level = c('*'=0.001), )+
  labs(y = "Mean log2(U-enrichment)")+
  coord_cartesian(ylim = c(-1.0,1.5))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 1))+
  theme_minimal()+
  scale_color_manual(values = c("red4","blue4","cyan4", "green4"))+
  facet_wrap(vars(category))

Dataset <- c(1,2,3)
Category <- c("RIF", "VAR", "stevor", "Others")

remove(dataset)

pairwise.t.test(dataset_2$mean, dataset_2$category, p.adjust.method = "BH")
pairwise.t.test(dataset_3$mean, dataset_3$category, p.adjust.method = "BH")
pairwise.t.test(dataset_4$mean, dataset_4$category, p.adjust.method = "BH")

pairwise.t.test(combined$mean, combined$dataset, pool.sd = FALSE,
                          subset = combined$category == 'VAR', 
                          p.adjust.method = "BH", alternative = "two.sided")
