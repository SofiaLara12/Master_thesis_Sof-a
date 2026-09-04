##Data analyses

busco_df <- read_csv("busco_results_summary.csv")

busco_df <- busco_df %>% 
  rename(sample_name = Species) %>%
  select(sample_name, `Complete_C%`)

captus_df <- summary_df %>% select(sample_name, n_genes_recovered)

datos_qc <- inner_join(captus_df, busco_df, by = "sample_name")

correlacion <- cor.test(datos_qc$n_genes_recovered, datos_qc$`Complete_C%`, method = "spearman")

print(correlacion)

library(ggplot2)

ggplot(datos_qc, aes(x = n_genes_recovered, y = `Complete_C%`)) +
  geom_point(aes(color = n_genes_recovered < 500), alpha = 0.7, size = 3) + 
  geom_smooth(method = "lm", color = "black", linetype = "dashed") +
  labs(
    x = "Number of Loci Recovered (Captus)",
    y = "BUSCO Completeness (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")



library(tidyverse)

ggplot(datos_qc, aes(x = n_genes_recovered, y = `Complete_C%`)) +
  geom_point(aes(color = n_genes_recovered >= 500), alpha = 0.7, size = 5) + 
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    xlim(400, 650) +
    scale_color_manual(
    name = "Number of loci recovered",
    values = c("TRUE" = "darkgrey", "FALSE" = "black"),
    labels = c("TRUE" = "≥ 500", "FALSE" = "< 500")
  ) +
    annotate(
    "text", 
    x = 406,                                    
    y = 60,                                     
    label = "R = 0.619\np < 0.05",
    size = 9,                                 
    hjust = 0                                  
  ) +
    labs(
    x = "Number of Loci Recovered",
    y = "BUSCO Completeness (%)",
    size = 7
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 25, face = "bold"),
    legend.text = element_text(size = 25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "inside",
    legend.position.inside = c(0.22, 0.85),   
    legend.background = element_rect(fill = "white", color = "NA", linewidth = 0) 
)


library(tidyverse)

summary_df <- summary_df %>%
  mutate(Genero = word(sample_name, 1)) %>%
  mutate(Familia = if_else(Genero == "Cycas", "Cycadaceae", "Zamiaceae"))

library(car)

leveneTest(n_genes_recovered ~ Familia, data = summary_df)

model <- lm(n_genes_recovered ~ Familia, data = summary_df)
summary_df$residuals <- residuals(model)
shapiro.test(summary_df$residuals)
qqnorm(residuals(model), main = "Normal Q-Q Plot of Model Residuals")
qqline(residuals(model), col = "red", lwd = 2)

modelGenero <- lm(n_genes_recovered ~ Genero, data = summary_df)
summary_df$residualsGenero <- residuals(modelGenero)
shapiro.test(summary_df$residualsGenero)


leveneTest(n_genes_recovered ~ Familia, data = summary_df)
leveneTest(n_genes_recovered ~ Genero, data = summary_df)

# Mann-Whitney U (Wilcoxon rank-sum test)
wilcox_result <- wilcox.test(n_genes_recovered ~ Familia, data = summary_df)
print(wilcox_result)

kruskal.test(n_genes_recovered ~ Genero, data = summary_df)

print(wilcox_result)

library(rstatix)
dunn_results <- summary_df %>%
  dunn_test(n_genes_recovered ~ Genero, p.adjust.method = "fdr")
sig_pairs <- dunn_results %>%
  filter(p.adj < 0.05)

print(sig_pairs)


library(ggplot2)

ggplot(summary_df, aes(x = Genero, y = n_genes_recovered, fill = Familia)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  theme_minimal() +
  labs(
    x = "Genus",
    y = "Number of Recovered Loci",
    title = "Gene Recovery Across Cycad Genera"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
library(tidyverse)

table_dunn <- sig_pairs %>%
  select(
    `Genus 1` = group1,
    `Genus 2` = group2,
    `n1` = n1,
    `n2` = n2,
    `Z statistic` = statistic,
    `p-value (FDR)` = p.adj,
    `Significance` = p.adj.signif
  ) %>%
  mutate(`p-value (FDR)` = formatC(`p-value (FDR)`, format = "e", digits = 2))

write.csv(table_dunn, "Dunn_posthoc_results.csv", row.names = FALSE)






