## Genetic structure analyses: PCA and ancestry coefficient
library(adegenet)
library(ggplot2)
library(dplyr)

X_mat <- tab(gen_obj)

pca_res <- dudi.pca(X_mat, scannf = FALSE, nf = 3)

variance_exp <- (pca_res$eig / sum(pca_res$eig)) * 100
pc1_var <- round(variance_exp[1], 2)
pc2_var <- round(variance_exp[2], 2)

pca_df <- as.data.frame(pca_res$li)

pca_df$Sample <- rownames(pca_df)

pca_df <- pca_df %>%
  mutate(
    Genus = sub("_.*", "", Sample)
  )

num_generos <- length(unique(pca_df$Genus))

pca_plot_genus <- ggplot(pca_df, aes(x = Axis1, y = Axis2, color = Genus, fill = Genus)) +
  geom_point(size = 3.5, alpha = 0.85, shape = 21, stroke = 0.8, color = "black") +
  labs(
    x = sprintf("PC1 (%.2f%%)", pc1_var),
    y = sprintf("PC2 (%.2f%%)", pc2_var),
    title = "Principal Component Analysis (PCA)",
    subtitle = sprintf("Based on 387 complete microhaplotype loci (0%% missing data, %d genera)", num_generos),
    color = "Genus",
    fill = "Genus"
  ) +
  scale_fill_viridis_d(option = "turbo") +
  scale_color_viridis_d(option = "turbo") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(color = "black", fill = NA, linewidth = 0.3),
    plot.title = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(color = "grey30", size = 11)
  )

# Mostrar el gráfico
print(pca_plot_genus)

# 7. Guardar las figuras en alta calidad
ggsave("PCA_Microhaplotypes_by_Genus_CompleteData.png", pca_plot_genus, width = 9, height = 6.5, dpi = 300)
ggsave("PCA_Microhaplotypes_by_Genus_CompleteData.pdf", pca_plot_genus, width = 9, height = 6.5)


### LEA
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("LEA")

library(LEA)
library(adegenet)
library(dplyr)
library(ggplot2)
library(tidyr)

matrix_to_export <- X_mat

if (nrow(matrix_to_export) == 320) {
  matrix_to_export <- t(matrix_to_export)
}

matrix_to_export[is.na(matrix_to_export)] <- 9

lines_geno <- vector("character", length = nrow(matrix_to_export))
for (i in 1:nrow(matrix_to_export)) {
  lines_geno[i] <- paste0(matrix_to_export[i, ], collapse = "")
}

writeLines(lines_geno, con = "microhaplotypes_transposed.geno")

cat(sprintf("-> Archivo exportado con %d líneas de %d caracteres (muestras para LEA).\n", 
            length(lines_geno), nchar(lines_geno[1])))

project_snmf <- snmf(
  "microhaplotypes_transposed.geno",
  K = 1:20,
  entropy = TRUE,
  repetitions = 5,
  project = "new"
)

mean_ce <- sapply(1:20, function(k) {
  mean(cross.entropy(project_snmf, K = k))
})

best_k_auto <- which.min(mean_ce)
K_selected <- best_k_auto
cat(sprintf("El K óptimo en el rango 1-20 es K = %d\n", best_k_auto))
# 6. Extraer y verificar la Matriz Q
best_run <- which.min(cross.entropy(project_snmf, K = K_selected))
Q_matrix <- Q(project_snmf, K = K_selected, run = best_run)

cat(sprintf("\n[VERIFICACIÓN FINAL]: Q_matrix tiene %d filas (muestras) y %d columnas (Ancestrías K).\n", 
            nrow(Q_matrix), ncol(Q_matrix)))

Q_df <- as.data.frame(Q_matrix)
colnames(Q_df) <- paste0("Ancestry_", 1:K_selected)

if (nrow(X_mat) == 320) {
  Q_df$Sample <- rownames(X_mat)
} else {
  Q_df$Sample <- colnames(X_mat)
}

Q_df$Genus <- sub("_.*", "", Q_df$Sample)

Q_long <- Q_df %>%
  pivot_longer(
    cols = starts_with("Ancestry_"),
    names_to = "Ancestry_Cluster",
    values_to = "Proportion"
  )

structure_plot <- ggplot(Q_long, aes(x = Sample, y = Proportion, fill = Ancestry_Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(~ Genus, scales = "free_x", space = "free_x") +
  labs(
    x = "Individual Samples grouped by Genus",
    y = "Ancestry Coefficients (Q)",
    title = sprintf("Genetic Population Structure (LEA snmf, K = %d)", K_selected),
    subtitle = "Proportion of individual genome assigned to each ancestral cluster",
    fill = "Ancestry Cluster"
  ) +
  scale_fill_viridis_d(option = "turbo") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.15, "lines"),
    strip.text = element_text(face = "bold.italic", size = 9, angle = 90),
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(structure_plot)

ggsave(sprintf("STRUCTURE_LEA20_K%d.png", K_selected), structure_plot, width = 12, height = 6, dpi = 300)
ggsave(sprintf("STRUCTURE_LEA20_K%d.pdf", K_selected), structure_plot, width = 12, height = 6)
