library(data.table) 
library(ggplot2)
library(ggridges)
library(tidyverse)
library(patchwork)
library(ggpubr)
rm(list = ls())
setwd("...")
out_tables="..."
# Load additive GWAS results
f_load_gwas<- function(file_01percent, file_05percent, file_1percent, file_5percent, group){
  gwas_01p <- fread(file_01percent, sep = "\t", fill = TRUE)
  gwas_01p$MAF <- "0.1% MAF"
  gwas_05p <- fread(file_05percent, sep = "\t", fill = TRUE)
  gwas_05p$MAF <- "0.5% MAF"
  gwas_1p <- fread(file_1percent, sep = "\t", fill = TRUE)
  gwas_1p$MAF <- "1% MAF"
  gwas_5p <- fread(file_5percent, sep = "\t", fill = TRUE)
  gwas_5p$MAF <- "5% MAF"
  gwas_m <- na.omit(rbind(gwas_01p, gwas_05p, gwas_1p, gwas_5p))
  gwas_m
  gwas_m$NegLogP <- -log10(gwas_m$V10)
  gwas_m$Group = group
  gwas_m
}
f_load_gwas
# Additive normal
dt_additive_normal <- f_load_gwas(
  file_01percent = "out_phenoall_additive_normal_min_pval_01percent.txt",
  file_05percent = "out_phenoall_additive_normal_min_pval_05percent.txt",
  file_1percent = "out_phenoall_additive_normal_min_pval_1percent.txt",
  file_5percent = "out_phenoall_additive_normal_min_pval_5percent.txt",
                             group = "Additive_normal")
dt_additive_normal

dt_additive_bmi <- f_load_gwas(
  file_01percent = "out_phenoall_additive_bmi_min_pval_01percent.txt",
  file_05percent = "out_phenoall_additive_bmi_min_pval_05percent.txt",
  file_1percent = "out_phenoall_additive_bmi_min_pval_1percent.txt",
  file_5percent = "out_phenoall_additive_bmi_min_pval_5percent.txt",
                          group = "Additive_bmi")
dt_additive_bmi
dt_additive_c <- rbind(dt_additive_normal, dt_additive_bmi)
dt_additive_c$MAF_group <- paste(dt_additive_c$Group, dt_additive_c$MAF)
dt_additive_c
dt_additive_c[MAF_group == "Additive_bmi 0.1% MAF",]

# Calculate True false significance rate: additive
additive_FP <- table(dt_additive_c$MAF_group, dt_additive_c$V10 < 5.0e-08)
additive_FP
additive_FPR <- data.frame(MAF_group = rownames(additive_FP), 
                            NonSig = additive_FP[, "FALSE"], 
                            GenomeWideSig = additive_FP[, "TRUE"])
additive_FPR$FPR <- (additive_FPR$GenomeWideSig/1000)*100
additive_FPR

dt_additive_percentile_95 <- aggregate(NegLogP ~ MAF_group, dt_additive_c[,], function(x) quantile(x, 0.95))
dt_additive_percentile_95
dt_additive_percentile_95$Pval <- 10^(-dt_additive_percentile_95$NegLogP)
dt_additive_percentile_95_bmi <- dt_additive_percentile_95[dt_additive_percentile_95$MAF_group %like% "Additive_bmi",]
dt_additive_percentile_95_bmi
dt_additive_percentile_95_normal <- dt_additive_percentile_95[dt_additive_percentile_95$MAF_group %like% "Additive_normal",]
dt_additive_percentile_95_normal

# Perfrom bootstrapping to get the CI
# Function to calculate the 95th percentile
percentile_95 <- function(x) quantile(x, 0.95)
dt_additive_c
# Create an empty data frame to store the results
result_df <- data.frame(MAF_group = unique(dt_additive_c$MAF_group), 
                        lower = rep(NA, length(unique(dt_additive_c$MAF_group))), 
                      upper = rep(NA, length(unique(dt_additive_c$MAF_group))))
result_df
# Bootstrap procedure within each group
for (group in unique(dt_additive_c$MAF_group)) {
  dt_group <- dt_additive_c[MAF_group == group, ]
  n_bootstrap <- 1000  # you can adjust this based on your needs
  bootstrap_results <- replicate(n_bootstrap, {
    sampled_data <- dt_group[sample(nrow(dt_group), replace = TRUE), ]
    aggregate(NegLogP ~ MAF_group, sampled_data, percentile_95)
  }, simplify = "data.frame")
  result_strings <- bootstrap_results[seq(1, n_bootstrap * 2, by = 2)]
  result_values <- bootstrap_results[seq(2, n_bootstrap * 2, by = 2)]
  result_numeric <- as.numeric(result_values)
  ci_percentile <- quantile(result_numeric, c(0.025, 0.975))
  result_df[result_df$MAF_group == group, "lower"] <- 10^(-ci_percentile[1])
  result_df[result_df$MAF_group == group, "upper"] <- 10^(-ci_percentile[2])
}
# Print the final results
result_df
dt_additive_percentile_95
m_CI_additive <- merge(dt_additive_percentile_95, result_df, by = "MAF_group")
m_CI_additive
m_CI_additive <- format(m_CI_additive, format = "e", digits = 3)
m_CI_additive$Meff <- as.integer(0.05/as.numeric(m_CI_additive$Pval))
m_CI_additive
# Save results
#fwrite(m_CI_additive, paste0(out_tables, "/001b_emprical_GWA_cut-off_n_CI_additive.tsv"), quote = F, col.names = T, row.names = F, sep = "\t")

# plot additive for the 95% cut-off
dt_additive_c$MAF_group <- factor(dt_additive_c$MAF_group, 
                                  levels = c( "Additive_bmi 0.1% MAF",
                                              "Additive_bmi 0.5% MAF",
                                              "Additive_bmi 1% MAF",
                                              "Additive_bmi 5% MAF", 
                                              "Additive_normal 0.1% MAF",
                                              "Additive_normal 0.5% MAF",
                                              "Additive_normal 1% MAF", 
                                              "Additive_normal 5% MAF") )
# Plot
custom_colors_normal <- c("Additive_normal 5% MAF" = "blue", 
                   "Additive_normal 1% MAF" = "green", 
                   "Additive_normal 0.5% MAF" = "orange",
                   "Additive_normal 0.1% MAF" = "grey"
                   )
## a) normal
plt_normal <- ggplot(dt_additive_c[Group == "Additive_normal", ], 
       aes(x = NegLogP, y = MAF_group, fill = MAF_group, color = MAF_group)) + 
  theme_classic() + 
  geom_density_ridges(scale = 2, alpha = 0.3, bandwidth = 0.2)  + 
  ylab("Density") + xlim(4, 13) +
  xlab(expression(log[10](paste(P[min])))) +
  geom_point(data = dt_additive_percentile_95_normal, 
             aes(x = NegLogP),color = "grey", size = 0.5, shape = 16, show.legend = FALSE) +
  geom_vline(xintercept = -log10(5.0e-8), linetype = "dashed", color = "black") + # gwas cut-off
  scale_fill_manual(values = custom_colors_normal) +
  scale_color_manual(values = custom_colors_normal) +
  theme(text = element_text(size = 18),
        legend.spacing.y  = unit(0.8, 'cm'),
        axis.text.y =element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none"
  ) + guides(fill = guide_legend(byrow = TRUE))
plt_normal 
## b) bmi
custom_colors_bmi <- c("Additive_bmi 5% MAF" = "blue", 
                          "Additive_bmi 1% MAF" = "green", 
                          "Additive_bmi 0.5% MAF" = "orange",
                          "Additive_bmi 0.1% MAF" = "grey"
)
plt_bmi <- ggplot(dt_additive_c[Group == "Additive_bmi", ], 
                     aes(x = NegLogP, y = MAF_group, fill = MAF_group, color = MAF_group)) + 
  theme_classic() + 
  geom_density_ridges(scale = 2, alpha = 0.3, bandwidth = 0.2)  + 
  ylab("") + xlim(4, 13) +
  xlab(expression(log[10](paste(P[min])))) + 
  geom_point(data = dt_additive_percentile_95_bmi, 
             aes(x = NegLogP),color = "grey", size = 0.5, shape = 16, show.legend = FALSE) +
  geom_vline(xintercept = -log10(5.0e-8), linetype = "dashed", color = "black") + # gwas cut-off
  scale_fill_manual(values = custom_colors_bmi) +
  scale_color_manual(values = custom_colors_bmi) +
  theme(text = element_text(size = 18),
        legend.spacing.y  = unit(0.8, 'cm'),
        axis.text.y =element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none"
  ) + guides(fill = guide_legend(byrow = TRUE))
plt_bmi
## combine plots
c.plts <- ggarrange(plt_normal, plt_bmi, labels = "AUTO")
c.plts
#  True false significance at GWAS cut-off
TFP_c <- rbind(additive_FPR)
TFP_c
#fwrite(TFP_c[, c(1, 3, 4)], paste0(out_tables, "/001a_True_false_positve_rate.tsv"), quote = F, col.names = T, row.names = F, sep = "\t")

