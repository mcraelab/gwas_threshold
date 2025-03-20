library(data.table) 
library(ggplot2)
library(ggridges)
library(tidyverse)
library(patchwork)
library(ggpubr)
rm(list = ls())
## ++++++++++++++++++++++++++++++++++++++++++++++++++
#
# 2: Visualise dominance GWAS threshold analysis
#
## +++++++++++++++++++++++++++++++++++++++++++++++++++
# Load dominance GWAS results
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
  gwas_m$NegLogP <- -log10(gwas_m$V12)
  gwas_m$Group = group
  gwas_m
}
f_load_gwas
# Dominance normalized traits
dt_dominance_normal <- f_load_gwas(
  file_01percent = "out_phenoall_dominance_normal_min_pval_01percent.txt",
  file_05percent = "out_phenoall_dominance_normal_min_pval_05percent.txt",
  file_1percent = "out_phenoall_dominance_normal_min_pval_1percent.txt",
  file_5percent = "out_phenoall_dominance_normal_min_pval_5percent.txt",
  group = "Dominance_normal")
dt_dominance_normal

# Load dominace gwas bmi
dt_dominance_bmi <- f_load_gwas(
  file_01percent = "out_phenoall_dominance_bmi_min_pval_01percent.txt",
  file_05percent = "out_phenoall_dominance_bmi_min_pval_05percent.txt",
  file_1percent = "out_phenoall_dominance_bmi_min_pval_1percent.txt",
  file_5percent = "out_phenoall_dominance_bmi_min_pval_5percent.txt",
  group = "Dominance_bmi")
dt_dominance_bmi
dt_dominance_c <- rbind(dt_dominance_normal, dt_dominance_bmi)
dt_dominance_c$MAF_group <- paste(dt_dominance_c$Group, dt_dominance_c$MAF)
dt_dominance_c

#  Calculate True significance positive rate: dominance
dominance_FP <- table(dt_dominance_c$MAF_group, dt_dominance_c$V12 < 5.0e-08)
dominance_FP
dominance_FPR <- data.frame(MAF_group = rownames(dominance_FP), 
                            NonSig = dominance_FP[, "FALSE"], 
                            GenomeWideSig = dominance_FP[, "TRUE"])
dominance_FPR$FPR <- (dominance_FPR$GenomeWideSig/1000)*100
dominance_FPR
dt_dominance_percentile_95 <- aggregate(NegLogP ~ MAF_group, dt_dominance_c[,], function(x) quantile(x, 0.95))
dt_dominance_percentile_95
dt_dominance_percentile_95$Pval <- 10^(-dt_dominance_percentile_95$NegLogP)
dt_dominance_percentile_95_bmi <- dt_dominance_percentile_95[dt_dominance_percentile_95$MAF_group %like% "Dominance_bmi",]
dt_dominance_percentile_95_normal <- dt_dominance_percentile_95[dt_dominance_percentile_95$MAF_group %like% "Dominance_normal",]
dt_dominance_percentile_95_normal

# Perform bootstrapping to get the CI
# Function to calculate the 95th percentile
percentile_95 <- function(x) quantile(x, 0.95)
dt_dominance_c
result_df <- data.frame(MAF_group = unique(dt_dominance_c$MAF_group), 
                        lower = rep(NA, length(unique(dt_dominance_c$MAF_group))), 
                        upper = rep(NA, length(unique(dt_dominance_c$MAF_group))))
result_df
# Bootstrap procedure within each group
for (group in unique(dt_dominance_c$MAF_group)) {
  dt_group <- dt_dominance_c[MAF_group == group, ]
  n_bootstrap <- 1000
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
result_df
dt_dominance_percentile_95
m_CI_dominance <- merge(dt_dominance_percentile_95, result_df, by = "MAF_group")
m_CI_dominance
m_CI_dominance <- format(m_CI_dominance, format = "e", digits = 3)
m_CI_dominance$Meff <- as.integer(0.05/as.numeric(m_CI_dominance$Pval))
m_CI_dominance
# Save results
#fwrite(m_CI_dominance, paste0(out_tables, "/001b_emprical_GWA_cut-off_n_CI_dominance.tsv"), quote = F, col.names = T, row.names = F, sep = "\t")
# Plot
dt_dominance_c$MAF_group <- factor(dt_dominance_c$MAF_group, 
                                   levels = c( "Dominance_bmi 0.1% MAF",
                                               "Dominance_bmi 0.5% MAF",
                                               "Dominance_bmi 1% MAF",
                                               "Dominance_bmi 5% MAF", 
                                               "Dominance_normal 0.1% MAF",
                                               "Dominance_normal 0.5% MAF",
                                               "Dominance_normal 1% MAF", 
                                               "Dominance_normal 5% MAF") )
# Plot
custom_colors_normal <- c("Dominance_normal 5% MAF" = "blue", 
                          "Dominance_normal 1% MAF" = "green", 
                          "Dominance_normal 0.5% MAF" = "orange",
                          "Dominance_normal 0.1% MAF" = "grey"
)
# a) normal
plt_normal <- ggplot(dt_dominance_c[Group == "Dominance_normal", ], 
                     aes(x = NegLogP, y = MAF_group, fill = MAF_group, color = MAF_group)) + 
  theme_classic() + 
  geom_density_ridges(scale = 2, alpha = 0.3, bandwidth = 0.2)  + 
  ylab("Density") + xlim(4, 13) +
  xlab(expression(log[10](paste(P[min])))) +
  geom_point(data = dt_dominance_percentile_95_normal, 
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
# b) bmi
custom_colors_bmi <- c("Dominance_bmi 5% MAF" = "blue", 
                       "Dominance_bmi 1% MAF" = "green", 
                       "Dominance_bmi 0.5% MAF" = "orange",
                       "Dominance_bmi 0.1% MAF" = "grey"
)
plt_bmi <- ggplot(dt_dominance_c[Group == "Dominance_bmi", ], 
                  aes(x = NegLogP, y = MAF_group, fill = MAF_group, color = MAF_group)) + 
  theme_classic() + 
  geom_density_ridges(scale = 2, alpha = 0.3, bandwidth = 0.2)  + 
  ylab("") + xlim(4, 13) +
  xlab(expression(log[10](paste(P[min])))) + 
  geom_point(data = dt_dominance_percentile_95_bmi, 
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
# Combine plots
c.plts <- ggarrange(plt_normal, plt_bmi, labels = "AUTO")
c.plts
#  true false signficance at GWAS cut-off
TFP_c <- rbind(dominance_FPR, dominance_FPR)
#fwrite(TFP_c[, c(1, 3, 4)], paste0(out_tables, "/001a_True_false_positve_rate.tsv"), quote = F, col.names = T, row.names = F, sep = "\t")
