##+++++++++++++++++++++++++++++
#+
#+ Bar plots for the FPR
#+
##+++++++++++++++++++++++++++++
rm(list = ls())
library(data.table)
library(ggplot2)
setwd("...")
# load data
dt.irnt <- fread("002a_res_cojo_additive_maf_0.001_irnt_traits.tsv")
dt.irnt[Trait_Abbreviation == "Microalbumin", ]
dt.raw <- fread("002a_res_cojo_additive_maf_0.001_raw_traits.tsv")
dt.raw
dt.raw[Trait_Abbreviation == "Microalbumin", ]

dt <- rbind(dt.irnt, dt.raw)
dt
# split the id column
dt[, c("Trait_ID", "Trait_type") := tstrsplit(Trait_id, "_", fixed = TRUE)]
dt
# mean FPR for irnt
mean_fpr_irnt <- mean(dt[Trait_type == "irnt", ]$FPR)
mean_fpr_irnt
range_fpr_irnt <- range(dt[Trait_type == "irnt", ]$FPR)
range_fpr_irnt
# mean FPR raw
mean_fpr_raw <- mean(dt[Trait_type == "raw", ]$FPR)
mean_fpr_raw
range_fpr_raw <- range(dt[Trait_type == "raw", ]$FPR)
range_fpr_raw
dt[Trait_type == "raw", ][which.min(FPR), ]
dt[Trait_type == "raw", ][which.max(FPR), ]

# subset traits with the highest fpr based in irnt trait
dt.top.fpr.irnt <- dt.irnt[order(FPR, decreasing = T)][1:31, ][order(FPR), ]
dt.top.fpr.irnt
# subset raw & irnt
dt.top.fpr.irnt.s <- dt[Trait_Abbreviation %in% dt.top.fpr.irnt$Trait_Abbreviation, ]
dt.top.fpr.irnt.s
# order based on 'irnt'
tr.order <- dt.top.fpr.irnt.s[Trait_type=="irnt", ][order(FPR, Trait_type),]
tr.order
dt.top.fpr.irnt.s$Trait_Description <- factor(dt.top.fpr.irnt.s$Trait_Description, levels = c(tr.order$Trait_Description))
dt.top.fpr.irnt.s$Trait_type <- factor(dt.top.fpr.irnt.s$Trait_type, levels = c('raw', 'irnt'))
dt.top.fpr.irnt.s
# bar plot
p <- ggplot(data = dt.top.fpr.irnt.s, aes(x = Trait_Description, y = FPR,
                                          fill = Trait_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "False Significance Rate (%)", x = "") + 
  coord_flip() +
  theme_bw() + 
  scale_fill_manual(values = c("grey", "black")) + 
  scale_color_discrete(labels = c("irnt", "raw")) +
  theme(text = element_text(size = 16, face = "bold", colour = "black"), 
        legend.title = element_blank()
        ) 
# mean fpr
p + geom_hline(yintercept = mean_fpr_irnt, linetype = "dashed", color = "red")
