library(spant)       # MRS tools
library(ggplot2)     # graphs 'n that
library(patchwork)   # plot glue
library(tidyverse)   # dataframe glue
library(kableExtra)  # nice table output

theme_set(theme_bw())

# load fit results
fit_list <- readRDS("all_fits_simple_basis.rds")

# remove fit 11 (MP35267_37) and 23 (MP39939_30) due to poor l/w and spurious
# echo respectively
fit_list <- fit_list[-c(11, 23)]

fit_result_table <- comb_fit_list_result_tables(fit_list)

fit_result_table$subj <- as.factor(fit_result_table$subj)

fit_result_table$subject <- fit_result_table$subj
levels(fit_result_table$subject) <- c("1", "2", "3", "4", "5", "6")
fit_result_table$participant <- fit_result_table$subject

SNR <- ggplot(fit_result_table, aes(x = subject, y = SNR)) +
  geom_point(alpha = 0.5, position = position_nudge(-0.1, 0)) +
  geom_boxplot(width = 0.1, position = position_nudge(0.1, 0),
               outlier.shape = NA) + ylab("SNR")

lw <- ggplot(fit_result_table, aes(x = subject, y = tNAA_lw)) +
  geom_point(alpha = 0.5, position = position_nudge(-0.1, 0)) +
  geom_boxplot(width = 0.1, position = position_nudge(0.1, 0),
               outlier.shape = NA) + ylab("Linewidth (PPM)")

qa_plot <- SNR + lw
qa_plot + plot_annotation(tag_levels = 'A')
ggsave("IMAGES/qa_plot.png", type="cairo-png", width = 7,
       height = 3.5)

# mean SNR
mean(fit_result_table$SNR)

# mean LW
mean(fit_result_table$tNAA_lw)

table_out <- fit_result_table %>%
  select(participant, GM_frac, SNR, tNAA_lw) %>%
  group_by(participant) %>% summarise(N = n(),
                                  mean_GM_frac = round(mean(GM_frac) * 100),
                                  sd_GM_frac = round(sd(GM_frac) * 100, 1),
                                  mean_SNR = round(mean(SNR), 0),
                                  sd_SNR = round(sd(SNR), 0),
                                  mean_lw = round(mean(tNAA_lw), 3),
                                  sd_lw = round(sd(tNAA_lw), 4))

kbl(table_out, col.names = c("participant", "N", "mean", "sd", "mean", "sd",
                             "mean", "sd")) %>% kable_classic() %>%
  add_header_above(c(" " = 1, " " = 1, "GM fraction (%)" = 2, 
                     "spectral SNR" = 2, "linewidth (ppm)" = 2)) %>%
  save_kable("TABLES/participant_table.html")

aov(GM_frac ~ subject, fit_result_table) %>% summary()
aov(SNR     ~ subject, fit_result_table) %>% summary()
aov(tNAA_lw ~ subject, fit_result_table) %>% summary()