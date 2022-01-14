library(spant)       # MRS tools
library(MASS)        # LDA
library(tidyverse)   # data wrangling + ggplot2
library(ggfortify)   # convenient plotting for PCA results
library(broom)       # tidy stat results
library(patchwork)   # plotting glue

# load fit results
fit_list <- readRDS("all_fits_simple_basis.rds")

# remove fit 11 (MP35267_37) and 23 (MP39939_30) due to poor l/w and spurious
# echo respectively
fit_list <- fit_list[-c(11, 23)]

# get all the fit results into a single table
full_fit_tab <- comb_fit_list_fit_tables(fit_list, inc_indices = FALSE,
                                         add_res_id = FALSE)

full_fit_tab$subj <- as.factor(full_fit_tab$subj)

full_fit_tab$subj_nice <- full_fit_tab$subj
levels(full_fit_tab$subj_nice) <- c("participant 1", "participant 2",
                                    "participant 3", "participant 4",
                                    "participant 5", "participant 6")

# add baseline corrected column
full_fit_tab <- full_fit_tab %>% mutate(DataBC = Data - Baseline)

# find the spectral indices spanning the tCr resonance
ppm_sc <- fit_list[[1]]$fits[[1]]$PPMScale
ind    <- which(ppm_sc > 2.9 & ppm_sc < 3.1)

# scale all spectra by the height of the tCr resonance
full_fit_tab <- full_fit_tab %>% group_by(scan_id) %>%
                mutate(DataSC = DataBC / max(DataBC[ind])) %>% ungroup

# plot all spectra
ggplot(full_fit_tab) +
  geom_line(aes(x = PPMScale, y = DataSC, group = scan_id), alpha = 0.05) +
  scale_x_reverse()

# all spectra facet
ggplot(full_fit_tab) + theme_bw() +
  geom_line(aes(x = PPMScale, y = DataSC, group = scan_id, alpha = subj)) +
  scale_x_reverse() + facet_wrap(~ subj_nice) + 
  scale_alpha_manual(values = c(0.2, 0.15, 0.15, 0.15, 0.15, 0.15)) +
  theme(legend.position = "none") + xlab("Chemical shift (ppm)") +
  ylab("Intensity (a.u.)")

ggsave("IMAGES/all_spectra_faceted.png", type="cairo-png", width = 7,
       height = 4.5)

# mean spectra across groups
full_fit_tab_subj_means <- full_fit_tab %>% group_by(subj, PPMScale) %>%
                           summarise(mean_spec = mean(DataSC))

# facet
ggplot(full_fit_tab_subj_means) +
  geom_line(aes(x = PPMScale, y = mean_spec)) +
  scale_x_reverse() + facet_wrap(~ subj)

# coloured
ggplot(full_fit_tab_subj_means) +
  geom_line(aes(x = PPMScale, y = mean_spec, col = subj)) +
  scale_x_reverse()

# create a list column
spec_list_col <- full_fit_tab %>% group_by(PPMScale) %>% nest()
aov_fun <- function(x) aov(DataSC ~ subj, data = x)
mean_fun <- function(x) mean(x$DataSC)

# perform ANOVA
spec_list_col <- spec_list_col %>%
                 mutate(model = purrr::map(data, aov_fun)) %>%
                 mutate(model_tidy = purrr::map(model, tidy))  %>%
                 mutate(F_stat = map_dbl(model_tidy,
                                         function(x) x$statistic[1]))

# calc mean spectrum
spec_list_col <- spec_list_col %>%
  mutate(mean_spec = purrr::map_dbl(data, mean_fun))

f_spec <- ggplot(spec_list_col) + geom_line(aes(x = PPMScale, y = F_stat)) +
  scale_x_reverse() + xlab("Chemical shift (ppm)") + ylab("F-statistic") +
  theme_bw() + 
  annotate("text", 3.00, 50, label = "tCho") +
  annotate("text", 3.38, 48, label = "sIns") +
  annotate("text", 2.00, 20, label = "tNAA") +
  annotate("text", 3.57, 20.5, label = "Ins")

mean_spec <- ggplot(spec_list_col) + geom_line(aes(x = PPMScale, y = mean_spec)) +
  scale_x_reverse() + xlab("Chemical shift (ppm)") + ylab("Intensity (a.u.)") +
  theme_bw() + annotate("text", 1.85, 1.45, label = "tNAA") + 
  annotate("text", 2.92, 1.05, label = "tCr") + 
  annotate("text", 4.04, 0.50, label = "tCr") + 
  annotate("text", 3.22, 0.92, label = "tCho") +
  annotate("text", 3.35, 0.50, label = "sIns") + 
  annotate("text", 3.62, 0.50, label = "Ins") + 
  annotate("text", 3.48, 0.32, label = "Tau") + 
  annotate("text", 3.75, 0.65, label = "Glx+GSH") + 
  annotate("text", 1.30, 0.15, label = "Lac") + 
  annotate("text", 2.20, 0.35, label = "Glx") + 
  annotate("text", 2.45, 0.35, label = "Glu") + 
  annotate("text", 2.75, 0.35, label = "tNAA") + 
  annotate("segment", x = 3.35, xend = 3.34, y = 0.4, yend = 0.15) +
  annotate("segment", x = 3.47, xend = 3.44, y = 0.26, yend = 0.07) +
  annotate("segment", x = 3.76, xend = 3.76, y = 0.55, yend = 0.23) +
  annotate("segment", x = 3.65, xend = 3.55, y = 0.22, yend = 0.10) +
  annotate("segment", x = 3.60, xend = 3.60, y = 0.42, yend = 0.16) +
  annotate("segment", x = 2.20, xend = 2.20, y = 0.28, yend = 0.14) +
  annotate("segment", x = 2.40, xend = 2.36, y = 0.28, yend = 0.22) +
  annotate("segment", x = 2.70, xend = 2.62, y = 0.28, yend = 0.16) +
  annotate("segment", x = 2.70, xend = 2.55, y = 0.10, yend = 0.22)

# f_spec / mean_spec + plot_annotation(tag_levels = 'A')

f_spec # / mean_spec + plot_annotation(tag_levels = 'A')

ggsave("IMAGES/spectral_anova.png", type="cairo-png", width = 5, height = 4)

#ggsave("IMAGES/spectral_anova.png", type="cairo-png", width = 7, height = 8)
