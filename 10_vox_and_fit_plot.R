library(spant)        # MRS tools
library(tidyverse)    # data wrangling + ggplot2
library(broom)        # data wrangling + ggplot2
library(corrplot)     # corrplot
library(RColorBrewer) # cols for corrplot
library(patchwork)    # plotting glue
library(cowplot)      # plotting glue
require(gridGraphics) # plotting glue

theme_set(theme_bw())

# load fit results
fit_list <- readRDS("all_fits_simple_basis.rds")
fit_list[[10]]$extra # print the scan id

mrs <- fit_list[[10]]$data
mri <- readNifti("./MRS_MRI_CLEANED/MP35267_03/MP35267_36/36_T1_mprage_sag_MD_1mm_20131030081038_2.nii.gz")

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

# create a list column
spec_list_col <- full_fit_tab %>% group_by(PPMScale) %>% nest()
aov_fun <- function(x) aov(DataSC ~ subj, data = x)
mean_fun <- function(x) mean(x$DataSC)

# calc mean spectrum
spec_list_col <- spec_list_col %>%
  mutate(mean_spec = purrr::map_dbl(data, mean_fun))

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


png("IMAGES/fit_vox_plot.png", width = 600, height = 320, type = "cairo")
par(mfrow = c(1, 3))
plot(fit_list[[10]], restore_def_par = FALSE, mar = c(4, 2, 2, 0.5))
mtext("A", line = 0.5, col = "black", cex = 1.5, adj = 0)
plot_voi_overlay(mri, mrs, bg = "white", mar = c(4, 0.5, 2, 2))
mtext("B", line = 0.5, col = "black", cex = 1.5, adj = 0.115)
mean_spec
dev.off()

plot_fit <- function() {
  par(cex = 0.9)
  plot(fit_list[[10]], restore_def_par = TRUE, mar = c(4, 2, 2, 0.5))
}

plot_overlay <- function() {
  #plot_voi_overlay(mri, mrs, bg = "white", mar = c(4, 0.5, 2, 2))
  plot_voi_overlay(mri, mrs, bg = "white", mar = c(2.0, 1.5, 1.5, 1.5))
}

top_row <- plot_grid(plot_overlay, plot_fit, labels = c("A", "B"))

png("IMAGES/fit_vox_plot.png", width = 500, height = 500, type = "cairo")
plot_grid(top_row, mean_spec, ncol = 1, labels = c("", "C"))
dev.off()
