library(spant)        # MRS tools
library(tidyverse)    # data wrangling + ggplot2
library(broom)        # data wrangling + ggplot2
library(kableExtra)   # nice table output
library(mlr3)         # machine learning
library(mlr3verse)    # machine learning
library(mlr3viz)      # machine learning
library(corrplot)     # corrplot
library(RColorBrewer) # cols for corrplot
library(patchwork)    # plotting glue

theme_set(theme_bw())

# load fit results
fit_list <- readRDS("all_fits_simple_basis.rds")

# remove fit 11 (MP35267_37) and 23 (MP39939_30) due to poor l/w and spurious
# echo respectively
fit_list <- fit_list[-c(11, 23)]

# scale metab levels to tCr
fit_list <- lapply(fit_list, scale_amp_ratio, "tCr")

# get a single fit results table
res_tab <- comb_fit_list_result_tables(fit_list)

# convert to factors
res_tab$subj <- as.factor(res_tab$subj)

res_tab$subject <- res_tab$subj
levels(res_tab$subject) <- c("subject 1", "subject 2", "subject 3",
                             "subject 4", "subject 5", "subject 6")

res_tab$subject_short <- res_tab$subj
levels(res_tab$subject_short) <- c("1", "2", "3",
                                   "4", "5", "6")

# write to file
write.table(res_tab, "full_res_tab.txt")

# ANOVA table
aov_fun <- function(x) aov(value ~ subj, data = x)
metab_anova <- res_tab %>% select(Ala, Asp, GABA, Glc, Gln, GSH, Glu, GPC,
                                  Ins, Lac, NAA, NAAG, PCh, sIns, Tau, tNAA,
                                  tCho, Glx, subj) %>%
  pivot_longer(-subj, names_to = "metab") %>%
  group_by(metab) %>% nest %>% 
  mutate(model = purrr::map(data, aov_fun)) %>%
  mutate(model_tidy = purrr::map(model, tidy)) %>%
  mutate(F_stat = map_dbl(model_tidy, function(x) x$statistic[1])) %>%
  mutate(p_value = map_dbl(model_tidy, function(x) x$p.value[1])) %>% arrange(desc(F_stat))

metab_anova_summary <- metab_anova %>% select(metab, F_stat, p_value)

# 0.05/18 = 2.78e-3

metab_anova_summary$p_value <- formatC(metab_anova_summary$p_value,
                                       format = "e", digits = 2)

metab_anova_summary$F_stat <- round(metab_anova_summary$F_stat, 1)

kable(metab_anova_summary, "html",
      col.names = c("metab. / tCr", "F-statistic", "p-value")) %>%
  kable_styling("striped") %>% kable_styling(full_width = F) %>%
  save_kable("TABLES/anova_table.png")

kable(metab_anova_summary, "html",
      col.names = c("metab. / tCr", "F-statistic", "p-value")) %>%
  kable_styling("striped") %>% kable_styling(full_width = F) %>%
  save_kable("TABLES/anova_table.html")

# library(gridExtra)
# tt <- ttheme_default(core=list(fg_params = list(hjust=1, x=0.9)),
#                      rowhead=list(fg_params = list(hjust=1, x=0.95)))
# png(filename = "IMAGES/output.png", width=480, height=480, bg = "white")
# grid.table(metab_anova_summary, theme = tt)
# dev.off()

# plot sIns vs tCho
ggplot(res_tab) + geom_point(aes(x = sIns, y = tCho, col = subject)) + 
  xlab("sIns / tCr") + ylab("tCho / tCr")
ggsave("IMAGES/tCho_vs_sIns.png")

# classifier simple
data_classify <- res_tab[c("subject_short", "tCho", "sIns")]
data_classify$subject_short <- as.factor(data_classify$subject_short)

task <- TaskClassif$new(id = "day2day", backend = data_classify,
                        target = "subject_short")

learner <- lrn("classif.svm", predict_type = "prob", kernel = "linear")

measure <- msr("classif.acc")
             
set.seed(100)
resampling <- rsmp("loo")
resampling$instantiate(task)

rr <- resample(task, learner, resampling)

print(rr$aggregate(measure))

plot_learner_prediction(learner, task) +
  xlab("sIns / tCr") + ylab("tCho / tCr") +
  guides(shape = FALSE, alpha = FALSE) + labs(fill = "participant")
ggsave("IMAGES/linear_svm.png", type="cairo-png", width = 5.5, height = 4)


# classifier more complex
data_classify <- res_tab[c("subject", "tCho", "sIns", "tNAA", "Ins")] # 97%
# data_classify <- res_tab[c("subj", "tCho", "sIns", "tNAA", "Tau")] # 97%
data_classify$subject <- as.factor(data_classify$subject)

task <- TaskClassif$new(id = "day2day", backend = data_classify,
                        target = "subject")

learner <- lrn("classif.svm", predict_type = "prob", kernel = "linear")

measure <- msr("classif.acc")

set.seed(100)
resampling <- rsmp("loo")
resampling$instantiate(task)

rr <- resample(task, learner, resampling)

print(rr$aggregate(measure))

# metab. level vs GM fraction
sins_gm <- ggplot(res_tab) +
  geom_point(aes(x = GM_frac * 100, y = sIns, col = subject_short)) + 
  ylab("sIns / tCr") + xlab("GM fraction (%)") +
  geom_smooth(aes(x = GM_frac * 100, y = sIns), method = "lm", col = "black") +
  annotate(geom = "text", x = 53, y = 0.072, 
    label = "r = 0.02\np = 0.87", hjust = 0, vjust = 1, size = 4) +
  labs(col = "participant")
  

# ggsave("IMAGES/GM_frac_vs_sIns.png")

cor.test(res_tab$sIns, res_tab$GM_frac, method = "pearson")

tcho_gm <- ggplot(res_tab) +
  geom_point(aes(x = GM_frac * 100, y = tCho, col = subject_short)) + 
  ylab("tCho / tCr") + xlab("GM fraction (%)") +
  geom_smooth(aes(x = GM_frac * 100, y = tCho), method = "lm", col = "black") +
  annotate(geom = "text", x = 58, y = 0.32, 
    label = "r = -0.15\np = 0.21", hjust = 0, vjust = 1, size = 4) +
  labs(col = "participant")

# ggsave("IMAGES/GM_frac_vs_tCho.png")

cor.test(res_tab$tCho, res_tab$GM_frac, method = "pearson")

gm_frac_plot <- sins_gm + tcho_gm
gm_frac_plot + plot_annotation(tag_levels = 'A')
ggsave("IMAGES/gm_frac.png", type="cairo-png", width = 9,
       height = 3.2)


# PCA plot

pca_tab <- res_tab %>% select(Ala, Asp, GABA, Glc, GSH, Ins, Lac, sIns, Tau,
                              tNAA, tCho, Glx)

pca_res <- prcomp(pca_tab, scale. = TRUE)

scores_1_2 <- autoplot(pca_res, data = res_tab, colour = 'subject_short',
                   loadings = TRUE, loadings.label = TRUE,
                   loadings.label.colour = "black",
                   loadings.colour = 'gray', x = 1, y = 2,
                   loadings.label.size = 3, loadings.label.hjust = 0.5,
                   loadings.label.vjust = 0.5) +
          labs(color='participant') 

scores_1_3 <- autoplot(pca_res, data = res_tab, colour = 'subject_short',
                   loadings = TRUE, loadings.label = TRUE,
                   loadings.label.colour = "black",
                   loadings.colour = 'gray', x = 1, y = 3,
                   loadings.label.size = 3, loadings.label.hjust = 0.5,
                   loadings.label.vjust = 0.5) +
          labs(color='participant') 

scree <- data.frame(sd = pca_res$sdev) %>%
  mutate(pct = 100 * (sd ^ 2 / sum(sd ^ 2))) %>%
  ggplot(aes(as.factor(1:12), pct)) + geom_col() +
  ylab("Explained varience (%)") + xlab("Principal component")

(scores_1_2 + scores_1_3) / scree + plot_annotation(tag_levels = 'A')
ggsave("IMAGES/pca.png", type="cairo-png", width = 9.5, height = 7)