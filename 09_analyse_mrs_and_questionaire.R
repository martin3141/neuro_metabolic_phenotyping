library(tidyverse)
library(corrplot)
library(RColorBrewer)  # cols for corrplot

theme_set(theme_bw())

# read the table
all_results <- read_csv("metab_and_q_data.csv")

all_results$subject <- as.factor(all_results$subj)
levels(all_results$subject) <- c("1", "2", "3", "4", "5", "6")
all_results$participant <- all_results$subject

# scale all metabolites by tCr
all_results <- all_results %>% mutate(across(c(Ala:Lac, NAA:Glx), ~ .x / tCr))

ggplot(all_results, aes(x = scanTime, y = Glx)) + 
  geom_point() + geom_smooth(method = "lm") + facet_wrap(~ subj_id)

ggplot(all_results, aes(x = daysFromFirstScan, y = sIns)) + 
  geom_point() + geom_smooth(method = "lm") + facet_wrap(~ subj_id)

ggplot(all_results, aes(x = daysFromFirstScan, y = tCho)) + 
  geom_point() + geom_smooth(method = "lm") + facet_wrap(~ subj_id)

ggplot(all_results, aes(x = weight_kg, y = tNAA)) + 
  geom_point() + geom_smooth(method = "lm") + facet_wrap(~ subj_id, scale = "free_x")

ggplot(all_results) + 
  geom_point(aes(x = weight_kg, y = tNAA, col = subj_id)) +
  geom_smooth(aes(x = weight_kg, y = tNAA), method = "lm")

lm(weight_kg ~ tNAA, all_results) %>% summary()

mean_weight <- all_results %>% select(tNAA, weight_kg, subj_id) %>%
  group_by(subj_id) %>%
  summarise(mean_tnaa = mean(tNAA),
            mean_weight = mean(weight_kg, na.rm = TRUE))

lm(mean_weight ~ mean_tnaa, mean_weight) %>% summary()

ggplot(mean_weight, aes(x = mean_weight, y = mean_tnaa)) +
  geom_point(aes(col = subj_id)) +
  geom_smooth(method = "lm", se = FALSE)

# interesting one >
ggplot(all_results, aes(y = all_results$bloodPressure_diastolic_mmHg, x = Glx)) + 
  geom_point() + geom_smooth(method = "lm") + facet_wrap(~ subj_id) 

ggplot(all_results, aes(y = all_results$bloodPressure_diastolic_mmHg, x = Glx, col = subj_id)) + 
  geom_point() + geom_smooth(method = "lm", se = FALSE)

library(lme4)
library(lmerTest) # to get p-value estimations that are not part of the standard lme4 packages

mixed_model <- lmer(formula = Glx ~ 1 + bloodPressure_diastolic_mmHg + (1 | subj_id),
     data    = all_results) #to run the model

summary(mixed_model)


ggplot(all_results, aes(x = tNAA, y = Glx, col = subj_id)) + 
  geom_point()

ggplot(all_results, aes(x = sIns, y = Ins, col = subj_id)) + 
  geom_point()

metab_only <- all_results %>% select(Ala:Lac, NAA:Glx, subj_id)

# not much drinking going on here
ggplot(all_results, aes(x = sIns, y = alcohol_last24hs_drinks, col = subj_id)) + 
       geom_point()

# time of day effects?
ggplot(all_results, aes(x = scanTime, y = Glx))  + geom_point()  + facet_wrap(~ subj_id)
ggplot(all_results, aes(x = scanTime, y = tNAA)) + geom_point() + facet_wrap(~ subj_id)
ggplot(all_results, aes(x = scanTime, y = sIns)) + geom_point() + facet_wrap(~ subj_id)

cor.test(all_results$Glx,  all_results$scanTime)
cor.test(all_results$tNAA, all_results$scanTime)
cor.test(all_results$sIns, all_results$scanTime)

ggplot(all_results, aes(x = scanTime, y = sIns, col = subj)) + geom_point()

# mixed effects model for time of day
lmer(formula = scanTime ~ 1 + Glx + (1 | subj_id), data = all_results) %>% summary
lmer(formula = scanTime ~ 1 + tNAA + (1 | subj_id), data = all_results) %>% summary


# this is a useful sanity check as weight should be stable
ggplot(all_results, aes(x = subj_id, y = weight_kg, color = subj_id)) + geom_point()




mean(all_results$sIns / all_results$Ins)

# subtract subject means
metab_only_mc <- metab_only %>% group_by(subj_id) %>%
  mutate(across(everything(), ~ .x - mean(.x))) %>% ungroup()

cor_tab <- metab_only_mc %>% select(Ala, Asp, GABA, Glc, Gln, GSH, Glu, GPC, Ins, Lac,
                              NAA, NAAG, PCh, sIns, Tau, tNAA, tCho, Glx)

M <- cor(cor_tab)
corrplot(M, type = "upper", order = "hclust",
         col = brewer.pal(n = 8, name = "RdYlBu"))

ggplot(metab_only_mc, aes(x = tNAA, y = tCho, col = subj_id)) + geom_point()

ggplot(metab_only_mc, aes(x = Glx, y = GABA, col = subj_id)) + geom_point()

ggplot(metab_only_mc, aes(x = tNAA, y = tCho)) + geom_point() +
  geom_smooth(method = "lm") + facet_wrap(~ subj_id)

ggplot(metab_only_mc, aes(x = tNAA, y = Ins)) + geom_point() +
  geom_smooth(method = "lm") + facet_wrap(~ subj_id)

ggplot(metab_only_mc, aes(x = Glx, y = GABA)) + geom_point() +
  geom_smooth(method = "lm") + facet_wrap(~ subj_id)

gm_frac_tab <- all_results %>% select(GM_frac)
metab_only_mc_gm_frac <- bind_cols(metab_only_mc, gm_frac_tab)

lm(GABA ~ Glx, metab_only_mc_gm_frac) %>% summary
lm(tNAA ~ tCho, metab_only_mc_gm_frac) %>% summary



lm(tNAA ~ GM_frac, metab_only_mc_gm_frac) %>% summary
lm(Glx ~ GM_frac, metab_only_mc_gm_frac) %>% summary

mixed_model <- lmer(formula = GABA ~ 1 + Glx + (1 | subj_id),
     data    = all_results) # to run the model

lmer(formula = GABA ~ 1 + Glx + (1 | subj_id), data = all_results) %>% summary
lmer(formula = tNAA ~ 1 + tCho + (1 | subj_id), data = all_results) %>% summary

# plot of scan dates
ggplot(all_results, aes(x = scanDate, y = participant)) + geom_point() + 
  xlab("Scan date") + ylab("Participant") + scale_x_date(date_labels = "%b %Y")

ggsave("IMAGES/scan_dates.png", type="cairo-png", width = 4, height = 4)

first_scan <- all_results %>% select(scanDate, participant) %>%
  group_by(participant) %>% filter(scanDate == min(scanDate))

last_scan <- all_results %>% select(scanDate, participant) %>%
  group_by(participant) %>% filter(scanDate == max(scanDate))

# mean time between first and last MRS scan
(last_scan$scanDate - first_scan$scanDate) %>% mean


# Simone suggestions

library(GGally)
ggplot(all_results, aes(x = participant, y = GM_frac * 100)) + geom_boxplot() +
  ylim(0, 100)

ggpairs(all_results %>% select(MR_HeliumLevel, MR_Room_Humidity,
                               MR_Room_Temperature, tNAA, tCho, sIns, Ins, Glx,
                               participant),
        aes(colour = participant))
