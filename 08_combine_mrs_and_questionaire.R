library(jsonlite)
library(spant)
library(dplyr)
library(readr)

# import the questionaire data
subj1_qresp <- as.data.frame(fromJSON("./QUESTIONNAIRE_DATA/QuestionnaireResponses_Subj01_clean.txt"))
subj1_qresp$subj_id <- "01"
subj1_qresp$session_id <- 1:nrow(subj1_qresp)
subj2_qresp <- as.data.frame(fromJSON("./QUESTIONNAIRE_DATA/QuestionnaireResponses_Subj02_clean.txt"))
subj2_qresp$subj_id <- "02"
subj2_qresp$session_id <- 1:nrow(subj2_qresp)
subj3_qresp <- as.data.frame(fromJSON("./QUESTIONNAIRE_DATA/QuestionnaireResponses_Subj03_clean.txt"))
subj3_qresp$subj_id <- "03"
subj3_qresp$session_id <- 1:nrow(subj3_qresp)
subj4_qresp <- as.data.frame(fromJSON("./QUESTIONNAIRE_DATA/QuestionnaireResponses_Subj04_clean.txt"))
subj4_qresp$subj_id <- "04"
subj4_qresp$session_id <- 1:nrow(subj4_qresp)
subj5_qresp <- as.data.frame(fromJSON("./QUESTIONNAIRE_DATA/QuestionnaireResponses_Subj05_clean.txt"))
subj5_qresp$subj_id <- "05"
subj5_qresp$session_id <- c(1:42, 44:(1 + nrow(subj5_qresp))) # 43 is missing (see day2day_dates_from_simone.csv)
subj6_qresp <- as.data.frame(fromJSON("./QUESTIONNAIRE_DATA/QuestionnaireResponses_Subj06_clean.txt"))
subj6_qresp$subj_id <- "06"
subj6_qresp$session_id <- 1:nrow(subj6_qresp)
subj7_qresp <- as.data.frame(fromJSON("./QUESTIONNAIRE_DATA/QuestionnaireResponses_Subj07_clean.txt"))
subj7_qresp$subj_id <- "07"
subj7_qresp$session_id <- c(1:7, 9:31, 33:40, 42:(3 + nrow(subj7_qresp))) # 8 and 32 and 41 are missing (see day2day_dates_from_simone.csv)
subj8_qresp <- as.data.frame(fromJSON("./QUESTIONNAIRE_DATA/QuestionnaireResponses_Subj08_clean.txt"))
subj8_qresp$subj_id <- "08"
subj8_qresp$session_id <- c(1:32, 34:(1 + nrow(subj8_qresp))) # 33 is missing (see day2day_dates_from_simone.csv)

# combine questionnaire data
q_data_full <- rbind(subj1_qresp, subj2_qresp, subj3_qresp, subj4_qresp,
                     subj5_qresp, subj6_qresp, subj7_qresp, subj8_qresp)

# convert keys to char and int
q_data_full$subj_id <- as.character(q_data_full$subj_id)

# load the MRS fit results
fit_list <- readRDS("all_fits_simple_basis.rds")

# remove fit 11 (MP35267_37) and 23 (MP39939_30) due to poor l/w and spurious
# echo respectively
fit_list <- fit_list[-c(11, 23)]

# get a table
fit_result_table <- comb_fit_list_result_tables(fit_list)

# convert keys to char and int
fit_result_table$subj_id <- as.character(fit_result_table$subj)
fit_result_table$session_id <- as.integer(unlist(lapply(
  strsplit(fit_result_table$scan_id, "_"), "[[", 2)))

# fix a naming error
fit_result_table[fit_result_table$subj_id == "07" & fit_result_table$session_id == "41",]$session_id <- 42

# join the tables
metab_and_q_data <- tibble(left_join(fit_result_table, q_data_full))

# save the table
write_csv(metab_and_q_data, "metab_and_q_data.csv")