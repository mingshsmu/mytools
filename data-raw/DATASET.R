## code to prepare `DATASET` dataset goes here

# out_path
out_path <- "/plots/"
usethis::use_data(out_path, overwrite = TRUE)

# candidate_drugs
candidate_drugs <- colnames(readRDS("D:/BaiduSyncdisk/Public_Database/R_functions/oncoPredict/Training Data/GDSC2_Res.rds"))
usethis::use_data(candidate_drugs, overwrite = TRUE)

