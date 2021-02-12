library(biomartr)
library(tidyverse)
library(vcfR)

args = commandArgs(trailingOnly = TRUE)

repmask_out <- args[1]
annot_vcf <- args[2]
annot_file <- args[3]

rep_mask <- read_rm(repmask_out) %>%
  mutate(match_len = qry_end - qry_start) %>%
  group_by(qry_id) %>%
  filter(match_len == max(match_len)) %>%
  filter(sw_score == max(sw_score)) %>%
  select(matching_class, repeat_id)

vcf <- read.vcfR(annot_vcf)
vcf_df <- tibble(CHR = getCHROM(vcf),
                 POS = getPOS(vcf),
                 qry_id = paste(seq_len(nrow(vcf)), getID(vcf), sep = "_"))

annot <- left_join(vcf_df, rep_mask, by = "qry_id") %>%
  replace_na(list(matching_class = "None", repeat_id = "None")) %>%
  mutate(TE=grepl("(SINE|LINE|LTR|DNA|SVA)", matching_class)) %>%
  select(-qry_id)

write_tsv(annot, file = annot_file)
