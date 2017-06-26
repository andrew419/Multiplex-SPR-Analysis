# Kinetic Workflow
# Andrew Crowley
# May - June, 2017

# PRE-AUTOMATION TASKS

#   visually inspect sensorgram and cull obvious non-interactors
#   create "user_valid" column and assign pass/fail (0/1) observations to it
#   assign title "spot_ID" to first column and reformat "res stdev" with underscore

# LIBRARIES

library(tidyverse)
library(readxl)
library(stringr)

# WRANGLING

# by habit, Scrubber2 results are saved as an Excel workbook with one analyte per worksheet

sheets <- excel_sheets("<FILE_NAME>.xlsx")

##### BEGIN TEMPORARY #####

for(i in 1:length(sheets)) {
  assign(as.character(sheets[i]), read_excel("<FILE_NAME>.xlsx", i))
}

R2A_R131$analyte <- "R2A_R131"
R2A_H131$analyte <- "R2A_H131"
R2B$analyte <- "R2B"
R3A_F158$analyte <- "R3A_F158"
R3A_V158$analyte <- "R3A_V158"
R3B_NA1$analyte <- "R3B_NA1"
R3B_SH$analyte <- "R3B_SH"

full <- do.call(rbind, mget(sheets))

##### END TEMPORARY #####

full$spot_ID <- tolower(full$spot_ID)
full$Rmax <- as.numeric(full$Rmax)

full$analyte <- as.factor(full$analyte)

# historically, the protein concentrations used to print spots generate very similar results; by extracting the spot concentration, 
#   results can be grouped by higher-order feature (eg. subclass or Fc point mutation for antibodies)

name_split <- strsplit(full$spot_ID, "_")

full$conc <- NA

for (i in 1:length(name_split)) {
  full$conc[i] <- last(name_split[[i]])
}

full$conc[str_detect(full$spot_ID, "blank") == T] <- "_"
full$conc[str_detect(full$spot_ID, "bsa") == T] <- "_"

full$conc <- as.factor(full$conc)

# kinetic analysis in Scrubber outputs KD values using molarity units taking the form "_M"; unable to properly compare without 
#   anchoring to a uniform baseline of 10e00 (ie. 1 M)

full$KD_value <- substr(full$KD, 1, nchar(full$KD) - 2)
full$KD_value <- as.numeric(full$KD_value)

full$KD_unit <- substr(full$KD, nchar(full$KD) - 1, nchar(full$KD))
full$KD_unit <- tolower(full$KD_unit)

scale_to_molar <- tibble(
  KD_unit = c("mm", "um", "nm", "pm"),
  scaling = c(1e-03, 1e-06, 1e-09, 1e-012)
)

full <- left_join(full, scale_to_molar, by = "KD_unit") %>% 
  mutate(KD = KD_value * scaling)

# first part of the spot name typically contains the subclass or specificity of an antibody and is therefore a major category of sample

full$major <- NA

for (i in 1:length(name_split)) {
  full$major[i] <- first(name_split[[i]])
}

# have not encountered samples that use > 4 parts to their name, so it is probably overkill to count the elements of name_split
#   and split the strings accordingly; the following should therefore be sufficient

full$minor_1 <- NA
full$minor_2 <- NA

# because Asp (N) -> Ala (A) is a somewhat common Fc mutation, a different default_missing value than NA is supplied

for (i in 1:length(name_split)) {
  full$minor_1[i] <- nth(name_split[[i]], 2, default = "_")
}

for (i in 1:length(name_split)) {
  full$minor_2[i] <- nth(name_split[[i]], 3, default = "_")
}

# ends up duplicating the concentrations into the sample name variables for shorter observations, so these are removed

full$minor_1[full$minor_1 == full$conc] <- "_"
full$minor_2[full$minor_2 == full$conc] <- "_"

# per the recommendation of the US distributor for the SPR, the hallmarks of genuine receptor-ligand interactions are:
#   Rmax value > 50 response units (RU)
#   residual stdev < 10-15 (but using a fraction of Rmax for the stdev threshold may be appropriate when signal is high)

full$auto_valid <- 0
full$auto_valid[full$Rmax >= 50 & full$Res_sd <= 15] <- 1
full$auto_valid[is.na(full$KD)] <- 0
full$auto_valid[full$KD == 0] <- 0
full$auto_valid[full$KD > 1] <- 0
full$auto_valid <- as.factor(full$auto_valid)

full <- select(full, major, minor_1, minor_2, conc, analyte, KD, user_valid, auto_valid) %>%
  arrange(analyte, major, minor_1, minor_2, conc)

# VISUALIZATION
