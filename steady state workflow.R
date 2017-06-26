# Complete Steady State Workflow
# Andrew Crowley
# May - June, 2017

# PRE-AUTOMATION TASKS

#   visually inspect sensorgram and cull obvious non-interactors
#   create "user_valid" column and assign pass/fail (0/1) observations to it

# LIBRARIES

library(tidyverse)
library(stringr)
library(readxl)

# WRANGLING

# by habit, Scrubber2 results are saved as an Excel workbook with one analyte per worksheet

sheets <- excel_sheets("<FILE_NAME.xlsx")

##### BEGIN TEMPORARY #####

for(i in 1:length(sheets)) {
  assign(as.character(sheets[i]), read_excel("<FILE_NAME>.xlsx", i))
}

# bypass process to insert analyte column for the moment, manually write in
R2A_R131$analyte <- "R2A_R131"
R2A_H131$analyte <- "R2A_H131"
R2B$analyte <- "R2B"
R3A_F158$analyte <- "R3A_F158"
R3A_V158$analyte <- "R3A_V158"
R3B_NA1$analyte <- "R3B_NA1"
R3B_SH$analyte <- "R3B_SH"

full <- do.call(rbind, mget(sheets))

##### END TEMPORARY #####

full <- rename(full, spot_ID = Name, Rmax_global = Rmax, KD_global = KD, Rmax_local = Rmax__1, KD_local = KD__1)
full <- select(full, spot_ID, analyte, Rmax_global, KD_global, Rmax_local, KD_local, user_valid)

full <- fill(full, Rmax_global)

full$spot_ID <- tolower(full$spot_ID)
full$Rmax_local <- as.numeric(full$Rmax_local)
full$Rmax_global <- as.numeric(full$Rmax_global)
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

# local results are already exported as properly scaled, numeric-only values, but global results need to be scaled to molar baseline

full$KD_global_value <- substr(full$KD_global, 1, nchar(full$KD_global) - 2)
full$KD_global_value <- as.numeric(full$KD_global_value)

full$KD_global_unit <- substr(full$KD_global, nchar(full$KD_global) - 1, nchar(full$KD_global))
full$KD_global_unit <- tolower(full$KD_global_unit)

scale_to_molar <- tibble(
  KD_global_unit = c("mm", "um", "nm", "pm"),
  scaling = c(1e-03, 1e-06, 1e-09, 1e-012)
)

full <- left_join(full, scale_to_molar, by = "KD_global_unit") %>% 
  mutate(KD_global = KD_global_value * scaling)

full <- select(full, major, minor_1, minor_2, conc, analyte, Rmax_global, KD_global, Rmax_local, KD_local, user_valid)

global_SS <- select(full, analyte, major, minor_1, minor_2, conc, Rmax_global, KD_global, user_valid)
global_SS <- rename(global_SS, Rmax = Rmax_global, KD = KD_global)

local_SS <- select(full, analyte, major, minor_1, minor_2, conc, Rmax_local, KD_local, user_valid)
local_SS <- rename(local_SS, Rmax = Rmax_local, KD = KD_local)

# per the recommendation of the US distributor for the SPR, the hallmarks of genuine receptor-ligand interactions are:
#   Rmax value > 50 response units (RU)
#   residual stdev < 10-15 (but using a fraction of Rmax for the stdev threshold may be appropriate when signal is high)

# steady state does not produce residual stdev, so only Rmax is used

global_SS$auto_valid <- 0
global_SS$auto_valid[global_SS$Rmax >= 50] <- 1
global_SS$auto_valid[is.na(global_SS$KD)] <- 0
global_SS$auto_valid[global_SS$KD == 0] <- 0
global_SS$auto_valid[global_SS$KD > 1] <- 0
global_SS$auto_valid <- as.factor(global_SS$auto_valid)

global_SS$user_valid <- as.factor(global_SS$user_valid)

global_SS <- select(global_SS, -Rmax) %>%
  arrange(analyte, major, minor_1, minor_2, conc)


local_SS$auto_valid <- 0
local_SS$auto_valid[local_SS$Rmax >= 50] <- 1
local_SS$auto_valid[is.na(local_SS$KD)] <- 0
local_SS$auto_valid[local_SS$KD == 0] <- 0
local_SS$auto_valid[local_SS$KD > 1] <- 0
local_SS$auto_valid <- as.factor(local_SS$auto_valid)

local_SS$user_valid <- as.factor(local_SS$user_valid)

local_SS <- select(local_SS, -Rmax) %>%
  arrange(analyte, major, minor_1, minor_2, conc)

# VISUALIZATION