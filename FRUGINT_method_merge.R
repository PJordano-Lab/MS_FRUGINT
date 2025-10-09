library(tidyverse)
library(readr)
library(tibble)
library(here)

# The workflow for integrating interaction data from multiple methods includes several steps:
# 1. Integrating datasets at different resolutions (record- and species-level) and times
# 2. Merging and standardizing by-method matrix (Grand Total Standardization)

#--------------------------------------------------------------

# GTS FUNCTION
# Standardize pairwise interaction number by the grand total
gts <- function(df) {
  total_sum <- sum(df$total.int, na.rm = TRUE)

  df %>%
    mutate(gts = total.int / total_sum,
           gts = ifelse(gts == 0, NA, round(gts, 8)))
}

# Standardize proportion data by the grand total
gts.prop <- function(df) {
  df %>%
    ungroup() %>%
    mutate(gts = gts0 / sum(gts0)) %>%
    mutate(gts = ifelse(gts == 0, NA, round(gts, 8)))
}

#--------------------------------------------------------------

## INPUT DATA
# INPUT DATA: RECORD LEVEL (1-5)

# 1.Cameras
CAM <- read_csv(here("FRUGINT/INPUT_DATA/CAM.csv"))

# 2.Barcoding-visits
BC_visit <- read_csv(here("FRUGINT/INPUT_DATA/BC_visit.csv"))

# 3.Barcoding-seeds
BC_seed <- read_csv(here("FRUGINT/INPUT_DATA/BC_seed.csv"))

# 4.Tracks
TR <- read_csv(here("FRUGINT/INPUT_DATA/TRACK.csv"))

# 5.Mist-netting (2024)
MN_24 <- read_csv(here("FRUGINT/INPUT_DATA/MN_2024.csv"))


# INPUT DATA: SPECIES LEVEL (6-7)

# 6.Mist-netting (1983)
MN_83 <- read_csv(here("FRUGINT/INPUT_DATA/MN_1983.csv"))

# 7.Direct observation data from two locations with different vegetation types
OBS_HR <- read_csv(here("FRUGINT/INPUT_DATA/OBS_HR.csv"))
OBS_JP <- read_csv(here("FRUGINT/INPUT_DATA/OBS_JP.csv"))


# INPUT DATA: FRUGIVORE DIET
bird_trait <- read_csv2(here("FRUGINT/TRAITS/FRUGINT_birdTraits.csv")) %>% select(1,6)
mamm_trait <- read_csv2(here("FRUGINT/TRAITS/FRUGINT_mammalTraits.csv")) %>% select(1,6)

  # Add missing species
new_spp <- tibble(
  Frug_species = c("Psammodromus algirus", "Sturnus vulgaris/unicolor"),
  Diet_fruit_prop = c(0.1,0.8))

frug_diet <- bind_rows(bird_trait,mamm_trait,new_spp)


#--------------------------------------------------------------

## AGGREGATE DATA: From record to species level

# 1.Cameras and tracks (pooled)
CAM_TRACK <- bind_rows(CAM,TR) %>%
  group_by(plantSp,animalSp) %>%
  summarize(total.int = n()) %>%
  mutate(method = "CAM",
         method = if_else(plantSp == "Chamaerops humilis", "TRACK", method))

# 2.Barcoding
BC_visit <- BC_visit  %>%
  # Estimate probability of fruit consumption during a visit as fruit proportion in the diet
  left_join(frug_diet, by = c("animalSp" = "Frug_species")) %>%
  group_by(plantSp,animalSp) %>%
  summarize(total.int = sum(Diet_fruit_prop))

BC_seed <- BC_seed %>%
  group_by(plantSp,animalSp) %>%
  summarize(total.int = n())

BC <- bind_rows(BC_visit, BC_seed) %>%
  # Sum intensities across both datasets for each plant-animal pair
  group_by(plantSp, animalSp) %>%
  summarise(total.int = sum(total.int), .groups = "drop") %>%
  mutate(method = "BC")

# 3.Mist-netting (2024)
MN_24 <- MN_24 %>%
  group_by(vegetation, plantSp, animalSp) %>%
  summarize(total.int = n(), .groups = "drop")

#--------------------------------------------------------------

## AGGREGATE DATA: different times and locations

# 1.Mist-netting
  # GTS by time period and average
MN_83_gts <- gts(MN_83)
MN_24_gts <- gts(MN_24)
MN_gts <- bind_rows(MN_24_gts, MN_83_gts) %>%
  group_by(plantSp, animalSp) %>%
  summarize(gts = mean(gts, na.rm = TRUE),
            total.int = sum(total.int), .groups = "drop") %>%
  rename(gts0 = gts) %>%
  mutate(method = "MIST-NET")

# 2.Observation
  # GTS by vegetation and average
OBS_HR_gts <- gts(OBS_HR)
OBS_JP_gts <- gts(OBS_JP)
OBS_gts <- bind_rows(OBS_HR_gts, OBS_JP_gts) %>%
  group_by(plantSp, animalSp) %>%
  summarize(gts = mean(gts, na.rm = TRUE),
            total.int = sum(total.int), .groups = "drop") %>%
  rename(gts0 = gts) %>%
  mutate(method = "OBS")

#--------------------------------------------------------------

## GRAND TOTAL STANDARDIZATION BY METHOD

# Standardization of interaction number
CAM_TRACK_gts <- CAM_TRACK %>%
  gts() %>%
  rename(PIE = gts)
BC_gts <- BC %>%
  gts() %>%
  rename(PIE = gts)

# Standardization of proportions
MN_gts <- MN_gts %>%
  gts.prop() %>%
  select(-gts0) %>%
  rename(PIE = gts)
OBS_gts <- OBS_gts %>%
  gts.prop() %>%
  select(-gts0) %>%
  rename(PIE = gts)

## Save file "FRUGINT_byMethod_gts.csv"
all_gts <- bind_rows(CAM_TRACK_gts,BC_gts,MN_gts,OBS_gts) %>%
  select(plantSp, animalSp, method, total.int, PIE) %>%
  # Taxonomic harmonization to AviList
  mutate(animalSp = if_else(animalSp == "Corvus monedula", "Coloeus monedula", animalSp))

write_csv(all_gts, file = here("FRUGINT", "INTERACTIONS/FRUGINT_byMethod_gts.csv"))


#--------------------------------------------------------------

## FINAL MERGE BY MEAN

# Calculate the mean
final_matrix <- all_gts %>%
  group_by(plantSp, animalSp) %>%
  summarize(meanPIE = round(mean(PIE, na.rm = TRUE),8),
            maxPIE = max(PIE, na.rm = TRUE),
            n_methods = n_distinct(method),
            .groups = "drop") %>%
  arrange(desc(maxPIE)) %>%
  mutate(meanPIE = meanPIE / sum(meanPIE),
         maxPIE  = maxPIE  / sum(maxPIE))

n_distinct(final_matrix$plantSp)
n_distinct(final_matrix$animalSp)

# Save file "FRUGINT_finalMatrix.csv"
write_csv(final_matrix, file = here("FRUGINT", "INTERACTIONS/FRUGINT_finalMatrix.csv"))



