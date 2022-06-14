

# Load raw data and format for breeding status analysis


rm(list=ls(all=T))
getwd()

# Load required package
library(dplyr)
library(tidyr)
library(data.table)


counts <- read.csv("0_data/raw/songcounts_22March2021.csv")

length(table(as.factor(counts$ARU_ID))) # The number of territories

counts$SS <- paste0(counts$ARU_ID, "-", counts$year)
length(table(as.factor(counts$SS))) # The number of territory-years 

counts$standardized_bs <- factor(counts$standardized_bs, levels = c("Feeding Young", "Paired", "Single"))

counts$pkey <- paste0(counts$ARU_ID, "-", counts$year, "-", counts$JDate)


# Calculate N, number of individual bird-days
daily_sc <- counts %>% 
  group_by(SS,JDate,pkey) %>%
  summarize(n_samples_per_day = n(),
            total_songs_per_day = sum(two_min_sc))

daily_sc <- daily_sc %>% 
  mutate(detected = ifelse(total_songs_per_day>0,1,0))

detection_days <- daily_sc %>% 
  group_by(SS) %>%
  summarize(n_samples = n(),
            days_detected = sum(detected))

# Determine individuals to remove (birds detected on <30% of sample dates)
detection_days <- detection_days %>% 
  mutate(percent_detected = (round(days_detected/n_samples*100,1)))

# Remove low detection birds from data set
birds_low_det <- detection_days %>% 
  filter(percent_detected<30)

cut <- unique(birds_low_det$SS)

counts <- counts %>% 
  filter(!SS %in% cut)

# Update bird-day summary values for trimmed data set (still called counts)
daily_sc <- counts %>% 
  group_by(SS,JDate,pkey,standardized_bs) %>%
  summarize(n_samples_per_day = n(),
            total_songs_per_day = sum(two_min_sc))

daily_sc <- daily_sc %>% 
  mutate(detected = ifelse(total_songs_per_day>0,1,0))

detection_days <- daily_sc %>% 
  group_by(SS) %>%
  summarize(n_samples = n(),
            days_detected = sum(detected))

# Calculate time relative to sunrise (negative numbers mean 2 minute recording started before sunrise, positive means after)
counts$sec_from_sunrise <- counts$BinStarts_sec_past_midn-counts$sunrise_sec_past_midn

# Check range of recording start times (should be between 45 minutes before sunrise and 15 minutes after)
min(counts$sec_from_sunrise)/60
max(counts$sec_from_sunrise)/60
# It seems that I have actually sampled from 47 minutes before sunrise until 13 minutes after sunrise...
# aka sampled 2-minute periods where the end of the period lies within 45 mins before to 15 mins after sunrise

# Create a counts data table (matrix where rows are bird-day, aka counts$pkey, columns are sub-samples for that bird-day, and values are total number of counts for each sample, using NA when sub-sample does not exist). Include the pkey as the unique identifier
counts.survey <- rbindlist(lapply(1:nlevels(as.factor(counts$pkey)), function(x){
  p <- levels(as.factor(counts$pkey))[x]
  d <- data.table(counts[which(counts$pkey == p), ])
  d$survey <- 1:nrow(d)
  return(d)
}))

y <- spread(counts.survey[, c("pkey", "survey", "two_min_sc")], survey, two_min_sc)

# Create the same matrix with time relative to sunrise
TimeRel2Sun <- spread(counts.survey[, c("pkey", "survey", "sec_from_sunrise")], survey, sec_from_sunrise)

# Get Julian day for each survey
Jday <- daily_sc[, c("pkey", "JDate")]
Jday$jday <- Jday$JDate/365
Jday$jday2 <- Jday$jday^2

bs <- daily_sc$standardized_bs

all.equal(y$pkey, TimeRel2Sun$pkey, Jday$pkey, daily_sc$pkey)

save(y, TimeRel2Sun, Jday, daily_sc, file = "0_data/processed/osfl_data_package.Rdata")
