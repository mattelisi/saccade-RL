rm(list=ls())
hablar::set_wd_to_script_path()
library(tidyverse)

d <- read_csv("../task/analysis/analysis_all_data.csv")
str(d)

# function that put data of a single participant into stan format
d_i <- d[d$participant_id==unique(d$participant_id)[1],]
str(d_i)

