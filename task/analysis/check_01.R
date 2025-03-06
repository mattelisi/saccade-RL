rm(list=ls())
hablar::set_wd_to_script_path()

library(tidyverse)
# d <- read_delim("../data/S3", col_names=F)
d <- read_delim("../data/S3-edit", col_names=F)
str(d)
