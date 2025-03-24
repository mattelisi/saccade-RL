rm(list=ls())
hablar::set_wd_to_script_path()
library(tidyverse)
library(mlisi)

d <- read_csv("../task/analysis/analysis_all_data.csv")
str(d)

# # remove this until you manage to solve issue from matlab:
# d <- d %>%
#   filter(str_detect(participant_id, "^[A-Za-z]"))

# general prep
d <- d %>%
  mutate(PID = str_sub(edf_file, 1, 4)) %>%
  rename(session_id = participant_id) %>%
  mutate(win=ifelse(win==-1, NA, win),
         P=ifelse(win==-1, NA, P)) %>%
  filter(!is.na(win))

d_i <- d[d$PID==unique(d$PID)[6],]
str(d_i)

# cart2pol(x, y)

par(mfrow=c(2,2))
hist(d_i$sacRT)
hist(d_i$sacAmp)
hist(d_i$sacVPeak)
with(d_i, plot(sacAmp, log10(sacVPeak)))

with(d_i, cor.test(sacAmp, sacVPeak))

which(is.na(d_i$sacAmp))
which(is.na(d_i$sacVPeak))

d_i$velres <- NA
d_i$velres[!is.na(d_i$sacAmp)] <- residuals(lm(sacVPeak~sacAmp, data=d_i))
#d_i$velres[!is.na(d_i$sacAmp)] <- residuals(lm(log(sacVPeak)~sacAmp, data=d_i))

d_i$hpTar <- ifelse(d_i$prob_2>0.5,2,1)
d_i$sacHP <- ifelse(d_i$sacChoice==d_i$hpTar, 1,0)

d_i$svelres <- ifelse(d_i$sacHP==1, abs(d_i$velres),
                      -abs(d_i$velres))

d_i$absvelres <- abs(d_i$velres)

agd <- d_i %>%
  filter(!is.na(sacAmp)) %>%
  group_by(trial_n, sacHP) %>%
  summarise(choice = mean(sacHP),
            sacVPeak = mean(abs(sacVPeak)),
            absvelres = mean(absvelres),
            sacAmp = mean(sacAmp),
            sacRT = mean(sacRT))

agd %>%
  mutate(sacHP = factor(sacHP)) %>%
  ggplot(aes(x=trial_n, y=absvelres, color=sacHP, group=sacHP))+
  #geom_line()+
  geom_point()+geom_smooth(method="lm")

library(lme4)
mm0 <- lmer(absvelres ~ trial_n*sacHP +(trial_n*sacHP|| PID), d_i)
summary(mm0)


