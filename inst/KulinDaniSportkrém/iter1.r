library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)

ojj <-  readxl::read_xlsx(here::here("inst","KulinDaniSportkrém","example.xlsx"))
bar <- ojj[,1][[1]]

hist( bar , breaks= 30)


# ------

stiff <-  readxl::read_xlsx(here::here("inst","KulinDaniSportkrém",
                                       "example2.xlsx"))

stiff$PatientId  <- as.factor(stiff$PatientId)
stiff$MyStudyDay <- as.factor(stiff$MyStudyDay)
stiff$Csoport <- as.factor(stiff$Csoport)
stiff$y <- stiff$`Stiffness index`

# TODO: long to wide baseline inclusion
stiff_wide <- stiff
stiff_wide$base <- NA

for (i in 1:nrow(stiff)) {
  act_base <- stiff %>%
                filter(PatientId == stiff$PatientId[i],
                       MeasurementNumberInDay == 1,
                       MyStudyDay == stiff$MyStudyDay[i]) %>%
                .$y
  if (length(act_base) != 0) {
    stiff_wide$base[i] <- act_base
  }

}
stiff_wide <- stiff_wide %>% filter(MeasurementNumberInDay == 2)

stiff_wide$change <- stiff_wide$y - stiff_wide$base
hist(stiff_wide$change, breaks = 30)


mod <- lmer( y ~
               Csoport +
               base +
               (1|MyStudyDay) +
               (1|PatientId)
             ,
             data = stiff_wide)

summary(mod)

plot(ranef(mod))
plot(mod)

mod %>% effects::predictorEffects() %>% plot()

emm <- emmeans(mod, specs = list("Csoport"))
emm_con <- pairs(contrast(emm))

(emm_con)

plot(emm_con)


