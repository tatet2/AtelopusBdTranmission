library(ggplot2)
library(survival)
library(dplyr)

source('ggsurv.R')

theme_set(theme_bw())
theme_update(legend.title = element_text(size=rel(1.2)), legend.text = element_text(size=13), axis.title = element_text(size=15), axis.text = element_text(size=15))

## survival plot
death <- read.csv("death.csv")
death$treatment <- relevel(death$treatment,"s")
deathsub <- droplevels(subset(death, days <= 62))

deathsub$comb <- with(deathsub, paste0(treatment,species))
surv.res <- with(deathsub, survfit( Surv(days, cens)~species ) )
coxph( Surv(days, cens)~species , data=deathsub) 
psplot <- ggsurv(surv.res)
psplot + guides(linetype = F) + scale_colour_discrete(name = 'Species', breaks = c("a","c"), labels=c('Atelopus', 'Espadarana')) + ggtitle('Atelopus and Espadarana Survival')
ggsave("survSpecies.tiff")


deathAt <- droplevels(subset(deathsub, species=="a"&treatment!="n"))
surv.atel <- with(deathAt, survfit( Surv(days, cens)~treatment ) )
atel.cph <- coxph( Surv(days, cens)~treatment-1 , data=deathAt) 
psplotA <- ggsurv(surv.atel)
psplotA + guides(linetype = F) + scale_colour_discrete(name = 'Treatment', breaks = c('c', 'm', 'n', 's'), labels=c('Conspecific', 'Mixed', 'Control', 'Single')) + ggtitle('Atelopus Survival')
ggsave("survAtel.tiff")

## load at death
dloads <- datanc %>%
    filter(Genus=="Atelopus") %>%
        group_by(Individual, Genus) %>%
            summarise(DeathLoad = tail(log,1)) 

write.csv(dloads, "deathLoad.csv", row.names=F)

medianDeath <- dloads %>%
    group_by(Genus) %>%
        summarise(load = median(DeathLoad), sd = sd(DeathLoad))

meand <- death %>%
    filter(species=="a") %>%
        group_by(treatment) %>%
            summarise(mean = mean(days))



