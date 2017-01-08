library(dplyr)
library(ggplot2)
library(lme4)
library(xtable)
library(rstan)

## ggplot settings
theme_set(theme_bw())
theme_update(legend.title = element_text(size=rel(1.2)), legend.text = element_text(size=13), axis.title = element_text(size=15), axis.text = element_text(size=15))

## read in data
data <- read.csv('dataForAnalysis.csv')

## just data with paired animals
datap <- droplevels(data[!is.na(data$pair),])

## just mixed speices pairs
datam <- droplevels(data[data$trt=='m',])

datam$combined <- with(datam, c(logfilt[Genus=="Atelopus"], log[Genus=='Espadarana']))

## remove controls
datanc <- droplevels(subset(data, trt!='n'))

## stan
datas <- datanc %>%
    select(samp, Week, trt, log, Genus) %>%
        na.omit()

standat <- list(id = as.numeric(datas$samp),
                Nfrog = length(unique(datas$samp)),
                N = nrow(datas),
                genus = as.numeric(datas$Genus) - 1,
                week = datas$Week,
                mixed = as.numeric(datas$trt=='m'),
                cons = as.numeric(datas$trt=='c'),
                lzsp = datas$log)

fit <- stan('regression.stan',
            data = standat,
            iter=1500,
            chains=4,
            control = list(adapt_delta = 0.9999, stepsize = 0.001, max_treedepth = 30) )


print(fit, c("mu_alpha", "mu_beta", "sigma_alpha", "sigma_beta", "alpha_g", "beta_g","alpha_c", "beta_c", "alpha_m", "beta_m", "alpha_gc", "beta_gc", "alpha_gm", "beta_gm", "sigma_y"))

pairs(fit, pars=c("mu_alpha", "mu_beta", "sigma_alpha", "sigma_beta", "alpha_g", "beta_g","alpha_c", "beta_c", "alpha_m", "beta_m", "sigma_y"))


modsimsb <- extract(fit)
newdatb <- data.frame( expand.grid( Week = seq(range(datas$Week)[1],range(datas$Week)[2], by = 1), Genus = c('Atelopus', 'Espadarana'), Treatment=c('Single','Consp','Hetero') ) )

Xmatb <- model.matrix(~Week*Genus*Treatment, data=newdatb)

parsb <- with(modsimsb, cbind(mu_alpha,
                              mu_beta,
                              alpha_g,
                              alpha_c,
                              alpha_m,
                              beta_g,
                              beta_c,
                              beta_m,
                              alpha_gc,
                              alpha_gm,
                              beta_gc,
                              beta_gm
                              ))

phat <- apply(parsb, 2, mean)
newdatb$mean <- Xmatb%*%phat

nsimb <- length(modsimsb$lp_)
lmatb <- matrix(ncol=nsimb, nrow=nrow(newdatb))
for(i in 1:nsimb){
    lmatb[,i] <- Xmatb%*%parsb[i,]
}

newdatb$lwr <- apply(lmatb, 1, quantile, prob=0.05)
newdatb$upr <- apply(lmatb, 1, quantile, prob=0.95)
newdatb$CIgroup <- with(newdatb, paste(Genus, Treatment))

ggplot(data=newdatb, aes(x=Week, y=mean)) +
    geom_line(aes(color=Genus, linetype=Treatment), size=1.2) +
        geom_ribbon(aes(ymin = lwr, ymax = upr, group=CIgroup), alpha=.2) +
            geom_point(data=datas, aes(x=Week, y=log, colour=Genus)) + 
                ylab("Log zoospores") +
                    scale_colour_discrete(labels=c(expression(italic("Atelopus")),expression(italic("Espadarana"))))

ggsave("Fig1.tiff")

## mixed with Atelopus filter and Ep swab

dataf <- datam %>%
    select(samp, Week, combined, Genus) %>%
        na.omit()

standatf <- list(id = as.numeric(dataf$samp),
                Nfrog = length(unique(dataf$samp)),
                N = nrow(dataf),
                genus = as.numeric(dataf$Genus) - 1,
                week = dataf$Week,
                lzsp = dataf$combined)

fitf <- stan('filter.stan',
             data = standatf,
             iter=1500,
             chains=4,
             control = list(adapt_delta = 0.9999, stepsize = 0.0001, max_treedepth = 40) )

print(fitf, c("mu_alpha", "mu_beta", "sigma_alpha", "sigma_beta", "alpha_g", "beta_g","sigma_y"))

## extract model
modsimsf <- extract(fitf)
newdatf <- data.frame( expand.grid( Week = seq(range(dataf$Week)[1],range(dataf$Week)[2], by = 1), Genus = c('Atelopus', 'Espadarana') ) )

Xmatf <- model.matrix(~Week*Genus, data=newdatf)

parsf <- with(modsimsf, cbind(mu_alpha,
                              mu_beta,
                              alpha_g,
                              beta_g
                              ))

phatf <- apply(parsf, 2, mean)
phatfl <- apply(parsf, 2, quantile, prob=0.05)
phatfu <- apply(parsf, 2, quantile, prob=0.95)
newdatf$mean <- Xmatf%*%phatf

nsimf <- length(modsimsf$lp_)
lmatf <- matrix(ncol=nsimf, nrow=nrow(newdatf))
for(i in 1:nsimf){
    lmatf[,i] <- Xmatf%*%parsf[i,]
}

newdatf$lwr <- apply(lmatf, 1, quantile, prob=0.05)
newdatf$upr <- apply(lmatf, 1, quantile, prob=0.95)

ggplot(data=newdatf, aes(x=Week, y=mean)) +
    geom_line(aes(linetype=Genus), size=1.2) +
        geom_ribbon(aes(ymin = lwr, ymax = upr, group=Genus), alpha=.2) +
            geom_point(data=dataf, aes(x=Week, y=combined, shape=Genus), size=2) + 
                ylab("Log zoospores") +
                    scale_linetype_discrete(labels=c(expression(italic("Atelopus")),expression(italic("Espadarana")))) +
                        scale_shape_discrete(guide=FALSE)

ggsave("Fig2.tiff")

## make parameters table

transcoefs <- data.frame(Asing = parsb[,"mu_alpha"],
                         Acons = parsb[,"mu_alpha"] + parsb[,"alpha_c"],
                         Ahet = parsb[,"mu_alpha"] + parsb[,"alpha_m"],
                         Asingslope = parsb[,"mu_beta"],
                         Aconsslope = parsb[,"mu_beta"] + parsb[,"beta_c"],
                         Ahetslope = parsb[,"mu_beta"] + parsb[,"beta_m"],
                         Esing = parsb[,"mu_alpha"] + parsb[,"alpha_g"],
                         Econs = parsb[,"mu_alpha"] + parsb[,"alpha_c"] + parsb[,"alpha_g"] + parsb[,"alpha_gc"],
                         Ehet = parsb[,"mu_alpha"] + parsb[,"alpha_m"] + parsb[,"alpha_g"] + parsb[,"alpha_gm"],
                         Esingslope = parsb[,"mu_beta"] + parsb[,"beta_g"],
                         Econsslope = parsb[,"mu_beta"] + parsb[,"beta_c"] + parsb[,"beta_g"] + parsb[,"beta_gc"],
                         Econsslope = parsb[,"mu_beta"] + parsb[,"beta_m"] + parsb[,"beta_g"] + parsb[,"beta_gm"])

tmean <- apply(transcoefs, 2, mean)
tlwr <- apply(transcoefs, 2, quantile, prob=0.05)
tupr <- apply(transcoefs, 2, quantile, prob=0.95)

restable <- data.frame(Mean=tmean, Lower=tlwr, Upper=tupr)
names(restable) <- c("Mean", expression("%5"), expression("%95"))
rownames(restable) <- c("Atelopus",
                        "Atelopus Cons",
                        "Atelopus Het",
                        "Atelopus x Week",
                        "Atelopus Cons x Week",
                        "Atelopus Het x Week",
                        "Espadarana",
                        "Espadarana Cons",
                        "Espadarana Het",
                        "Espadarana x Week",
                        "Espadarana Cons x Week",
                        "Espadarana Het x Week"
                        )

write.csv(round(restable,2), "paramTable.csv")




