rm(list=ls())

##load packages
library(dplyr)
library(reshape2)
library(reshape)
library(ggplot2)
library(lme4)
library(nlme)
library(tidyverse)
library(brms)
library(viridis)
library(dplyr)
library(gmodels)
library(ggthemes)

##read in files
dat6 = read.csv("1Master/final files/international_lambdas.csv")
head(dat6)

##metadata explanation##
# wyear - the winter year when the data were collected, e.g. Nov 2015 - March 2016 is winter year 2016
# species - four letter code indicating bat species
# country - country where sampling occured
# ysw - year since WNS arrived; an arbitrarily large number is displayed for non-US sites. This isn't used for analyses but helps to maintain consistency in column types. 
# WNS - 1 = WNS present, 0 = WNS absent
# years - the year of sampling, generally identically to wyear unless there is a December count
# x_lags - the number of years between counts
# lambda = (count/previous count) ^(1/x_lags)
# country2 = country name repeated
# ysw2 - ysw repeated
# country3 - cleaned up properly capitalized country name


prior=set_prior("normal(0,10)", class = "b")

bmod_ilam=brm(lambda~species+(1|site),data=dat6, prior = prior,chains = 4,
             control = list(adapt_delta = 0.99),cores=2,
             save_all_pars = TRUE);summary(bmod_ilam, waic = TRUE)

m<-marginal_effects(bmod_ilam,"species", probs = c(0.025, 0.975))
m

dat6$species = droplevels(dat6$species)
newdat_i=expand.grid(unique(dat6$species))
names(newdat_i)=c("species")
head(newdat_i)
fit= as.data.frame(fitted(bmod_ilam,newdata=newdat_i,re_formula=NA,type="response"))
newdat_i=bind_cols(newdat_i,fit)

head(newdat_i)
#newdat_i$lower=newdat_i$`2.5%ile`
#newdat_i$upper=newdat_i$`97.5%ile`
#depending on the version of brms, the naming of the CIs is different
#if running the most recent version of brms, the naming is as below
newdat_i$lower=newdat_i$Q2.5
newdat_i$upper=newdat_i$Q97.5

##clean up the species names
dat6$species2[dat6$species=="MYNA"]="Myotis nattereri"
dat6$species2[dat6$species=="MUHI"]="Murina hilgendorfi"
dat6$species2[dat6$species=="RHCO"]="Rhinolophus cornutus"
dat6$species2[dat6$species=="BABA"]="Barbastella barbastellus"
dat6$species2[dat6$species=="RHEU"]="Rhinolophus euryale"
dat6$species2[dat6$species=="RHHI"]="Rhinolophus hipposideros"
dat6$species2[dat6$species=="MIFU"]="Miniopterus fuliginosus"
dat6$species2[dat6$species=="PLAU"]="Plecotus auritus"
dat6$species2[dat6$species=="MYMBA"]="Myotis mys/be/al"
dat6$species2[dat6$species=="MYEM"]="Myotis emarginatus"
dat6$species2[dat6$species=="MYDAU"]="Myotis daubentonii"
dat6$species2[dat6$species=="MYBLM"]="Myotis blythii/myotis"
dat6$species2[dat6$species=="RHFE"]="Rhinolophus ferrimequinum"
dat6$species2[dat6$species=="HYAL"]="Hypsugo alaschanicus"
dat6$species2[dat6$species=="MULE"]="Murina leucogaster"
dat6$species2[dat6$species=="MYPE"]="Myotis petax"

newdat_i$species2=dat6$species2[match(newdat_i$species,dat6$species)]


#make some nicer colors for plotting
flatdesign=c("#EFC94C","#334D5C","#E27A3F","red4","#45B29D")#
dat6$country3=NA
dat6$country3[dat6$country=="china"]="China"
dat6$country3[dat6$country=="hungary"]="Hungary"
dat6$country3[dat6$country=="japan"]="Japan"
dat6$country3[dat6$country=="mongolia"]="Mongolia"
dat6$country3[dat6$country=="uk"]="United Kingdom"

library(viridis)
cols=viridis(10)
names(cols) <- levels(dat6$country3)
cols
cols["China"] = "#112F41"
cols["Japan"] = "#547980"#068587
cols["Hungary"] = "#CFD11A"#E5FCC2
cols["United Kingdom"] = "#6FB07F"#9DE0AD
cols["United States"] = "#FC5B3F"
cols["Mongolia"] = "cyan4"
cols["Israel"] = "darkolivegreen3"
cols["Georgia"] = "paleturquoise3"
colScale <- scale_colour_manual(name = "country3",values = cols)


p=ggplot(data=dat6, aes(x=species2, y=lambda))+#fill=species, 
  geom_abline(intercept = 1,slope=0,linetype="dashed", color="red")+
  geom_point(aes(color=country3),alpha=.6)+
  guides(size=FALSE)+
  geom_point(data = newdat_i,aes(x=species2, y=Estimate),color="black")+
  geom_pointrange(data=newdat_i,aes(x=species2, y=Estimate, ymin = lower, ymax = upper), 
                color="black", size=.5)+
  ylab("Population growth rate")+
  xlab("")+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  colScale+
  guides(color = guide_legend(nrow = 1, override.aes = list(size=4, alpha=1)))+
  theme_tufte()+
  theme(axis.title=element_text(size=20),axis.text.x=element_text(angle=45,hjust=1,size=8,face = "italic"),
        axis.text=element_text(size=15),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.text=element_text(size=15),legend.title=element_blank())
p

dat4 = read.csv("1Master/final files/north_america_lambdas.csv")
head(dat4)

##metadata explanation##
# wyear - the winter year when the data were collected, e.g. Nov 2015 - March 2016 is winter year 2016
# species - four letter code indicating bat species
# country - country where sampling occured
# ysw - years since WNS arrived with 0 as the winter year that WNS was first detected
# WNS - 1 = WNS present, 0 = WNS absent
# years - the year of sampling, identically to wyear (only late winter counts)
# x_lags - the number of years between counts
# lambda = (count/previous count) ^(1/x_lags)
# country2 = country name pasted together with years since WNS, e.g. us 0 is the year WNS arrived
# mlambda - this is the mean lambda of a species; it is just used for sorting the dataframe
# ysw2 - this is ysw repeated, but groups all the >=(-1) data into -1, e.g. pre-WNS arrival
# lambdaL - this is log10(lambda)
# country3 - cleaned up properly capitalized country name
# ysw3 - sampe as ysw2


prior=set_prior("normal(0,10)", class = "b")
dat4$ysw3 = as.factor(dat4$ysw3)

bmod_lam=brm(lambda~species+ysw3+(1|site),data=dat4, prior = prior,chains = 4,
             control = list(adapt_delta = 0.98),cores=2,
             save_all_pars = TRUE);summary(bmod_lam, waic = TRUE)#control = list(max_treedepth = <x>)

#####
newdat_l=expand.grid(factor(c("MYLU","MYSE","PESU","EPFU", "MYSO")),factor(c("-1","0","1","2","3")))
names(newdat_l)=c("species","ysw3")
head(newdat_l)
fit= as.data.frame(fitted(bmod_lam,newdata=newdat_l,re_formula=NA,type="response"))
newdat_l=bind_cols(newdat_l,fit)

head(newdat_l)

###aggregate data for plotting
dat4$Estimate=dat4$lambda
newdat_l1<-data.frame(newdat_l)

dat4$new_spec=1
dat4$new_spec[dat4$species=="MYSE"]="Myotis septentrionalis"
dat4$new_spec[dat4$species=="MYLU"]="Myotis lucifugus"
dat4$new_spec[dat4$species=="PESU"]="Perimyotis subflavus"
dat4$new_spec[dat4$species=="EPFU"]="Eptesicus fuscus"
dat4$new_spec[dat4$species=="MYSO"]="Myotis sodalis"

newdat_l1$new_spec[newdat_l1$species=="MYSE"]="Myotis septentrionalis"
newdat_l1$new_spec[newdat_l1$species=="MYLU"]="Myotis lucifugus"
newdat_l1$new_spec[newdat_l1$species=="PESU"]="Perimyotis subflavus"
newdat_l1$new_spec[newdat_l1$species=="EPFU"]="Eptesicus fuscus"
newdat_l1$new_spec[newdat_l1$species=="MYSO"]="Myotis sodalis"


newdat_l1$spec_ysw=paste0(newdat_l1$species,newdat_l1$ysw3)

#there is not enough data for MYSO in these years so remove them
newdat_l1=subset(newdat_l1, spec_ysw!="MYSO-1")
newdat_l1=subset(newdat_l1, spec_ysw!="MYSO0")
head(newdat_l1)

ysw_names <- c(
  "-1" = "Pre-Pd Invasion",
  "0" = "Pd Invasion Year 1",
  "1" = "Year 2 since Pd",
  "2" = "Year 3 since Pd",
  "3" = "Year 4 since Pd")

ylab=expression("Population growth rate" (lambda))
newdat_l1=newdat_l1[!(is.na(newdat_l1$ysw3)),]

names(newdat_l1)

###THIS MAKES FIGURE 1B###
p1=ggplot(data=dat4, aes(x=new_spec, y=Estimate, color=ysw3))+#fill=species, 
  facet_wrap(~ysw3,ncol=5, scales = "free_x",labeller = as_labeller(ysw_names))+#scales = "free_x",
  geom_abline(intercept = 1,slope=0,linetype="dashed", color="red")+
  geom_point(data=dat4,aes(color=ysw2),position = position_dodge(1),alpha=.5)+#size=2.5\
  geom_vline(data=filter(dat4, ysw2=="-1"|ysw2=="0"), aes(xintercept=4.5)) + 
  geom_vline(data=filter(dat4, ysw2=="1"|ysw2=="2"|ysw2=="3"), aes(xintercept=5.5)) + 
  geom_point(data=newdat_l1, aes(x=new_spec,y=Estimate),position = position_dodge(1), size=4, color="black")+#size=2.5
  geom_pointrange(data=newdat_l1,aes(x=new_spec, y=Estimate, ymin = Q2.5, ymax = Q97.5), 
                  color="black", size=.5)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  scale_colour_gradient2(low = "orange",mid="red",midpoint = 2, high = "darkorange", breaks=c(-.5,1.8),labels=c("Stable", "Declining"),
                         limits=c(-1,3))+
  ylab("Population growth rate")+
  xlab("Species")+
  theme_tufte()+
  guides(size=FALSE, colour = guide_colourbar(override.aes = list(alpha = 0.2)))+
  theme(strip.background = element_blank(),legend.key.width = unit(.4, "cm"), axis.title=element_text(size=20),strip.text = element_text(size=10),
        axis.text.x=element_text(angle=45,hjust=1,size=8,face = "italic"),axis.text=element_text(size=15),
        panel.grid = element_blank(), axis.line=element_line(),legend.position = c(0.92, 0.9),
        legend.text=element_text(size=7),legend.title=element_blank(),legend.box = "horizontal",legend.direction = "horizontal")
p1



library(gridExtra)
library(cowplot)

j2<-plot_grid(p,p1, ncol = 1, labels=c('A', 'B'),hjust = 0, label_x =.07,rel_heights = c(1/2, 1/2))
j2

ggsave(file="Fig 1_final.pdf",j2,width=9,height=10,units="in",limitsize=FALSE,useDingbats=FALSE)



