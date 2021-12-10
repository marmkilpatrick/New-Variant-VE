lapply(c("ggplot2","ggthemes","tidyverse","dplyr","lme4","mgcv"),require,character.only=T) #load packages
setwd("") #set working directory
FC=read.csv("fold_figure.csv") #load VE, nAbs data

FC2=aggregate(fold_red~Variant+Study,data=FC,FUN=mean) #average across estimates w/in a study
#remove studies w/ <3 estimates
x=as.data.frame.matrix(table(FC2$Study,FC2$Variant));x$varN=rowSums(x);#rownames(x[x$varN>2,]) 
FC2=FC2[FC2$Study%in%rownames(x[x$varN>2,]),]
FC2$Variant=factor(FC2$Variant,levels=c("Beta","Gamma","Alpha","Delta"));FC2$Study=as.factor(FC2$Study) #reordering variants
FC2$figure="A";FC2$upper=NA;FC2$lower=NA

r1a=lmer(log(fold_red)~Variant+(1|Study),data = FC2);summary(r1a)
r1=lm(log(fold_red)~Variant,data = FC2);summary(r1);AIC(r1,r1a) #simpler model is preferred

mean_var=data.frame(Variant=c("Alpha","Beta","Delta","Gamma"))
mean_var$logmean=unlist(predict(r1,newdata=mean_var,se.fit = TRUE)[1])
mean_var$logse=unlist(predict(r1,newdata=mean_var,se.fit = TRUE)[2])
mean_var$logupper_bound = (mean_var$logmean + (1.96 * mean_var$logse))
mean_var$loglower_bound = (mean_var$logmean- (1.96 * mean_var$logse))
mean_var$fold_red=exp(mean_var$logmean)
mean_var$upper=exp(mean_var$logupper_bound)
mean_var$lower=exp(mean_var$loglower_bound)
mean_var$figure="B";mean_var$Study=NA
FC3=mean_var%>%
  select(Study,Variant,fold_red,upper,lower,figure)
FC_fig=(rbind(FC2,FC3))

ggplot(data=FC_fig[FC_fig$figure=="A",],aes(y=fold_red, x=Variant,color=Study,group=Study))+
  geom_point(size=3)+ geom_line()+
  scale_y_continuous(trans='log2',breaks = c(1,2,4,8,16,32))+theme_few()+
  geom_point(data=FC_fig[FC_fig$figure=="B",],aes(y=fold_red, x=Variant),position=position_nudge(x = 0.05),col="black",size=3,alpha=1.,shape=15)+
  geom_errorbar(data=FC_fig[FC_fig$figure=="B",],aes(x=Variant,ymin=lower,ymax=upper),position=position_nudge(x = 0.05),width=0,alpha=1.)+
  xlab("Variant")+  ylab("Ratio of reduction in neutralizing antibody titer to WT")+
  theme(axis.title=element_text(size=25),legend.position = c(.75, .720),#legend.position = "top",
        axis.text=element_text(size=25),legend.text=element_text(size=20),
        legend.title=element_text(size=20))

#standard errors and CVs for text
(FC3$upper-FC3$fold_red)/1.96
((FC3$upper-FC3$fold_red)/1.96)/FC3$fold_red

FC3[5,]=c("","D614G",1,1,1,"B") #added for later merge

#Hospital graph
m_sev_all=read.csv("severe_disease.csv") #load VE, nAbs data
m_sev_all=m_sev_all[m_sev_all$all_time=="1",] #removing limited follow-up windows
m_sev_all$Variant=m_sev_all$dom_var 
m_sev_all=merge(m_sev_all,FC3,by="Variant") #merging VE data w/ Nabs data
m_sev_all$fold_red=as.numeric(m_sev_all$fold_red) 
m_sev_all$rel_Nabs=m_sev_all$Nabs_ratio/m_sev_all$fold_red #adjusting Nabs for variants

#remove time window and age specific data from main data frame
m_hosp=m_sev_all[m_sev_all$dom_var%in%c("Alpha","Beta","Delta","Gamma")&m_sev_all$all_time==1&m_sev_all$population=="gen pop"&m_sev_all$endpoint%in%c("Hospitalization","Hospitalization or Death"),]
m_hosp_p=m_sev_all[m_sev_all$vaccine=="Pfizer"&m_sev_all$dom_var%in%c("Alpha","Beta","Delta","Gamma")&m_sev_all$all_time==1&m_sev_all$population=="gen pop"&m_sev_all$endpoint%in%c("Hospitalization","Hospitalization or Death"),]

v0=glm(VE~log2(rel_Nabs),data=m_hosp,family=binomial,weights=eff_sample);summary(v0)
m_hosp$fold_red_jit=m_hosp$fold_red+rnorm(nrow(m_hosp),0,0.03) #jitter xvals
m_hosp$Nabs_ratio_jit=m_hosp$Nabs_ratio+rnorm(nrow(m_hosp),0,0.03) #jitter xvals
m_hosp$rel_nabs_jit=m_hosp$rel_Nabs+rnorm(nrow(m_hosp),0,0.03) #jitter xvals

ggplot()+
  geom_point(data=m_hosp,aes(x=fold_red_jit,y=VE,color=Variant),size=3)+
  geom_errorbar(data=m_hosp,aes(x=fold_red_jit,ymin=lower.x,ymax=upper.x,color=Variant))+
  # geom_line(data=mL_hosp,aes(x=fold_red,y=VE))+
  # geom_ribbon(data=mL_hosp,aes(x=fold_red,ymin=lower,ymax=upper),fill="gray",alpha=0.5)+
  scale_x_continuous(trans="log2",breaks=c(1:8,16,32,40),limits=c(1.5,40))+
  theme_few()+scale_shape_manual(values = c(15:18))+ 
  xlab("Fold reduction in neutralizing antibody titers")+
  ylab("VE against hospitalization")+  labs(col="Variant")+
  theme(axis.title=element_text(size=20),#legend.position = c(.8, .750),#legend.position = "top",
        axis.text=element_text(size=20),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  facet_wrap(.~vaccine)+theme(strip.text.x = element_text(size = 15))

#Pfizer waning ratios relative to 1 mo from dose 2 from Falsey et al 2021 Fig 1A WT, age weighted (71%, 29%)
wane_ratio=1/(.71*497/83+.29*538/41)
boost_ratio=.71*1754/497+.29*1318/538
#Moderna boosting @28d relative to post 2nd dose from Chu et al medRxiv 0,181, 209: 1268, 150.2, 1951.7
Mod_boost_ratio=1951.7/1268
Mod_wane_ratio=150.2/1268
Nabs_Pfiz_WT=2.3723404;Nabs_Mod_WT=4.1332912

#Predicted hospitalized omicron
p_h_o=data.frame(rel_Nabs=c(Nabs_Pfiz_WT/40*wane_ratio,Nabs_Pfiz_WT/40,Nabs_Pfiz_WT/40*boost_ratio,Nabs_Pfiz_WT/2.744457*wane_ratio,Nabs_Pfiz_WT/2.744457,Nabs_Pfiz_WT/2.744457*boost_ratio,Nabs_Pfiz_WT*wane_ratio,Nabs_Pfiz_WT,Nabs_Pfiz_WT*boost_ratio,
                            4.1332912/40*Mod_wane_ratio,4.1332912/40,4.1332912/40*Mod_boost_ratio,4.1332912/2.744457*Mod_wane_ratio,4.1332912/2.744457,4.1332912/2.744457*Mod_boost_ratio,4.1332912*Mod_wane_ratio,4.1332912,4.1332912*Mod_boost_ratio),
                 VE=predict(v0,newdata=data.frame(rel_Nabs=c(Nabs_Pfiz_WT/40*wane_ratio,Nabs_Pfiz_WT/40,Nabs_Pfiz_WT/40*boost_ratio,Nabs_Pfiz_WT/2.744457*wane_ratio,Nabs_Pfiz_WT/2.744457,Nabs_Pfiz_WT/2.744457*boost_ratio,Nabs_Pfiz_WT*wane_ratio,Nabs_Pfiz_WT,Nabs_Pfiz_WT*boost_ratio,
                                                             4.1332912/40*Mod_wane_ratio,4.1332912/40,4.1332912/40*Mod_boost_ratio,4.1332912/2.744457*Mod_wane_ratio,4.1332912/2.744457,4.1332912/2.744457*Mod_boost_ratio,4.1332912*Mod_wane_ratio,4.1332912,4.1332912*Mod_boost_ratio)),type="response")
                 ,Vaccine=rep(c("Pfizer","Moderna"),each=9),Variant=rep(c("Omicron","Omicron","Omicron","Delta","Delta","Delta","Wild type","Wild type","Wild type"),by=2),Dose=rep(c("Waned","Two dose","Boosted"),by=6))
p_h_o$VE_logit=unlist(predict(v0,newdata=p_h_o, type="link", se.fit = TRUE)[1])
p_h_o$VE_se=unlist(predict(v0,newdata=p_h_o, type="link", se.fit = TRUE)[2])
p_h_o$upper = plogis(p_h_o$VE_logit + (1.96 * p_h_o$VE_se))
p_h_o$lower = plogis(p_h_o$VE_logit - (1.96 * p_h_o$VE_se))

#max fold Nabs range
max(m_hosp$rel_Nabs)/min(m_hosp$rel_Nabs)

#Symptomatic disease
m_symp=read.csv("symp_disease.csv") #load VE, nAbs data
m_symp=m_symp[m_symp$all_time=="1",]
m_symp=m_symp[m_symp$vaccine!="Johnson and Johnson",] #remove J&J - only 1 estimate
m_symp$lower[m_symp$vaccine=="Astrazeneca"&m_symp$Variant=="Beta"]=0.01 #reset lower CI of AZ-Beta to plot it
m_symp=merge(m_symp,FC3,by="Variant")
m_symp$fold_red=as.numeric(m_symp$fold_red)
m_symp$rel_Nabs=m_symp$Nabs_ratio/m_symp$fold_red

m_symp_p=m_symp[m_symp$vaccine=="Pfizer",]

v1=glm(VE~log2(rel_Nabs),data=m_symp,family=binomial,weights=eff_sample);summary(v1)

#predicted values for symptomatic disease for Omicron
p_s_o=data.frame(rel_Nabs=c(Nabs_Pfiz_WT/40*wane_ratio,Nabs_Pfiz_WT/40,Nabs_Pfiz_WT/40*boost_ratio,Nabs_Pfiz_WT/2.744457*wane_ratio,Nabs_Pfiz_WT/2.744457,Nabs_Pfiz_WT/2.744457*boost_ratio,Nabs_Pfiz_WT*wane_ratio,Nabs_Pfiz_WT,Nabs_Pfiz_WT*boost_ratio,
                            4.1332912/40*Mod_wane_ratio,4.1332912/40,4.1332912/40*Mod_boost_ratio,4.1332912/2.744457*Mod_wane_ratio,4.1332912/2.744457,4.1332912/2.744457*Mod_boost_ratio,4.1332912*Mod_wane_ratio,4.1332912,4.1332912*Mod_boost_ratio),
                 VE=predict(v1,newdata=data.frame(rel_Nabs=c(Nabs_Pfiz_WT/40*wane_ratio,Nabs_Pfiz_WT/40,Nabs_Pfiz_WT/40*boost_ratio,Nabs_Pfiz_WT/2.744457*wane_ratio,Nabs_Pfiz_WT/2.744457,Nabs_Pfiz_WT/2.744457*boost_ratio,Nabs_Pfiz_WT*wane_ratio,Nabs_Pfiz_WT,Nabs_Pfiz_WT*boost_ratio,
                                                             4.1332912/40*Mod_wane_ratio,4.1332912/40,4.1332912/40*Mod_boost_ratio,4.1332912/2.744457*Mod_wane_ratio,4.1332912/2.744457,4.1332912/2.744457*Mod_boost_ratio,4.1332912*Mod_wane_ratio,4.1332912,4.1332912*Mod_boost_ratio)),type="response")
                 ,Vaccine=rep(c("Pfizer","Moderna"),each=9),Variant=rep(c("Omicron","Omicron","Omicron","Delta","Delta","Delta","Wild type","Wild type","Wild type"),by=2),Dose=rep(c("Waned","Two dose","Boosted"),by=6))
p_s_o$VE_logit=unlist(predict(v1,newdata=p_s_o, type="link", se.fit = TRUE)[1])
p_s_o$VE_se=unlist(predict(v1,newdata=p_s_o, type="link", se.fit = TRUE)[2])
p_s_o$upper = plogis(p_s_o$VE_logit + (1.96 * p_s_o$VE_se))
p_s_o$lower = plogis(p_s_o$VE_logit - (1.96 * p_s_o$VE_se))

v3=glm(VE~log2(fold_red)*vaccine,data=m_symp[!m_symp$vaccine%in%c("Sinovac","Sputnik V "),],
       family=binomial,weights=eff_sample);summary(v3)
fr=seq(1.5,40,by=0.1)
mL_symp=data.frame(fold_red=rep(fr,length(unique(m_symp[!m_symp$vaccine%in%c("Sinovac","Sputnik V "),]$vaccine))),
                   vaccine=rep(unique(m_symp[!m_symp$vaccine%in%c("Sinovac","Sputnik V "),]$vaccine),each=length(fr)))
mL_symp$VE=predict(v3,newdata=mL_symp,type="response")
mL_symp$VE_logit=unlist(predict(v3,newdata=mL_symp, type="link", se.fit = TRUE)[1])
mL_symp$VE_se=unlist(predict(v3,newdata=mL_symp, type="link", se.fit = TRUE)[2])
mL_symp$upper = plogis(mL_symp$VE_logit + (1.96 * mL_symp$VE_se))
mL_symp$lower = plogis(mL_symp$VE_logit- (1.96 * mL_symp$VE_se))

m_symp$fold_red_jit=m_symp$fold_red+runif(nrow(m_symp),-.03,0.03) #jitter xvals
m_symp$Nabs_ratio_jit=m_symp$Nabs_ratio+rnorm(nrow(m_symp),0,0.03) #jitter xvals
m_symp$rel_nabs_jit=m_symp$rel_Nabs+rnorm(nrow(m_symp),0,0.03) #jitter xvals

ggplot(data=m_symp[!m_symp$vaccine%in%c("Sinovac","Sputnik V "),])+
  geom_point(aes(x=fold_red_jit,y=VE,color=Variant),size=3)+
  geom_errorbar(aes(x=fold_red_jit,ymin=lower.x,ymax=upper.x,color=Variant))+
  # geom_line(data=mL_symp,aes(x=fold_red,y=VE))+
  # geom_ribbon(data=mL_symp,aes(x=fold_red,ymin=lower,ymax=upper),fill="gray",alpha=0.5)+
  scale_x_continuous(trans="log2",breaks=c(1:8,16,32,40),limits=c(.95,40))+
  scale_y_continuous(breaks=c(seq(0,1,by=.2)),limits=c(0,1))+
  theme_few()+scale_shape_manual(values = c(15:18))+ 
  xlab("Reduction in neutralizing antibody titers")+
  ylab("VE against symptomatic disease")+  labs(col="Variant")+
  facet_wrap(.~vaccine)+theme(strip.text.x = element_text(size = 15))+
  theme(axis.title=element_text(size=20),#legend.position = c(.8, .750),#legend.position = "top",
        axis.text=element_text(size=20),legend.text=element_text(size=15),
        legend.title=element_text(size=15))

#Documented infection
m_inf=read.csv("doc_inf.csv") #load VE, nAbs data
m_inf$eff_sample[m_inf$eff_sample>500]=500

m_inf=merge(m_inf,FC3,by="Variant")
m_inf$fold_red=as.numeric(m_inf$fold_red)
m_inf$rel_Nabs=m_inf$Nabs_ratio/m_inf$fold_red

m_inf_p=m_inf[m_inf$vaccine=="Pfizer",]

m_inf$fold_red_jit=m_inf$fold_red+rnorm(nrow(m_inf),0,0.03) #jitter xvals
m_inf$Nabs_ratio_jit=m_inf$Nabs_ratio+rnorm(nrow(m_inf),0,0.03) #jitter xvals
m_inf$rel_nabs_jit=m_inf$rel_Nabs+rnorm(nrow(m_inf),0,0.03) #jitter xvals

v6=glm(VE~log2(rel_Nabs),data=m_inf,family=binomial,weights=eff_sample);summary(v6)

ggplot(data=m_inf)+
  geom_point(aes(x=fold_red_jit,y=VE,color=Variant),size=3)+
  geom_errorbar(aes(x=fold_red_jit,ymin=lower.x,ymax=upper.x,color=Variant))+
  # geom_line(data=mL_inf,aes(x=fold_red,y=VE))+
  # geom_ribbon(data=mL_inf,aes(x=fold_red,ymin=lower,ymax=upper),fill="gray",alpha=0.5)+
  theme_few()+labs(col="Variant")+
  scale_x_continuous(trans="log2",breaks=c(1:8,16,32,40),limits=c(1.5,40))+
  theme_few()+scale_shape_manual(values = c(15:18))+ 
  xlab("Reduction in neutralizing antibody titers")+
  ylab("VE against documented infection")+  labs(col="Variant")+
  facet_wrap(.~vaccine,nrow=2)+theme(strip.text.x = element_text(size = 20))+
  theme(axis.title=element_text(size=20),legend.position = c(.9, .150),#legend.position = "top",
        axis.text=element_text(size=20),legend.text=element_text(size=20),#legend.position = "right",
        legend.title=element_text(size=20))

v2=glm(VE~log2(rel_Nabs),data=m_inf,family=binomial,weights=eff_sample);summary(v2)

p_i_o=data.frame(rel_Nabs=c(Nabs_Pfiz_WT/40*wane_ratio,Nabs_Pfiz_WT/40,Nabs_Pfiz_WT/40*boost_ratio,Nabs_Pfiz_WT/2.744457*wane_ratio,Nabs_Pfiz_WT/2.744457,Nabs_Pfiz_WT/2.744457*boost_ratio,Nabs_Pfiz_WT*wane_ratio,Nabs_Pfiz_WT,Nabs_Pfiz_WT*boost_ratio,
                            4.1332912/40*Mod_wane_ratio,4.1332912/40,4.1332912/40*Mod_boost_ratio,4.1332912/2.744457*Mod_wane_ratio,4.1332912/2.744457,4.1332912/2.744457*Mod_boost_ratio,4.1332912*Mod_wane_ratio,4.1332912,4.1332912*Mod_boost_ratio),
                 VE=predict(v2,newdata=data.frame(rel_Nabs=c(Nabs_Pfiz_WT/40*wane_ratio,Nabs_Pfiz_WT/40,Nabs_Pfiz_WT/40*boost_ratio,Nabs_Pfiz_WT/2.744457*wane_ratio,Nabs_Pfiz_WT/2.744457,Nabs_Pfiz_WT/2.744457*boost_ratio,Nabs_Pfiz_WT*wane_ratio,Nabs_Pfiz_WT,Nabs_Pfiz_WT*boost_ratio,
                                                             4.1332912/40*Mod_wane_ratio,4.1332912/40,4.1332912/40*Mod_boost_ratio,4.1332912/2.744457*Mod_wane_ratio,4.1332912/2.744457,4.1332912/2.744457*Mod_boost_ratio,4.1332912*Mod_wane_ratio,4.1332912,4.1332912*Mod_boost_ratio)),type="response")
                 ,Vaccine=rep(c("Pfizer","Moderna"),each=9),Variant=rep(c("Omicron","Omicron","Omicron","Delta","Delta","Delta","Wild type","Wild type","Wild type"),by=2),Dose=rep(c("Waned","Two dose","Boosted"),by=6))
p_i_o$VE_logit=unlist(predict(v2,newdata=p_i_o, type="link", se.fit = TRUE)[1])
p_i_o$VE_se=unlist(predict(v2,newdata=p_i_o, type="link", se.fit = TRUE)[2])
p_i_o$upper = plogis(p_i_o$VE_logit + (1.96 * p_i_o$VE_se))
p_i_o$lower = plogis(p_i_o$VE_logit - (1.96 * p_i_o$VE_se))

m_hosp$endpoint="Hospitalization"
comp_hi=rbind(m_hosp%>%select("rel_Nabs","endpoint","VE","eff_sample"),m_inf%>%select("rel_Nabs","endpoint","VE","eff_sample"))

#All infetions
mAll=read.csv("VE_nAbs_data.csv") #load VE, nAbs data
m1=read.csv("Nabs ratios Delta WT.csv") #load VE, nAbs data

#Create separate data frame for transmission data; #remove transmission data from main data frame
mAllt=mAll[mAll$Groups=="Delta - transmission",]; mAll=mAll[mAll$Groups!="Delta - transmission",]

#Pfizer waning ratios from Falsey et al 2021 Fig 1A WT, age weighted (71%, 29%)
wane_ratio=1/(.71*497/83+.29*538/41)
boost_ratio=.71*1754/497+.29*1318/538
boost_ratio/wane_ratio
hybrid_ratio=6300/800 #Goel et al
#Moderna boosting @28d relative to post 2nd dose from Chu et al medRxiv 0,181, 209: 1268, 150.2, 1951.7
Mod_boost_ratio=1951.7/1268
Mod_wane_ratio=150.2/1268

Nabs_Pfiz_delta_inf=mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"]


#Table S3: Fitted model for protection vs nAbs
f2=glm(VE~log2(Nabs_Ratio),data=mAll[mAll$Vacc_avg==1,],family=binomial,weights=N_eff);summary(f2)

p_a_o=data.frame(Nabs_Ratio=c(Nabs_Pfiz_WT/40*wane_ratio,Nabs_Pfiz_WT/40,Nabs_Pfiz_WT/40*boost_ratio,Nabs_Pfiz_WT/2.744457*wane_ratio,Nabs_Pfiz_WT/2.744457,Nabs_Pfiz_WT/2.744457*boost_ratio,Nabs_Pfiz_WT*wane_ratio,Nabs_Pfiz_WT,Nabs_Pfiz_WT*boost_ratio,
                              4.1332912/40*Mod_wane_ratio,4.1332912/40,4.1332912/40*Mod_boost_ratio,4.1332912/2.744457*Mod_wane_ratio,4.1332912/2.744457,4.1332912/2.744457*Mod_boost_ratio,4.1332912*Mod_wane_ratio,4.1332912,4.1332912*Mod_boost_ratio),
                 VE=predict(f2,newdata=data.frame(Nabs_Ratio=c(Nabs_Pfiz_WT/40*wane_ratio,Nabs_Pfiz_WT/40,Nabs_Pfiz_WT/40*boost_ratio,Nabs_Pfiz_WT/2.744457*wane_ratio,Nabs_Pfiz_WT/2.744457,Nabs_Pfiz_WT/2.744457*boost_ratio,Nabs_Pfiz_WT*wane_ratio,Nabs_Pfiz_WT,Nabs_Pfiz_WT*boost_ratio,
                                                               4.1332912/40*Mod_wane_ratio,4.1332912/40,4.1332912/40*Mod_boost_ratio,4.1332912/2.744457*Mod_wane_ratio,4.1332912/2.744457,4.1332912/2.744457*Mod_boost_ratio,4.1332912*Mod_wane_ratio,4.1332912,4.1332912*Mod_boost_ratio)),type="response")
                 ,Vaccine=rep(c("Pfizer","Moderna"),each=9),Variant=rep(c("Omicron","Omicron","Omicron","Delta","Delta","Delta","Wild type","Wild type","Wild type"),by=2),Dose=rep(c("Waned","Two dose","Boosted"),by=6))
p_a_o$VE_logit=unlist(predict(f2,newdata=p_a_o, type="link", se.fit = TRUE)[1])
p_a_o$VE_se=unlist(predict(f2,newdata=p_a_o, type="link", se.fit = TRUE)[2])
p_a_o$upper = plogis(p_a_o$VE_logit + (1.96 * p_a_o$VE_se))
p_a_o$lower = plogis(p_a_o$VE_logit - (1.96 * p_a_o$VE_se))

#Transmission model
f3=glm(VE~log2(Nabs_Ratio),data=mAllt,family=binomial,weights=N_eff);summary(f3)

p_t_o=data.frame(Nabs_Ratio=c(Nabs_Pfiz_WT/40*wane_ratio,Nabs_Pfiz_WT/40,Nabs_Pfiz_WT/40*boost_ratio,Nabs_Pfiz_WT/2.744457*wane_ratio,Nabs_Pfiz_WT/2.744457,Nabs_Pfiz_WT/2.744457*boost_ratio,Nabs_Pfiz_WT*wane_ratio,Nabs_Pfiz_WT,Nabs_Pfiz_WT*boost_ratio,
                              4.1332912/40*Mod_wane_ratio,4.1332912/40,4.1332912/40*Mod_boost_ratio,4.1332912/2.744457*Mod_wane_ratio,4.1332912/2.744457,4.1332912/2.744457*Mod_boost_ratio,4.1332912*Mod_wane_ratio,4.1332912,4.1332912*Mod_boost_ratio),
                 VE=predict(f3,newdata=data.frame(Nabs_Ratio=c(Nabs_Pfiz_WT/40*wane_ratio,Nabs_Pfiz_WT/40,Nabs_Pfiz_WT/40*boost_ratio,Nabs_Pfiz_WT/2.744457*wane_ratio,Nabs_Pfiz_WT/2.744457,Nabs_Pfiz_WT/2.744457*boost_ratio,Nabs_Pfiz_WT*wane_ratio,Nabs_Pfiz_WT,Nabs_Pfiz_WT*boost_ratio,
                                                               4.1332912/40*Mod_wane_ratio,4.1332912/40,4.1332912/40*Mod_boost_ratio,4.1332912/2.744457*Mod_wane_ratio,4.1332912/2.744457,4.1332912/2.744457*Mod_boost_ratio,4.1332912*Mod_wane_ratio,4.1332912,4.1332912*Mod_boost_ratio)),type="response")
                 ,Vaccine=rep(c("Pfizer","Moderna"),each=9),Variant=rep(c("Omicron","Omicron","Omicron","Delta","Delta","Delta","Wild type","Wild type","Wild type"),by=2),Dose=rep(c("Waned","Two dose","Boosted"),by=6))
p_t_o$VE_logit=unlist(predict(f3,newdata=p_t_o, type="link", se.fit = TRUE)[1])
p_t_o$VE_se=unlist(predict(f3,newdata=p_t_o, type="link", se.fit = TRUE)[2])
p_t_o$upper = plogis(p_t_o$VE_logit + (1.96 * p_t_o$VE_se))
p_t_o$lower = plogis(p_t_o$VE_logit - (1.96 * p_t_o$VE_se))

nAbr=seq(.04,4.5,by=0.1)
hos_symp1=data.frame(rel_Nabs=nAbr,figure="A")
hos_symp1$VE=predict(v0,newdata=data.frame(rel_Nabs=nAbr),type="response")
hos_symp1$VE_logit=unlist(predict(v0,newdata=hos_symp1, type="link", se.fit = TRUE)[1])
hos_symp1$VE_se=unlist(predict(v0,newdata=hos_symp1, type="link", se.fit = TRUE)[2])
hos_symp1$upper = plogis(hos_symp1$VE_logit + (1.96 * hos_symp1$VE_se))
hos_symp1$lower = plogis(hos_symp1$VE_logit- (1.96 * hos_symp1$VE_se))

hos_symp2=data.frame(rel_Nabs=nAbr,figure="B")
hos_symp2$VE=predict(v1,newdata=data.frame(rel_Nabs=nAbr),type="response")
hos_symp2$VE_logit=unlist(predict(v1,newdata=hos_symp2, type="link", se.fit = TRUE)[1])
hos_symp2$VE_se=unlist(predict(v1,newdata=hos_symp2, type="link", se.fit = TRUE)[2])
hos_symp2$upper = plogis(hos_symp2$VE_logit + (1.96 * hos_symp2$VE_se))
hos_symp2$lower = plogis(hos_symp2$VE_logit- (1.96 * hos_symp2$VE_se))

hos_symp3=data.frame(rel_Nabs=nAbr,figure="C")
hos_symp3$VE=predict(v6,newdata=data.frame(rel_Nabs=nAbr),type="response")
hos_symp3$VE_logit=unlist(predict(v6,newdata=hos_symp3, type="link", se.fit = TRUE)[1])
hos_symp3$VE_se=unlist(predict(v6,newdata=hos_symp3, type="link", se.fit = TRUE)[2])
hos_symp3$upper = plogis(hos_symp3$VE_logit + (1.96 * hos_symp3$VE_se))
hos_symp3$lower = plogis(hos_symp3$VE_logit- (1.96 * hos_symp3$VE_se))

hos_symp=rbind(hos_symp1,hos_symp2,hos_symp3)

pred_fig=rbind(p_h_o[c(2,3,11,12),],p_s_o[c(2,3,11,12),],p_i_o[c(2,3,11,12),])
pred_fig$figure=rep(c("A","B","C"),each=4)
pred_fig=rbind(p_h_o[c(2,3,11,12),],p_s_o[c(2,3,11,12),],p_i_o[c(2,3,11,12),])
pred_fig$figure=rep(c("A","B","C"),each=4)
pred_fig$vaccine=pred_fig$Vaccine

m_inf$figure="C";m_symp$figure="B";m_hosp$figure="A"
tlab2=data.frame(figure=c("A","B","C"),text1=c("A: Hospitalization","B: Symptomatic disease","C: Documented infection"),rel_Nabs=c(1.8,1.8,1.8),VE=c(0.4,0.12,.42))

m_hosp$rel_nabs_jit=m_hosp$rel_Nabs+rnorm(nrow(m_hosp),0,0.03) #jitter xvals

ggplot(data=hos_symp)+
  geom_point(data=pred_fig,aes(x=rel_Nabs,y=VE,color=Variant,shape=vaccine),size=3)+
  geom_errorbar(data=pred_fig,aes(x=rel_Nabs,ymin=lower,ymax=upper,color=Variant),width=0)+
  geom_point(data=m_symp,aes(x=rel_nabs_jit,y=VE,color=Variant,shape=vaccine),size=3)+
  geom_errorbar(data=m_symp,aes(x=rel_nabs_jit,ymin=lower.x,ymax=upper.x,color=Variant),width=0)+
  geom_point(data=m_inf,aes(x=rel_nabs_jit,y=VE,color=Variant,shape=vaccine),size=3)+
  geom_errorbar(data=m_inf,aes(x=rel_nabs_jit,ymin=lower.x,ymax=upper.x,color=Variant),width=0)+
  geom_point(data=m_hosp,aes(x=rel_nabs_jit,y=VE,color=Variant,shape=vaccine),size=3)+
  geom_errorbar(data=m_hosp,aes(x=rel_nabs_jit,ymin=lower.x,ymax=upper.x,color=Variant),width=0)+
  geom_line(data=hos_symp[hos_symp$figure=="A",],aes(x=rel_Nabs,y=VE))+
  geom_ribbon(data=hos_symp[hos_symp$figure=="A",],aes(x=rel_Nabs,ymin=lower,ymax=upper),fill="gray",alpha=0.5)+
  geom_line(data=hos_symp[hos_symp$figure=="B",],aes(x=rel_Nabs,y=VE))+
  geom_ribbon(data=hos_symp[hos_symp$figure=="B",],aes(x=rel_Nabs,ymin=lower,ymax=upper),fill="gray",alpha=0.5)+
  geom_line(data=hos_symp[hos_symp$figure=="C",],aes(x=rel_Nabs,y=VE))+
  geom_ribbon(data=hos_symp[hos_symp$figure=="C",],aes(x=rel_Nabs,ymin=lower,ymax=upper),fill="gray",alpha=0.5)+
  scale_x_continuous(trans="log2",limits=c(0.04,4.5),breaks=c(0.0625*2^c(0:7)))+
  theme_few()+scale_shape_manual(values = c(15:21))+
  xlab("Ratio of antibody neutralization relative to convalescent/WT")+
  ylab("Vaccine effectiveness")+  labs(col="Variant", shape="Vaccine")+
  theme(axis.title=element_text(size=25),#legend.position = c(.8, .750),#legend.position = "top",
        axis.text=element_text(size=25),legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  facet_wrap(.~figure,nrow=3, scales="free_y")+theme(strip.text.x = element_blank())+
  geom_text(data=tlab2,aes(x=rel_Nabs,y=VE,label=text1),size=6,hjust=0)

