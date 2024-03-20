#install.packages(c("gtools","merDeriv","brms","mvtnorm","brglm2","ggplot2","ggthemes","tidyverse","dplyr","lme4","mgcv","bbmle","RColorBrewer","viridis")) #install packages
lapply(c("lmerTest","scales","gtools","merDeriv","brms","mvtnorm","brglm2","ggplot2","ggthemes","tidyverse","dplyr","lme4","mgcv","bbmle","RColorBrewer","viridis"),require,character.only=T) #load packages

#load VE, nAbs data
m=read.csv("Omicron_data_long2.csv") 
m$ID=as.factor(m$ID)

#add vaccine codes
vac_comp=data.frame(vaccine=c("Astrazeneca","Johnson and Johnson","Moderna","Novavax","Pfizer","Sinovac","Sputnik V"),
                    code=c("ChAdOx1 nCoV-19","Ad26.COV2.S","mRNA-1273","NVX-CoV2373","BNT162b2","CoronaVac","Sputnik V"))

m=merge(m,vac_comp,by="vaccine")

#Begin original fold change figure: Load data
FC=read.csv("fold_figure_all.csv") #load VE, nAbs data
FC$Variant[FC$Variant=="Omicron (initial)"]="BA.1 Dec. 2021"

FC2=aggregate(fold_red~Variant+Study,data=FC,FUN=mean) #average across estimates w/in a study
#remove studies w/ <3 estimates
x=as.data.frame.matrix(table(FC2$Study,FC2$Variant));
x$`BA.1 Dec. 2021`=3*x$`BA.1 Dec. 2021`
x$varN=rowSums(x);#rownames(x[x$varN>2,]) 
FC2=FC2[FC2$Study%in%rownames(x[x$varN>1,]),]
FC2$Variant=factor(FC2$Variant,levels=c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021"));FC2$Study=as.factor(FC2$Study) #reordering variants
FC2$figure="A";FC2$upper=NA;FC2$lower=NA

r1a=lmer(log(fold_red)~Variant+(1|Study),data = FC2);summary(r1a)
r1=lm(log(fold_red)~Variant,data = FC2);summary(r1);AIC(r1,r1a) #random effects model is preferred

#Means and SEs for each variant
mean_varO=data.frame(Variant=c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021"),logmean=NA,logse=NA,
                     logupper_bound=NA,loglower_bound=NA,fold_red=NA,figure="B",upper=NA,lower=NA,Study=NA)

coefs = fixef(r1a)
V = vcov(r1a)

DF.O = data.frame(Variant = factor("Alpha", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.O = model.matrix(~ Variant, data = DF.O)
mean_varO$logmean[1]=c(X.O %*% coefs)
mean_varO$logse[1]=sqrt(diag(X.O %*% V %*% t(X.O)))

DF.1 = data.frame(Variant = factor("Gamma", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.1 = model.matrix(~ Variant, data = DF.1)
mean_varO$logmean[2]=c(X.1 %*% coefs)
mean_varO$logse[2]=sqrt(diag(X.1 %*% V %*% t(X.1)))

DF.2 = data.frame(Variant = factor("Delta", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.2 = model.matrix(~ Variant, data = DF.2)
mean_varO$logmean[3]=c(X.2 %*% coefs)
mean_varO$logse[3]=sqrt(diag(X.2 %*% V %*% t(X.2)))

DF.3 = data.frame(Variant = factor("Beta", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.3 = model.matrix(~ Variant, data = DF.3)
mean_varO$logmean[4]=c(X.3 %*% coefs)
mean_varO$logse[4]=sqrt(diag(X.3 %*% V %*% t(X.3)))

DF.4 = data.frame(Variant = factor("BA.1", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.4 = model.matrix(~ Variant, data = DF.4)
mean_varO$logmean[5]=c(X.4 %*% coefs)
mean_varO$logse[5]=sqrt(diag(X.4 %*% V %*% t(X.4)))

DF.5 = data.frame(Variant = factor("BA.2", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.5 = model.matrix(~ Variant, data = DF.5)
mean_varO$logmean[6]=c(X.5 %*% coefs)
mean_varO$logse[6]=sqrt(diag(X.5 %*% V %*% t(X.5)))

DF.6 = data.frame(Variant = factor("BA.4/5", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.6 = model.matrix(~ Variant, data = DF.6)
mean_varO$logmean[7]=c(X.6 %*% coefs)
mean_varO$logse[7]=sqrt(diag(X.6 %*% V %*% t(X.6)))

DF.7 = data.frame(Variant = factor("BA.1 Dec. 2021", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.7 = model.matrix(~ Variant, data = DF.7)
mean_varO$logmean[8]=c(X.7 %*% coefs)
mean_varO$logse[8]=sqrt(diag(X.7 %*% V %*% t(X.7)))

mean_varO$logupper_bound = (mean_varO$logmean + (1.96 * mean_varO$logse))
mean_varO$loglower_bound = (mean_varO$logmean- (1.96 * mean_varO$logse))
mean_varO$fold_red=exp(mean_varO$logmean)
mean_varO$upper=exp(mean_varO$logupper_bound)
mean_varO$lower=exp(mean_varO$loglower_bound)

FC3=mean_varO[c("Study","Variant","fold_red","upper","lower","figure")]
FC_fig=(rbind(FC2,FC3))

FC_fig$Variant=factor(FC_fig$Variant,levels=c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021"));FC_fig$Study=as.factor(FC_fig$Study)

#standard errors and CVs for text
(FC3$upper-FC3$fold_red)/1.96
((FC3$upper-FC3$fold_red)/1.96)/FC3$fold_red

FC3[9,]=c("","D614G",1,1,1,"B") #added for later merge
FC3$fold_red=as.numeric(FC3$fold_red);FC3$lower=as.numeric(FC3$lower);FC3$upper=as.numeric(FC3$upper);

#Merge fold change data with main dataframe
m=merge(m,FC3[,c("Variant","fold_red")],by="Variant") #merging VE data w/ Nabs data
m$fold_red=as.numeric(m$fold_red) 
m$rel_Nabs=m$Nabs_ratio/m$fold_red #adjusting Nabs for variants

#Jitter values for plots
m$fold_red_jit=m$fold_red+rnorm(nrow(m),0,0.03) #jitter xvals
m$Nabs_ratio_jit=m$Nabs_ratio+rnorm(nrow(m),0,0.03) #jitter xvals
m$rel_nabs_jit=m$rel_Nabs+rnorm(nrow(m),0,0.03) #jitter xvals

#Fig 1 
ggplot(data=FC_fig[FC_fig$figure=="A",],aes(y=fold_red, x=Variant,color=Study,group=Study))+
  geom_point(size=3,alpha=0.8)+ geom_line(alpha=0.8)+
  scale_y_continuous(trans='log2',breaks = c(1,2,4,8,16,32,64),labels=c("1","1/2","1/4","1/8","1/16","1/32","1/64"))+theme_few()+
  geom_point(data=FC_fig[FC_fig$figure=="B",],aes(y=fold_red, x=Variant),position=position_nudge(x = 0.05),col="black",size=5,alpha=1.,shape=15)+
  geom_errorbar(data=FC_fig[FC_fig$figure=="B",],aes(x=Variant,ymin=lower,ymax=upper),position=position_nudge(x = 0.05),col="black",width=0,alpha=1.)+
  xlab("Variant")+  ylab(expression("Neutralizing antibody titer ratio ("*NATR[var]*")"))+
  theme(axis.title=element_text(size=17),legend.position = "top",#legend.position = c(.75, .720),
        axis.text=element_text(size=15),legend.text=element_text(size=13),
        legend.title=element_text(size=16))

#Fig S1
ggplot()+
  geom_point(data=m[m$endpoint=="Hospitalization",],aes(x=fold_red_jit,y=VE,color=Variant),size=3)+
  geom_errorbar(data=m[m$endpoint=="Hospitalization",],aes(x=fold_red_jit,ymin=lower,ymax=upper,color=Variant))+
  # geom_line(data=mL_hosp,aes(x=fold_red,y=VE))+
  # geom_ribbon(data=mL_hosp,aes(x=fold_red,ymin=lower,ymax=upper),fill="gray",alpha=0.5)+
  scale_x_continuous(trans="log2",breaks=c(1:8),limits=c(0.7,9))+
  theme_few()+scale_shape_manual(values = c(15:18))+ 
  xlab("Fold reduction in neutralizing antibody titers")+
  ylab("VE against hospitalization")+  labs(col="Variant")+
  theme(axis.title=element_text(size=20),#legend.position = c(.8, .750),#legend.position = "top",
        axis.text=element_text(size=20),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  facet_wrap(.~code)+theme(strip.text.x = element_text(size = 15))

#Fig S2
ggplot()+
  geom_point(data=m[m$endpoint=="Symptomatic disease"&m$vaccine!="Johnson and Johnson"&m$vaccine!="Sinovac"&m$vaccine!="Sputnik V",],aes(x=fold_red_jit,y=VE,color=Variant),size=3)+
  geom_errorbar(data=m[m$endpoint=="Symptomatic disease"&m$vaccine!="Johnson and Johnson"&m$vaccine!="Sinovac"&m$vaccine!="Sputnik V",],aes(x=fold_red_jit,ymin=lower,ymax=upper,color=Variant))+
  # geom_line(data=mL_hosp,aes(x=fold_red,y=VE))+
  # geom_ribbon(data=mL_hosp,aes(x=fold_red,ymin=lower,ymax=upper),fill="gray",alpha=0.5)+
  scale_x_continuous(trans="log2",breaks=c(1:8),limits=c(0.7,9))+
  theme_few()+scale_shape_manual(values = c(15:18))+ 
  xlab("Fold reduction in neutralizing antibody titers")+
  ylab("VE against symptomatic disease")+  labs(col="Variant")+
  theme(axis.title=element_text(size=20),#legend.position = c(.8, .750),#legend.position = "top",
        axis.text=element_text(size=20),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  facet_wrap(.~code)+theme(strip.text.x = element_text(size = 15))


#calculate effective events in each arm
VE_ci=function (par,mrow) {
  #  print(mrow)
  if(mrow$VE>=1){
    E_rat=exp(par[1])
    E1s=0;N1s=1000000;N2s=1000000;E2s=100
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    VE=1
    m1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m1,family=binomial,weights=N,method = "brglmFit");#summary(f2)
    m1$fit=predict(f2,type="link",se.fit=T)$fit
    m1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m1$fit[m1$Group=="Vaccine"],
                       m1$se.fit[m1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m1$fit[m1$Group=="Control"],
                   m1$se.fit[m1$Group=="Control"]))
    VE=data.frame(VE_m=NA,VE_up=NA,VE_low=NA)
    VE$VE_m=1
    VE$VE_up=1-quantile(xnas1,0.025)
    VE$VE_low=1-quantile(xnas1,0.975)
    err=(mrow$upper-VE$VE_up)^2+(mrow$lower-VE$VE_low)^2
    print(paste(c(E_rat,VE,err),sep=" "))
    err
    
    
  }
  else{
    E_rat=exp(par[1])
    E1s=100;N1s=1000000;N2s=1000000;E2s=(E1s/N1s)*N2s/(1-mrow$VE)
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    VE=1-(E1/N1)/(E2/N2)
    m1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m1,family=binomial,weights=N);#summary(f2)
    m1$fit=predict(f2,type="link",se.fit=T)$fit
    m1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m1$fit[m1$Group=="Vaccine"],
                       m1$se.fit[m1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m1$fit[m1$Group=="Control"],
                   m1$se.fit[m1$Group=="Control"]))
    VE=data.frame(VE_m=NA,VE_up=NA,VE_low=NA)
    VE$VE_m=1-plogis(m1$fit[m1$Group=="Vaccine"])/plogis(m1$fit[m1$Group=="Control"])
    VE$VE_up=1-quantile(xnas1,0.025)
    VE$VE_low=1-quantile(xnas1,0.975)
    err=(mrow$upper-VE$VE_up)^2+(mrow$lower-VE$VE_low)^2
    #  print(paste(c(E_rat,VE,err),sep=" "))
    err
  }
}


for (i in 1:nrow(m)) {
  mrow=m[i,c("Study","VE","upper","lower")]
  print(mrow)
  if(mrow$VE>=1){
    o1=optim(par=c(1),fn=VE_ci,control=list(trace=T,reltol=1e-4),method="Nelder-Mead",mrow=mrow);o1
    E1s=0;N1s=1000000;N2s=1000000;E2s=100
    E_rat=exp(o1$par[1])
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    1
    m$I_v[i]=E1;m$I_c[i]=E2;m$N_v[i]=N1;m$N_c[i]=N2
    m1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m1,family=binomial,weights=N,method = "brglmFit");#summary(f2)
    m1$fit=predict(f2,type="link",se.fit=T)$fit
    m1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m1$fit[m1$Group=="Vaccine"],
                       m1$se.fit[m1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m1$fit[m1$Group=="Control"],
                   m1$se.fit[m1$Group=="Control"]))
    m$VE_m[i]=1
    m$VE_up[i]=1-quantile(xnas1,0.025)
    m$VE_low[i]=1-quantile(xnas1,0.975)
  }else{
    o1=optim(par=c(1),fn=VE_ci,control=list(trace=T,reltol=1e-4),method="Nelder-Mead",mrow=mrow);o1
    E1s=100;N1s=1000000;N2s=1000000;E2s=(E1s/N1s)*N2s/(1-mrow$VE)
    E_rat=exp(o1$par[1])
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    1-(E1/N1)/(E2/N2)
    m$I_v[i]=E1;m$I_c[i]=E2;m$N_v[i]=N1;m$N_c[i]=N2
    m1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m1,family=binomial,weights=N);#summary(f2)
    m1$fit=predict(f2,type="link",se.fit=T)$fit
    m1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m1$fit[m1$Group=="Vaccine"],
                       m1$se.fit[m1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m1$fit[m1$Group=="Control"],
                   m1$se.fit[m1$Group=="Control"]))
    m$VE_m[i]=1-plogis(m1$fit[m1$Group=="Vaccine"])/plogis(m1$fit[m1$Group=="Control"])
    m$VE_up[i]=1-quantile(xnas1,0.025)
    m$VE_low[i]=1-quantile(xnas1,0.975)
  }
}

#Cap events to account for very large studies having heavy influence
#Create separate dataframes to compare caps of 10^5, 10^4, 10^3
#Use 10^3 = 1000 in paper
#Cap events at 10^5
cap=100000
for (i in 1:nrow(m)) {
  
  if(m$I_c[i]>cap){
    m$I_v[i]=m$I_v[i]*cap/m$I_c[i]
    m$I_c[i]=cap
  }
}
m_10e5=m

#prevalence
m_10e5$b_s=m_10e5$I_c/m_10e5$N_c

#Cap events at 10^4
cap=10000
for (i in 1:nrow(m)) {
  
  if(m$I_c[i]>cap){
    m$I_v[i]=m$I_v[i]*cap/m$I_c[i]
    m$I_c[i]=cap
  }
}

m_10e4=m

#prevalence
m_10e4$b_s=m_10e4$I_c/m_10e4$N_c


cap=1000
for (i in 1:nrow(m)) {
  
  if(m$I_c[i]>cap){
    m$I_v[i]=m$I_v[i]*cap/m$I_c[i]
    m$I_c[i]=cap
  }
}

#prevalence
m$b_s=m$I_c/m$N_c


lkf=function(c0,c1) {   #This function estimates the neg log likelihood using data in data.frame m1
  # ud=unique(m1$week)#unique time points
  lk=matrix(0,(nrow(mm1)),1)
  llk=matrix(0,(nrow(mm1)),1)
  for (i in 1:nrow(mm1)) {
    I_c=mm1$I_c[i] #total cases control
    N_c=mm1$N_c[i] #total control
    I_v=mm1$I_v[i] #total cases vax
    N_v=mm1$N_v[i] #total vax
    b_s=mm1$b_s[i] #baseline risk (prev control)
    rel_Nabs=mm1$rel_Nabs[i] #relative nabs vax
    #lk[i]=(b_s^I_c)*((1-b_s)^(N_c-I_c))*(b_s*(1-(1/(1+exp(-(c0+c1*log2(rel_Nabs))))))^I_v)*((1-b_s*(1-(1/(1+exp(-(c0+c1*log2(rel_Nabs)))))))^(N_v-I_v))#likelihood of data
    llk[i]=I_c*log(b_s)+(N_c-I_c)*log(1-b_s)+(I_v)*log(b_s*((1/(1+exp(-(c0+c1*log2(rel_Nabs)))))))+(N_v-I_v)*log(1-b_s*((1/(1+exp(-(c0+c1*log2(rel_Nabs)))))))
    
  }
  tkl=-sum(llk) #sum likelihood 
  return(tkl)	}

#number of draws 
nds=10000
qdraws=runif(nds)

#create dataframes for model outputs 
y=matrix(data=NA,nrow =nds, ncol = 21 )
z=data.frame(endpoint=NA,nabs=rep(0.0078125*2^((0:20)/2),4),mean=NA)
CI=data.frame(endpoint=NA,nabs=rep(0.0078125*2^((0:20)/2),4),mean=NA)
VE_l=matrix(data=NA,nrow =nds, ncol = 21 )

y_10e4=matrix(data=NA,nrow =nds, ncol = 21 )
z_10e4=data.frame(endpoint=NA,nabs=rep(0.0078125*2^((0:20)/2),4),mean=NA)
CI_10e4=data.frame(endpoint=NA,nabs=rep(0.0078125*2^((0:20)/2),4),mean=NA)
VE_l_10e4=matrix(data=NA,nrow =nds, ncol = 21 )

y_10e5=matrix(data=NA,nrow =nds, ncol = 21 )
z_10e5=data.frame(endpoint=NA,nabs=rep(0.0078125*2^((0:20)/2),4),mean=NA)
CI_10e5=data.frame(endpoint=NA,nabs=rep(0.0078125*2^((0:20)/2),4),mean=NA)
VE_l_10e5=matrix(data=NA,nrow =nds, ncol = 21 )

base_prev=115/10000 #baseline prev from original UK report Andrews et al.

#baseline nabs
pfizer_nabs=2.3723404
moderna_nabs=4.1332912
#Pfizer waning ratios relative to 1 mo from dose 2 from Falsey et al 2021 Fig 1A WT, age weighted (71%, 29%)
wane_ratio=1/(.71*497/83+.29*538/41)
boost_ratio=.71*1754/497+.29*1318/538
#Moderna boosting @28d relative to post 2nd dose from Chu et al medRxiv 0,181, 209: 1268, 150.2, 1951.7
Mod_boost_ratio=1951.7/1268
Mod_wane_ratio=150.2/1268

#create dataframe for predictions
pred_om=data.frame(Variant=rep(c(rep(c("Omicron","Omicron","Delta"),each=6)),4),vaccine=rep(c(rep(c("Pfizer","Moderna"),each=3)),12),status=rep(c("waned","full","boosted"),24),
                   endpoint=c(rep(c("documented infection","All infections","Hospitalization","Symptomatic disease"),each=18)),
                   Prediction=rep(c(rep(c("Dec. 2021","Updated","Delta"),each=6)),4),
                   baseline_nabs=NA,nabs=NA,x_lower=NA,x_upper=NA, mean=NA,
                   lowerPI=NA,upperPI=NA,lowerCI=NA,upperCI=NA)

#baseline nabs for Pfizer/immune status
pred_om$baseline_nabs[pred_om$vaccine=="Pfizer"&pred_om$status=="waned"]=pfizer_nabs*wane_ratio
pred_om$baseline_nabs[pred_om$vaccine=="Pfizer"&pred_om$status=="full"]=pfizer_nabs
pred_om$baseline_nabs[pred_om$vaccine=="Pfizer"&pred_om$status=="boosted"]=pfizer_nabs*boost_ratio

#baseline nabs for Moderna/immune status
pred_om$baseline_nabs[pred_om$vaccine=="Moderna"&pred_om$status=="waned"]=moderna_nabs*Mod_wane_ratio
pred_om$baseline_nabs[pred_om$vaccine=="Moderna"&pred_om$status=="full"]=moderna_nabs
pred_om$baseline_nabs[pred_om$vaccine=="Moderna"&pred_om$status=="boosted"]=moderna_nabs*Mod_boost_ratio

#baseline nabs for Updated BA.1 prediction including upper/lower CIs
pred_om$nabs[pred_om$Prediction=="Updated"]=pred_om$baseline_nabs[pred_om$Prediction=="Updated"]/FC3$fold_red[FC3$Variant=="BA.1"]
pred_om$x_lower[pred_om$Prediction=="Updated"]=pred_om$baseline_nabs[pred_om$Prediction=="Updated"]/FC3$upper[FC3$Variant=="BA.1"]
pred_om$x_upper[pred_om$Prediction=="Updated"]=pred_om$baseline_nabs[pred_om$Prediction=="Updated"]/FC3$lower[FC3$Variant=="BA.1"]

#baseline nabs for Dec. 2021 BA.1 prediction including upper/lower CIs
pred_om$nabs[pred_om$Prediction=="Dec. 2021"]=pred_om$baseline_nabs[pred_om$Prediction=="Dec. 2021"]/FC3$fold_red[FC3$Variant=="BA.1 Dec. 2021"]
pred_om$x_lower[pred_om$Prediction=="Dec. 2021"]=pred_om$baseline_nabs[pred_om$Prediction=="Dec. 2021"]/FC3$upper[FC3$Variant=="BA.1 Dec. 2021"]
pred_om$x_upper[pred_om$Prediction=="Dec. 2021"]=pred_om$baseline_nabs[pred_om$Prediction=="Dec. 2021"]/FC3$lower[FC3$Variant=="BA.1 Dec. 2021"]

#baseline nabs for Delta prediction including upper/lower CIs
pred_om$nabs[pred_om$Prediction=="Delta"]=pred_om$baseline_nabs[pred_om$Prediction=="Delta"]/FC3$fold_red[FC3$Variant=="Delta"]
pred_om$x_lower[pred_om$Prediction=="Delta"]=pred_om$baseline_nabs[pred_om$Prediction=="Delta"]/FC3$upper[FC3$Variant=="Delta"]
pred_om$x_upper[pred_om$Prediction=="Delta"]=pred_om$baseline_nabs[pred_om$Prediction=="Delta"]/FC3$lower[FC3$Variant=="Delta"]

#Calculate coefficient values for model for each endpoint, fitted line/CIs
ep=unique(m$endpoint[m$endpoint!="Transmission"])

cvalues=data.frame(endpoint=c("documented infection","All infections","Hospitalization","Symptomatic disease"),c0=NA,c1=NA,sigma11=NA,sigma12=NA,sigma21=NA,sigma22=NA)
cvalues_10e4=data.frame(endpoint=c("documented infection","All infections","Hospitalization","Symptomatic disease"),c0=NA,c1=NA,sigma11=NA,sigma12=NA,sigma21=NA,sigma22=NA)
cvalues_10e5=data.frame(endpoint=c("documented infection","All infections","Hospitalization","Symptomatic disease"),c0=NA,c1=NA,sigma11=NA,sigma12=NA,sigma21=NA,sigma22=NA)

for (i in 1:length(ep)) {
  mm1=m[m$endpoint==ep[i],]
  fp2 <- mle2(lkf,start=list(c0=0,c1=0),
              fixed=list(),control=list(trace=3))
  h1=fp2@details$hessian
  mean=c(coef(fp2)[1],coef(fp2)[2])
  sigma=solve(h1)
  x=rmvnorm(n=nds, mean=mean, sigma=sigma)
  
  cvalues$endpoint[i]=ep[i]
  cvalues$c0[i]=coef(fp2)[1]
  cvalues$c1[i]=coef(fp2)[2]
  cvalues$sigma11[i]=sigma[1,1]
  cvalues$sigma12[i]=sigma[1,2]
  cvalues$sigma21[i]=sigma[2,1]
  cvalues$sigma22[i]=sigma[2,2]
  
  if(i==3)
  {hosp_sum=summary(fp2)}  
  if(i==4)
  {symp_sum=summary(fp2)}
  for (j in 0:20) {
    y[,j+1]=1-(1/(1+exp(-(x[,1]+x[,2]*log2(0.0078125*2^(j/2))))))
    
  }
  
  for (k in 1:ncol(y)) {
    CI$endpoint[k+ncol(y)*(i-1)]=ep[i]
    CI$lower[k+ncol(y)*(i-1)]=quantile(y[,k],probs=0.025)
    CI$upper[k+ncol(y)*(i-1)]=quantile(y[,k],probs=0.975)
  }
  for(l in 1:nrow(y)){
    for(n in 1:ncol(y)){
      VE_l[l,n]=1-rbinom(1,10000,prob=base_prev*(1-y[l,n]))/rbinom(1,10000,prob=base_prev)
    }
  }
  for (p in 1:ncol(y)) {
    z$endpoint[p+ncol(y)*(i-1)]=ep[i]
    z$mean[p+ncol(y)*(i-1)]=quantile(VE_l[,p],probs=0.5)
    z$lower[p+ncol(y)*(i-1)]=quantile(VE_l[,p],probs=0.025)
    z$upper[p+ncol(y)*(i-1)]=quantile(VE_l[,p],probs=0.975)
  }
}

for (i in 1:length(ep)) {
  mm1=m_10e4[m_10e4$endpoint==ep[i],]
  fp2 <- mle2(lkf,start=list(c0=0,c1=0),
              fixed=list(),control=list(trace=3))
  h1=fp2@details$hessian
  mean=c(coef(fp2)[1],coef(fp2)[2])
  sigma=solve(h1)
  x=rmvnorm(n=nds, mean=mean, sigma=sigma)
  
  cvalues_10e4$endpoint[i]=ep[i]
  cvalues_10e4$c0[i]=coef(fp2)[1]
  cvalues_10e4$c1[i]=coef(fp2)[2]
  cvalues_10e4$sigma11[i]=sigma[1,1]
  cvalues_10e4$sigma12[i]=sigma[1,2]
  cvalues_10e4$sigma21[i]=sigma[2,1]
  cvalues_10e4$sigma22[i]=sigma[2,2]
  
  for (j in 0:20) {
    y_10e4[,j+1]=1-(1/(1+exp(-(x[,1]+x[,2]*log2(0.0078125*2^(j/2))))))
    
  }
  
  for (k in 1:ncol(y_10e4)) {
    CI_10e4$endpoint[k+ncol(y_10e4)*(i-1)]=ep[i]
    CI_10e4$lower[k+ncol(y_10e4)*(i-1)]=quantile(y_10e4[,k],probs=0.025)
    CI_10e4$upper[k+ncol(y_10e4)*(i-1)]=quantile(y_10e4[,k],probs=0.975)
  }
  for(l in 1:nrow(y_10e4)){
    for(n in 1:ncol(y_10e4)){
      VE_l_10e4[l,n]=1-rbinom(1,10000,prob=base_prev*(1-y_10e4[l,n]))/rbinom(1,10000,prob=base_prev)
    }
  }
  for (p in 1:ncol(y_10e4)) {
    z_10e4$endpoint[p+ncol(y_10e4)*(i-1)]=ep[i]
    z_10e4$mean[p+ncol(y_10e4)*(i-1)]=quantile(VE_l_10e4[,p],probs=0.5)
    z_10e4$lower[p+ncol(y_10e4)*(i-1)]=quantile(VE_l_10e4[,p],probs=0.025)
    z_10e4$upper[p+ncol(y_10e4)*(i-1)]=quantile(VE_l_10e4[,p],probs=0.975)
  }
}

for (i in 1:length(ep)) {
  mm1=m_10e5[m_10e5$endpoint==ep[i],]
  fp2 <- mle2(lkf,start=list(c0=0,c1=0),
              fixed=list(),control=list(trace=3))
  h1=fp2@details$hessian
  mean=c(coef(fp2)[1],coef(fp2)[2])
  sigma=solve(h1)
  x=rmvnorm(n=nds, mean=mean, sigma=sigma)
  
  cvalues_10e5$endpoint[i]=ep[i]
  cvalues_10e5$c0[i]=coef(fp2)[1]
  cvalues_10e5$c1[i]=coef(fp2)[2]
  cvalues_10e5$sigma11[i]=sigma[1,1]
  cvalues_10e5$sigma12[i]=sigma[1,2]
  cvalues_10e5$sigma21[i]=sigma[2,1]
  cvalues_10e5$sigma22[i]=sigma[2,2]
  
  for (j in 0:20) {
    y_10e5[,j+1]=1-(1/(1+exp(-(x[,1]+x[,2]*log2(0.0078125*2^(j/2))))))
    
  }
  
  for (k in 1:ncol(y_10e5)) {
    CI_10e5$endpoint[k+ncol(y_10e5)*(i-1)]=ep[i]
    CI_10e5$lower[k+ncol(y_10e5)*(i-1)]=quantile(y_10e5[,k],probs=0.025)
    CI_10e5$upper[k+ncol(y_10e5)*(i-1)]=quantile(y_10e5[,k],probs=0.975)
  }
  for(l in 1:nrow(y_10e5)){
    for(n in 1:ncol(y_10e5)){
      VE_l_10e5[l,n]=1-rbinom(1,10000,prob=base_prev*(1-y_10e5[l,n]))/rbinom(1,10000,prob=base_prev)
    }
  }
  for (p in 1:ncol(y_10e5)) {
    z_10e5$endpoint[p+ncol(y_10e5)*(i-1)]=ep[i]
    z_10e5$mean[p+ncol(y_10e5)*(i-1)]=quantile(VE_l_10e5[,p],probs=0.5)
    z_10e5$lower[p+ncol(y_10e5)*(i-1)]=quantile(VE_l_10e5[,p],probs=0.025)
    z_10e5$upper[p+ncol(y_10e5)*(i-1)]=quantile(VE_l_10e5[,p],probs=0.975)
  }
}

#Create dataframe for prediction outputs from model
y1=matrix(data=NA,nrow =nds, ncol = 18 )
VE_l1=matrix(data=NA,nrow =nds, ncol = 18 )

#Predictions
for (i in 1:length(ep)) {
  mm1=pred_om[pred_om$endpoint==ep[i],]
  mean=c(cvalues$c0[cvalues$endpoint==ep[i]],cvalues$c1[cvalues$endpoint==ep[i]])
  sigma[1,1]=cvalues[i,4];sigma[1,2]=cvalues[i,5];sigma[2,1]=cvalues[i,6];sigma[2,2]=cvalues[i,7]
  x=rmvnorm(n=nds, mean=mean, sigma=sigma)
  
  for (j in 1:nrow(mm1)) {
    y1[,j]=1-(1/(1+exp(-(x[,1]+x[,2]*log2(mm1$nabs[j])))))
    # yl[,j+1]=1-(1/(1+exp(-(x[,1]+x[,2]*log2(mm1$xlower[j])))))
    # yu[,j+1]=1-(1/(1+exp(-(x[,1]+x[,2]*log2(mm1$xupper[j])))))
    
  }
  
  for(l in 1:nrow(y1)){
    for(n in 1:ncol(y1)){
      VE_l1[l,n]=1-rbinom(1,10000,prob=base_prev*(1-y1[l,n]))/rbinom(1,10000,prob=base_prev)
    }
  }
  for (k in 1:nrow(mm1)) {
    # z$endpoint[k+19*(i-1)]=ep[i]
    pred_om$mean[k+nrow(mm1)*(i-1)]=quantile(y1[,k],probs=0.5)
    pred_om$lowerPI[k+nrow(mm1)*(i-1)]=quantile(VE_l1[,k],probs=0.025)
    pred_om$upperPI[k+nrow(mm1)*(i-1)]=quantile(VE_l1[,k],probs=0.975)
    pred_om$lowerCI[k+nrow(mm1)*(i-1)]=quantile(y1[,k],probs=0.025)
    pred_om$upperCI[k+nrow(mm1)*(i-1)]=quantile(y1[,k],probs=0.975)
  }
}

# #Relative risk - will need later
pred_om$RR=1-pred_om$mean

#Label for fig 2
tlab4=data.frame(endpoint=c("Hospitalization","Symptomatic disease")
                 ,text1=c("B: Hospitalization","A: Symptomatic disease")
                 ,nabs=c(0.0105,0.0105),mean=c(0.95,0.938))
tlab4$endpoint=factor(tlab4$endpoint,levels = c("Symptomatic disease","Hospitalization"))

#Split data by type for ggplot2
z$data="mean"
CI$data="CI"
m$data="data"


#Modify prediction data for plot
pred_om_fig2=pred_om[pred_om$Prediction=="Updated",]
pred_om_fig2$Variant="Omicron (predicted)"
pred_om_fig2$data="pred"
pred_om_fig2=merge(pred_om_fig2,vac_comp,by="vaccine")

#Data into single datafrae for plot
fig2=bind_rows(z,CI,m,pred_om_fig2)
fig2=fig2[fig2$endpoint=="Symptomatic disease"|fig2$endpoint=="Hospitalization",]
fig2$endpoint=factor(fig2$endpoint,levels = c("Symptomatic disease","Hospitalization"))

#modify lower bound of study to fit on plot
fig2$lower[fig2$Study=="Madhi"]=-0.2

#Fig 2
ggplot(data=fig2)+
  geom_line(data=fig2[fig2$data=="mean",],aes(x=nabs,y=mean),color="black")+
  geom_ribbon(data=fig2[fig2$data=="mean",],aes(x=nabs,ymin=lower,ymax=upper),fill="gray",alpha=0.3)+
  geom_ribbon(data=fig2[fig2$data=="CI",],aes(x=nabs,ymin=lower,ymax=upper),fill="gray",alpha=0.5)+
  geom_point(data=fig2[fig2$data=="data",],aes(x=rel_nabs_jit,y=VE,color=Variant,shape=code),size=3)+
  geom_errorbar(data=fig2[fig2$data=="data",],aes(x=rel_nabs_jit,ymin=lower,ymax=upper,color=Variant),width=0)+
  #geom_point(data=fig2[fig2$data=="cap",],aes(x = rel_nabs_jit,y = VE),shape=8,color="black",size=1)+
  geom_point(data=fig2[fig2$data=="pred",],aes(x=nabs,y=mean,color=Variant,shape=code),size=3)+
  geom_errorbar(data=fig2[fig2$data=="pred",],aes(x=nabs,ymin=lower,ymax=upper,color=Variant),width=0)+
  geom_errorbar(data=fig2[fig2$data=="pred",],aes(y=mean,xmin=x_lower,xmax=x_upper,color=Variant),width=0)+
  # geom_segment(data = arrow_df,aes(x = rel_nabs_jit,
  #                                 y = VE-0.1,
  #                                 xend = rel_nabs_jit,
  #                                 yend = VE),
  #              arrow = arrow(angle=30,length = unit(0.15, "cm")))+
  facet_wrap(.~endpoint,nrow=2, scales="free_y")+
  scale_x_continuous(expand = c(0, 0),trans="log2",breaks=c(0.0078125*2^c(1:10)),labels = label_number(accuracy = 0.001, scale = 1, 
                                                                                                       prefix = "", suffix = "",
                                                                                                       big.mark = " ", decimal.mark = "."))+
  scale_y_continuous(expand = c(0, 0))+
  scale_color_manual(breaks=c("D614G","Alpha","Beta","Delta","Gamma","Omicron (predicted)"),values=c( "#3e049c","#8a09a5","#c8437b","#e97257","#fba238","black"))+
  coord_cartesian(ylim = c(NA, 1),xlim = c(0.01,5)) +
  theme_few()+scale_shape_manual(values = c(15,18,17,3,19,4,9))+
  xlab(expression("Neutralizing antibody titer ratio ("*NATR[tot]*")"))+
  ylab("Vaccine effectiveness (VE)")+  labs(col="Variant", shape="Vaccine")+
  theme(axis.title=element_text(size=20),#legend.position = c(.8, .750),#legend.position = "top",
        axis.text=element_text(size=18),legend.text=element_text(size=20),
        legend.title=element_text(size=20))+theme(strip.text.x = element_blank())+
  geom_text(data=tlab4,aes(x=nabs,y=mean,label=text1),size=6,hjust=0)

#Compile dataframe for figure S4
z_10e3=z
z_10e3$data="Cap = 1,000"
z_10e4$data="Cap = 10,000"
z_10e5$data="Cap = 100,000"

cap_fig=bind_rows(z_10e3,z_10e4,z_10e5)
cap_fig=cap_fig[cap_fig$endpoint=="Symptomatic disease"|cap_fig$endpoint=="Hospitalization",]
cap_fig$endpoint=factor(cap_fig$endpoint,levels = c("Symptomatic disease","Hospitalization"))

#Figure S4 - comparing
ggplot(data=cap_fig)+
  geom_line(aes(x=nabs,y=mean,color=data))+
  facet_wrap(.~endpoint,nrow=2, scales="free_y")+
  geom_point(data=fig2[fig2$data=="data",],aes(x=rel_nabs_jit,y=VE),size=3)+
  geom_errorbar(data=fig2[fig2$data=="data",],aes(x=rel_nabs_jit,ymin=lower,ymax=upper),width=0)+
  scale_x_continuous(expand = c(0, 0),trans="log2",breaks=c(0.0078125*2^c(1:10)),labels = label_number(accuracy = 0.001, scale = 1, 
                                                                                                       prefix = "", suffix = "",
                                                                                                       big.mark = " ", decimal.mark = "."))+
  scale_y_continuous(expand = c(0, 0))+
  #scale_color_manual(breaks=c("D614G","Alpha","Beta","Delta","Gamma","Omicron (predicted)"),values=c( "#3e049c","#8a09a5","#c8437b","#e97257","#fba238","black"))+
  coord_cartesian(ylim = c(NA, 1),xlim = c(0.01,5)) +
  theme_few()+scale_shape_manual(values = c(15,18,17,3,19,4,9))+
  xlab(expression("Neutralizing antibody titer ratio ("*NATR[tot]*")"))+
  ylab("Vaccine effectiveness (VE)")+  labs(col="Data Cap")+
  theme(axis.title=element_text(size=20),#legend.position = c(.8, .750),#legend.position = "top",
        axis.text=element_text(size=18),legend.text=element_text(size=20),
        legend.title=element_text(size=20))+theme(strip.text.x = element_blank())+
  geom_text(data=tlab4,aes(x=nabs,y=mean,label=text1),size=6,hjust=0)

#load/clean validation data and merge with NATR
vdf_all=read.csv("validation.csv") 
vdf=vdf_all[vdf_all$first.author!="Altarawneh"&vdf_all$endpoint!="documented infection",]
vdf$fold_red=FC3$fold_red[FC3$Variant=="BA.1"]
vdf$fold_red_upper=FC3$upper[FC3$Variant=="BA.1"]
vdf$fold_red_lower=FC3$lower[FC3$Variant=="BA.1"]

#modify nabs by timing of dose and NATRs
vdf$nabs=(vdf$timing_x*vdf$vax_nabs/vdf$fold_red)
vdf$x_lower=(vdf$timing_x*vdf$vax_nabs/vdf$fold_red_upper)
vdf$x_upper=(vdf$timing_x*vdf$vax_nabs/vdf$fold_red_lower)

#Calculate effective events
vdf=vdf[order(vdf$endpoint),]
ep2=unique(vdf$endpoint)

VE_ci_vdf=function (par,m_vdfrow) {
  #  print(mrow)
  if(m_vdfrow$VE>=1){
    E_rat=exp(par[1])
    E1s=0;N1s=10000000;N2s=10000000;E2s=100
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    VE=1
    m_vdf1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m_vdf1,family=binomial,weights=N,method = "brglmFit");#summary(f2)
    m_vdf1$fit=predict(f2,type="link",se.fit=T)$fit
    m_vdf1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m_vdf1$fit[m_vdf1$Group=="Vaccine"],
                       m_vdf1$se.fit[m_vdf1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m_vdf1$fit[m_vdf1$Group=="Control"],
                   m_vdf1$se.fit[m_vdf1$Group=="Control"]))
    VE=data.frame(VE_m=NA,VE_up=NA,VE_low=NA)
    VE$VE_m=1
    VE$VE_up=1-quantile(xnas1,0.025)
    VE$VE_low=1-quantile(xnas1,0.975)
    err=(m_vdfrow$upper-VE$VE_up)^2+(m_vdfrow$lower-VE$VE_low)^2
    print(paste(c(E_rat,VE,err),sep=" "))
    err
    
  }
  else{
    E_rat=exp(par[1])
    E1s=100;N1s=10000000;N2s=10000000;E2s=(E1s/N1s)*N2s/(1-m_vdfrow$VE)
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    VE=1-(E1/N1)/(E2/N2)
    m_vdf1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m_vdf1,family=binomial,weights=N);#summary(f2)
    m_vdf1$fit=predict(f2,type="link",se.fit=T)$fit
    m_vdf1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m_vdf1$fit[m_vdf1$Group=="Vaccine"],
                       m_vdf1$se.fit[m_vdf1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m_vdf1$fit[m_vdf1$Group=="Control"],
                   m_vdf1$se.fit[m_vdf1$Group=="Control"]))
    VE=data.frame(VE_m=NA,VE_up=NA,VE_low=NA)
    VE$VE_m=1-plogis(m_vdf1$fit[m_vdf1$Group=="Vaccine"])/plogis(m_vdf1$fit[m_vdf1$Group=="Control"])
    VE$VE_up=1-quantile(xnas1,0.025)
    VE$VE_low=1-quantile(xnas1,0.975)
    err=(m_vdfrow$upper-VE$VE_up)^2+(m_vdfrow$lower-VE$VE_low)^2
    #  print(paste(c(E_rat,VE,err),sep=" "))
    err
  }
}

m_vdf=vdf
m_vdf$VE=m_vdf$mean


for (i in 1:nrow(m_vdf)) {
  m_vdfrow=m_vdf[i,c("VE","upper","lower")]
  print(m_vdfrow)
  if(m_vdfrow$VE>=1){
    o1=optim(par=c(1),fn=VE_ci_vdf,control=list(trace=T,reltol=1e-4),method="Nelder-Mead",m_vdfrow=m_vdfrow);o1
    E1s=0;N1s=1000000;N2s=1000000;E2s=100
    E_rat=exp(o1$par[1])
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    1
    m_vdf$I_v[i]=E1;m_vdf$I_c[i]=E2;m_vdf$N_v[i]=N1;m_vdf$N_c[i]=N2
    m_vdf1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m_vdf1,family=binomial,weights=N,method = "brglmFit");#summary(f2)
    m_vdf1$fit=predict(f2,type="link",se.fit=T)$fit
    m_vdf1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m_vdf1$fit[m_vdf1$Group=="Vaccine"],
                       m_vdf1$se.fit[m_vdf1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m_vdf1$fit[m_vdf1$Group=="Control"],
                   m_vdf1$se.fit[m_vdf1$Group=="Control"]))
    m_vdf$VE_m[i]=1
    m_vdf$VE_up[i]=1-quantile(xnas1,0.025)
    m_vdf$VE_low[i]=1-quantile(xnas1,0.975)
  }else{
    o1=optim(par=c(1),fn=VE_ci_vdf,control=list(trace=T,reltol=1e-4),method="Nelder-Mead",m_vdfrow=m_vdfrow);o1
    E1s=100;N1s=1000000;N2s=1000000;E2s=(E1s/N1s)*N2s/(1-m_vdfrow$VE)
    E_rat=exp(o1$par[1])
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    1-(E1/N1)/(E2/N2)
    m_vdf$I_v[i]=E1;m_vdf$I_c[i]=E2;m_vdf$N_v[i]=N1;m_vdf$N_c[i]=N2
    m_vdf1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m_vdf1,family=binomial,weights=N);#summary(f2)
    m_vdf1$fit=predict(f2,type="link",se.fit=T)$fit
    m_vdf1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m_vdf1$fit[m_vdf1$Group=="Vaccine"],
                       m_vdf1$se.fit[m_vdf1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m_vdf1$fit[m_vdf1$Group=="Control"],
                   m_vdf1$se.fit[m_vdf1$Group=="Control"]))
    m_vdf$VE_m[i]=1-plogis(m_vdf1$fit[m_vdf1$Group=="Vaccine"])/plogis(m_vdf1$fit[m_vdf1$Group=="Control"])
    m_vdf$VE_up[i]=1-quantile(xnas1,0.025)
    m_vdf$VE_low[i]=1-quantile(xnas1,0.975)
  }
}

#Merge with vaccine codes and modify dataframe
vdf_fig4=merge(m_vdf,vac_comp,by="vaccine")
vdf_fig4$Variant="Omicron (validation)"
vdf_fig4$data="validation"
vdf_fig4$rel_nabs_jit=vdf_fig4$nabs+rnorm(nrow(vdf_fig4),0,0.01)
vdf_fig4=vdf_fig4[vdf_fig4$first.author!="Andrews (preprint)"&vdf_fig4$Dose!="full",]


vdf$Data="Validation Data"
pred_om$Dose=pred_om$status
pred_om$Data[pred_om$Prediction=="Dec. 2021"]="Prediction Dec. 2021"
pred_om$Data[pred_om$Prediction=="Updated"]="Updated Prediction"

#Prepare to calculate variance weighted means for validation data
vdf$logitmean=qlogis(vdf$mean)
vdf$logitSE=(qlogis(vdf$upper)-qlogis(vdf$mean))/1.96
vdf$precision=1/(vdf$logitSE^2)
vdf$logitVar=vdf$logitSE^2

vdf_select=vdf%>%
  unite(c("study.ID","endpoint"),col="Study_End",remove=F)

vdf_GM=vdf_select[vdf_select$first.author!="Andrews (preprint)"&(vdf_select$Study_End!="130_Symptomatic disease"),]
vdf_GM$weight=vdf_GM$precision*vdf_GM$mean

#Subset dataframe
grand_test1=vdf_GM%>%
  unite(c("Dose","vaccine","endpoint"),col="Unique")%>%
  select(c("Unique","mean","precision","lower","upper","logitVar","logitSE","weight"))

#Continue preparing for grand mean calcultion
x0=as.data.frame(table(grand_test1$Unique));
x0=x0%>%
  rename(Unique=Var1)
grand_test2=merge(grand_test1,x0,all.x=T)
grand_test=grand_test2[grand_test2$Freq>=2,]
grand_test_hold=grand_test2[grand_test2$Freq==1,]

#Calculate frand means for endpoint/vaccine/immune status combination with multiple validation points
GM_test=data.frame(aggregate(precision~Unique,data = grand_test,FUN=sum))
GM_test[,3]=aggregate(weight~Unique,data = grand_test,FUN=sum)[2]
GM_test$mean=GM_test$weight/GM_test$precision
GM_test[,5]=aggregate(logitVar~Unique,data = grand_test,FUN=sum)[2]
GM_test$lowerVar=qlogis(GM_test$mean)-1.96*sqrt(GM_test$logitVar)
GM_test$upperVar=qlogis(GM_test$mean)+1.96*sqrt(GM_test$logitVar)
GM_test$lower=plogis(GM_test$lowerVar)
GM_test$upper=plogis(GM_test$upperVar)

#Subset grand mean data and put into single dataframe
GM_test_hold=grand_test_hold%>%
  select(c("Unique","mean","lower","upper"))%>%
  separate(col="Unique",sep = "_",remove=T,into = c("Dose","vaccine","endpoint"))

GM_test1=GM_test%>%
  select(c("Unique","mean","lower","upper"))%>%
  separate(col="Unique",sep = "_",remove=T,into = c("Dose","vaccine","endpoint"))

GM_fig3=rbind(GM_test1,GM_test_hold)

GM_fig3$Data="Grand Mean"

#Calculate absolute error, relative error, fits within PIs, etc. 
MAEpred=pred_om[pred_om$Data=="Updated Prediction"&pred_om$Variant=="Omicron",]%>%
  unite(c("Dose","vaccine","endpoint"),col="Unique")%>%
  select(c("Unique","mean","lowerPI","upperPI"))

MAEvdf=vdf[vdf$first.author!="Andrews (preprint)",]
MAEvdf$validmean=MAEvdf$mean
MAEvdf=MAEvdf%>%
  unite(c("Dose","vaccine","endpoint"),col="Unique")%>%
  select(c("Unique","validmean","lower","upper"))


GM_MAE=GM_fig3
GM_MAE$grandmean=GM_MAE$mean
GM_MAE$grandlower=GM_MAE$lower
GM_MAE$grandupper=GM_MAE$upper

GM_MAE=GM_MAE%>%  
  unite(c("Dose","vaccine","endpoint"),col="Unique")%>%
  select(c("Unique","grandmean","grandlower","grandupper"))

MAEpred2=merge(MAEpred,GM_MAE,all.x = T)
MAEpred3=na.omit(MAEpred2)
MAEpred3$AE=abs(MAEpred3$mean-MAEpred3$grandmean)
MAEpred3$Err=(MAEpred3$mean-MAEpred3$grandmean)
MAEpred3$lowerErr=(MAEpred3$grandmean-MAEpred3$lowerPI)
MAEpred3$upperErr=(MAEpred3$upperPI-MAEpred3$grandmean)
MAEpred3$RR=1-MAEpred3$mean
MAEpred3$relE=MAEpred3$AE/MAEpred3$RR

MAEpred3=MAEpred3%>%
  separate(col="Unique",into = c("Dose","vaccine","endpoint"),sep="_",remove=F)

MAEpred3a=MAEpred3[MAEpred3$Dose!="full",]

MAEpredall2=merge(MAEpred,MAEvdf,all.x = T)

MAEpredall3=na.omit(MAEpredall2)
MAEpredall3$AE=abs(MAEpredall3$mean-MAEpredall3$validmean)
MAEpredall3$Err=(MAEpredall3$mean-MAEpredall3$validmean)


#Means and vars of absolute errors
mean(MAEpred3a$AE[MAEpred3a$endpoint=="Hospitalization"]) #Hosp
sqrt(var(x=MAEpred3a$AE[MAEpred3a$endpoint=="Hospitalization"],y=NULL,na.rm=FALSE))

mean(MAEpred3a$AE[MAEpred3a$endpoint=="Symptomatic disease"]) #Symp
sqrt(var(x=MAEpred3a$AE[MAEpred3a$endpoint=="Symptomatic disease"],y=NULL,na.rm=FALSE))

mean(MAEpred3a$AE[MAEpred3a$vaccine=="Pfizer"]) #Pfizer
sqrt(var(x=MAEpred3a$AE[MAEpred3a$vaccine=="Pfizer"],y=NULL,na.rm=FALSE))

mean(MAEpred3a$AE[MAEpred3a$vaccine=="Moderna"]) #Moderna
sqrt(var(x=MAEpred3a$AE[MAEpred3a$vaccine=="Moderna"],y=NULL,na.rm=FALSE))


#Final dataframe for fig 3
GM_fig3=merge(GM_fig3,vac_comp,by="vaccine")
fig3=bind_rows(pred_om[pred_om$endpoint!="All infections"&pred_om$Variant=="Omicron",],vdf_select[vdf_select$endpoint!="documented infection"&vdf_select$first.author!="Andrews (preprint)"&(vdf_select$Study_End!="130_Symptomatic disease"),])
fig3=merge(fig3,vac_comp,by="vaccine")
fig3a=bind_rows(fig3,GM_fig3)

#Manually create x-values 
fig3a$x1=NA;fig3a$x2=NA;fig3a$x3=NA;fig3a$x4=NA

fig3a$x1[fig3a$Dose=="boosted"]=4.8
fig3a$x1[fig3a$Dose=="waned"]=3
fig3a$x2[fig3a$Data=="Prediction Dec. 2021"]=-0.15
fig3a$x2[fig3a$Data=="Updated Prediction"]=-0.05
fig3a$x2[fig3a$Data=="Validation Data"]=0.05
fig3a$x2[fig3a$Data=="Grand Mean"]=0.15
fig3a$x3[fig3a$vaccine=="Pfizer"]=-0.35
fig3a$x3[fig3a$vaccine=="Moderna"]=0.35

fig3a$x4=fig3a$x1+fig3a$x2+fig3a$x3


#Jitter values for plots
fig3a$x4=fig3a$x4+rnorm(nrow(fig3a),0,0.05) #jitter xvals

#Clean/modify data for final figure
fig3a2=fig3a[fig3a$endpoint!="documented infection"&fig3a$Dose!="full",]
fig3a2$Data[fig3a2$Data=="Grand Mean"]="Weighted Mean"
fig3a2$Data[fig3a2$Data=="Prediction Dec. 2021"]="Prediction Dec. 11 2021"

fig3a2$endpoint=factor(fig3a2$endpoint,levels = c("Symptomatic disease","Hospitalization"))

#Figure 3
ggplot(data=fig3a2)+
  geom_point(data=fig3a2[fig3a2$Data=="Validation Data",],aes(x=x4,y=mean,color=Data,shape=code),size=2.5,alpha=1)+
  geom_errorbar(data=fig3a2[fig3a2$Data=="Validation Data",],aes(x=x4,ymin=lower,ymax=upper,color=Data),linewidth=1.5,width=0,alpha=0.5)+
  geom_errorbar(data=fig3a2[fig3a2$Data=="Prediction Dec. 11 2021",],aes(x=x4,ymin=lowerCI,ymax=upperCI,color=Data),linewidth=1.5,width=0,alpha=0.5)+
  geom_point(data=fig3a2[fig3a2$Data=="Prediction Dec. 11 2021",],aes(x=x4,y=mean,color=Data,shape=code),size=2.5,alpha=0.5)+
  geom_errorbar(data=fig3a2[fig3a2$Data=="Prediction Dec. 11 2021",],aes(x=x4,ymin=lowerPI,ymax=upperPI,color=Data),width=0,alpha=0.5)+
  geom_errorbar(data=fig3a2[fig3a2$Data=="Updated Prediction",],aes(x=x4,ymin=lowerCI,ymax=upperCI,color=Data),linewidth=1.5,width=0,alpha=0.5)+
  geom_point(data=fig3a2[fig3a2$Data=="Updated Prediction",],aes(x=x4,y=mean,color=Data,shape=code),size=2.5,alpha=0.5)+
  geom_errorbar(data=fig3a2[fig3a2$Data=="Updated Prediction",],aes(x=x4,ymin=lowerPI,ymax=upperPI,color=Data),width=0,alpha=0.5)+
  geom_point(data=fig3a2[fig3a2$Data=="Weighted Mean",],aes(x=x4,y=mean,color=Data,shape=code),size=2.5,alpha=1)+
  geom_errorbar(data=fig3a2[fig3a2$Data=="Weighted Mean",],aes(x=x4,ymin=lower,ymax=upper,color=Data),linewidth=1.5,width=0,alpha=0.5)+
  scale_x_continuous(limits=c(2.4,5.4),labels=c("Two-dose waned","Three-dose boosted"),breaks=c(3,4.8))+ 
  scale_color_manual(values=c( "Blue", "Red","Dark Gray","Black"))+
  theme_few()+
  # theme(strip.text.x = element_blank(),axis.title=element_text(size=30),#legend.position = c(.8, .750),legend.position = "top",
  #                                                     axis.text=element_text(size=30),legend.text=element_text(size=30),
  #                                                     legend.title=element_text(size=30))+
  theme(text=element_text(size=25),axis.title=element_text(size=23),#legend.position = c(.8, .750),#legend.position = "top",
        axis.text=element_text(size=20),legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  ylab("Vaccine effectiveness (VE)")+ xlab("Immune status")+ labs(col="Prediction/Data", shape="Vaccine")+
  
  facet_wrap(.~endpoint,nrow=2,scales="free_y")



#Figure S3
pred_om_now=pred_om[pred_om$Prediction=="Updated",]
pred_om_delta=pred_om[pred_om$Prediction=="Delta",]

compare_delt_om=data.frame(vaccine=pred_om_now$vaccine,status=pred_om_now$status,
                           endpoint=pred_om_now$endpoint,
                           RR=pred_om_now$RR/pred_om_delta$RR)

pred_om_compare=rbind(pred_om_delta,pred_om_now)


pred_om_compare$x1=NA;pred_om_compare$x2=NA;pred_om_compare$x3=NA;pred_om_compare$x4=NA

pred_om_compare$x1[pred_om_compare$Dose=="boosted"]=2.8
pred_om_compare$x1[pred_om_compare$Dose=="full"]=2
pred_om_compare$x1[pred_om_compare$Dose=="waned"]=1.2
pred_om_compare$x2[pred_om_compare$Prediction=="Delta"]=0.1
pred_om_compare$x2[pred_om_compare$Prediction=="Updated"]=-0.1
pred_om_compare$x3[pred_om_compare$vaccine=="Pfizer"]=-0.05
pred_om_compare$x3[pred_om_compare$vaccine=="Moderna"]=0

pred_om_compare$x4=pred_om_compare$x1+pred_om_compare$x2+pred_om_compare$x3

pred_om_compare$endpoint[pred_om_compare$endpoint=="documented infection"]="Documented infection"

pred_om_compare$endpoint=factor(pred_om_compare$endpoint,levels=c("All infections","Documented infection","Symptomatic disease","Hospitalization"))

pred_om_compare$endpoint=factor(pred_om_compare$endpoint,levels=c("Symptomatic disease","Hospitalization","Documented infection","All infections"))

pred_om_compare=merge(pred_om_compare,vac_comp,by="vaccine")

ggplot(data=pred_om_compare[pred_om_compare$endpoint!="All infections"&pred_om_compare$endpoint!="Documented infection"&pred_om_compare$status!="full",])+
  geom_point(aes(x=x4,y=mean,color=Variant,shape=code),size=3)+
  geom_errorbar(aes(x=x4,ymin=lowerCI,ymax=upperCI,color=Variant),width=0)+theme_few()+#theme(strip.text.x = element_blank())+
  scale_x_continuous(limits=c(0.5,3.5),labels=c("Two-dose waned","Three-dose boosted"),breaks=c(1.2,2.8))+
  ylab("Vaccine effectiveness")+ xlab("Immune status")+ labs(col="Variant", shape="Vaccine")+
  scale_color_manual(values=c("#00BFC4","#F564E3"))+
  #theme(axis.title.x = element_blank())+
  theme(text=element_text(size=25),axis.title=element_text(size=20),#legend.position = c(.8, .750),#legend.position = "top",
        axis.text=element_text(size=15),legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  #geom_text(data=tlab2,aes(x=Variant,y=VE,label=text1),size=5,hjust=0)+
  facet_wrap(.~endpoint,nrow=1,scales="free_y")+
  coord_cartesian(ylim = c(NA, 1))

#Figure S5
m_supp=m[m$endpoint=="Hospitalization"|m$endpoint=="Symptomatic disease",]
m_supp=m_supp[order(m_supp$VE, decreasing = TRUE), ]
m_supp$est_ID=c(LETTERS,paste("A",LETTERS,sep=""),paste("B",LETTERS[1:(length(m_supp$Variant)-52)],sep=""))
m_supp$dummy_x1=seq(0.8,63.8,by=1)
m_supp$dummy_x2=seq(1.2,64.2,by=1)


labels_m_supp=c(LETTERS,paste("A",LETTERS,sep=""),paste("B",LETTERS[1:(length(m_supp$Variant)-52)],sep=""))

ggplot(data=m_supp)+
  geom_point(aes(x=dummy_x1,y=VE),shape=19)+
  geom_errorbar(aes(x=dummy_x1,ymax=upper,ymin=VE),color="red",width=0)+
  geom_errorbar(aes(x=dummy_x1,ymin=lower,ymax=VE),color="blue",width=0)+
  geom_point(aes(x=dummy_x2,y=VE_m),shape=17)+
  geom_errorbar(aes(x=dummy_x2,ymax=VE_up,ymin=VE_m),color="red",width=0)+
  geom_errorbar(aes(x=dummy_x2,ymin=VE_low,ymax=VE_m),color="blue",width=0)+
  scale_x_continuous(expand = c(0, 0),limits=c(0.6,64.4),labels=m_supp$est_ID,breaks=seq(1,64,by=1))+
  ylab("Vaccine effectiveness")+ xlab("Study")+ 
  theme_few() + 
  theme(text=element_text(size=25),axis.title=element_text(size=20),#legend.position = c(.8, .750),#legend.position = "top",
        axis.text=element_text(size=10),legend.text=element_text(size=15),
        legend.title=element_text(size=15))

