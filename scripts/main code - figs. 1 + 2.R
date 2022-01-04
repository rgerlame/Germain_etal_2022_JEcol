
{
  
require(dplyr);require(ggplot2);require(lmerTest);require(visreg);require(pscl);require(car):require(doBy);require(scales)
require(grid); require(cowplot); library(effects); require(glmmTMB); require(ggExtra); require(ggeffects)

}

###############
#data prep
###############
  
#site 6
df1<-read.csv("data/site 6 - seed data.csv")
df2<-read.csv("data/site 6 - neighbourhood data.csv")
df2<-df2[,-2]

df<-merge(df1,df2,by="id")
df$seed.number<-as.numeric(as.character(df$seed.number))
df$any.seed<-ifelse(df$seed.number>0,1,0)
df$hard.seed<-df$seed.number*(1-df$prop.mushy)

df.6save<-df[,c(1,2,3,7,8,9)]

df$density<-rowSums(df[,10:44],na.rm=TRUE)

df$source.pop<-df$population
df$population<-as.character(df$population)

df$population<-ifelse(df$population=="grey","sym",df$population)
df$population<-ifelse(df$population=="red","allo1",df$population)
df$population<-ifelse(df$population=="yellow","allo2",df$population)
df$population<-ifelse(df$population=="green","allo3",df$population)

df$focal<-factor(df$pop,levels=c("sym","allo1","allo2","allo3"))
df$competition<-as.factor(ifelse(df$type=="alpha","yes","no"))

#df<-subset(df,germ=="y"&competition=="yes")
df$history<-ifelse(df$focal=="sym","local","foreign")

df6.full<-df
df6<-df6.full
colnames(df6)[3]<-c("seeds")
df6$seeds<-ceiling(df6$seeds)
df6$density.inter<-rowSums(df[,c(10:11,13:44)],na.rm=TRUE)

df6<-df6[,c("seeds","focal","germ","brohor","vulmic","clagra","plaere","hemcon","lolmul","lotwra","density","competition")]
df6<-na.omit(df6)
df6$history<-ifelse(df6$focal=="sym","local","foreign")


#site 18
df1<-read.csv("data/site 18 - seed data.csv")
df2<-read.csv("data/site 18 - neighbourhood data.csv")
df2<-df2[,-2]

df<-merge(df1,df2,by="id")
df$seed.number<-as.numeric(as.character(df$seed.number))
df$any.seed<-ifelse(df$seed.number>0,1,0)
df$hard.seed<-df$seed.number*(1-df$prop.mushy)

df.18save<-df[,c(1,2,3,8,9,10)]

df$density<-rowSums(df[,11:59],na.rm=TRUE)

df$source.pop<-df$population
df$population<-as.character(df$population)

df$population<-ifelse(df$population=="grey","allo1",df$population)
df$population<-ifelse(df$population=="red","sym",df$population)
df$population<-ifelse(df$population=="yellow","allo2",df$population)
df$population<-ifelse(df$population=="green","allo3",df$population)


df$focal<-factor(df$pop,levels=c("sym","allo1","allo2","allo3"))
df$competition<-as.factor(ifelse(df$type=="alpha","yes","no"))

#df<-subset(df,germ=="y"&competition=="yes")
df$history<-ifelse(df$focal=="sym","local","foreign")

df18.full<-df
df18<-df18.full
colnames(df18)[3]<-c("seeds")
df18$seeds<-ceiling(df18$seeds)

df18<-df18[,c("seeds","focal","germ","cenmel","lotpur","brodia","hemcon","avefat","brohor","plaere","density","competition")]
df18<-na.omit(df18)
df18$history<-ifelse(df18$focal=="sym","local","foreign")


#comparing neighbour density among sites
mean.6<-mean(df6$seeds, na.rm=TRUE)
mean.18<-mean(df18$seeds, na.rm=TRUE)
mean.18/mean.6

mean(rowSums(df6[,4:10])/df6[,11],na.rm=TRUE) #% of community made up by 7 species, 0.9559
mean(rowSums(df18[,4:10])/df18[,11],na.rm=TRUE) #% of community made up by 7 species, 0.9364058


######################################################################
#main analysis (fig. 1a,b) - fitness vs competition, local adaptation
#####################################################################

#statistical analysis
lm6.nb1z<-glmmTMB(seeds~competition*history + (1|focal), data=subset(df6,germ=="y"), ziformula=~., family="nbinom1")
lm6.nb2z<-glmmTMB(seeds~competition*history + (1|focal), data=subset(df6,germ=="y"), ziformula=~., family="nbinom2")
anova(lm6.nb1z,lm6.nb2z)

lm6<-lm6.nb2z #models similar, going with nbinom2 since it matters for site 18, i.e., consistency

lm18.nb1z<-glmmTMB(seeds~competition*history + (1|focal), data=df18, ziformula=~., family="nbinom1")
lm18.nb2z<-glmmTMB(seeds~competition*history + (1|focal), data=df18, ziformula=~., family="nbinom2")
anova(lm18.nb1z,lm18.nb2z)

lm18<-lm18.nb2z #best model

Anova(lm6); Anova(lm18)
summary(lm6)
summary(lm18)

vis6<-visreg(lm6,xvar="competition",by="history",overlay=TRUE, scale="response")
vis18<-visreg(lm18,xvar="competition",by="history",overlay=TRUE, scale="response")


#log ratio and standard error
mydf6.la <- ggpredict(lm6, terms = c("competition","history"), type="fe.zi")
mydf18.la <- ggpredict(lm18, terms = c("competition","history"), type="fe.zi")

plot(mydf6.la)
plot(mydf18.la)

mydf6.la <- ggpredict(lm6, terms = c("history","competition"), type="fe.zi")
mydf18.la <- ggpredict(lm18, terms = c("history","competition"), type="fe.zi")

lrr6.no<-log(mydf6.la[3,2]/mydf6.la[1,2])
lrr6.yes<-log(mydf6.la[4,2]/mydf6.la[2,2])

se6.no<-sqrt(((mydf6.la$std.error[1]^2)/(mydf6.la[1,2]^2)) + ((mydf6.la$std.error[3]^2)/(mydf6.la[3,2]^2)))
se6.yes<-sqrt(((mydf6.la$std.error[2]^2)/(mydf6.la[2,2]^2)) + ((mydf6.la$std.error[4]^2)/(mydf6.la[4,2]^2)))

lrr18.no<-log(mydf18.la[3,2]/mydf18.la[1,2])
lrr18.yes<-log(mydf18.la[4,2]/mydf18.la[2,2])

se18.no<-sqrt(((mydf18.la$std.error[1]^2)/(mydf18.la[1,2]^2)) + ((mydf18.la$std.error[3]^2)/(mydf18.la[3,2]^2)))
se18.yes<-sqrt(((mydf18.la$std.error[2]^2)/(mydf18.la[2,2]^2)) + ((mydf18.la$std.error[4]^2)/(mydf18.la[4,2]^2)))

la.df6<-data.frame(comp=c("no","yes"), lrr=c(lrr6.no,lrr6.yes),se=c(se6.no,se6.yes))
la.df18<-data.frame(comp=c("no","yes"), lrr=c(lrr18.no,lrr18.yes),se=c(se18.no,se18.yes))


###############################################
#main analysis (fig. 1c,d) - fitness vs density
###############################################

hist(df6$density,xlim=c(0,150))
hist(df18$density,xlim=c(0,150))

#stats - density as continuous variable
lm6.nb1z<-glmmTMB(seeds~poly(density,2)*history + (1|focal), data=subset(df6,germ=="y"), ziformula=~., family="nbinom1",control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")))
lm6.nb2z<-glmmTMB(seeds~poly(density,2)*history + (1|focal), data=subset(df6,germ=="y"), ziformula=~., family="nbinom2")
anova(lm6.nb1z,lm6.nb2z)

lm6<-lm6.nb2z

lm18.nb1z<-glmmTMB(seeds~poly(density,2)*history + (1|focal), data=subset(df18,germ=="y"), ziformula=~., family="nbinom1")
lm18.nb2z<-glmmTMB(seeds~poly(density,2)*history + (1|focal), data=subset(df18,germ=="y"), ziformula=~., family="nbinom2")
anova(lm18.nb1z,lm18.nb2z)

lm18<-lm18.nb2z 

#site 6: nonzero part = poly 1 + 2, poly1:interaction, zero part = poly 1 + 2, history (overall marg. sig history x density)
#site 18: density only
Anova(lm6); Anova(lm18)
summary(lm6)
summary(lm18)


#fig1c,d materials
mydf6 <- ggpredict(lm6, terms = c("density [all]","history"), type="fe.zi")
mydf18 <- ggpredict(lm18, terms = c("density [all]","history"), type="fe.zi")


#################
#figure 1a,b,c,d
#################

par(mfrow=c(1,2))

fig1.a<-ggplot(la.df6,aes(x=comp,y=lrr, group=comp)) + 
  geom_errorbar(aes(ymin=lrr-se, ymax=lrr+se), colour="black", width=.05)+
  geom_hline(yintercept=0,color="grey",linetype=2) + geom_point(size=3,pch=21,fill="cornflowerblue") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Neighbours present",y="Strength of LA (log ratio)") +
  scale_x_discrete(labels= c("no","yes"),limits = levels(c("no","yes")))  + ylim(-2.8,2.5) +
  geom_text(x=0.9, y=0.1, label="local adaptation", col="grey") + geom_text(x=1, y=-0.1, label="local maladaptation", col="grey") 

fig1.b<-ggplot(la.df18,aes(x=comp,y=lrr, group=comp)) + 
  geom_errorbar(aes(ymin=lrr-se, ymax=lrr+se), colour="black", width=.05)+
  geom_hline(yintercept=0,color="grey",linetype=2) + geom_point(size=3,pch=21,fill="cornflowerblue") +
  theme_classic() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x="Neighbours present",y="Strength of LA (log ratio)") +
  scale_x_discrete(labels= c("no","yes"),limits = levels(c("no","yes")))  + ylim(-2.8,2.5) 



tmp2<-subset(df6,germ=="y")
colnames(tmp2)[c(1,11,13)]<-c("predicted","x","group")
tmp2$type<-ifelse(tmp2$predicted>0,"closed","open")

fig1.c<-ggplot(data=mydf6, aes(x, predicted, color=group, group=group, fill=group)) +
  geom_hline(yintercept=1,color="grey",linetype=2) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high,color=NULL, fill=group,alpha=group)) +
  geom_line(size=1, aes(linetype=group))+
  xlab("Neighbour density")+
  ylab("Individual fitness") +
  scale_alpha_discrete(range = c(0.15, 0.4)) +
  scale_fill_manual(values=c("deeppink4","cornflowerblue")) +
  scale_color_manual(values=c("deeppink4","cornflowerblue")) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  geom_point(data=tmp2,aes(shape=type),alpha=0.3,size=2) +
  theme_classic() +  theme(legend.position="none"); fig1.c


tmp2<-subset(df18,germ=="y")
colnames(tmp2)[c(1,11,13)]<-c("predicted","x","group")
tmp2$type<-ifelse(tmp2$predicted>0,"closed","open")

fig1.d<-ggplot(mydf18, aes(x, predicted, color=group, group=group, fill=group)) +
  geom_hline(yintercept=1,color="grey",linetype=2) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high,color=NULL, fill=group,alpha=group)) +
  geom_line(size=1, aes(linetype=group))+
  xlab("Neighbour density")+
  ylab("Individual fitness") +
  scale_alpha_discrete(range = c(0.2, 0.4)) +
  scale_fill_manual(values=c("deeppink4","cornflowerblue")) +
  scale_color_manual(values=c("deeppink4","cornflowerblue")) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  labs(fill = "history", color="history",alpha="history",linetype="history") +
  geom_point(data=tmp2,aes(shape=type),alpha=0.3,size=2) +
  theme_classic() +
  scale_shape(guide = 'none') +
  theme(legend.position = c(0.75, 0.5)); fig1.d

plot_grid(fig1.b,fig1.a, labels = "auto")
ggsave("plots/figure 1ab.pdf", w=8,h=4)

plot_grid(fig1.d,fig1.c, labels = "auto")
ggsave("plots/figure 1cd.pdf", w=8,h=4)

plot_grid(fig1.b, fig1.a, fig1.d, fig1.c, labels = "auto", scale=0.95)
ggsave("plots/figure 1.pdf", w=7,h=7)


###############################################
#main analysis for fig. 2 - species specific
###############################################

#############
#site 6

  lm6.0.nb2<-glmmTMB(seeds~
                     brohor*history + 
                     vulmic*history +
                     lotwra*history +
                     hemcon*history +
                     lolmul*history +
                     clagra*history +
                     plaere*history + (1|focal),
                   family="nbinom2",  ziformula=~., data=subset(df6,germ=="y"))

  lm6.1.nb2<-glmmTMB(seeds~density*history +
                     brohor*history + 
                     vulmic*history +
                     lotwra*history +
                     hemcon*history +
                     lolmul*history +
                     clagra*history +
                     plaere*history + (1|focal),
                   family="nbinom2",  ziformula=~., data=subset(df6,germ=="y"))
   
  lm6.2.nb2<-glmmTMB(seeds~poly(density,2)*history +
                     brohor*history + 
                     vulmic*history +
                     lotwra*history +
                     hemcon*history +
                     lolmul*history +
                     clagra*history +
                     plaere*history + (1|focal),
                   family="nbinom2",  ziformula=~., data=subset(df6,germ=="y"))

   
   anova(lm6.0.nb2, lm6.1.nb2, lm6.2.nb2) #most complex model has sig lowest deviance

   lm6<-lm6.2.nb2; Anova(lm6); summary(lm6)

   # par(mfrow=c(2,4))
   # 
   # visreg(lm6, xvar="density", by="history", overlay=TRUE)  
   # visreg(lm6, xvar="brohor", by="history", overlay=TRUE)
   # visreg(lm6, xvar="vulmic", by="history", overlay=TRUE)
   # visreg(lm6, xvar="lolmul", by="history", overlay=TRUE)
   # visreg(lm6, xvar="hemcon", by="history", overlay=TRUE)
   # visreg(lm6, xvar="lotwra", by="history", overlay=TRUE)
   # visreg(lm6, xvar="clagra", by="history", overlay=TRUE)
   # visreg(lm6, xvar="plaere", by="history", overlay=TRUE)

   #multiple comps to extract interaction coefficients
   lm6.l<-glmmTMB(seeds~poly(density,2) +
                   brohor + 
                   vulmic+
                   lotwra +
                   hemcon +
                   lolmul +
                   clagra +
                   plaere,
                 family="nbinom2", ziformula=~.,  data=subset(df6,history=="local"&germ=="y"))
  
   lm6.f<-glmmTMB(seeds~poly(density,2) +
                   brohor + 
                   vulmic+
                   lotwra +
                   hemcon +
                   lolmul +
                   clagra +
                   plaere + (1|focal),
                 family="nbinom2", ziformula=~., data=subset(df6,history=="foreign"&germ=="y"))
  
  
  # #######################################
  # #result without accounting for density
  # #######################################
  # 
  #  lm6<-lm6.0.nb2; Anova(lm6); summary(lm6)
  #  
  #  anova(lm6.0.nb2,lm6.2.nb2)
  #  
  #  lm6.3.nb2<-glmmTMB(seeds~
  #                       poly(brohor,2)*history + 
  #                       poly(vulmic,2)*history +
  #                       poly(lotwra,2)*history +
  #                       poly(hemcon,2)*history +
  #                       poly(lolmul,2)*history +
  #                       poly(clagra,2)*history +
  #                       poly(plaere,2)*history + (1|focal),
  #                     family="nbinom2",  ziformula=~., data=subset(df6,germ=="y"))
  #  
  #  lm6<-lm6.3.nb2; Anova(lm6); summary(lm6)
  #  
  #  anova(lm6.2.nb2,lm6.3.nb2) #density is a better model than poly for each species by AIC
   

  #############
  #site 18 

  lm18.0.nb2<-glmmTMB(seeds~
                    brohor*history + 
                    brodia*history +
                    hemcon*history +
                    lotpur*history +
                    cenmel*history +
                    avefat*history +
                    plaere*history + (1|focal),
                  family="nbinom2", ziformula=~., subset(df18,germ=="y"))
  
  lm18.1.nb2<-glmmTMB(seeds~density*history +
                   brohor*history + 
                   brodia*history +
                   hemcon*history +
                   lotpur*history +
                   cenmel*history +
                   avefat*history +
                   plaere*history + (1|focal),
                  family="nbinom2", ziformula=~., subset(df18,germ=="y"))
  
  lm18.2.nb2<-glmmTMB(seeds~poly(density,2)*history +
                        brohor*history + 
                        brodia*history +
                        hemcon*history +
                        lotpur*history +
                        cenmel*history +
                        avefat*history +
                        plaere*history + (1|focal),
                  family="nbinom2", ziformula=~., subset(df18,germ=="y"),
                  control = glmmTMBControl(optCtrl=list(iter.max=1e3,eval.max=1e3))) #convergence issue
  

  #model comparison
  anova(lm18.0.nb2, lm18.1.nb2, lm18.2.nb2) #no density model best, 3rd model overfit?
  
  #common response of just being there and generically taking up space and light
  lm18<-lm18.0.nb2; Anova(lm18); summary(lm18)
  
  # par(mfrow=c(2,4))
  # 
  # #visreg(lm18, xvar="density", by="history", overlay=TRUE, scale = "response")  
  # visreg(lm18, xvar="brohor", by="history", overlay=TRUE, scale = "response")   
  # visreg(lm18, xvar="brodia", by="history", overlay=TRUE, scale = "response")     
  # visreg(lm18, xvar="hemcon", by="history", overlay=TRUE, scale = "response")   
  # visreg(lm18, xvar="lotpur", by="history", overlay=TRUE, scale = "response")   
  # visreg(lm18, xvar="cenmel", by="history", overlay=TRUE, scale = "response")    
  # visreg(lm18, xvar="avefat", by="history", overlay=TRUE, scale = "response")   
  # visreg(lm18, xvar="plaere", by="history", overlay=TRUE, scale = "response")   
  

  #multiple comps to extract interaction coefficients
  lm18.l<-glmmTMB(seeds~
                   brohor + 
                   brodia+
                   hemcon +
                   lotpur +
                   cenmel +
                   avefat +
                   plaere,
                 family="nbinom2", ziformula=~.,data=subset(df18,history=="local"&germ=="y")); Anova(lm18.l)
  
  lm18.f<-glmmTMB(seeds~
                   brohor + 
                   brodia +
                   hemcon +
                   lotpur +
                   cenmel +
                   avefat +
                   plaere + (1|focal),
                 family="nbinom2", ziformula=~., data=subset(df18,history=="foreign"&germ=="y")); Anova(lm18.f)

  
  #diagnostics #looks nice :-)
  require(DHARMa) 
  plot(simulateResiduals(lm6))
  plot(simulateResiduals(lm18))
  
  ##########################
  #figure 2
  ##########################
 
  require(egg)
  
  #changes to make in Affinity:
  #remove border around facets
  #add line drawings
  
  #coefficient df
  #site 6
  site6.local<-as.data.frame(confint(lm6.l))[-c(1:3,11:13,22),]
  site6.local$history<-"local"
  site6.local$species<-c("brohor","vulmic","lotwra","hemcon","lolmul","clagra","plaere","brohor","vulmic","lotwra","hemcon","lolmul","clagra","plaere")
  site6.local$type<-c(rep("cond",7),rep("zi",7))
  
  site6.foreign<-as.data.frame(confint(lm6.f))[-c(1:3,11:14,22),]
  site6.foreign$history<-"foreign"
  site6.foreign$species<-c("brohor","vulmic","lotwra","hemcon","lolmul","clagra","plaere","brohor","vulmic","lotwra","hemcon","lolmul","clagra","plaere")
  site6.foreign$type<-c(rep("cond",7),rep("zi",7))
  
  site6.all<-rbind(site6.local,site6.foreign)
  colnames(site6.all)<-c("CI.low","CI.high","Estimate","history","species","type")
  
  #coefficient df
  #site 18
  site18.local<-as.data.frame(confint(lm18.l))[-c(1,9),]
  site18.local$history<-"local"
  site18.local$species<-c("brohor","brodia","hemcon","lotpur","cenmel","avefat","plaere","brohor","brodia","hemcon","lotpur","cenmel","avefat","plaere")
  site18.local$type<-c(rep("cond",7),rep("zi",7))
  
  site18.foreign<-as.data.frame(confint(lm18.f))[-c(1,9:10,18),]
  site18.foreign$history<-"foreign"
  site18.foreign$species<-c("brohor","brodia","hemcon","lotpur","cenmel","avefat","plaere","brohor","brodia","hemcon","lotpur","cenmel","avefat","plaere")
  site18.foreign$type<-c(rep("cond",7),rep("zi",7))
  
  site18.all<-rbind(site18.local,site18.foreign)
  colnames(site18.all)<-c("CI.low","CI.high","Estimate","history","species","type")
  
  
  
  ##site 6########
  site6.all$species = factor(site6.all$species, levels=c("brohor","vulmic","clagra","plaere","hemcon","lotwra","lolmul"))
  
  my_tag <- c(".","","..","","..","",".")
  
  fig2c<-ggplot(data=subset(site6.all,type=="cond"),aes(x=history,y=Estimate, group=species, color=species)) + 
    scale_x_discrete(labels= c("F","L")) + 
    geom_hline(yintercept=0, color="grey", size=0.5, linetype="dashed") + 
    geom_errorbar(aes(ymin=CI.low, ymax=CI.high),color="black" ,width=.1) +
    geom_line(show.legend = FALSE) +
    geom_point(size=4,pch=21,color="black",aes(fill=species), show.legend = FALSE)+ 
    facet_wrap(~species,nrow=1) + theme_bw() +
    theme(strip.background = element_blank(), strip.text.x = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
    scale_fill_manual(values=c("cornflowerblue","seagreen3","yellow3","orange2","firebrick2","violetred","darkorchid")) +
    scale_color_manual(values=c("cornflowerblue","seagreen3","yellow3","orange2","firebrick2","violetred","darkorchid")) +
    xlab("Population history") + ylab("Interaction strength") +
    scale_y_continuous(limits = c(-0.95, 0.3), breaks =  c(-0.9,-0.6,-0.3,0,0.30)); fig2c
  
  fig2c<-tag_facet(fig2c,
            x = 1.5, y = 0.25, 
            vjust = 0, hjust = 0.25,
            open = "", close = "",
            fontface = 4,
            size = 7,
            family = "serif",
            tag_pool = my_tag)
  
  x <- ggpredict(lm6, terms = c("brohor[all]","history"), type="fe")
  
  fig2d<-ggplot(data=x,aes(x=x,y=predicted, group=group, color=NULL, fill=group)) +
    labs(y="Individual fitness",x=expression(paste("Density of ",italic("Bromus hordeaceus")))) +
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
    scale_color_manual(values=c("deeppink4","cornflowerblue")) +
    scale_fill_manual(values=c("deeppink4","cornflowerblue")) +
    scale_linetype_manual(values=c("dashed","solid")) +
    scale_alpha_discrete(range = c(0.15, 0.4)) +
    geom_line(size=1, show.legend = FALSE, aes(linetype=group,color=group))+
    geom_ribbon(data=x$fit,aes(ymin=conf.low, ymax=conf.high,fill=group, alpha=group), show.legend = FALSE) +
    coord_cartesian(ylim = c(0, 40))
  
  
  ##site 18########
  site18.all$species = factor(site18.all$species, levels=c("brohor","lotpur","avefat","plaere","hemcon","brodia","cenmel"))
  
  my_tag <- c("","","",".","","..","")

  fig2a<-ggplot(data=subset(site18.all,type=="cond"),aes(x=history,y=Estimate, group=species, color=species)) + 
    scale_x_discrete(labels= c("F","L")) + 
    geom_hline(yintercept=0, color="grey", size=0.5, linetype="dashed") + 
    geom_errorbar(aes(ymin=CI.low, ymax=CI.high),color="black" ,width=.1) +
    geom_line(show.legend = FALSE) +
    geom_point(size=4,pch=21,color="black",aes(fill=species), show.legend = FALSE)+ 
    facet_wrap(~species,nrow=1) + theme_bw() +
    theme(strip.background = element_blank(), strip.text.x = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
    scale_fill_manual(values=c("cornflowerblue","seagreen3","yellow3","orange2","firebrick2","violetred","darkorchid")) +
    scale_color_manual(values=c("cornflowerblue","seagreen3","yellow3","orange2","firebrick2","violetred","darkorchid")) +
    xlab(NULL) + ylab("Interaction strength") + 
    scale_y_continuous(limits = c(-0.95, 0.25), breaks =  c(-0.9,-0.6,-0.3,0,0.30)); fig2a
    
  fig2a<-tag_facet(fig2a,
              x = 1.5, y = 0.25, 
              vjust = 0, hjust = 0.25,
              open = "", close = "",
              fontface = 4,
              size = 7,
              family = "serif",
              tag_pool = my_tag)
  

    y <- ggpredict(lm18, terms = c("brohor[all]","history"), type="fe")
  
    fig2b<-ggplot(data=y,aes(x=x,y=predicted, group=group, color=NULL, fill=group)) +
      labs(y="Individual fitness",x=expression(paste("Density of ",italic("Bromus hordeaceus")))) +
      theme_classic() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
      scale_color_manual(values=c("deeppink4","cornflowerblue")) +
      scale_fill_manual(values=c("deeppink4","cornflowerblue")) +
      scale_linetype_manual(values=c("dashed","solid")) +
      scale_alpha_discrete(range = c(0.15, 0.4)) +
      geom_line(size=1, show.legend = FALSE, aes(linetype=group,color=group))+
      geom_ribbon(data=x$fit,aes(ymin=conf.low, ymax=conf.high,fill=group, alpha=group), show.legend = FALSE) +
      coord_cartesian(ylim = c(0, 40))
    
  plot_grid(fig2a,fig2b,fig2c,fig2d,labels=c("a","b","c","d"),rel_widths=c(3,1))
  ggsave("plots/figure 2.pdf",height=5.5,width=11)
  
#########################
#figure S4 
#########################
  
  plot(mydf6.la)
  plot(mydf18.la)
  
  mydf6.la<-as.data.frame(mydf6.la)
  mydf18.la<-as.data.frame(mydf18.la)

  lrr6.f<-log(mydf6.la[2,2]/mydf6.la[1,2])
  lrr6.l<-log(mydf6.la[4,2]/mydf6.la[3,2])
  
  se6.f<-sqrt(((mydf6.la[2,3]^2)/(mydf6.la[1,2]^2)) + ((mydf6.la[2,3]^2)/(mydf6.la[1,2]^2)))
  se6.l<-sqrt(((mydf6.la[4,3]^2)/(mydf6.la[3,2]^2)) + ((mydf6.la[4,3]^2)/(mydf6.la[3,2]^2)))

  lrr18.f<-log(mydf18.la[2,2]/mydf18.la[1,2])
  lrr18.l<-log(mydf18.la[4,2]/mydf18.la[3,2])
  
  se18.f<-sqrt(((mydf18.la[2,3]^2)/(mydf18.la[1,2]^2)) + ((mydf18.la[2,3]^2)/(mydf18.la[1,2]^2)))
  se18.l<-sqrt(((mydf18.la[4,3]^2)/(mydf18.la[3,2]^2)) + ((mydf18.la[4,3]^2)/(mydf18.la[3,2]^2)))
  
  la.df6.comp<-data.frame(history=c("f","l"), lrr=c(lrr6.f,lrr6.l),se=c(se6.f,se6.l))
  la.df18.comp<-data.frame(history=c("f","l"), lrr=c(lrr18.f,lrr18.l),se=c(se18.f,se18.l))
  

  #site 6

 figS4.a <- ggplot(mydf18.la,aes(x=x,y=predicted, group=group, fill=group)) + 
   geom_errorbar(aes(ymin=predicted-std.error, ymax=predicted+std.error), width=.05, col="black") + geom_line() + geom_point(pch=21,size=3) + 
   theme_classic() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   labs(x="Population history",y="Individual fitness") +
   scale_x_discrete(labels= c("Foreign","Local"))+
   scale_fill_manual(name = "neighbors present", values = c("grey","black")) +
   theme_classic() +  theme(legend.position="none"); figS4.a
  
  figS4.b <- ggplot(mydf6.la,aes(x=x,y=predicted, group=group, fill=group)) + 
    geom_errorbar(aes(ymin=predicted-std.error, ymax=predicted+std.error), width=.05, col="black") + geom_line() + geom_point(pch=21,size=3) + 
    theme_classic() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x="Population history",y="Individual fitness") +
    scale_x_discrete(labels= c("Foreign","Local"))+
    scale_fill_manual(name = "neighbors present", values = c("grey","black")) +
    theme_classic() +  theme(legend.position=c(0.7,0.7)); figS4.b
  
  figS4.c <- ggplot(la.df18.comp,aes(x=history,y=lrr)) + 
    geom_errorbar(aes(ymin=lrr-se, ymax=lrr+se), width=.05, col="black") + 
    geom_point(pch=21,size=3,fill="black") + 
    theme_classic() + theme(plot.title = element_text(hjust = 0.5),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x="Population history",y="Individual fitness") +
    scale_x_discrete(labels= c("Foreign","Local"))+
    geom_hline(yintercept=0,color="grey",linetype=2) +
    ylim(-2,2) +
    geom_text(x=0.8, y=0.1, label=" facilitation", col="grey") + 
    geom_text(x=0.85, y=-0.1, label="competition", col="grey") +
    theme_classic();   figS4.c

  figS4.d <- ggplot(la.df6.comp,aes(x=history,y=lrr, fill)) + 
    geom_errorbar(aes(ymin=lrr-se, ymax=lrr+se), width=.05, col="black") + 
    geom_point(pch=21,size=3,fill="black") + 
    theme_classic() + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x="Population history",y="Individual fitness") +
    scale_x_discrete(labels= c("Foreign","Local"))+
    theme_classic()  +
    geom_hline(yintercept=0,color="grey",linetype=2) +
    ylim(-2,2) +
    geom_text(x=0.8, y=0.1, label=" facilitation", col="grey") + 
    geom_text(x=0.85, y=-0.1, label="competition", col="grey"); figS4.d
  

  plot_grid(figS4.a, figS4.b, figS4.c, figS4.d, labels = "auto", scale=0.95)
  ggsave("plots/Fig S4.pdf", w=8,h=8)
  

  
  ###############################################
  #including all transplants (table S2)
  ###############################################

  #stats - density as continuous variable
  lm6.nb2z<-glmmTMB(seeds~poly(density,2)*history + (1|focal), data=df6, ziformula=~., family="nbinom2")
  lm6<-lm6.nb2z
  
  lm18.nb2z<-glmmTMB(seeds~poly(density,2)*history + (1|focal), data=df18, ziformula=~., family="nbinom2")
  lm18<-lm18.nb2z 
  
  #site 6: nonzero part = poly 1 + 2, poly1:interaction, zero part = poly 1 + 2, history (overall marg. sig history x density)
  #site 18: density only
  Anova(lm6); Anova(lm18)
  summary(lm6)
  summary(lm18)
  

  
  ######################
  #figure S5
  ######################

  tmp<-subset(df18,germ=="y")
  tmp$fits<-fitted.values(lm18.nb2z)
  tmp2<-subset(df18,germ=="y")
  colnames(tmp2)[1]<-"fits"
  tmp2$type<-ifelse(tmp2$fits>0,"closed","open")
  
  a<-ggplot(data=tmp, aes(density, fits, color=history, group=focal, fill=history)) +
    geom_hline(yintercept=1,color="grey",linetype=2) +
    #geom_ribbon(aes(ymin=conf.low, ymax=conf.high,color=NULL, fill=group), alpha=0.3) +
    geom_line(size=1, aes(linetype=history))+
    xlab("Neighbour density")+
    ylab("Individual fitness") +
    scale_fill_manual(values=c("deeppink4","cornflowerblue")) +
    scale_color_manual(values=c("deeppink4","cornflowerblue")) +
    scale_linetype_manual(values=c("dashed", "solid")) +
    geom_point(data=tmp2,aes(shape=type),alpha=0.3,size=2) +
    theme_classic() +  theme(legend.position="none") +
    facet_wrap(~focal);b
  
  
  tmp<-subset(df6,germ=="y")
  tmp$fits<-fitted.values(lm6.nb2z)
  tmp2<-subset(df6,germ=="y")
  colnames(tmp2)[1]<-"fits"
  tmp2$type<-ifelse(tmp2$fits>0,"closed","open")
  
  b<-ggplot(data=tmp, aes(density, fits, color=history, group=focal, fill=history)) +
    geom_hline(yintercept=1,color="grey",linetype=2) +
    #geom_ribbon(aes(ymin=conf.low, ymax=conf.high,color=NULL, fill=group), alpha=0.3) +
    geom_line(size=1, aes(linetype=history))+
    xlab("Neighbour density")+
    ylab("Individual fitness") +
    scale_fill_manual(values=c("deeppink4","cornflowerblue")) +
    scale_color_manual(values=c("deeppink4","cornflowerblue")) +
    scale_linetype_manual(values=c("dashed", "solid")) +
    geom_point(data=tmp2,aes(shape=type),alpha=0.3,size=2) +
    theme_classic() +  theme(legend.position="none") +
    facet_wrap(~focal);a
  
  plot_grid(a,b,labels=c("a","b"))
  
  #################################
  #figure S6 - fitness components 
  #################################
  
  par(mfrow=c(2,2))
  
  #germination
  #site 6
  df6$germ.bin<-ifelse(df6$germ=="y",1,0)
  
  lm<-glmmTMB(germ.bin~history*poly(density,2)+(1|focal), data=df6, family="binomial")
  Anova(lm)
  summary(lm)
  vis6.g<-visreg(lm,xvar="density", by="history", overlay=TRUE,scale="response")
  
  #site 18
  df18$germ.bin<-ifelse(df18$germ=="y",1,0)
  
  lm<-glmmTMB(germ.bin~history*poly(density,2)+(1|focal), data=df18, family="binomial")
  Anova(lm)
  summary(lm)
  vis18.g<-visreg(lm,xvar="density", by="history", overlay=TRUE,scale="response")
  
  #plots
  fig.s5a<-ggplot(data=vis6.g$fit,aes(x = density, y = visregFit, color = history, group=history)) + 
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr,color=NULL, fill=history),alpha=0.3) +
    geom_line(size=1,aes(linetype=history))+
    xlab("Neighbour density")+
    ylab("Germination rate") +
    scale_fill_manual(values=c("deeppink4","cornflowerblue")) +
    scale_color_manual(values=c("deeppink4","cornflowerblue")) +
    scale_linetype_manual(values=c("dashed", "solid")) +
    theme_classic() +  theme(legend.position="none") + ylim(0,1) +
    theme(legend.position= c(0.5,0.5)) 
  
  fig.s5c<-ggplot(data=vis18.g$fit,aes(x = density, y = visregFit, color = history, group=history)) + 
    geom_ribbon(aes(ymin=visregLwr, ymax=visregUpr,color=NULL, fill=history),alpha=0.3) +
    geom_line(size=1,aes(linetype=history))+
    xlab("Neighbour density")+
    ylab("Germination rate") +
    scale_fill_manual(values=c("deeppink4","cornflowerblue")) +
    scale_color_manual(values=c("deeppink4","cornflowerblue")) +
    scale_linetype_manual(values=c("dashed", "solid")) +
    theme_classic() +  theme(legend.position="none") + ylim(0,1) 
  
  
  plot_grid(fig.s5c, fig.s5a, labels = "auto", scale=0.95)
  ggsave("plots/Fig S5.pdf", w=7,h=3.5)
  
  
  ###################################
  #figure S8 - supplement for Zi part
  ###################################
  
  ##site 6########
  site6.all$species = factor(site6.all$species, levels=c("brohor","vulmic","clagra","plaere","hemcon","lotwra","lolmul"))
  
  figs8b<-ggplot(data=subset(site6.all,type=="zi"),aes(x=history,y=Estimate, group=species, color=species)) + 
    scale_x_discrete(labels= c("F","L")) + 
    geom_hline(yintercept=0, color="grey", size=0.5, linetype="dashed") + 
    geom_errorbar(aes(ymin=CI.low, ymax=CI.high),color="black" ,width=.1) +
    geom_line(show.legend = FALSE) +
    geom_point(size=4,pch=21,color="black",aes(fill=species), show.legend = FALSE)+ 
    facet_wrap(~species,nrow=1) + theme_bw() +
    theme(strip.background = element_blank(), strip.text.x = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
    scale_fill_manual(values=c("cornflowerblue","seagreen3","yellow3","orange2","firebrick2","violetred","darkorchid")) +
    scale_color_manual(values=c("cornflowerblue","seagreen3","yellow3","orange2","firebrick2","violetred","darkorchid")) +
    xlab("Population history") + ylab("Interaction strength") +
    scale_y_continuous(limits = c(-3.4, 1)); figs8b
  
  
  ##site 18########
  site18.all$species = factor(site18.all$species, levels=c("brohor","lotpur","avefat","plaere","hemcon","brodia","cenmel"))
  
  my_tag <- c(".","","","",".","","..")
  
  figs8a<-ggplot(data=subset(site18.all,type=="zi"),aes(x=history,y=Estimate, group=species, color=species)) + 
    scale_x_discrete(labels= c("F","L")) + 
    geom_hline(yintercept=0, color="grey", size=0.5, linetype="dashed") + 
    geom_errorbar(aes(ymin=CI.low, ymax=CI.high),color="black" ,width=.1) +
    geom_line(show.legend = FALSE) +
    geom_point(size=4,pch=21,color="black",aes(fill=species), show.legend = FALSE)+ 
    facet_wrap(~species,nrow=1) + theme_bw() +
    theme(strip.background = element_blank(), strip.text.x = element_blank())+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
    scale_fill_manual(values=c("cornflowerblue","seagreen3","yellow3","orange2","firebrick2","violetred","darkorchid")) +
    scale_color_manual(values=c("cornflowerblue","seagreen3","yellow3","orange2","firebrick2","violetred","darkorchid")) +
    xlab(NULL) + ylab("Interaction strength") +
    scale_y_continuous(limits = c(-3.4, 1))
  
  figs8a<-tag_facet(figs8a,
                   x = 1.5, y = 1, 
                   vjust = 0, hjust = 0.25,
                   open = "", close = "",
                   fontface = 4,
                   size = 7,
                   family = "serif",
                   tag_pool = my_tag)
  
  
  plot_grid(figs8a,figs8b,labels=c("a","b"),nrow=2)
  ggsave("plots/Fig s8zi.pdf",height=5,width=10)
  