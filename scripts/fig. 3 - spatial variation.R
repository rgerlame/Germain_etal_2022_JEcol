
{
require(dplyr);require(ggplot2);require(lmerTest);require(visreg);require(pscl);require(car):require(doBy);require(scales)
require(grid); require(cowplot); library(effects); require(buildmer)
}

###############
#data prep

df1<-read.csv("data/site 6 - seed data.csv")
df2<-read.csv("data/site 6 - neighbourhood data.csv")
df2<-df2[,-2]

df<-merge(df1,df2,by="id")
df$seed.number<-as.numeric(as.character(df$seed.number))
df$any.seed<-ifelse(df$seed.number>0,1,0)
df$hard.seed<-df$seed.number*(1-df$prop.mushy)

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
df1<-df
colnames(df1)[3]<-c("seeds")

df1$density.inter<-rowSums(df[,c(10:11,13:44)],na.rm=TRUE)

df1<-df1[,c("seeds","focal","brohor","vulmic","clagra","plaere","hemcon","lolmul","lotwra","competition","density","id")]


df1<-na.omit(df1)
df1$history<-ifelse(df1$focal=="sym","local","foreign")

#df1$history<-df1$focal

fig1.1<-ggplot(data=subset(df1,competition=="yes"),aes(y=brohor,x=as.numeric(id))) + geom_line(col="cornflowerblue") +
  geom_smooth(formula = y~poly(x,26),method="glm", se=FALSE,method.args = list(family = "poisson"),col="cornflowerblue") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="brohor") + theme(axis.text.x=element_blank())

fig1.2<-ggplot(data=subset(df1,competition=="yes"),aes(y=vulmic,x=as.numeric(id))) + geom_line(col="seagreen3") +
  geom_smooth(formula = y~poly(x,26),method="glm", se=FALSE,method.args = list(family = "poisson"),col="seagreen3") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="vulmic") + theme(axis.text.x=element_blank())

fig1.3<-ggplot(data=subset(df1,competition=="yes"),aes(y=clagra,x=as.numeric(id))) + geom_line(col="yellow3") +
  geom_smooth(formula = y~poly(x,26),method="glm", se=FALSE,method.args = list(family = "poisson"),col="yellow3") + 
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="clagra") + theme(axis.text.x=element_blank())

fig1.4<-ggplot(data=subset(df1,competition=="yes"),aes(y=plaere,x=as.numeric(id))) + geom_line(col="orange2") +
  geom_smooth(formula = y~poly(x,26),method="glm", se=FALSE,method.args = list(family = "poisson"),col="orange2") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="plaere") + theme(axis.text.x=element_blank())

fig1.5<-ggplot(data=subset(df1,competition=="yes"),aes(y=hemcon,x=as.numeric(id))) + geom_line(col="firebrick2") +
  geom_smooth(formula = y~poly(x,26),method="glm", se=FALSE,method.args = list(family = "poisson"),col="firebrick2") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="hemcon") + theme(axis.text.x=element_blank())

fig1.6<-ggplot(data=subset(df1,competition=="yes"),aes(y=lotwra,x=as.numeric(id))) + geom_line(col="violetred") +
  geom_smooth(formula = y~poly(x,26),method="glm", se=FALSE,method.args = list(family = "poisson"),col="violetred") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="lotwra") + theme(axis.text.x=element_blank())

fig1.7<-ggplot(data=subset(df1,competition=="yes"),aes(y=lolmul,x=as.numeric(id))) + geom_line(col="darkorchid") +
  geom_smooth(formula = y~poly(x,8),method="glm", se=FALSE,method.args = list(family = "poisson"),col="darkorchid") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="lolmul") + theme(axis.text.x=element_blank())

site6<-plot_grid(fig1.1,fig1.2,fig1.3,fig1.4,fig1.5,fig1.6,fig1.7, ncol=1,align="v")

mean(df1$brohor); median(df1$brohor)
mean(df1$vulmic); median(df1$vulmic)
mean(df1$clagra); median(df1$clagra)
mean(df1$plaere); median(df1$plaere)
mean(df1$hemcon); median(df1$hemcon)
mean(df1$lotwra); median(df1$lotwra)
mean(df1$lolmul); median(df1$lolmul)

######################
#site 18

###############
#data prep

df1<-read.csv("data/site 18 - seed data.csv")
df2<-read.csv("data/site 18 - neighbourhood data.csv")
df2<-df2[,-2]

df<-merge(df1,df2,by="id")
df$seed.number<-as.numeric(as.character(df$seed.number))
df$any.seed<-ifelse(df$seed.number>0,1,0)
df$hard.seed<-df$seed.number*(1-df$prop.mushy)

df$source.pop<-df$population
df$population<-as.character(df$population)

df$population<-ifelse(df$population=="grey","sym",df$population)
df$population<-ifelse(df$population=="red","allo1",df$population)
df$population<-ifelse(df$population=="yellow","allo2",df$population)
df$population<-ifelse(df$population=="green","allo3",df$population)

df$focal<-factor(df$pop,levels=c("sym","allo1","allo2","allo3"))
df$competition<-as.factor(ifelse(df$type=="alpha","yes","no"))

#df<-subset(df,germ=="y"&competition=="yes")
df1<-df
colnames(df1)[3]<-c("seeds")

df1<-df1[,c("seeds","focal","brohor","brodia","cenmel","plaere","hemcon","avefat","lotpur","competition","id")]

df1$id<-as.numeric(as.character(df1$id))


df1<-na.omit(df1)
df1$history<-ifelse(df1$focal=="sym","local","foreign")

df1<-subset(df1,competition=="yes")

#df1$history<-df1$focal

fig1.1<-ggplot(data=subset(df1,competition=="yes"),aes(y=brohor,x=as.numeric(id))) + geom_line(col="cornflowerblue") +
  geom_smooth(formula = y~poly(x,25),method="glm", se=FALSE,method.args = list(family = "poisson"),col="cornflowerblue") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="brohor") + theme(axis.text.x=element_blank())

fig1.2<-ggplot(data=subset(df1,competition=="yes"),aes(y=lotpur,x=as.numeric(id))) + geom_line(col="seagreen3") +
  geom_smooth(formula = y~poly(x,25),method="glm", se=FALSE,method.args = list(family = "poisson"),col="seagreen3") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="lotpur") + theme(axis.text.x=element_blank())

fig1.3<-ggplot(data=subset(df1,competition=="yes"),aes(y=avefat,x=as.numeric(id))) + geom_line(col="yellow3") +
  geom_smooth(formula = y~poly(x,25),method="glm", se=FALSE,method.args = list(family = "poisson"),col="yellow3") + 
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="avefat") + theme(axis.text.x=element_blank())

fig1.4<-ggplot(data=subset(df1,competition=="yes"),aes(y=plaere,x=as.numeric(id))) + geom_line(col="orange2") +
  geom_smooth(formula = y~poly(x,25),method="glm", se=FALSE,method.args = list(family = "poisson"),col="orange2") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="plaere") + theme(axis.text.x=element_blank())

fig1.5<-ggplot(data=subset(df1,competition=="yes"),aes(y=hemcon,x=as.numeric(id))) + geom_line(col="firebrick2") +
  geom_smooth(formula = y~poly(x,25),method="glm", se=FALSE,method.args = list(family = "poisson"),col="firebrick2") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="he,con") + theme(axis.text.x=element_blank())

fig1.6<-ggplot(data=subset(df1,competition=="yes"),aes(y=brodia,x=as.numeric(id))) + geom_line(col="violetred") +
  geom_smooth(formula = y~poly(x,25),method="glm", se=FALSE,method.args = list(family = "poisson"),col="violetred") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="brodia") + theme(axis.text.x=element_blank())

fig1.7<-ggplot(data=subset(df1,competition=="yes"),aes(y=cenmel,x=as.numeric(id))) + geom_line(col="darkorchid") +
  geom_smooth(formula = y~poly(x,8),method="glm", se=FALSE,method.args = list(family = "poisson"),col="darkorchid") +
  ylim(0,150) + theme_classic() +   labs(x=NULL,y="cenmel") + theme(axis.text.x=element_blank())

site18<-plot_grid(fig1.1,fig1.2,fig1.3,fig1.4,fig1.5,fig1.6,fig1.7, ncol=1,align="v")

plot_grid(site18,site6, labels = c("a","b"))
ggsave("plots/community turnover all - wo loess.pdf",h=10,w=10)




mean(df1$brohor); median(df1$brohor)
mean(df1$lotpur); median(df1$lotpur)
mean(df1$avefat); median(df1$avefat)
mean(df1$plaere); median(df1$plaere)
mean(df1$hemcon); median(df1$hemcon)
mean(df1$brodia); median(df1$brodia)
mean(df1$cenmel); median(df1$cenmel)

##########
#pie chart site 18

#get other!!!

df <- data.frame(
  group = c("bh", "lp", "af","pe","hc","bd","cm"),
  value = c(43.5, 8.8, 5.2,1.3,1.6,1.6,1))
head(df)

sum(df$value)

df$group<-factor(df$group,levels=c("bh", "lp", "af","pe","hc","bd","cm"))

bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
bp

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("cornflowerblue", "seagreen3", "yellow2","orange2","firebrick2","violetred","darkorchid"))

#pie chart site 6

df <- data.frame(
  group = c("bh", "lp", "af","pe","hc","bd","cm"),
  value = c(21, 9.2, 4.0,3.9,2.2,1,0.8))
head(df)

sum(df$value)

df$group<-factor(df$group,levels=c("bh", "lp", "af","pe","hc","bd","cm"))

bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
bp

pie <- bp + coord_polar("y", start=0)

pie + scale_fill_manual(values=c("cornflowerblue", "seagreen3", "yellow2","orange2","firebrick2","violetred","darkorchid"))
