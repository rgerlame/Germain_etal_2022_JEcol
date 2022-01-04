
library(tidyverse);require(grid); require(cowplot); require(car); require(glmmTMB)
require(ggplot2); require(magick); require(visreg)
library(DHARMa)


df <- read.csv("data/seed mass.csv")
df1 <- read.csv("data/gh traits.csv")



#plot aesthetics
ng1 <- theme(aspect.ratio=1.0,panel.background = element_blank(), 
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             panel.border=element_blank(),
             axis.line = element_line(size=1), 
             axis.line.x = element_line(color="black", size = 1),
             axis.line.y = element_line(color="black", size = 1),
             axis.ticks=element_line(color="black"), 
             axis.text=element_text(color="black"), 
             axis.title=element_text(color="black"), 
             axis.title.y=element_text(vjust=0.2, size=12),
             axis.title.x=element_text(vjust=0.1,size=12),
             axis.text.x=element_text(size=10),
             axis.text.y=element_text(size=10),
             #legend.position="none",
             #legend.title=element_blank(),
             legend.key = element_rect(fill = "white"),
             plot.title = element_blank())


###################
#total seed mass

df1$tot.seed.mass.perplant <- df1$total.seed.mass / df1$num.plants

df2 <- df1 %>% 
  mutate(`origin` = stringr::str_replace(`origin`, "G", "Greenhouse")) %>% 
  mutate(`origin` = stringr::str_replace(`origin`, "F", "Field"))

df2$pop<-as.factor(df2$pop)
df2$origin<-as.factor(df2$origin)

figS1c<- ggplot(df2, aes(x =  pop, y = tot.seed.mass.perplant, fill = origin)) +
  geom_boxplot() +  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  #scale_x_discrete(labels=1:4)+
  xlab("Population") +
  ylab("Fecundity (g seed/plant)") + theme(legend.position=c(1.5,0.5)) +
  ng1; figS1c

lm1 <- glmmTMB(tot.seed.mass.perplant ~ pop * origin, data = df1, family=gaussian)
Anova(lm1)
visreg(lm1,xvar="pop",by="origin", overlay=TRUE)


###################
# mass per seed figure

df <- read.csv("data/seed mass.csv")
df_sm <- df %>% 
  mutate(`origin` = stringr::str_replace(`origin`, "G", "Greenhouse")) %>% 
  mutate(`origin` = stringr::str_replace(`origin`, "F", "Field"))

df_sm$origin<-as.factor(df_sm$origin)

figS1a <- ggplot(df_sm, aes(x = pop, y = weight, fill = origin)) +
  geom_boxplot() +  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  #scale_x_discrete(labels=1:4)+
  xlab("Population") +
  ylab("Seed mass (g)") + theme(legend.position="none") +
  ng1; figS1a


lm1 <- glmmTMB(weight ~ pop * origin, data = df, family=gaussian)
Anova(lm1)
visreg(lm1,xvar="pop",by="origin", overlay=TRUE)


###################
# est.height.avg figure

df1 <- read.csv("data/gh traits.csv")
df_gh <- df1 %>% 
  mutate(`origin` = stringr::str_replace(`origin`, "G", "Greenhouse")) %>% 
  mutate(`origin` = stringr::str_replace(`origin`, "F", "Field"))

df_gh$origin<-as.factor(df_gh$origin)
df_gh<-as.data.frame(df_gh)

figS1b <- ggplot(df_gh, aes(x = pop, y = est.height.avg, fill = origin)) +
  geom_boxplot() +  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  #scale_x_discrete(labels=1:4)+
  xlab("Population") +
  ylab("Plant height (cm)") + theme(legend.position="none") +
  ng1; figS1b 

lm1 <- glmmTMB(est.height.avg ~ pop * origin, data = df1, family=gaussian)
Anova(lm1)
summary(lm1)
visreg(lm1,xvar="pop",by="origin", overlay=TRUE)



#############
#all together

pdf('plots/Fig S2.pdf',width=7,height=7)

plot_grid(figS1a,figS1b,figS1c, labels = c("a","b","c"),align="hv",label_fontface="bold")

dev.off()
