#### HOI SUPP. GRAPHS ####

#code started by M. Raymundo, finished by R. Germain

require(vegan); library(ggplot2); require(grid); require(cowplot)

#data
site6Diversity <- read.csv('data/site 6 - neighbourhood data.csv', header=T, sep=',', strip.white=T)
site18Diversity <- read.csv('data/site 18 - neighbourhood data.csv', header=T, sep=',', strip.white=T)


## site 6 ##
str(site6Diversity)
site6Diversity$site <- as.factor(site6Diversity$site)

site6Diversity<-subset(site6Diversity,type=="alpha")
site.6.Matrix <- site6Diversity[,6:40] ### Create abundance matrix for calculation (1:36, because last 2 columns are site and id)

Indices6 <- site6Diversity[,c("site", "id","type")] 
Indices6$Diversity <- diversity(site.6.Matrix)
Indices6$Evenness <- (diversity(site.6.Matrix))/log(specnumber(site.6.Matrix)) 
Indices6$SppRichness <- specnumber(site.6.Matrix)
Indices6$Density <- rowSums(site.6.Matrix)


## site 18 ##
str(site18Diversity)
site18Diversity$site <- as.factor(site18Diversity$site)

site18Diversity<-subset(site18Diversity,type=="alpha")
site.18.Matrix <- site18Diversity[,c(6:54)] 

Indices18 <- site18Diversity[,c("site", "id","type")] 
Indices18$Diversity <- diversity(site.18.Matrix)
Indices18$Evenness <- (diversity(site.18.Matrix))/log(specnumber(site.18.Matrix))
Indices18$SppRichness <- specnumber(site.18.Matrix)
Indices18$Density <- rowSums(site.18.Matrix)

#combine site 6 and 18
Indices<-rbind(Indices6,Indices18)

################ Diversity and Density ################
Diversity <- Indices

DensityStats <- lm(Density ~ site, data = Diversity)
anova(DensityStats)

DiversityStats <- lm(Diversity ~ site, data = Diversity)
anova(DiversityStats)

EvennessStats <- lm(Evenness ~ site, data = Diversity)
anova(EvennessStats)

#plotting
fig.s4c<-ggplot(Diversity, aes(x=site, y=SppRichness)) + geom_violin(fill="dodgerblue",color="white",alpha=0.5) +
         geom_boxplot(fill="dodgerblue",width = 0.25)+
         #geom_jitter(alpha=0.3,pch=19,color="dodgerblue") + 
         theme_classic()  + labs(x="",y="Species richness") + 
         scale_x_discrete(labels=c("Low fitness","High fitness"))

fig.s4d<-ggplot(Diversity, aes(x=site, y=Density)) + geom_violin(fill="dodgerblue",color="white",alpha=0.5) +
         geom_boxplot(fill="dodgerblue", width = 0.25)+
         #geom_jitter(alpha=0.3,pch=19,fill="dodgerblue") + 
         theme_classic()  +  labs(x="Common garden",y="Plant density") + 
         scale_x_discrete(labels=c("Low fitness","High fitness"))

fig.s4e<-ggplot(Diversity, aes(x=site, y=Evenness)) + geom_violin(fill="dodgerblue",color="white",alpha=0.5) +
         geom_boxplot(fill="dodgerblue",width = 0.25)+
         #geom_jitter(alpha=0.3,pch=19,color="dodgerblue") + 
         theme_classic()  + labs(x="Common garden",y="Community evenness") + 
         scale_x_discrete(labels=c("Low fitness","High fitness"))


colSums(site.6.Matrix, na.rm=TRUE)


ab.6<-as.data.frame(cbind(colSums(site.6.Matrix, na.rm=TRUE),colnames(site.6.Matrix)))
colnames(ab.6)<-c("abundance","species")
ab.6$abundance<-as.numeric(as.character(ab.6$abundance))
ab.6<-ab.6[order(-ab.6$abundance),]
ab.6$rank<-1:nrow(ab.6)
ab.6<-subset(ab.6,abundance>0)
ab.6$color<-ifelse(ab.6$rank<8,"dodgerblue","black")

ab.18<-as.data.frame(cbind(colSums(site.18.Matrix, na.rm=TRUE),colnames(site.18.Matrix)))
colnames(ab.18)<-c("abundance","species")
ab.18$abundance<-as.numeric(as.character(ab.18$abundance))
ab.18<-ab.18[order(-ab.18$abundance),]
ab.18$rank<-1:nrow(ab.18)
ab.18<-subset(ab.18,abundance>0)
ab.18$color<-ifelse(ab.18$rank<8,"dodgerblue","black")


fig.s4a<-ggplot(ab.6, aes(x=rank, y=abundance, fill=color)) + geom_vline(xintercept=7.5, color = "grey", linetype="dotted") +
         geom_line() +  
         geom_line(data=ab.18, aes(x=rank, y=abundance),linetype="dashed") +
         geom_point(data=ab.18, aes(x=rank, y=abundance),pch=24, size=3) +  
         geom_point(pch=25, size=3) + 
         theme_classic()  + labs(x="Rank",y="Total abundance in region") +
         scale_fill_manual(values=c("black","dodgerblue")) +
         scale_y_log10() + 
         theme(legend.position="none")


row2<-plot_grid(fig.s4c,fig.s4d,fig.s4e, align="hv", nrow=2,labels=c("(b)","(c)","(d)"), label_x = -0.1)
plot_grid(fig.s4a,row2, labels=c("(a)",""))

ggsave("plots/Fig S4.pdf",w=10,h=5)

