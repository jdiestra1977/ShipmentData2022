################ Para generar redes donde haya ocurrido n brotes ##########
library("dplyr")
library("ggplot2")
library("plyr")
library("stringr")
library("igraph")
library("tnet")
library("Matrix")
library("grid")
library("gridExtra")
library("cowplot")
library("tidygraph")
library("tidyverse")
library("univariateML")
library("sp")
library("ggnet")
library("network")
library("sna")
library("GGally")
library("intergraph")
library("ggraph")

remove(list = ls())

####### Data ####

TUR0<-readRDS("~/Projects/FMD-Surveillance/ManuscriptPlosCompBio/DataAndFigsCode/gadm36_TUR_0_sp.rds")
TUR1<-readRDS("~/Projects/FMD-Surveillance/ManuscriptPlosCompBio/DataAndFigsCode/gadm36_TUR_1_sp.rds")
load("~/Projects/FMD-Surveillance/ManuscriptPlosCompBio/DataAndFigsCode/temporalData.RData")
load("~/Projects/FMD-Surveillance/ManuscriptPlosCompBio/DataAndFigsCode/todosStatsComplete.RData")
load("ListaConResultados.RData")
load("~/Projects/FMD-Surveillance/ManuscriptPlosCompBio/DataAndFigsCode/OddRatioAllCases.RData")

####### Figure 1 ########

TUR_for<-fortify(TUR0)
epiunitsLocation<-todosStatsComplete %>% filter(Strain=="All") %>% 
  select(Epiunit=Node,Latitude=Latitude,Longitude=Longitud)

mapa<-ggplot(TUR_for) +
  geom_point(data=epiunitsLocation,aes(x=Longitude,y=Latitude),size=0.5,color="purple")+
  geom_polygon( aes(x = long, y = lat, group = group),
                color = "black", fill = "NA") +
  theme_bw()+
  theme(text=element_text(size=20))+
  xlab("")+ ylab("")

tr <- make_tree(21, children = 3, mode = "undirected")
V(tr)$group<-c("I","I",
               "Hr","Hr","Hr","Hr","Hr",
               "Lr","Lr","Lr","Lr","Lr","Lr","Lr","Lr","Lr","Lr","Lr","Lr","Lr","Lr")

tr1<-tr
grafo1<-ggnet2(tr1, color = "group", palette = c("Hr" = "#0072B2", "I" = "#D55E00",
                                                 "Lr"="#009E73"),
               label="group",edge.size=1.5) +
  theme(legend.position = "none",plot.margin=margin(t=0.5,r=0.5,b=0.5,l=0.5,unit="cm"))+
  geom_point(aes(color = color), size = 17, color = "white") +
  geom_point(aes(color = color), size = 17, alpha = 0.5) +
  geom_point(aes(color = color), size = 15) +
  geom_text(aes(label = color), color = "white", fontface = "bold",size=8) +
  guides(color = TRUE)

cowplot::plot_grid(mapa,grafo1,rel_widths = c(2,1),rel_heights = c(2,1.7),labels = "auto",label_size = 25)

####### Figure 2 #######

archivo <-todosStatsComplete %>%
  filter(Strain=="All")
head(archivo)

#data base
grados<-archivo %>% select(Node,InDegree,OutDegree) %>%
  pivot_longer(cols = contains("Degree"),names_to = "Mode",values_to = "Value") %>%
  mutate(Measure="Degree")
grados$Mode<-str_replace_all(grados$Mode, 
                             c("InDegree" = "In", 
                               "OutDegree" = "Out"))
capas<-archivo %>% select(Node,InCoreness,OutCoreness) %>%
  pivot_longer(cols = contains("Coreness"),names_to = "Mode",values_to = "Value")%>%
  mutate(Measure="Coreness")
capas$Mode<-str_replace_all(capas$Mode, 
                            c("InCoreness" = "In", 
                              "OutCoreness" = "Out"))
pesos<-archivo %>% select(Node,InStrength,OutStrength) %>%
  pivot_longer(cols = contains("Strength"),names_to = "Mode",values_to = "Value")%>%
  mutate(Measure="Strength")
pesos$Mode<-str_replace_all(pesos$Mode, 
                            c("InStrength" = "In", 
                              "OutStrength" = "Out"))

statsOfNodesCompleteNet<-rbind(grados,capas,pesos)

#Frequency k-shell (c)
g1f<-ggplot(subset(statsOfNodesCompleteNet,Measure=="Coreness"),
            aes(Value,color=factor(Mode),fill=factor(Mode))) + 
  theme_bw()+
  geom_histogram(position = "identity",alpha=0.6,bins=50) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.title = element_blank()) +
  xlab(expression(italic("k-shell"))) + ylab("Frequency") +
  theme(legend.position = c(0.85,0.75)) +
  theme(text=element_text(size=20),plot.margin = margin(10, 15, 10, 10))

#Average strength versus degree (b)

pesosVersusGrados<-rbind(archivo %>% select(Degree="InDegree",Strength="InStrength") %>% mutate(Mode="In"),
                         archivo %>% select(Degree="OutDegree",Strength="OutStrength") %>% mutate(Mode="Out")) %>%
  group_by(Degree,Mode) %>%summarize_each(funs(Mean=mean,SD=sd))

pesosGrados1<-ggplot(subset(pesosVersusGrados,Degree>0),aes(x=Degree,y=Mean,color=Mode,fill=Mode)) +
  geom_point(size=4,shape=1)+
  theme_bw()+
  facet_wrap(~Mode,ncol=1)+
  geom_abline(slope = 1,intercept = 0,linetype=2) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  #  geom_smooth(method="lm",color="purple",se=F)+
  scale_colour_brewer(palette = "Set1") +
  xlab("Epiunit's degree") + ylab("Average strength of epiunits by degree")+
  #  xlab(expression(italic(k))) + ylab(expression(italic(bar(s)(k))))+
  theme(text = element_text(size=20),legend.position = "none") 

#cummulative degree (a)

degreeIn<-subset(statsOfNodesCompleteNet,Measure=="Degree" & Mode=="In")
degreeOut<-subset(statsOfNodesCompleteNet,Measure=="Degree" & Mode=="Out")
mean(degreeIn$Value)
# Let's count the frequencies of each degree
G.degree.histogram.In <- as.data.frame(table(degreeIn$Value))
G.degree.histogram.Out <- as.data.frame(table(degreeOut$Value))

G.degree.histogram.In$Var1<-as.numeric(as.character(G.degree.histogram.In$Var1))
G.degree.histogram.Out$Var1<-as.numeric(as.character(G.degree.histogram.Out$Var1))

G.degree.histogram.In$CumSum<-1-cumsum(G.degree.histogram.In$Freq)/sum(G.degree.histogram.In$Freq)
G.degree.histogram.In$LogCumSum<-log(G.degree.histogram.In$CumSum)
G.degree.histogram.In$Mode<-"In"
G.degree.histogram.In[G.degree.histogram.In==-Inf]<-0
G.degree.histogram.Out$CumSum<-1-cumsum(G.degree.histogram.Out$Freq)/sum(G.degree.histogram.Out$Freq)
G.degree.histogram.Out$LogCumSum<-log(G.degree.histogram.Out$CumSum)
G.degree.histogram.Out$Mode<-"Out"
G.degree.histogram.Out[G.degree.histogram.Out==-Inf]<-0

G.degree.Both<-rbind(G.degree.histogram.In,G.degree.histogram.Out)
lm_fit <- lm(CumSum ~ Var1, data=subset(G.degree.histogram.In,Mode="In" & Var1>1500 & Var1<5000))
summary(lm_fit)

main.plot <- 
  ggplot()+
  geom_point(data=subset(subset(G.degree.histogram.Out,CumSum>0),Var1>0),
             aes(x=Var1,y=CumSum),color="#377EB8",size=4,shape=1)+
  geom_point(data=subset(subset(G.degree.histogram.In,CumSum>0),Var1>0),
             aes(x=Var1,y=CumSum),color="#E41A1C",size=4,shape=1) +
  theme_bw() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_colour_brewer(palette = "Set1") +
  geom_smooth(data=subset(G.degree.histogram.In,Var1>200 & Var1<2500),
              aes(x=Var1,y=CumSum),method = "lm",color="yellow",se=F,linetype=2)+
  geom_smooth(data=subset(G.degree.histogram.Out,Var1>310 & Var1<750),
              aes(x=Var1,y=CumSum),method = "lm",color="yellow",linetype=2) +
  xlab("Epiunit's Degree") + 
  ylab("Cummulative degree distribution" ) +
  theme(text=element_text(size=20),plot.margin = margin(10, 20, 10, 10))

inset.plot <- ggplot(subset(G.degree.histogram.In,CumSum>0),aes(x=Var1,y=CumSum)) +
  geom_point(color="#E41A1C") +
  theme_bw()+
  scale_y_log10(limits=c(1e-5,1e-2),
                #        breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  #  scale_x_log10()+
  geom_smooth(data=subset(G.degree.histogram.In,Var1>2900 & Var1<10000),
              method = "lm",color="black",se=F,linetype=2) +
  geom_point(data=subset(G.degree.histogram.Out,CumSum>0),aes(x=Var1,y=CumSum),color="#377EB8")+
  geom_smooth(data=subset(G.degree.histogram.Out,Var1>1500 & Var1<10000),
              method = "lm",color="black",se=F) +
  xlim(1000,10000) +xlab("")+ylab("")+ 
  theme(text=element_text(size=10))

plot.with.inset <-
  ggdraw() +
  draw_plot(main.plot) +
  draw_plot(inset.plot, x = 0.23, y = .145, width = .5, height = .34)

cowplot::plot_grid(plot.with.inset,pesosGrados1,g1f,labels = c('a','b','c'),label_size = 22,ncol=3)

####### Figure 3 ######

temporalData %>% glimpse()

colors <- c("Serotype 1"="green","Serotype A"="blue","Serotype O"="yellow")
shipmentRange<-data.frame (xmin="2007-01-01", xmax="2012-08-03", ymin=0, ymax=Inf)

temporalData %>%
  ggplot(aes(x=Month,y=freq,fill=Serotype)) +
  geom_rect(data=shipmentRange, aes(xmin=as.Date(xmin), xmax=as.Date(xmax), ymin=ymin, ymax=ymax), 
            fill="gray", alpha=0.7, inherit.aes = FALSE)+
  geom_bar(stat="identity",position = "identity",alpha=0.6)+
  theme_bw()+
  xlab("") + ylab("Infected farms") +
  scale_fill_manual(values = colors) +
  theme(legend.title = element_blank(),legend.position = c(0.1,0.6),
        text=element_text(size=20))
  
# mapas con serotypes (b)

todosStatsComplete %>% head()

TUR1_for<-fortify(TUR1)

mapaAllOutbreaks<-ggplot(TUR_for) +
  theme_bw()+
  geom_polygon( aes(x = long, y = lat, group = group),
                color = "gray", fill = "NA") +
  geom_polygon( data=TUR1_for,aes(x = long, y = lat, group = group),
                color = "black", fill = "NA") +
  geom_point(data=subset(todosStatsComplete,Condicion=="Infected" & Strain!="All"),
             aes(x=Longitud,y=Latitude),size=1,color="red") + 
  theme(legend.position = "none",text = element_text(size=20)) +
  facet_wrap(~Strain) +
  xlab("") + ylab("")

#Fraction of farms (c)

todosStatsCompleteWithoutAll<-
  todosStatsComplete %>% filter(Strain!="All") %>%
  select(Strain,Condicion) %>%
  group_by(Strain,Condicion) %>% count() %>%
  mutate(Fractions=freq/49580)
todosStatsCompleteWithoutAll$CondNueva<-str_replace_all(todosStatsCompleteWithoutAll$Condicion, 
                                                        c("At risk"="High risk", 
                                                          "Safe"="Low risk",
                                                          "Infected"="Infected"))
todosStatsCompleteWithoutAll$Strain<-str_replace_all(todosStatsCompleteWithoutAll$Strain, 
                                                     c("Outbreak 1"="Serotype 1", 
                                                       "Outbreak A"="Serotype A",
                                                       "Outbreak O"="Serotype O"))

todosStatsCompleteWithoutAll <- todosStatsCompleteWithoutAll %>%
  mutate(Strain = Strain %>% fct_relevel(c("Serotype O","Serotype A","Serotype 1")))
fracFarms<-ggplot(todosStatsCompleteWithoutAll,
                  aes(x=Strain,y=Fractions,fill=CondNueva,
                      color=CondNueva)) +
  theme_bw()+
  scale_fill_manual(breaks = c("Infected", "High risk", "Low risk"), 
                    values=c("#D55E00", "#0072B2", "#009E73")) +
  scale_color_manual(breaks = c("Infected", "High risk", "Low risk"),
                     values=c("#D55E00", "#0072B2", "#009E73")) +
  scale_x_discrete(labels=c("Outbreak 1"="Serotype 1",
                            "Outbreak A"="Serotype A",
                            "Outbreak O"="Serotype O"))+
  geom_point(size=5) +
  theme(legend.position = c(0.8,0.8),legend.title = element_blank(),
        text=element_text(size=17)) +
  xlab("") +ylab("Frac. of farms") +
  ylim(0,1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

abajo<-plot_grid(mapaAllOutbreaks,fracFarms,rel_widths = c(3, 1),labels = c("b","c"),label_size = 25)
plot_grid(infFarmsTime,abajo,ncol=1,rel_heights = c(1,0.67),labels = c("a"," "),label_size = 25)

####### Figure 4 ####
netStatsAll<-todosStatsComplete %>% filter(Strain!="All")
head(netStatsAll)
netStatsAll$Strain<-ifelse(netStatsAll$Strain=="Outbreak 1","Serotype 1",
                           ifelse(netStatsAll$Strain=="Outbreak A","Serotype A","Serotype O"))
netStatsAll$Strain<-
  factor(netStatsAll$Strain, levels = c("Serotype O","Serotype A", "Serotype 1"))

netStatsAll$Condicion1<-ifelse(netStatsAll$Condicion=="Infected","Infected",
                               ifelse(netStatsAll$Condicion=="Safe","Low risk","High risk"))

netStatsAll$Condicion1<-
  factor(netStatsAll$Condicion1, levels = c("Low risk","High risk", "Infected"))

aveEigenTotal<-mean(subset(netStatsAll,Strain=="Serotype O")$EigenNormal,rm.na=T)
aveInDegreeTotal<-mean(subset(netStatsAll,Strain=="Serotype O")$InDegree)
aveInCorenessTotal<-mean(subset(netStatsAll,Strain=="Serotype O")$InCoreness)
aveOutDegreeTotal<-mean(subset(netStatsAll,Strain=="Serotype O")$OutDegree)
aveOutCorenessTotal<-mean(subset(netStatsAll,Strain=="Serotype O")$OutCoreness)
aveBetweenTotal<-mean(subset(netStatsAll,Strain=="Serotype O")$RelativeBet)

uno<-
  ggplot(subset(netStatsAll,EigenNormal>0),
         #       aes(x=InDegree/max(InDegree),y=EigenNormal/max(EigenNormal),color=factor(WithOutbreak))) +
         aes(y=InDegree/aveInDegreeTotal,x=EigenNormal/aveEigenTotal,color=factor(Condicion1))) +  
  geom_point(shape=19,size=3,alpha=0.2) +
  theme_bw()+
  #  geom_rug(sides="b",length = unit(0.05, "npc"))+
  #  annotate("text", label = "Test", size = 4, x = 15, y = 5)+
  facet_grid(Strain~Condicion1) +
  scale_fill_manual(breaks = c("Infected", "High risk", "Low risk"), 
                    values=c("#D55E00", "#0072B2", "#009E73")) +
  scale_color_manual(breaks = c("Infected", "High risk", "Low risk"),
                     values=c("#D55E00", "#0072B2", "#009E73")) +
  geom_vline(xintercept = 1,color="black",linetype=1) +
  geom_hline(yintercept = 1,color="black",linetype=1)+
  #  scale_x_log10(name="eigenvector centrality (scaled)") + #for poster
  #  scale_y_log10(name="In-degree (scaled)",   #for poster
  #                limits=c(1E-2,1E2)) +   #for poster
  scale_x_log10(name=expression(italic("ec (scaled)"))) +
  scale_y_log10(name=expression(paste(italic(k["i"]^{"in"}),italic("(scaled)"))),
                limits=c(1E-2,1E2)) +
  theme(legend.position = "none") +
  theme(text=element_text(size=22),
        axis.text.x = element_text( angle=45,hjust = 1))

dos <- ggplot(netStatsAll, 
              aes(x = InDegree/aveInDegreeTotal, fill = Condicion1)) + 
  geom_density(alpha = 0.4) + theme_void() + 
  facet_wrap(~Strain,ncol=1,scales = "free_y") +
  geom_vline(xintercept = 1)+
  scale_x_log10(limits=c(1E-2,1E2))+
  scale_fill_manual(breaks = c("Infected", "High risk", "Low risk"), 
                    values=c("#D55E00", "#0072B2", "#009E73")) +
  theme(legend.position = "none",plot.margin = margin(29,0,87, 0),
        text=element_text(size=22),strip.text = element_blank(),strip.background = element_blank()) + 
  coord_flip()

unoDos<-cowplot::plot_grid(uno,dos,rel_widths = c(6,1))

tres<-ggplot(netStatsAll,
             aes(x=EigenNormal/aveEigenTotal,fill=factor(Condicion1))) +
  geom_density(alpha=0.4)+
  theme_bw()+
  ylab("")+
  facet_wrap(~Strain,ncol=3,scales = "free_y") +
  #  scale_x_log10(name="eigenvector centrality (scaled)")+ #Poster
  scale_x_log10(name=expression(italic("ec (scaled)")))+
  scale_fill_manual(breaks = c("Infected", "High risk", "Low risk"), 
                    values=c("#D55E00", "#0072B2", "#009E73")) +
  geom_vline(xintercept = 1) +
  theme(legend.position = "none",
        text=element_text(size=22),plot.margin = margin(t=0,b=0,l=-2,r=17))

cuatro<-
  ggplot(subset(netStatsAll,EigenNormal>0),
         aes(y=InCoreness/aveInCorenessTotal,
             x=EigenNormal/aveEigenTotal,color=factor(Condicion1))) +  
  geom_point(shape=19,size=3,alpha=0.2) +
  theme_bw()+
  facet_grid(Strain~Condicion1) +
  scale_fill_manual(breaks = c("Infected", "High risk", "Low risk"), 
                    values=c("#D55E00", "#0072B2", "#009E73")) +
  scale_color_manual(breaks = c("Infected", "High risk", "Low risk"),
                     values=c("#D55E00", "#0072B2", "#009E73")) +
  geom_vline(xintercept = 1,color="black",linetype=1) +
  geom_hline(yintercept = 1,color="black",linetype=1)+
  scale_x_log10(name=expression(italic("ec (scaled)"))) +
  scale_y_log10(name=expression(italic(paste(italic("kC"["i"]^{"in"}),italic("(scaled)"))))) +
  theme(legend.position = "none") +
  theme(text=element_text(size=22),
        axis.text.x = element_text( angle=45,hjust = 1))

cinco <- 
  ggplot(netStatsAll,
         aes(x = InCoreness/aveInCorenessTotal,fill = Condicion1)) + 
  geom_density(alpha = 0.4) + theme_void() + 
  facet_wrap(~Strain,ncol=1,scales = "free_y") +
  geom_vline(xintercept = 1)+
  scale_x_log10()+
  scale_fill_manual(breaks = c("Infected", "High risk", "Low risk"), 
                    values=c("#D55E00", "#0072B2", "#009E73")) +
  theme(legend.position = "none",plot.margin = margin(29,0,84, 0),
        text=element_text(size=22),strip.text = element_blank()
        ,strip.background = element_blank()) + 
  coord_flip()

cuatroCinco<-cowplot::plot_grid(cuatro,cinco,rel_widths = c(6,1))

arriba<-cowplot::plot_grid(unoDos,cuatroCinco,ncol=2)

cowplot::plot_grid(arriba,tres,ncol=1,rel_heights = c(3,1),labels = "auto",label_size = 25)

####### Figure 5 #########

OddRatioAllCases %>% glimpse()

OddRatioAllCases$Region<-factor(OddRatioAllCases$Region,
                                levels=c("Complete","Endemic","Epidemic","Epi-partial"))
colorBlind<-c("#F0E442", "#0072B2", "#D55E00", "#009E73")

dat_text1 <- data.frame(
  label = c("***", "***","",""),
  Source   = c("Infected","Low-risk","Infected","Low-risk"),
  x=c(2,1,0,0),
  y=c(2.1,0.4,0,0),
  Model=c("Univariate","Univariate","Multivariate","Multivariate")
)
dat_text1$Model<-factor(dat_text1$Model,levels = c("Univariate","Multivariate"))

labels1 <- data.frame(
  label = c("a", "b","c","d"),
  Source   = c("Infected","Low-risk","Infected","Low-risk"),
  x=c(0.3,0.3,0.3,0.3),
  y=c(8,8,8,8),
  Model=c("Univariate","Univariate","Multivariate","Multivariate")
)
labels1$Model<-factor(labels1$Model,levels = c("Univariate","Multivariate"))

ggplot(OddRatioAllCases,aes(x=Variable,y=OR,color=Region,shape=Model))+
  theme_bw()+
  geom_point(position=position_dodge(width = .6),size=5)+
  geom_errorbar(aes(ymin=LowConf,ymax=HighConf),position=position_dodge(width = .6),size=1) +
  geom_text(data=dat_text1,aes(x=x,y=y,label=label),size=12,color="purple")+
  geom_text(data=labels1,aes(x=x,y=y,label=label),size=10,color="black",fontface="bold")+
  geom_hline(yintercept = 1,lty=2)+
  xlab("")+ ylab("Odds ratio")+
  #  annotate(geom="text", x=-2, y=3, label="a",size=13) +
  facet_grid(Model~Source)+
  scale_color_manual(values = colorBlind)+
  theme(legend.title = element_blank(),legend.position = c(0.9,0.85),
        text=element_text(size=25),axis.text.x = element_text(angle=45,hjust = 1))

################ Global functions #########################

#Function to create directed weighted networks. Need a file with 3 columns"
#Source, Target and Weight
networkCreationDW<-function(x){
  el1 <-as.data.frame(x) #Verify that it's a data.frame
#  el1<-subset(el1,source_epiunit_id!=destination_epiunit_id) #Avoid loops, self-connections
#  el1<-subset(el1,i!=j) #Avoid loops, self-connections
#  el1<-subset(el1,Source!=Target) #Avoid loops, self-connections
  ell3<-as.matrix(el1)
  ell3[,1]=as.character(ell3[,1]) #Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
  ell3[,2]=as.character(ell3[,2])
  ggg3=graph.edgelist(ell3[,1:2],directed = TRUE) #Create the directed network
  E(ggg3)$weight=as.numeric(ell3[,3]) #weights using the number of cattle made
#  E(ggg3)$weight=as.numeric(ell3[,4]) #weights using the number of shipments made
  ggg3<-simplify(ggg3) #Verify that there are not multiple connection or remaining loops
  compsG<-clusters(ggg3) #Check if the network is broken
  g<-induced.subgraph(ggg3,compsG$membership==which.max(compsG$csize)) #Only select the largest CC
  #If the network is not broken, keeps the networks as is.
}

#Function to given a network in igraph format
#creates a file with some measures
netsRealStats<-function(x){
  pi<-transitivity(x,type="weighted",weights = E(x)$weight)
  is.na(pi)<-sapply(pi,is.infinite)
  reci=paste(vcount(x),ecount(x),mean(pi,na.rm=T),
             eigen_centrality(x,directed = T,weights = E(x)$weight)$value,
             reciprocity(x),
             mean_distance(x, directed = TRUE),
             assortativity(x, types1 = graph.strength(x), directed = T),
             diameter(x,directed = TRUE,weights =NULL),
             sep=" ")
  statsGlobalReal<-reci
  #  write.table(reci,file="~/Projects/FMD-Surveillance/raw-data/FrequencyRelated/RealCompleteNetworkFrequency.csv",append = TRUE,
  #              sep = ",",col.names=FALSE,row.names = F)
  statsGlobalReal
}

#This function creates random networks according to the tutorial in
#https://toreopsahl.com/tnet/weighted-networks/random-networks/
#where, links i1 -> j1  and i2 -> j2 are rewired as i1 -> j2 and i2 -> j1,
# keeping the weight from the outgoing link.
randomNetStatistics <- function(listaEnlaces){
  #Select source, destination and the number of shipments made
  prueba<-edgeListCompleteNetwork %>% select(source_epiunit_id,
                                             destination_epiunit_id,Number)
  #  prueba<-listaEnlaces %>% select(source_epiunit_id,destination_epiunit_id,Number)
  #Delete loops, where source and destination are equal
  prueba<-subset(prueba,source_epiunit_id!=destination_epiunit_id)
  #create a random network using algorithm explained above and conver to igraph object
  prueba1<-rg_reshuffling_w(prueba,option = "links",directed = T) 
  esto<-rg_reshuffling_w(prueba1,option = "weights",directed = T) %>%
    tnet_igraph(.)
  #calculate network measures
  pi<-transitivity(esto,type="weighted",weights = E(esto)$weight)
  is.na(pi)<-sapply(pi,is.infinite)
  return(c(reciprocity(esto),
           eigen_centrality(esto,directed = T,weights = E(esto)$weight)$value,
           assortativity(esto, types1 = graph.strength(esto), directed = T),
           mean(pi,na.rm=T),
           mean_distance(esto, directed = TRUE),
           diameter(esto,directed = TRUE)))
}

# These are for individual characterization directed-weighted network (all but clustering)
nodeStatsDW<-function(x){
  V(x)$Eig<-eigen_centrality(x,directed = T,weights = E(x)$weight)$vector
  V(x)$InDegree<-degree(x,mode = "in")
  V(x)$OutDegree<-degree(x,mode = "out")
  V(x)$InStrength<-strength(x,mode="in",weights = E(x)$weight)
  V(x)$OutStrength<-strength(x,mode="out",weights = E(x)$weight)
  V(x)$Betweenness<-betweenness(x, directed=T, weights=E(x)$weight,nobigint = TRUE, normalized = FALSE)
  V(x)$InCoreness<-coreness(x,mode="in")
  V(x)$OutCoreness<-coreness(x,mode="out")
  y<-as.undirected(x)
  V(x)$ClusterW<-transitivity(y,type="weighted",weights = E(y)$weight)
  centrality <- data.frame(Node = as.numeric(V(x)),
                           OutDegree    = V(x)$OutDegree,
                           InDegree    = V(x)$InDegree,
                           OutStrength = V(x)$OutStrength,
                           InStrength   = V(x)$InStrength,
                           OutCoreness = V(x)$OutCoreness,
                           InCoreness = V(x)$InCoreness,
                           ClusterW=V(x)$ClusterW,
                           Betweenness  = V(x)$Betweenness,
                           EigenCent = V(x)$Eig)
  centrality
}





