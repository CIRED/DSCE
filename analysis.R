# Draw graphs for Fillon, Guivarch, Taconet, 2023.

#setwd 
#setwd("C:/Users/")
rm(list=ls())
library(zoo)
library(akima)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(plotly)
library(scales)
library(cowplot)
library(gcookbook) 
library(readxl)
library(directlabels)
library(cowplot)
library(gcookbook) 
library(readxl)
library(directlabels)
library(tikzDevice)

matrix_add <- read_excel("run/SCC_matrix_runs_add_3105b.xlsx")
matrix_rs <- read_excel("run/SCC_matrix_runs_rs_3105b.xlsx")

colnames(matrix_add)<-c("id","run_no","rho1y","theta","rra","J","broken","tp","deterministic","SCC0","finalT","threshold_max","damage")
colnames(matrix_rs)<-c("id","run_no","rho1y","theta","epsilon","J","broken_b","tp","deterministic","SCC0_rs","finalT_b","threshold_max","damage")

#Revisions JEEM : J < 10%
#matrix_add=subset(matrix_add,matrix_add$deterministic==0 & matrix_add$SCC0!=0 & matrix_add$J<=0.1 & matrix_add$rho1y==0.015 & matrix_add$theta==1.5 & matrix_add$damage!=0)
matrix_add=subset(matrix_add,matrix_add$deterministic==0 & matrix_add$SCC0!=0 & matrix_add$J<=0.2 & matrix_add$rho1y==0.015 & matrix_add$theta==1.5 & matrix_add$threshold_max==5.7 & matrix_add$damage!=0)
matrix_add=subset(matrix_add, select=c("theta","J","SCC0"))

#Revisions JEEM : J < 10%
matrix_rs=subset(matrix_rs,matrix_rs$deterministic==0 & matrix_rs$SCC0_rs!=0 & matrix_rs$theta==1.5 & matrix_rs$J<=0.2 & matrix_rs$threshold_max==5.7 & matrix_rs$damage!=0)
#matrix_rs=subset(matrix_rs,matrix_rs$deterministic==0 & matrix_rs$SCC0_rs!=0 & matrix_rs$theta==1.5 & matrix_rs$J<=0.1 & matrix_rs$damage!=0)
matrix_rs$SWF='risk-sensitive'
matrix_add=matrix_add
matrix_add$epsilon=0
matrix_add$SWF='additive'
matrix_rs<-matrix_rs[ which(  matrix_rs$epsilon==0.001 | matrix_rs$epsilon == 0.133 | matrix_rs$epsilon == 0.3), ]
matrix_rs=subset(matrix_rs, select=c("theta","J","SCC0_rs","epsilon","SWF"))
colnames(matrix_rs)=c("theta","J","SCC0","epsilon","SWF")

Matrix_comp2=rbind(matrix_add,matrix_rs)
colnames(matrix_rs)=c("theta","J","SCC0_rs","epsilon_b","SWF")
Matrix_comp=merge(matrix_add,matrix_rs, by=c("theta","J"))
Matrix_comp=as.data.frame(Matrix_comp)
Matrix_comp<-Matrix_comp %>%
  mutate(ratio= SCC0_rs/SCC0)
Matrix_comp<-Matrix_comp[ which( Matrix_comp$epsilon_b==0.001 | Matrix_comp$epsilon_b == 0.133 | Matrix_comp$epsilon_b == 0.3), ]
Matrix_comp<-Matrix_comp[ which( Matrix_comp$J<=0.1), ]
Matrix_comp2<-Matrix_comp2[ which( Matrix_comp2$J<=0.1), ]

plot1<-ggplot(Matrix_comp, aes(x=J*100, y=ratio, group=epsilon_b, linetype=factor(epsilon_b))) +
  xlab("Irreversible shock to the output J (in %)") +
  geom_line(color='black')+
  ylab("Ratio of risk-sensitive to additive SCC") +
  scale_linetype_manual("Temporal risk aversion \u03B5", values = c("0.001" = 8, "0.133" = 6, "0.3" = 3))+
  scale_y_discrete(limits=c(1, 1.2,1.4,1.6, 1.8, 2))+
  theme_bw()+
  theme(legend.position="none")

plot2<-ggplot(Matrix_comp2, aes(x=J*100, y=SCC0, group=epsilon, linetype=factor(epsilon))) +
  xlab("Irreversible shock to the output J (in %)") +
  geom_line()+
  ylab("Absolute value of the SCC (in $/tC)") +
  scale_y_discrete(limits=c(75, 150, 225, 300, 375))+
  scale_linetype_manual("Temporal risk aversion \u03B5", values = c("0" = 1, "0.001" = 8, "0.133" = 6, "0.3" = 3))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12))

library(patchwork)
a=plot2+plot1
a

###EQUIVALENCE IN J / rho (FIGURE 2)
matrix_crra <- read_excel("run/SCC_matrix_runs_add_3105b_revisions.xlsx")
matrix_rs <- read_excel("run/SCC_matrix_runs_rs_3105b.xlsx")
colnames(matrix_crra)<-c("id","run_no","rho1y","theta","rra","J","broken","tp","deterministic","SCC0","finalT","threshold_max","damage")
colnames(matrix_rs)<-c("id","run_no","rho1y","theta","epsilon","J","broken_b","tp","deterministic","SCC0_rs","finalT_b","threshold_max","damage")

matrix_crra=subset(matrix_crra,matrix_crra$deterministic==0 & matrix_crra$SCC0!=0 & matrix_crra$rho1y==0.015 & matrix_crra$damage!=0 & matrix_crra$theta==1.5 & matrix_crra$threshold_max==5.7 & matrix_crra$J<=0.44)
matrix_rs=subset(matrix_rs,matrix_rs$deterministic==0 & matrix_rs$SCC0_rs!=0 & matrix_rs$damage!=0 & matrix_rs$theta==1.5 & matrix_rs$epsilon==0.133 & matrix_rs$threshold_max==5.7 & matrix_rs$J<=0.2)
matrix_rs=subset(matrix_rs, select=c("theta","J","SCC0_rs"))
colnames(matrix_rs)<-c("theta","J_rs","SCC0_rs")

#Merge based on closest values
ans15 <- vapply(matrix_crra$SCC0, function(x) x-matrix_rs$SCC0_rs, numeric(10))
indx15 <- apply(abs(ans15), 2, which.min)
nearest15<-cbind(matrix_crra, matrix_rs[indx15, ][-2])
nearest15<-subset(nearest15,select=c("theta","J","SCC0","SCC0_rs"))
nearest15<-merge(nearest15,matrix_rs, by= c("theta","SCC0_rs"))
#drop duplicated data
nearest<-nearest15[which(duplicated(nearest15$SCC0_rs)!="TRUE"),]

plotJ<-ggplot(nearest, aes(x=J_rs*100, y=J*100)) +
  geom_line(size=1, colour="black") +
  xlab("Irreversible shock J (%) under risk-sensitive preferences") +
  ylab("Irreversible shock J (%) under additive preferences") +
  labs(color = "\u03B8") +
  geom_vline(xintercept = 10, linetype="dashed", colour="black")+
  geom_hline(yintercept = 14, linetype="dashed",colour="black")+
  geom_vline(xintercept = 20, linetype="dashed", colour="black")+
  geom_hline(yintercept = 40, linetype="dashed", colour="black")+
  geom_segment(aes(x = 0, xend = 20, y = 0, yend = 20),
               colour = "black",size=0.4, linetype=3)+
  scale_y_continuous(breaks=c(0,10,14,20,30,40))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12))

SCC_matrix_runs_add <- read_excel("run/SCC_matrix_runs_add_3105b.xlsx")
colnames(SCC_matrix_runs_add)<-c("id","run_no","rho1y","theta","rra","J","broken","tp","deterministic","SCC0","finalT","threshold_max","damage")
Matrix_CRRA=subset(SCC_matrix_runs_add,SCC_matrix_runs_add$rho1y!=0.015 & SCC_matrix_runs_add$J != 0)
Matrix_CRRA=subset(Matrix_CRRA, select=c("theta","J","rho1y"))
origin=c(1.5, 0, 0.015)
Matrix_CRRA=rbind(Matrix_CRRA, origin)

plotrho<-ggplot(Matrix_CRRA, aes(x=J*100, y=rho1y*100)) +
  xlab("Irreversible shock to the output J (%)") +
  geom_line(size=1)+
  ylab("Pure time preference \u03C1 (% yearly)") +
  geom_hline(yintercept = 0.1, linetype="dotted")+
  geom_hline(yintercept = 1.5, linetype="dotted")+
  scale_y_continuous(breaks=c(0.1, 0.5, 1, 1.5), labels=c("Stern = 0.1", "0.5", "1", "Nordhaus = 1.5")) +
  theme_bw()+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12))
library(patchwork)
plotparam=plotrho+plotJ

#Graph decomposition (FIGURE 3 & Annex 3)
DWI_rs_1_1 <- read.csv("outputs/runs_revisions9rs_run0026-2023-03-08-09h16/complete_DWI_rs.csv")
MHE_rs_1_1 <- read.csv("outputs/runs_revisions9rs_run0026-2023-03-08-09h16/complete_MHE_rs.csv")
DWI_rs_10_1 <- read.csv("outputs/runs_revisions9rs_run0027-2023-03-08-15h27/complete_DWI_rs.csv")
MHE_rs_10_1 <- read.csv("outputs/runs_revisions9rs_run0027-2023-03-08-15h27/complete_MHE_rs.csv")
DWI_rs_1_2 <- read.csv("outputs/runs_revisions9rs_run0024-2023-03-06-20h53/complete_DWI_rs.csv")
MHE_rs_1_2 <- read.csv("outputs/runs_revisions9rs_run0024-2023-03-06-20h53/complete_MHE_rs.csv")
DWI_rs_10_2 <- read.csv("outputs/runs_revisions9rs_run0025-2023-03-06-11h39/complete_DWI_rs.csv")
MHE_rs_10_2 <- read.csv("outputs/runs_revisions9rs_run0025-2023-03-06-11h39/complete_MHE_rs.csv")
DWI_rs_1_3 <- read.csv("outputs/runs_revisions9rs_run0028-2023-03-08-22h55/complete_DWI_rs.csv")
MHE_rs_1_3 <- read.csv("outputs/runs_revisions9rs_run0028-2023-03-08-22h55/complete_MHE_rs.csv")
DWI_rs_10_3 <- read.csv("outputs/runs_revisions9rs_run0029-2023-03-09-08h55/complete_DWI_rs.csv")
MHE_rs_10_3 <- read.csv("outputs/runs_revisions9rs_run0029-2023-03-09-08h55/complete_MHE_rs.csv")

matrix_decomp <- rbind(DWI_rs_1_1,MHE_rs_1_1,DWI_rs_10_1,MHE_rs_10_1,DWI_rs_1_2,MHE_rs_1_2,DWI_rs_10_2,MHE_rs_10_2,DWI_rs_1_3,MHE_rs_1_3,DWI_rs_10_3,MHE_rs_10_3)

matrix_decomp$channel=rep(1:2, times = 6)
matrix_decomp$id=rep(c("0.001","0.133","0.3"), each = 4)
matrix_decomp$J=rep(c(1,10),each=2,times=3)
colnames(matrix_decomp)=c("value","channel","id","J")
matrix_decomp$value=abs(matrix_decomp$value)

matrix_decomp_mhe=subset(matrix_decomp, channel == 2) 

histo1<-ggplot(matrix_decomp_mhe, aes(x = id, y = value,fill=factor(J))) +
  theme(legend.position="none")+
  scale_fill_manual("legend", values = c("1" = "grey","10" = "black"), name="New")+
  geom_col(position = "dodge") +  
  ylab("MHE contribution to SCC (in $/tC)") +
  xlab("Temporal risk aversion \u03B5") 

matrix_decomp_dwi=subset(matrix_decomp, channel == 1) 

histo2<-ggplot(matrix_decomp_dwi, aes(x = id, y = value,fill=factor(J))) +
  scale_fill_manual("Shock J (in %)", values = c("1" = "grey","10" = "black"))+
  ylab("DWI contribution to SCC (in $/tC)") +
  xlab("Temporal risk aversion \u03B5") +
  geom_col(position = "dodge") 

library(patchwork)
histo=histo1+histo2


#Graph decomposition (FIGURE 3 & Annex 3)
DWI_rs_1_1 <- read.csv("outputs/runs_revisions9rs_run0026-2023-03-08-09h16/dwi_immediate.csv")
DWI_rs_1_1=mean(DWI_rs_1_1[,1])
MHE_rs_1_1 <- read.csv("outputs/runs_revisions9rs_run0026-2023-03-08-09h16/mhe_immediate.csv")
MHE_rs_1_1=mean(MHE_rs_1_1[,1])
DWI_rs_10_1 <- read.csv("outputs/runs_revisions9rs_run0027-2023-03-08-15h27/dwi_immediate.csv")
DWI_rs_10_1=mean(DWI_rs_10_1[,1])
MHE_rs_10_1 <- read.csv("outputs/runs_revisions9rs_run0027-2023-03-08-15h27/mhe_immediate.csv")
MHE_rs_10_1=mean(MHE_rs_10_1[,1])
DWI_rs_1_2 <- read.csv("outputs/runs_revisions9rs_run0024-2023-03-06-20h53/dwi_immediate.csv")
DWI_rs_1_2=mean(DWI_rs_1_2[,1])
MHE_rs_1_2 <- read.csv("outputs/runs_revisions9rs_run0024-2023-03-06-20h53/mhe_immediate.csv")
MHE_rs_1_2=mean(MHE_rs_1_2[,1])
DWI_rs_10_2 <- read.csv("outputs/runs_revisions9rs_run0025-2023-03-06-11h39/dwi_immediate.csv")
DWI_rs_10_2=mean(DWI_rs_10_2[,1])
MHE_rs_10_2 <- read.csv("outputs/runs_revisions9rs_run0025-2023-03-06-11h39/mhe_immediate.csv")
MHE_rs_10_2=mean(MHE_rs_10_2[,1])
DWI_rs_1_3 <- read.csv("outputs/runs_revisions9rs_run0028-2023-03-08-22h55/dwi_immediate.csv")
DWI_rs_1_3=mean(DWI_rs_1_3[,1])
MHE_rs_1_3 <- read.csv("outputs/runs_revisions9rs_run0028-2023-03-08-22h55/mhe_immediate.csv")
MHE_rs_1_3=mean(MHE_rs_1_3[,1])
DWI_rs_10_3 <- read.csv("outputs/runs_revisions9rs_run0029-2023-03-09-08h55/dwi_immediate.csv")
DWI_rs_10_3=mean(DWI_rs_10_3[,1])
MHE_rs_10_3 <- read.csv("outputs/runs_revisions9rs_run0029-2023-03-09-08h55/mhe_immediate.csv")
MHE_rs_10_3=mean(MHE_rs_10_3[,1])

matrix_decomp <- rbind(DWI_rs_1_1,MHE_rs_1_1,DWI_rs_10_1,MHE_rs_10_1,DWI_rs_1_2,MHE_rs_1_2,DWI_rs_10_2,MHE_rs_10_2,DWI_rs_1_3,MHE_rs_1_3,DWI_rs_10_3,MHE_rs_10_3)
#matrix_decomp <- read_excel("Decomp_matrix_runs_rs_3105b.xlsx")
matrix_decomp=as.data.frame(matrix_decomp)

matrix_decomp$channel=rep(1:2, times = 6)
matrix_decomp$id=rep(c("0.001","0.133","0.3"), each = 4)
matrix_decomp$J=rep(c(1,10),each=2,times=3)
colnames(matrix_decomp)=c("value","channel","id","J")
matrix_decomp$value=abs(matrix_decomp$value)

matrix_decomp_mhe=subset(matrix_decomp, channel == 2) 

histo3<-ggplot(matrix_decomp_mhe, aes(x = id, y = value,fill=factor(J))) +
  theme(legend.position="none")+
  scale_fill_manual("legend", values = c("1" = "grey","10" = "black"), name="New")+
  geom_col(position = "dodge") +  
  ylab("mhe contribution to SCC (in $/tC)") +
  xlab("Temporal risk aversion \u03B5") 

matrix_decomp_dwi=subset(matrix_decomp, channel == 1) 

histo4<-ggplot(matrix_decomp_dwi, aes(x = id, y = value,fill=factor(J))) +
  scale_fill_manual("Shock J (in %)", values = c("1" = "grey","10" = "black"))+
  ylab("dwi contribution to SCC (in $/tC)") +
  xlab("Temporal risk aversion \u03B5") +
  geom_col(position = "dodge") 

library(patchwork)
histob=histo3+histo4


#####Temperature intervals (Annex 1)
Matrix_add <- read_excel("run/SCC_matrix_runs_add_3105b.xlsx")
Matrix_rs <- read_excel("run/SCC_matrix_runs_rs_3105b.xlsx")
colnames(Matrix_add)<-c("id","run_no","rho1y","theta","rra","J","broken","tp","deterministic","SCC0","finalT","threshold_max","damage")
colnames(Matrix_rs)<-c("id","run_no","rho1y","theta","epsilon","J","broken_b","tp","deterministic","SCC0_rs","finalT_b","threshold_max","damage")

Matrix_add=subset(Matrix_add,Matrix_add$deterministic==0 & Matrix_add$rho1y==0.015 & Matrix_add$SCC0!=0 & Matrix_add$theta==1.5 & Matrix_add$J==0.2)
Matrix_add=subset(Matrix_add, select=c("theta","J","SCC0","threshold_max"))
Matrix_rs=subset(Matrix_rs,Matrix_rs$deterministic==0 & Matrix_rs$theta==1.5 & Matrix_rs$epsilon==0.133 & Matrix_rs$J==0.2)
Matrix_comp=merge(Matrix_add,Matrix_rs, by=c("theta","J","threshold_max"))
Matrix_comp=as.data.frame(Matrix_comp)
Matrix_comp<-Matrix_comp %>%
  mutate(ratio= SCC0_rs/SCC0)

#plot
ggplot(Matrix_comp, aes(x=threshold_max, y=ratio))+
  xlab("Higher temperature interval (°C)") +
  geom_line(size=1.2,color = "black")+
  ylab("Ratio of risk-sensitive to additive SCC") +
  scale_linetype_manual("Temporal risk aversion", values = c("0.001" = 1, "0.133" = 6, "0.3" = 3))+
  theme_bw()+
 # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12))


##### Sensitivity elasmu (Annex 2)
Matrix_add <- read_excel("run/SCC_matrix_runs_add_3105b.xlsx")
Matrix_rs <- read_excel("run/SCC_matrix_runs_rs_3105b.xlsx")
colnames(Matrix_add)<-c("id","run_no","rho1y","theta","rra","J","broken","tp","deterministic","SCC0","finalT","threshold_max","damage")
colnames(Matrix_rs)<-c("id","run_no","rho1y","theta","epsilon","J","broken_b","tp","deterministic","SCC0_rs","finalT_b","threshold_max","damage")

Matrix_add=subset(Matrix_add,Matrix_add$deterministic==0 & Matrix_add$SCC0!=0 & Matrix_add$rho1y==0.015 & Matrix_add$threshold_max==5.7)
Matrix_add=subset(Matrix_add, select=c("theta","J","SCC0"))
Matrix_rs=subset(Matrix_rs,Matrix_rs$deterministic==0 & Matrix_rs$epsilon==0.133 & Matrix_rs$threshold_max==5.7)
Matrix_comp2=merge(Matrix_add,Matrix_rs, by=c("theta","J"))
Matrix_comp2=subset(Matrix_comp2,Matrix_comp2$tp==1 & Matrix_comp2$epsilon==0.133)

Matrix_comp2<-Matrix_comp2[ which( Matrix_comp2$theta==0.5 | Matrix_comp2$theta==1.5  |Matrix_comp2$theta==2.5  ), ]
Matrix_comp2<-Matrix_comp2[ which( Matrix_comp2$J>=0.01), ]

Matrix_comp2<-as.data.frame(Matrix_comp2)
Matrix_comp2<-Matrix_comp2 %>%
  mutate(ratio= log(SCC0_rs/SCC0))

#plot
ggplot(Matrix_comp2, aes(x=J*100, y=ratio, group=theta, linetype=factor(theta))) +
  xlab("Irreversible increase J in the damage factor (in %)") +
  geom_line(size=1)+
  ylab("Ratio of risk-sensitive to additive SCC (log scale)") +
  labs(color = "\u03B8") +
  geom_hline(yintercept = 0, linetype="dotted")+
  scale_linetype_discrete("Inequality aversion \u03B7")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12))


#Revisions JEEM
#SCCed/SCCstoch for different rra under EZW preferences

matrix_EZW <- read_excel("run/SCC_matrix_runs_revisions7.xlsx")
colnames(matrix_EZW)<-c("id","run_no","rho1y","theta","epsilon_b","J","broken","tp","deterministic","SCC0","finalT","threshold_max","damage","gamma")

matrix_det10=subset(matrix_EZW,matrix_EZW$deterministic==1 & matrix_EZW$tp==1 & matrix_EZW$J==0.1 & matrix_EZW$theta==0.8)
matrix_det10$det=1
matrix_stoch10=subset(matrix_EZW,matrix_EZW$deterministic==0 & matrix_EZW$tp==1 & matrix_EZW$J==0.1 & matrix_EZW$theta==0.8)
matrix_stoch10$det=0

total <- merge(matrix_stoch10, matrix_det10, by=c("theta","gamma","J"))
total$value = (total$SCC0.y/total$SCC0.x)*100
total$value2 = (total$SCC0.y/total$SCC0.x)*100

total=subset(total, select=c("theta","J","value","value2","gamma"))
total_min=subset(total, total$gamma<=2)

plot_rev1<-ggplot(total_min, aes(x=gamma, y=value)) +
  xlab("Relative risk aversion \u03B3") +
  geom_line(color='black')+
  ylab("Share of SCC with expected damages in stochastic SCC (in %)") +
  theme_bw()+
  geom_hline(yintercept=100,linetype="dotted")

plot_rev2<-ggplot(total, aes(x=gamma, y=value2)) +
  xlab("Relative risk aversion \u03B3") +
  geom_line(color='black')+
  ylab("Share of SCC with expected damages in stochastic SCC (in %)") +
  theme_bw()+
  geom_hline(yintercept=100,linetype="dotted")
library(patchwork)
a=plot_rev1+plot_rev2


#SCC under EZW preferences
matrix_EZW=subset(matrix_EZW,matrix_EZW$deterministic==0 & matrix_EZW$J==0.1 & matrix_EZW$theta==0.8)
matrix_EZW$J=100*matrix_EZW$J

plot_SCCEZW<-ggplot(matrix_EZW, aes(x=gamma, y=SCC0)) +
  xlab("Relative risk aversion \u03B3") +
  geom_line(color='black')+
  ylim(0, 600000)+
  ylab("Stochastic SCC (in $/tC)") +
  theme_bw()

plot_logSCCEZW<-ggplot(matrix_EZW, aes(x=gamma, y=log(SCC0))) +
  xlab("Relative risk aversion \u03B3") +
  geom_line(color='black')+
  ylab("Stochastic SCC (log scale)") +
  ylim(4, 14)+
  theme_bw()

#SCC under RS preferences
matrix_RS <- read_excel("run/SCC_matrix_runs_revisions9rs.xlsx")
colnames(matrix_RS)<-c("id","run_no","rho1y","theta","epsilon_b","J","broken","tp","deterministic","SCC0","finalT","threshold_max","damage","gamma")
matrix_RS=subset(matrix_RS,matrix_RS$deterministic==0 & matrix_RS$J!=0)
matrix_RS$J=100*matrix_RS$J

plot_SCCRS<-ggplot(matrix_RS, aes(x=epsilon_b, y=SCC0)) +
  xlab("Temporal risk aversion \u03B5") +
  geom_line(color='black')+
  ylim(0, 600000)+
  ylab("Stochastic SCC (in $/tC)") +
  theme_bw()

plot_logSCCRS<-ggplot(matrix_RS, aes(x=epsilon_b, y=log(SCC0))) +
  xlab("Temporal risk aversion \u03B5") +
  geom_line(color='black')+
  ylab("Stochastic SCC (log scale)") +
  ylim(4, 14)+
  theme_bw()

library(patchwork)
SCC=plot_SCCEZW+plot_SCCRS
SCClog=plot_logSCCEZW+plot_logSCCRS

matrix_rs_ed <- read_excel("run/SCC_matrix_runs_revisions9rs.xlsx")
colnames(matrix_rs_ed)<-c("id","run_no","rho1y","theta","epsilon_b","J","broken","tp","deterministic","SCC0","finalT","threshold_max","damage","gamma")

matrix_rs_eddet=subset(matrix_rs_ed,matrix_rs_ed$deterministic==1 & matrix_rs_ed$tp==1 & matrix_rs_ed$J==0.1 & matrix_rs_ed$damage!=0)
matrix_rs_eddet$det=1
matrix_rs_edstoch=subset(matrix_rs_ed,matrix_rs_ed$deterministic==0 & matrix_rs_ed$tp==1 & matrix_rs_ed$J==0.1 & matrix_rs_ed$damage!=0)
matrix_rs_edstoch$det=0

total2 <- merge(matrix_rs_edstoch, matrix_rs_eddet,by=c("epsilon_b","J"))
total2$value = (total2$SCC0.y/total2$SCC0.x)*100
total2=subset(total2, select=c("epsilon_b","J","value"))
                
plot_rev3<-ggplot(total2, aes(x=epsilon_b, y=value)) +
  xlab("Temporal risk aversion \u03B5") +
  geom_line(color='black')+
  ylab("Share of SCC with expected damages in stochastic SCC (in %)") +
  theme_bw()+
  scale_x_continuous(breaks=c(0,0.3,1,2,3,4))+
  geom_hline(yintercept=50,linetype="dotted")+
  geom_vline(xintercept=0.3,linetype="dotted")

final<- merge(matrix_rs_edstoch, matrix_rs_eddet,by=c("epsilon_b","J"))
final$value = (final$SCC0.x-final$SCC0.y)
matrix_add <- read_excel("run/SCC_matrix_runs_revisions8.xlsx")
colnames(matrix_add)<-c("id","run_no","rho1y","theta","epsilon_b","J","broken","tp","deterministic","SCC0","finalT","threshold_max","damage","gamma")
matrix_add_det10=subset(matrix_add,matrix_add$deterministic==1 & matrix_add$tp==1 & matrix_add$theta==1.5  & matrix_add$J==0.1)
matrix_add_det10$det=1
matrix_add_stoch10=subset(matrix_add,matrix_add$deterministic==0 & matrix_add$tp==1 & matrix_add$theta==1.5 & matrix_add$J==0.1)
matrix_add_stoch10$det=0
diff=matrix_add_stoch10$SCC0-matrix_add_det10$SCC0

final$value=(diff/final$value)*100

final_min=subset(final, final$epsilon_b<=0.3)
#final_min=subset(final, final$epsilon_b<=1)

plot_rev4<-ggplot(final_min, aes(x=epsilon_b, y=value)) +
  xlab("Temporal risk aversion \u03B5") +
  geom_line(color='black')+
  ylab("Share of the additive premium in the risk-sensitive premium (in %)") +
  theme_bw()

library(patchwork)
b=plot_rev3+plot_rev4
b

#SCCed/SCCstoch for different theta under additive preferences
matrix_add <- read_excel("run/SCC_matrix_runs_revisions8.xlsx")
colnames(matrix_add)<-c("id","run_no","rho1y","theta","epsilon_b","J","broken","tp","deterministic","SCC0","finalT","threshold_max","damage","gamma")

matrix_add_det10=subset(matrix_add,matrix_add$deterministic==0 & matrix_add$tp==1)
matrix_add_det10$det=1
matrix_add_stoch10=subset(matrix_add,matrix_add$deterministic==1 & matrix_add$tp==1)
matrix_add_stoch10$det=0

total <- merge(matrix_add_det10,matrix_add_stoch10,by=c("theta","J"))
total$value = total$SCC0.x/total$SCC0.y
total=subset(total, select=c("theta","J","value"))

plot_rev2<-ggplot(total, aes(x=theta, y=value, group=J, linetype=factor(J))) +
  xlab("Elasticity of marginal utility") +
  geom_line(color='black')+
  ylab("Ratio of stochastic SCC to the SCC under expected damages") +
  scale_linetype_manual("Shock J (in %)", values = c("1" = 8, "10" = 3))+
  theme_bw()+
  geom_hline(yintercept=1)

##Stochastic time paths
library(matrixStats)

processing2 <- function(traj_SCC) {
  probs <- c(0.05, 0.5, 0.95)
  traj_SCC=data.matrix(traj_SCC)
  # Row quantiles
  q1 <- as.data.frame(rowQuantiles(traj_SCC, probs = 0.01))
  q1$years=rownames(q1)
  colnames(q1)=c("value_1%","years")
  q50 <- as.data.frame(rowMeans(traj_SCC))
  q50$years=rownames(q50)
  colnames(q50)=c("value_50%","years")
  q99 <- as.data.frame(rowQuantiles(traj_SCC, probs = 0.99))
  q99$years=rownames(q99)
  colnames(q99)=c("value_99%","years")
  #q=rbind(q1,q50,q99)
  q=merge(q1,q50, by="years")
  q=merge(q,q99, by="years")
  q=as.data.frame(q)
  q$years=as.numeric(q$years)
  qcheck=subset(q, q$years<=18)
  return(qcheck)
}

temp_add1 <- read.table("outputs/runs_revisions9add_run0004-2023-03-07-10h40/cumem.csv", sep=";", quote="\"")
add1=processing2(temp_add1)
add1$id="Additive 1%"
temp_add2 <- read.table("outputs/runs_revisions9add_run0005-2023-03-07-16h34/cumem.csv", sep=";", quote="\"")
add2=processing2(temp_add2)
add2$id="Additive 10%"
temp_rs1 <- read.table("outputs/runs_revisions9rs_run0024-2023-03-06-20h53/cumem.csv", sep=";", quote="\"")
rs1=processing2(temp_rs1)
rs1$id="Risk-sensitive 1%"
temp_rs2 <- read.table("outputs/runs_revisions9rs_run0025-2023-03-06-11h39/cumem.csv", sep=";", quote="\"")
rs2=processing2(temp_rs2)
rs2$id="Risk-sensitive 10%"

a=rbind(add1, add2,rs1,rs2)
a$years=2015 + (as.numeric(a$years)-1)*5

plottraj1<- ggplot(a, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  geom_line(color='black')+
  ylab("Stochastic mean temperature increase (in °C)") +
  scale_linetype_manual("Mean path",values = c("Additive 1%" = "longdash", "Risk-sensitive 1%" = "solid", "Additive 10%" = "dotted",  "Risk-sensitive 10%" = "dashed"))+
  theme_bw()+
  theme(legend.position="none")

b1=subset(a, a$id=="Additive 1%")

plottraj2<- ggplot(b1, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic temperature increase (in °C)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Additive 1%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+ 
  theme(legend.position="none")+
  #theme(legend.position="none")+
  ggtitle("Additive, J=1%")

b2=subset(a, a$id=="Additive 10%")

plottraj2b<- ggplot(b2, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic temperature increase (in °C)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Additive 10%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+
  theme(legend.position="none")+
  ggtitle("Additive, J=10%")

b3=subset(a, a$id=="Risk-sensitive 1%")

plottraj2c<- ggplot(b3, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic temperature increase (in °C)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Risk-sensitive 1%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+
  theme(legend.position="none")+
  ggtitle("Risk-sensitive, J=1%")

b4=subset(a, a$id=="Risk-sensitive 10%")

plottraj2d<- ggplot(b4, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic temperature increase (in °C)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Risk-sensitive 10%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+
  theme(legend.position="none")+
  ggtitle("Risk-sensitive, J=10%")

b1$ratio=b3$`value_50%`/b1$`value_50%`
b2$ratio=b4$`value_50%`/b2$`value_50%`
c=rbind(b1,b2)

plottrajc<- ggplot(c, aes(x=years, y=ratio, group=id, linetype=factor(id))) +
  xlab("Years") +
  geom_line(color='black')+
  ylab("Ratio median risk-sensitive SCC to median additive SCC") +
  scale_linetype_manual(values = c("Additive 1%" = "solid", "Additive 10%" = "dashed"))+
  geom_hline(yintercept=1, linetype='dotted')+
  theme_bw()

library(patchwork)
traj=plottraj2+plottraj2c+plottraj2b+plottraj2d

traj_SCC_add1 <- read.table("outputs/runs_revisions9add_run0004-2023-03-07-10h40/traj_SCC.csv", sep=";", quote="\"")
add1=processing2(traj_SCC_add1)
add1$id="Additive 1%"
traj_SCC_add2 <- read.table("outputs/runs_revisions9add_run0005-2023-03-07-16h34/traj_SCC.csv", sep=";", quote="\"")
add2=processing2(traj_SCC_add2)
add2$id="Additive 10%"
traj_SCC_rs1 <- read.table("outputs/runs_revisions9rs_run0024-2023-03-06-20h53/traj_SCC.csv", sep=";", quote="\"")
rs1=processing2(traj_SCC_rs1)
rs1$id="Risk-sensitive 1%"
traj_SCC_rs2 <- read.table("outputs/runs_revisions9rs_run0025-2023-03-06-11h39/traj_SCC.csv", sep=";", quote="\"")
rs2=processing2(traj_SCC_rs2)
rs2$id="Risk-sensitive 10%"

a=rbind(add1, add2,rs1,rs2)
a$years=2015 + (as.numeric(a$years)-1)*5

plottraj1b<- ggplot(a, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  geom_line(color='black')+
  ylab("Mean stochastic SCC (in $/tC)") +
  scale_linetype_manual("Mean path",values = c("Additive 1%" = "longdash", "Risk-sensitive 1%" = "solid", "Additive 10%" = "dotted",  "Risk-sensitive 10%" = "dashed"))+
  theme_bw()+
  theme(legend.position="none")

#temp for all paths

b1=subset(a, a$id=="Additive 1%")

plottraj2<- ggplot(b1, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic SCC (in $/tC)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Additive 1%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+ 
  ylim(0, 2500)+
  theme(legend.position="none")+
  ggtitle("Additive, J=1%")

b2=subset(a, a$id=="Additive 10%")

plottraj2b<- ggplot(b2, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic SCC (in $/tC)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Additive 10%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+
  ylim(0, 2500)+
  theme(legend.position="none")+
  ggtitle("Additive, J=10%")

b3=subset(a, a$id=="Risk-sensitive 1%")

plottraj2c<- ggplot(b3, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic SCC (in $/tC)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Risk-sensitive 1%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+
  ylim(0, 2500)+
  theme(legend.position="none")+
  ggtitle("Risk-sensitive, J=1%")

b4=subset(a, a$id=="Risk-sensitive 10%")

plottraj2d<- ggplot(b4, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic SCC (in $/tC)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Risk-sensitive 10%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+
  ylim(0, 2500)+
  theme(legend.position="none")+
  ggtitle("Risk-sensitive, J=10%")

b1$ratio=b3$`value_50%`/b1$`value_50%`
b2$ratio=b4$`value_50%`/b2$`value_50%`
c=rbind(b1,b2)

plottrajc<- ggplot(c, aes(x=years, y=ratio, group=id, linetype=factor(id))) +
  xlab("Years") +
  geom_line(color='black')+
  ylab("Ratio median risk-sensitive SCC to median additive SCC") +
    scale_linetype_manual(values = c("Additive 1%" = "solid", "Additive 10%" = "dashed"))+
  geom_hline(yintercept=1, linetype='dotted')+
  theme_bw()+
  theme(legend.position="none")

library(patchwork)
traj2=plottraj2+plottraj2c+plottraj2b+plottraj2d

##### abatmeent
abatement_add1 <- read.table("outputs/runs_revisions9add_run0004-2023-03-07-10h40/abat_save.csv", sep=";", quote="\"")
add1=processing2(abatement_add1)
add1$id="Additive 1%"
abatement_add2 <- read.table("outputs/runs_revisions9add_run0005-2023-03-07-16h34/abat_save.csv", sep=";", quote="\"")
add2=processing2(abatement_add2)
add2$id="Additive 10%"
abatement_rs1 <- read.table("outputs/runs_revisions9rs_run0024-2023-03-06-20h53/abat_save.csv", sep=";", quote="\"")
rs1=processing2(abatement_rs1)
rs1$id="Risk-sensitive 1%"
abatement_rs2 <- read.table("outputs/runs_revisions9rs_run0025-2023-03-06-11h39/abat_save.csv", sep=";", quote="\"")
rs2=processing2(abatement_rs2)
rs2$id="Risk-sensitive 10%"

a=rbind(add1, add2,rs1,rs2)
a$years=2015 + (as.numeric(a$years)-1)*5
a$`value_1%`=a$`value_1%`*100
a$`value_50%`=a$`value_50%`*100
a$`value_99%`=a$`value_99%`*100

plottraj1c<- ggplot(a, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  geom_line(color='black')+
  ylab("Mean abatement rate (in %)") +
  scale_linetype_manual("Mean path",values = c("Additive 1%" = "longdash", "Risk-sensitive 1%" = "solid", "Additive 10%" = "dotted",  "Risk-sensitive 10%" = "dashed"))+
  theme_bw()


#abatement for all paths

b1=subset(a, a$id=="Additive 1%")

plottraj2<- ggplot(b1, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic abatement rate (in %)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Additive 1%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+ 
  theme(legend.position="none")+
  ggtitle("Additive, J=1%")

b2=subset(a, a$id=="Additive 10%")

plottraj2b<- ggplot(b2, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic abatement rate (in %)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Additive 10%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+
  theme(legend.position="none")+
  ggtitle("Additive, J=10%")

b3=subset(a, a$id=="Risk-sensitive 1%")

plottraj2c<- ggplot(b3, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic abatement rate (in %)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Risk-sensitive 1%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+
  theme(legend.position="none")+
  ggtitle("Risk-sensitive, J=1%")
                
b4=subset(a, a$id=="Risk-sensitive 10%")

plottraj2d<- ggplot(b4, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  ylab("Stochastic abatement rate (in %)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Risk-sensitive 10%" = "solid"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`),fill = "grey70")+
  geom_line(color='black',size=1.2)+
  theme_bw()+
  theme(legend.position="none")+
  ggtitle("Risk-sensitive, J=10%")

b1$ratio=b3$`value_50%`/b1$`value_50%`
b2$ratio=b4$`value_50%`/b2$`value_50%`
c=rbind(b1,b2)

plottrajc<- ggplot(c, aes(x=years, y=ratio, group=id, linetype=factor(id))) +
  xlab("Years") +
  geom_line(color='black')+
  ylab("Ratio median risk-sensitive SCC to median additive SCC") +
  scale_linetype_manual(values = c("Additive 1%" = "solid", "Additive 10%" = "dashed"))+
  geom_hline(yintercept=1, linetype='dotted')+
  theme_bw()+
  theme(legend.position="none")

library(patchwork)
traj3=plottraj2+plottraj2c+plottraj2b+plottraj2d

c=subset(a,a$id=="Risk-sensitive 1%" | a$id=="Risk-sensitive 10%")

plottraj3<- ggplot(c, aes(x=years, y=`value_50%`, group=id, linetype=factor(id))) +
  xlab("Years") +
  geom_line(color='black')+
  ylab("Stochastic SCC (in $/tC)") +
  scale_linetype_manual("Quantiles of the distribution",values = c("Risk-sensitive 1%" = "longdash", "Risk-sensitive 10%" = "dotted"))+
  geom_ribbon(aes(ymin=`value_99%`, ymax=`value_1%`, fill = id))+
  theme_bw()

library(patchwork)
plot_traj_global=plottraj1+plottraj1b+plottraj1c

temp <- read.table("tempt.csv", sep=";", quote="\"")
library(matrixStats)
probs <- c(0.01,0.05, 0.1,0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
cume=data.matrix(temp)
# Row quantiles
q1 <- as.data.frame(rowQuantiles(cume, probs = 0.01))
q1$years=rownames(q1)
q1$id="1%"
colnames(q1)=c("value","years","id")
q50 <- as.data.frame(rowQuantiles(cume, probs = 0.5))
q50$years=rownames(q50)
q50$id="50%"
colnames(q50)=c("value","years","id")
q99 <- as.data.frame(rowQuantiles(cume, probs = 0.99))
q99$years=rownames(q99)
q99$id="99%"
colnames(q99)=c("value","years","id")
q=rbind(q1,q50,q99)
q=as.data.frame(q)
q$years=as.numeric(q$years)

qcheck=subset(q, q$years<=20)

ggplot(qcheck, aes(x=years, y=value, group=id, linetype=factor(id))) +
  xlab("Years") +
  geom_line(color='black')+
  ylab("Temperature (in °C)") +
  scale_linetype_manual("Quantiles of the distribution", values = c("1%" = 8, "99%" = 6, "50%" = 3))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

