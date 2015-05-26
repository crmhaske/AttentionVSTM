
#****CHANGE THIS TO THE LOCATION OF YOUR GIT REPOSITORY****
setwd("C:/Users/Christie/Documents/Research/gits/attnGuessNew")

#Required custom functions
source(paste(getwd(),"/rLibraries/funcs.R",sep=""))
source(paste(getwd(),"/rLibraries/readdata.R",sep=""))

#Required libraries
library(ggplot2)
library(plyr)
library(heplots)
library(e1071)
library(grid)
library(lme4)
library(lmerTest)
library(nlme)
library(multcomp) 

#Read in the data, calculate the corrected mean, angular difference, 
#and remove practice trials
#------------------

#100% cue, one sequential trial type
d1<-read.data(paste(getwd(),"/exp1/data",sep=""))
  #Corrected angdiff
  d1<-ddply(d1,.(userid,validCue,blockType),transform,meanAngdiff=mean(angdiff))
  d1$angdiff.cor<-d1$angdiff-d1$meanAngdiff

  #Remove practice trials
  dnp1<-subset(d1,d1$practice=="FALSE")

#80% cue, one sequential trial type
d2<-read.data(paste(getwd(),"/exp2/data",sep=""))
  #Corrected angdiff
  d2<-ddply(d2,.(userid,validCue,blockType),transform,meanAngdiff=mean(angdiff))
  d2$angdiff.cor<-d2$angdiff-d2$meanAngdiff

  #Remove practice trials
  dnp2<-subset(d2,d2$practice=="FALSE")
  
  #Invalid trial anomaly
  dnp2<-dnp2[!(dnp2$userid==17),]
  dnp2<-dnp2[!(dnp2$userid==7),]

#50% cue, one sequential trial type
d3<-read.data(paste(getwd(),"/exp3/data",sep=""))
  #Corrected angdiff
  d3<-ddply(d3,.(userid,validCue,blockType),transform,meanAngdiff=mean(angdiff))
  d3$angdiff.cor<-d3$angdiff-d3$meanAngdiff

  #Remove practice trials
  dnp3<-subset(d3,d3$practice=="FALSE")
  
#80% cue, two sequential trial types
d4<-read.data(paste(getwd(),"/exp4/data",sep=""))
  #Corrected angdiff
  d4<-ddply(d4,.(userid,validCue,blockType),transform,meanAngdiff=mean(angdiff))
  d4$angdiff.cor<-d4$angdiff-d4$meanAngdiff

  #Remove practice trials
  dnp4<-subset(d4,d4$practice=="False")
  
  #Remove subject 21 as they did not finish
  dnp4<-subset(dnp4,dnp4$userid!=21)

#50% cue, two sequential trial types
d5<-read.data(paste(getwd(),"/exp5/data",sep=""))
  #Corrected angdiff
  d5<-ddply(d5,.(userid,validCue,blockType),transform,meanAngdiff=mean(angdiff))
  d5$angdiff.cor<-d5$angdiff-d5$meanAngdiff

  #Remove practice trials
  dnp5<-subset(d5,d5$practice=="False")
  
#Remove participants based on the selection criteria
#---------------------------------------------------------------------
#Participants who rotated the grating in the wrong direction
#less often than chance were removed

removal<-function(d) {
  sel<-prop.cor(d)
  rm<-sel[which(sel$prop.cor <= 0.5),]$userid
  if (length(rm>0)) {
    for (i in 1:length(rm)) {
      d<-subset(d,d$userid!=rm[i])
    }
  }
  d$numRem<-rep(length(rm),length(d$angdiff))
  return(d)
}

dnp1<-removal(dnp1)
unique(dnp1$numRem)

dnp2<-removal(dnp2)
unique(dnp2$numRem)

dnp3<-removal(dnp3)
unique(dnp3$numRem)

dnp4<-removal(dnp4)
unique(dnp4$numRem)

dnp5<-removal(dnp5)
unique(dnp5$numRem)


#Experiment 1
#--------------------------------------------------

#Covert userid from numeric to factor
dnp1$userid<-as.factor(dnp1$userid)

#Dataframes for blockType and no cue
dsimul1<-subset(dnp1,dnp1$blockType=="simul")
dseq1<-subset(dnp1,dnp1$blockType=="seq")
dneut1<-subset(dnp1,dnp1$validCue=="neutral")

#Replication of Liu & Becker, uncued trials, variance
var.user.nC1<-ddply(dneut1,.(userid,blockType),summarise,var=var(angdiff))
var.nC.aov1<-aov(var ~ blockType + Error(userid/blockType),data = var.user.nC1)
summary(var.nC.aov1)
etasq(var.nC.aov1[[3]])

#Variance plot
dplot<-var.user.nC1
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq","simul"),
                           to=c("Sequential", "Simultaneous"))
p<-ggplot(dplot,aes(x=blockType,y=var))+geom_boxplot()+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Block Type")+ylab("Variance")+ylim(c(0,3000))
ggsave(p,file=paste(getwd(),"/Christie/Fig2.jpg",sep=""),dpi=225,width=6,height=5,units="cm")

#Cueing and block type, variance
var.user1<-ddply(dnp1,.(userid,blockType,validCue),summarise,var=var(angdiff))
var.aov1<-aov(var ~ blockType+validCue + Error(userid/(blockType+validCue)),data = var.user1)
summary(var.aov1)
etasq(var.aov1[[3]])

#Cueing and block type, median abs error
err.user1<-ddply(dnp1,.(userid,blockType,validCue),summarise,medErr=median(abs(angdiff.cor)))
err.aov1<-aov(medErr ~ blockType*validCue + Error(userid/(blockType*validCue)),data = err.user1)
summary(err.aov1)
etasq(err.aov1[[5]])

#Cueing on simultaneous trials
err.simul.user1<-ddply(dsimul1,.(userid,validCue),summarise,medErr=median(abs(angdiff.cor)))
err.simul.aov1<-aov(medErr ~ validCue + Error(userid/validCue),data = err.simul.user1)
summary(err.simul.aov1)
etasq(err.simul.aov1[[3]])

#Cueing on sequential trials
err.seq.user1<-ddply(dseq1,.(userid,validCue),summarise,medErr=median(abs(angdiff.cor)))
err.seq.aov1<-aov(medErr ~ validCue + Error(userid/validCue),data = err.seq.user1)
summary(err.seq.aov1)

#Median, absolute error plot
dplot<-err.user1
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq","simul"),
                           to=c("Sequential", "Simultaneous"))
dplot$validCue<-mapvalues(dplot$validCue,
                           from=c("neutral","TRUE"),
                           to=c("Uncued", "Validly Cued"))
p<-ggplot(dplot,aes(x=validCue,y=medErr))+geom_boxplot()+facet_grid(.~blockType)+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Cueing")+ylab("Med. Abs. Error (deg.)")+ylim(c(0,50))
ggsave(p,file=paste(getwd(),"/Christie/Fig3.jpg",sep=""),dpi=225,width=12,height=5,units="cm")

#Model fitting, g
fit.user1<-ddply(dnp1,.(userid,blockType,validCue),summarise,g=mixedfit(angdiff)$g,sigma=mixedfit(angdiff)$sigma,mu=mixedfit(angdiff)$mu)
fit.aov1.g<-aov(g ~ blockType*validCue + Error(userid/(blockType*validCue)),data = fit.user1)
summary(fit.aov1.g)
etasq(fit.aov1.g[[3]])

#Model fitting, sigma
fit.aov1.sigma<-aov(sigma ~ blockType*validCue + Error(userid/(blockType*validCue)),data = fit.user1)
summary(fit.aov1.sigma)
etasq(fit.aov1.sigma[[4]])
etasq(fit.aov1.sigma[[5]])

#Model fitting, sigma, simultaneous
fit.simul.aov1<-aov(sigma~validCue+Error(userid/validCue),data=fit.user1[which(fit.user1$blockType=="simul"),])
summary(fit.simul.aov1)
etasq(fit.simul.aov1[[3]])

#Model fitting, sigma, sequential
fit.seq.aov1<-aov(sigma~validCue+Error(userid/validCue),data=fit.user1[which(fit.user1$blockType=="seq"),])
summary(fit.seq.aov1)

#Guessing plot
dplot<-fit.user1
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq","simul"),
                           to=c("Sequential", "Simultaneous"))
dplot$validCue<-mapvalues(dplot$validCue,
                          from=c("neutral","TRUE"),
                          to=c("Uncued", "Validly Cued"))
p<-ggplot(dplot,aes(x=validCue,y=g))+geom_boxplot()+facet_grid(.~blockType)+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Cueing")+ylab("Guessing Prop.")+ylim(c(0,1))
ggsave(p,file=paste(getwd(),"/Christie/Fig4a.jpg",sep=""),dpi=225,width=12,height=5,units="cm")

#Precision plot
dplot<-fit.user1
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq","simul"),
                           to=c("Sequential", "Simultaneous"))
dplot$validCue<-mapvalues(dplot$validCue,
                          from=c("neutral","TRUE"),
                          to=c("Uncued", "Validly Cued"))
p<-ggplot(dplot,aes(x=validCue,y=sigma*180/pi))+geom_boxplot()+facet_grid(.~blockType)+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Cueing")+ylab("Precision (deg.)")+ylim(c(0,50))
ggsave(p,file=paste(getwd(),"/Christie/Fig4b.jpg",sep=""),dpi=225,width=12,height=5,units="cm")


#Experiment 2
#--------------------------------------------------

#Covert userid from numeric to factor
dnp2$userid<-as.factor(dnp2$userid)

#Dataframes for blockType and no cue
dsimul2<-subset(dnp2,dnp2$blockType=="simul")
dseq2<-subset(dnp2,dnp2$blockType=="seq")
dneut2<-subset(dnp2,dnp2$validCue=="neutral")

#Replication of Liu & Becker, uncued trials, variance
var.user.nC2<-ddply(dneut2,.(userid,blockType),summarise,var=var(angdiff))
var.nC.aov2<-aov(var ~ blockType + Error(userid/blockType),data = var.user.nC2)
summary(var.nC.aov2)
etasq(var.nC.aov2[[3]])

#Cueing and block type, median abs error
err.user2<-ddply(dnp2,.(userid,blockType,validCue),summarise,medErr=median(abs(angdiff.cor)))
err.aov2<-aov(medErr ~ blockType*validCue + Error(userid/(blockType*validCue)),data = err.user2)
summary(err.aov2)
etasq(err.aov2[[4]])
etasq(err.aov2[[5]])

#Median, absolute error plot
dplot<-err.user2
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq","simul"),
                           to=c("Sequential", "Simultaneous"))
dplot$validCue<-mapvalues(dplot$validCue,
                          from=c("neutral","TRUE","FALSE"),
                          to=c("Uncued", "Validly Cued","Invalidly Cued"))
p<-ggplot(dplot,aes(x=validCue,y=medErr))+geom_boxplot()+facet_grid(.~blockType)+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Cueing")+ylab("Med. Abs. Error (deg.)")+ylim(c(0,60))
ggsave(p,file=paste(getwd(),"/Christie/Fig5.jpg",sep=""),dpi=225,width=12,height=5,units="cm")

#Linear Tukey constrats between the levels of cueing, simultaneous
err.simul.user2<-err.user2[which(err.user2$blockType=="simul"),]
err.simul.user2$validCue<-factor(err.simul.user2$validCue)
cueing.simul.lme2<-lme(medErr ~ validCue, random = ~1 | userid/validCue , data = err.simul.user2)
cueing.simul.lme.glht2 <- glht(cueing.simul.lme2,linfct=mcp(validCue = "Tukey"))
summary(cueing.simul.lme.glht2)

#Linear Tukey constrats between the levels of cueing, sequential
err.seq.user2<-err.user2[which(err.user2$blockType=="seq"),]
err.seq.user2$validCue<-factor(err.seq.user2$validCue)
cueing.seq.lme2<-lme(medErr ~ validCue, random = ~1 | userid/validCue , data = err.seq.user2)
cueing.seq.lme.glht2 <- glht(cueing.seq.lme2,linfct=mcp(validCue = "Tukey"))
summary(cueing.seq.lme.glht2)

#Experiment 3
#--------------------------------------------------

#Covert userid from numeric to factor
dnp3$userid<-as.factor(dnp3$userid)

#Dataframes for blockType and no cue
dsimul3<-subset(dnp3,dnp3$blockType=="simul")
dseq3<-subset(dnp3,dnp3$blockType=="seq")
dneut3<-subset(dnp3,dnp3$validCue=="neutral")

#Replication of Liu & Becker, uncued trials, variance
var.user.nC3<-ddply(dneut3,.(userid,blockType),summarise,var=var(angdiff))
var.nC.aov3<-aov(var ~ blockType + Error(userid/blockType),data = var.user.nC3)
summary(var.nC.aov3)
etasq(var.nC.aov3[[3]])

#Cueing and block type, median abs error
err.user3<-ddply(dnp3,.(userid,blockType,validCue),summarise,medErr=median(abs(angdiff.cor)))
err.aov3<-aov(medErr ~ blockType*validCue + Error(userid/(blockType*validCue)),data = err.user3)
summary(err.aov3)
etasq(err.aov3[[4]])
etasq(err.aov3[[5]])

dplot<-err.user3
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq","simul"),
                           to=c("Sequential", "Simultaneous"))
dplot$validCue<-mapvalues(dplot$validCue,
                          from=c("neutral","TRUE","FALSE"),
                          to=c("Uncued", "Validly Cued","Invalidly Cued"))
p<-ggplot(dplot,aes(x=validCue,y=medErr))+geom_boxplot()+facet_grid(.~blockType)+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Cueing")+ylab("Med. Abs. Error (deg.)")+ylim(c(0,60))
ggsave(p,file=paste(getwd(),"/Christie/Fig6.jpg",sep=""),dpi=225,width=12,height=5,units="cm")


#Experiment 4
#--------------------------------------------------

#Covert userid from numeric to factor
dnp4$userid<-as.factor(dnp4$userid)

#Dataframes for blockType
dsimul4<-subset(dnp4,dnp4$blockType=="simul")
dseq4<-subset(dnp4,dnp4$blockType!="simul")

#Trials where the cue and target appeared temporally together
dseq4.1<-subset(dseq4,dseq4$blockType=="seq1" & dseq4$stimResp==1)
dseq4.2<-subset(dseq4,dseq4$blockType=="seq2" & dseq4$stimResp==2)
dseq4<-rbind(dseq4.1,dseq4.2)
dnp4<-rbind(dsimul4,dseq4)

#Dataframe for noe cue
dneut4<-subset(dnp4,dnp4$validCue=="neutral")

#Replication of Liu & Becker, uncued trials, variance
var.user.nC4<-ddply(dneut4,.(userid,blockType),summarise,var=var(angdiff))
var.nC.aov4<-aov(var ~ blockType + Error(userid/blockType),data = var.user.nC4)
summary(var.nC.aov4)
etasq(var.nC.aov4[[3]])

#Variance plot
dplot<-var.user.nC4
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq1","seq2","simul"),
                           to=c("Sequential:\n Cue First", "Sequential:\n Cue Second","Simultaneous"))
p<-ggplot(dplot,aes(x=blockType,y=var))+geom_boxplot()+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Block Type")+ylab("Variance")
ggsave(p,file=paste(getwd(),"/Christie/Fig7.jpg",sep=""),dpi=225,width=8,height=5,units="cm")


#Cueing and block type, median abs error
err.user4<-ddply(dnp4,.(userid,blockType,validCue),summarise,medErr=median(abs(angdiff.cor)))
err.aov4<-aov(medErr ~ blockType*validCue + Error(userid/(blockType*validCue)),data = err.user4)
summary(err.aov4)
etasq(err.aov4[[4]])
etasq(err.aov4[[5]])

#Median, absolute error plot
dplot<-err.user4
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq1","seq2","simul"),
                           to=c("Sequential:\n Cue First", "Sequential:\n Cue Second","Simultaneous"))
dplot$validCue<-mapvalues(dplot$validCue,
                          from=c("neutral","true","false"),
                          to=c("Uncued", "Validly Cued","Invalidly Cued"))
p<-ggplot(dplot,aes(x=validCue,y=medErr))+geom_boxplot()+facet_grid(.~blockType)+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Cueing")+ylab("Med. Abs. Error (deg.)")+ylim(c(0,60))
ggsave(p,file=paste(getwd(),"/Christie/Fig8.jpg",sep=""),dpi=225,width=18,height=5,units="cm")

#Linear Tukey constrats between the levels of block type, Invalid
err.inval.user4<-err.user4[which(err.user4$validCue=="false"),]
err.inval.user4$blockType<-factor(err.inval.user4$blockType)
cueing.inval.lme4<-lme(medErr ~ blockType, random = ~1 | userid/blockType , data = err.inval.user4)
cueing.inval.lme.glht4 <- glht(cueing.inval.lme4,linfct=mcp(blockType = "Tukey"))
summary(cueing.inval.lme.glht4)

#Linear Tukey constrats between the levels of block type, uncued
err.nC.user4<-err.user4[which(err.user4$validCue=="neutral"),]
err.nC.user4$blockType<-factor(err.nC.user4$blockType)
cueing.nC.lme4<-lme(medErr ~ blockType, random = ~1 | userid/blockType , data = err.nC.user4)
cueing.nC.lme.glht4 <- glht(cueing.nC.lme4,linfct=mcp(blockType = "Tukey"))
summary(cueing.nC.lme.glht4)

#Linear Tukey constrats between the levels of block type, valid
err.val.user4<-err.user4[which(err.user4$validCue=="true"),]
err.val.user4$blockType<-factor(err.val.user4$blockType)
cueing.val.lme4<-lme(medErr ~ blockType, random = ~1 | userid/blockType , data = err.val.user4)
cueing.val.lme.glht4 <- glht(cueing.val.lme4,linfct=mcp(blockType = "Tukey"))
summary(cueing.val.lme.glht4)

#Different invalid trials for sequential trials
dseq.inval4<-subset(dseq4,dseq4$validCue=="false")
err.seqInval.user<-ddply(dseq.inval4,.(userid,atStimNotCued,blockType),summarise,medErr=median(abs(angdiff.cor)))

#Invalid trial, median absolute error plot
dplot<-err.seqInval.user
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq1","seq2"),
                           to=c("Sequential:Cue First", "Sequential:Cue Second"))
dplot$atStimNotCued<-mapvalues(dplot$atStimNotCued,
                          from=c("true","false"),
                          to=c("Cued Non-Target","Empty Space"))
p<-ggplot(dplot,aes(x=atStimNotCued,y=medErr))+geom_boxplot()+facet_grid(.~blockType)+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Invalid Cue Type")+ylab("Med. Abs. Error (deg.)")
ggsave(p,file=paste(getwd(),"/Christie/Fig9.jpg",sep=""),dpi=225,width=12,height=5,units="cm")


#Different invalid trials for sequential trials, anova
err.seqInval.aov<-aov(medErr~atStimNotCued*blockType+Error(userid/(atStimNotCued*blockType)), data = err.seqInval.user)
summary(err.seqInval.aov)
etasq(err.seqInval.aov[[3]])

#Experiment 5
#--------------------------------------------------

#Covert userid from numeric to factor
dnp5$userid<-as.factor(dnp5$userid)

#Dataframes for blockType
dsimul5<-subset(dnp5,dnp5$blockType=="simul")
dseq5<-subset(dnp5,dnp5$blockType!="simul")

#Trials where the cue and target appeared temporally together
dseq5.1<-subset(dseq5,dseq5$blockType=="seq1" & dseq5$stimResp==1)
dseq5.2<-subset(dseq5,dseq5$blockType=="seq2" & dseq5$stimResp==2)
dseq5<-rbind(dseq5.1,dseq5.2)
dnp5<-rbind(dsimul5,dseq5)

#Dataframe for noe cue
dneut5<-subset(dnp5,dnp5$validCue=="neutral")

#Replication of Liu & Becker, uncued trials, variance
var.user.nC5<-ddply(dneut5,.(userid,blockType),summarise,var=var(angdiff))
var.nC.aov5<-aov(var ~ blockType + Error(userid/blockType),data = var.user.nC5)
summary(var.nC.aov5)
etasq(var.nC.aov5[[3]])

#Cueing and block type, median abs error
err.user5<-ddply(dnp5,.(userid,blockType,validCue),summarise,medErr=median(abs(angdiff.cor)))
err.aov5<-aov(medErr ~ blockType*validCue + Error(userid/(blockType*validCue)),data = err.user5)
summary(err.aov5)
etasq(err.aov5[[4]])
etasq(err.aov5[[5]])

#Median absolute erorr plot
dplot<-err.user5
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq1","seq2","simul"),
                           to=c("Sequential:\n Cue First", "Sequential:\n Cue Second","Simultaneous"))
dplot$validCue<-mapvalues(dplot$validCue,
                          from=c("neutral","true","false"),
                          to=c("Uncued", "Validly Cued","Invalidly Cued"))
p<-ggplot(dplot,aes(x=validCue,y=medErr))+geom_boxplot()+facet_grid(.~blockType)+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Cueing")+ylab("Med. Abs. Error (deg.)")+ylim(c(0,60))
ggsave(p,file=paste(getwd(),"/Christie/Fig10.jpg",sep=""),dpi=225,width=18,height=5,units="cm")

#Collapsing across all experiments
dnp4$validCue<-ifelse(dnp4$validCue=="false","FALSE",
                      ifelse(dnp4$validCue=="true","TRUE",
                             "neutral"))
dnp5$validCue<-ifelse(dnp5$validCue=="false","FALSE",
                      ifelse(dnp4$validCue=="true","TRUE",
                             "neutral"))
dall<-rbind(dnp1,dnp2,dnp3,dnp4,dnp5)
dall$userid<-as.factor(dall$userid)
dall$isSeq<-dall$blockType!="simul"
err.user.all<-ddply(dall,.(Experiment,userid,isSeq,validCue),summarise,medErr=median(abs(angdiff)))
err.aov.all<-aov(medErr~Experiment*validCue*as.factor(isSeq)+Error(userid/(validCue*as.factor(isSeq))),data=err.user.all)
summary(err.aov.all)

dplot<-err.user.all
dplot$isSeq<-mapvalues(dplot$isSeq,
                           from=c("TRUE","FALSE"),
                           to=c("Sequential", "Simultaneous"))
dplot$validCue<-mapvalues(dplot$validCue,
                          from=c("neutral","TRUE","FALSE"),
                          to=c("Uncued", "Validly Cued","Invalidly Cued"))
dplot$validCue<-factor(dplot$validCue,levels=c("Invalidly Cued","Uncued","Validly Cued"))
p<-ggplot(dplot,aes(x=validCue,y=medErr))+geom_boxplot()+facet_grid(.~isSeq)+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Cueing")+ylab("Med. Abs. Error (deg.)")
ggsave(p,file=paste(getwd(),"/Christie/Fig11.jpg",sep=""),dpi=225,width=12,height=5,units="cm")

#Different invalid trials for sequential trials
dseq.inval5<-subset(dseq5,dseq5$validCue=="false")
err.seqInval.user5<-ddply(dseq.inval5,.(userid,atStimNotCued,blockType),summarise,medErr=median(abs(angdiff.cor)))

dplot<-err.seqInval.user5
dplot$blockType<-mapvalues(dplot$blockType,
                           from=c("seq1","seq2"),
                           to=c("Sequential:Cue First", "Sequential:Cue Second"))
dplot$atStimNotCued<-mapvalues(dplot$atStimNotCued,
                               from=c("true","false"),
                               to=c("Cued Non-Target","Empty Space"))
p<-ggplot(dplot,aes(x=atStimNotCued,y=medErr))+geom_boxplot()+facet_grid(.~blockType)+
  theme_bw()+ theme(text=element_text(size=10),legend.position="none",plot.margin=unit(x=c(0,0,0,0),units="mm"))+
  xlab("Invalid Cue Type")+ylab("Med. Abs. Error (deg.)")
ggsave(p,file=paste(getwd(),"/Christie/Fig12.jpg",sep=""),dpi=225,width=12,height=5,units="cm")


#Different invalid trials for sequential trials, anova
err.seqInval.aov5<-aov(medErr~atStimNotCued*blockType+Error(userid/(atStimNotCued*blockType)), data = err.seqInval.user5)
summary(err.seqInval.aov5)
etasq(err.seqInval.aov[[3]])