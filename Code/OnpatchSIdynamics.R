library(ggplot2)
library(reshape)

setwd("C:\\Users\\lucyj\\Downloads")

Time<-0:400

onepatchinfect<-read.csv("SImodel_onepatch_infected_total_12_10_19.csv",header=FALSE)
onepatchRT<-read.csv("SImodel_onepatch_RT_total_12_10_19.csv",header=FALSE)
onepatchWT<-read.csv("SImodel_onepatch_WT_total_12_10_19.csv",header=FALSE)
onepatchinfectRT<-read.csv("SImodel_onepatch_infected_RT_12_10_19.csv",header=FALSE)
onepatchinfectWT<-read.csv("SImodel_onepatch_infected_WT_12_10_19.csv",header=FALSE)

tiff("onepatchSImodel_12_10_19_truncatedversion.tiff", width=10,height=10, units='in',res=600)
par(mfrow=c(2,2))
par(mar=c(2,4.5,1,2))



plot(Time,onepatchWT[,1]+onepatchRT[,1],type="l",lwd=2,ylab=expression(bold("Host_total: Robust+Wild")),xlab="",ylim=c(2250,2400),cex.lab=1.5)
points(Time,onepatchWT[,4]+onepatchRT[,4],type="l",lwd=2,lty=3,col="")
points(Time,onepatchWT[,7]+onepatchRT[,7],type="l",lwd=2,lty=5,col="")
#points(Time,onepatchWT[,4]+onepatchRT[,4],type="l",lwd=2,lty=8)
points(Time,onepatchWT[,9]+onepatchRT[,9],type="l",lwd=2,lty=10,col="")
legend(1000,2000,c("1","1.3","1.6","1.8"),lty=c(1,3,5,10),lwd=c(2,2,2,2),title="Mortality Ratio",bty="n")

# plot(Time,onepatchRT[,1],type="l",lwd=2,main="Robust type: S+I",xlab="",ylab="",ylim=c(200,2500))
# points(Time,onepatchRT[,4],type="l",lwd=2,lty=3)
# points(Time,onepatchRT[,7],type="l",lwd=2,lty=5)
# points(Time,onepatchRT[,9],type="l",lwd=2,lty=10)
# #points(Time,onepatchRT[,9],type="l",lwd=2,lty=10)
# #legend(1500,250,c("1","1.4","1.8"),lty=c(1,3,5),type=c(1,1,1),lwd=c(2,2,2),title="Mortality Ratio")

plot(Time,onepatchWT[,1]/onepatchRT[,1],type="l",lwd=2,ylab=expression(bold("Host_ratio: Wild/Robust")),xlab="",ylim=c(1.442,1.4435),cex.lab=1.5)
points(Time,onepatchWT[,4]/onepatchRT[,4],type="l",lwd=2,lty=3,col="")
points(Time,onepatchWT[,7]/onepatchRT[,7],type="l",lwd=2,lty=5,col="")
     #points(Time,onepatchWT[,4],type="l",lwd=2,lty=8)
points(Time,onepatchWT[,9]/onepatchRT[,9],type="l",lwd=2,lty=10,col="")
     # legend(1500,3000,c("1","1.4","1.8"),lty=c(1,3,5),lwd=c(2,2,2),title="Mortality Ratio",bty="n"),type="l",lwd=2,main="Wild/Robust",xlab="",ylab="")


plot(Time,onepatchinfect[,1],type="l",lwd=2,ylab=expression(bold("Infected_total: Robust+Wild")),xlab="",cex.lab=1.5)
points(Time,onepatchinfect[,4],type="l",lwd=2,lty=3)
points(Time,onepatchinfect[,7],type="l",lwd=2,lty=5)
#points(Time,onepatchinfect[,4],type="l",lwd=2,lty=8)
points(Time,onepatchinfect[,9],type="l",lwd=2,lty=10)
#legend(1500,250,c("1","1.4","1.8"),lty=c(1,3,5),type=c(1,1,1),lwd=c(2,2,2),title="Mortality Ratio")

plot(Time,onepatchinfectWT[,1]/onepatchinfectRT[,1],type="l",lwd=2,ylab=expression(bold("Infect_ratio: Wild/Robust")),xlab="",ylim=c(0.8,1.5),cex.lab=1.5)
points(Time,onepatchinfectWT[,4]/onepatchinfectRT[,4],type="l",lwd=2,lty=3)
points(Time,onepatchinfectWT[,7]/onepatchinfectRT[,7],type="l",lwd=2,lty=5)
#points(Time,onepatchWT[,4],type="l",lwd=2,lty=8)
points(Time,onepatchinfectWT[,9]/onepatchinfectRT[,9],type="l",lwd=2,lty=10)
# legend(1500,3000,c("1","1.4","1.8"),lty=c(1,3,5),lwd=c(2,2,2),title="Mortality Ratio",bty="n")




dev.off()










tiff("onepatchSImodel_12_10_19_logversion.tiff", width=10,height=10, units='in',res=600)
par(mfrow=c(2,2))
par(mar=c(2,4,3,2))



plot(Time,log(onepatchWT[,1]+onepatchRT[,1]),type="l",lwd=2,ylab="Host_total:Robust+Wild",xlab="",ylim=c(log(400),log(2400)))
points(Time,log(onepatchWT[,4]+onepatchRT[,4]),type="l",lwd=2,lty=3)
points(Time,log(onepatchWT[,7]+onepatchRT[,7]),type="l",lwd=2,lty=5)
#points(Time,onepatchWT[,4]+onepatchRT[,4],type="l",lwd=2,lty=8)
points(Time,log(onepatchWT[,9]+onepatchRT[,9]),type="l",lwd=2,lty=10)
legend(1000,2000,c("1","1.3","1.6","1.8"),lty=c(1,3,5,10),lwd=c(2,2,2,2),title="Mortality Ratio",bty="n")

# plot(Time,onepatchRT[,1],type="l",lwd=2,main="Robust type: S+I",xlab="",ylab="",ylim=c(200,2500))
# points(Time,onepatchRT[,4],type="l",lwd=2,lty=3)
# points(Time,onepatchRT[,7],type="l",lwd=2,lty=5)
# points(Time,onepatchRT[,9],type="l",lwd=2,lty=10)
# #points(Time,onepatchRT[,9],type="l",lwd=2,lty=10)
# #legend(1500,250,c("1","1.4","1.8"),lty=c(1,3,5),type=c(1,1,1),lwd=c(2,2,2),title="Mortality Ratio")

plot(Time,onepatchWT[,1]/onepatchRT[,1],type="l",lwd=2,ylab="Host_ratio: Wild/Robust",xlab="",ylim=c(1.3,1.45))
points(Time,onepatchWT[,4]/onepatchRT[,4],type="l",lwd=2,lty=3)
points(Time,onepatchWT[,7]/onepatchRT[,7],type="l",lwd=2,lty=5)
#points(Time,onepatchWT[,4],type="l",lwd=2,lty=8)
points(Time,onepatchWT[,9]/onepatchRT[,9],type="l",lwd=2,lty=10)
# legend(1500,3000,c("1","1.4","1.8"),lty=c(1,3,5),lwd=c(2,2,2),title="Mortality Ratio",bty="n"),type="l",lwd=2,main="Wild/Robust",xlab="",ylab="")


plot(Time,onepatchinfect[,1],type="l",lwd=2,ylab="Infected_total",xlab="")
points(Time,onepatchinfect[,4],type="l",lwd=2,lty=3)
points(Time,onepatchinfect[,7],type="l",lwd=2,lty=5)
#points(Time,onepatchinfect[,4],type="l",lwd=2,lty=8)
points(Time,onepatchinfect[,9],type="l",lwd=2,lty=10)
#legend(1500,250,c("1","1.4","1.8"),lty=c(1,3,5),type=c(1,1,1),lwd=c(2,2,2),title="Mortality Ratio")

plot(Time,onepatchinfectWT[,1]/onepatchinfectRT[,1],type="l",lwd=2,ylab="Infect_ratio: Wild/Robust",xlab="",ylim=c(0.8,1.5))
points(Time,onepatchinfectWT[,4]/onepatchinfectRT[,4],type="l",lwd=2,lty=3)
points(Time,onepatchinfectWT[,7]/onepatchinfectRT[,7],type="l",lwd=2,lty=5)
#points(Time,onepatchWT[,4],type="l",lwd=2,lty=8)
points(Time,onepatchinfectWT[,9]/onepatchinfectRT[,9],type="l",lwd=2,lty=10)
# legend(1500,3000,c("1","1.4","1.8"),lty=c(1,3,5),lwd=c(2,2,2),title="Mortality Ratio",bty="n")




dev.off()











library(ggplot2)
library(reshape)

setwd("C:\\Users\\lucyj\\Downloads")

Time<-0:2500

onepatchinfect<-read.csv("SImodel_onepatch_infected_total_long12_10_19.csv",header=FALSE)
onepatchRT<-read.csv("SImodel_onepatch_RT_total_long12_10_19.csv",header=FALSE)
onepatchWT<-read.csv("SImodel_onepatch_WT_total_long12_10_19.csv",header=FALSE)
onepatchinfectRT<-read.csv("SImodel_onepatch_infected_RT_long12_10_19.csv",header=FALSE)
onepatchinfectWT<-read.csv("SImodel_onepatch_infected_WT_long12_10_19.csv",header=FALSE)

tiff("onepatchSImodel_12_10_19_long250genetruncatedversion.tiff", width=10,height=10, units='in',res=600)
par(mfrow=c(2,2))
par(mar=c(2,4.5,1,2))



plot(Time,onepatchWT[,1]+onepatchRT[,1],type="l",lwd=2,ylab=expression(bold("Host_total: Robust+Wild")),xlab="",ylim=c(2200,2400),cex.lab=1.5)
points(Time,onepatchWT[,4]+onepatchRT[,4],type="l",lwd=2,lty=3)
points(Time,onepatchWT[,7]+onepatchRT[,7],type="l",lwd=2,lty=5)
#points(Time,onepatchWT[,4]+onepatchRT[,4],type="l",lwd=2,lty=8)
points(Time,onepatchWT[,9]+onepatchRT[,9],type="l",lwd=2,lty=10)
legend(1000,2000,c("1","1.3","1.6","1.8"),lty=c(1,3,5,10),lwd=c(2,2,2,2),title="Mortality Ratio",bty="n")

# plot(Time,onepatchRT[,1],type="l",lwd=2,main="Robust type: S+I",xlab="",ylab="",ylim=c(200,2500))
# points(Time,onepatchRT[,4],type="l",lwd=2,lty=3)
# points(Time,onepatchRT[,7],type="l",lwd=2,lty=5)
# points(Time,onepatchRT[,9],type="l",lwd=2,lty=10)
# #points(Time,onepatchRT[,9],type="l",lwd=2,lty=10)
# #legend(1500,250,c("1","1.4","1.8"),lty=c(1,3,5),type=c(1,1,1),lwd=c(2,2,2),title="Mortality Ratio")

plot(Time,onepatchWT[,1]/onepatchRT[,1],type="l",lwd=2,ylab=expression(bold("Host_ratio: Wild/Robust")),xlab="",ylim=c(0,1.45),cex.lab=1.5)
points(Time,onepatchWT[,4]/onepatchRT[,4],type="l",lwd=2,lty=3)
points(Time,onepatchWT[,7]/onepatchRT[,7],type="l",lwd=2,lty=5)
#points(Time,onepatchWT[,4],type="l",lwd=2,lty=8)
points(Time,onepatchWT[,9]/onepatchRT[,9],type="l",lwd=2,lty=10)
# legend(1500,3000,c("1","1.4","1.8"),lty=c(1,3,5),lwd=c(2,2,2),title="Mortality Ratio",bty="n"),type="l",lwd=2,main="Wild/Robust",xlab="",ylab="")


plot(Time,onepatchinfect[,1],type="l",lwd=2,ylab=expression(bold("Infected_total: Robust+Wild")),xlab="",cex.lab=1.5)
points(Time,onepatchinfect[,4],type="l",lwd=2,lty=3)
points(Time,onepatchinfect[,7],type="l",lwd=2,lty=5)
#points(Time,onepatchinfect[,4],type="l",lwd=2,lty=8)
points(Time,onepatchinfect[,9],type="l",lwd=2,lty=10)
#legend(1500,250,c("1","1.4","1.8"),lty=c(1,3,5),type=c(1,1,1),lwd=c(2,2,2),title="Mortality Ratio")

plot(Time,onepatchinfectWT[,1]/onepatchinfectRT[,1],type="l",lwd=2,ylab=expression(bold("Infect_ratio: Wild/Robust")),xlab="",ylim=c(0,1.5),cex.lab=1.5)
points(Time,onepatchinfectWT[,4]/onepatchinfectRT[,4],type="l",lwd=2,lty=3)
points(Time,onepatchinfectWT[,7]/onepatchinfectRT[,7],type="l",lwd=2,lty=5)
#points(Time,onepatchWT[,4],type="l",lwd=2,lty=8)
points(Time,onepatchinfectWT[,9]/onepatchinfectRT[,9],type="l",lwd=2,lty=10)
# legend(1500,3000,c("1","1.4","1.8"),lty=c(1,3,5),lwd=c(2,2,2),title="Mortality Ratio",bty="n")




dev.off()



