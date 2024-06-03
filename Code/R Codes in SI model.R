library(ggplot2)
library(reshape)

setwd("C:\\Users\\lucyj\\Downloads")


Time<-0:60

onepatchinfect<-read.csv("SImodel_onepatch_infected_total_15gene.csv",header=FALSE)
onepatchRT<-read.csv("SImodel_onepatch_RT_total_15gene.csv",header=FALSE)
onepatchWT<-read.csv("SImodel_onepatch_WT_total_15gene.csv",header=FALSE)
onepatchinfectRT<-read.csv("SImodel_onepatch_infected_RT_15gene.csv",header=FALSE)
onepatchinfectWT<-read.csv("SImodel_onepatch_infected_WT_15gene.csv",header=FALSE)

theme_set(theme_bw(20))



###FIGURE 1

tiff("onepatchSImodel_1_3_20.tiff", width=10,height=10, units='in',res=600)
par(mfrow=c(3,1))
par(mar=c(2,5,1,4.5))

plot(Time,onepatchWT[,1]/onepatchRT[,1],type="l",lwd=2,xlim=c(8,60),ylab=expression(bold("Host_ratio: Wild/Robust")),cex.axis=1.5,xlab="",cex.lab=1.5,xaxt="n",ylim=c(1.19,1.389)) #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,onepatchWT[,2]/onepatchRT[,2],type="l",lwd=2,lty=2,col="black")
points(Time,onepatchWT[,3]/onepatchRT[,3],type="l",lwd=2,lty=3,col="black")
points(Time,onepatchWT[,4]/onepatchRT[,4],type="l",lwd=2,lty=6,col="red")
legend(15,1.34,c("1","1.5","2"),lty=c(1,2,3),col=c("black","black","black"),lwd=c(2,2,2),title="Mortality Ratio",bty="n",cex=2)


plot(Time,onepatchWT[,1]+onepatchRT[,1],type="l",lwd=2,ylim=c(1490,1920),xlim=c(8,60),cex.axis=1.5,ylab=expression(bold("Host_total: Robust+Wild")),xlab="",cex.lab=1.5,xaxt="n") #ylim=c(1800,2200),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,onepatchWT[,2]+onepatchRT[,2],type="l",lwd=2,lty=2,col="black")
points(Time,onepatchWT[,3]+onepatchRT[,3],type="l",lwd=2,lty=3,col="black")
#points(Time,onepatchWT[,4]+onepatchRT[,4],type="l",lwd=2,lty=8)
points(Time,onepatchWT[,4]+onepatchRT[,4],type="l",lwd=2,lty=6,col="red")


plot(Time,onepatchinfect[,1],type="l",lwd=2,xlim=c(8,60),cex.axis=1.5,ylab=expression(bold("Infected_total: Robust+Wild")),xlab="",cex.lab=1.5,ylim=c(0.8,300),xaxt="n")
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,onepatchinfect[,2],type="l",lwd=2,lty=2,col="black")
points(Time,onepatchinfect[,3],type="l",lwd=2,lty=3,col="black")
#points(Time,onepatchinfect[,4],type="l",lwd=2,lty=8)
points(Time,onepatchinfect[,4],type="l",lwd=2,lty=6,col="red")
#legend(1500,250,c("1","1.4","1.8"),lty=c(1,3,5),type=c(1,1,1),lwd=c(2,2,2),title="Mortality Ratio")


dev.off()



###FIGURE S1
tiff("onepatchSImodel_wild_robust_1_3_20.tiff", width=10,height=10, units='in',res=600)
par(mfrow=c(2,1))
par(mar=c(2,5,1,4.5))

plot(Time,onepatchWT[,1],type="l",lwd=2,xlim=c(8,60),ylab=expression(bold("Wild-type Size")),cex.axis=1.5,xlab="",ylim=c(810,1130),cex.lab=1.5,xaxt="n") #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,onepatchWT[,2],type="l",lwd=2,lty=2,col="black")
points(Time,onepatchWT[,3],type="l",lwd=2,lty=3,col="black")
legend(8,1000,c("1","1.5","2"),lty=c(1,2,3),col=c("black","black","black"),lwd=c(2,2,2),title="Mortality Ratio",bty="n",cex=2)


plot(Time,onepatchRT[,1],type="l",lwd=2,xlim=c(8,60),ylab=expression(bold("Robust-type Size")),ylim=c(610,810),cex.axis=1.5,xlab="",cex.lab=1.5,xaxt="n") #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,onepatchRT[,2],type="l",lwd=2,lty=2,col="black")
points(Time,onepatchRT[,3],type="l",lwd=2,lty=3,col="black")

dev.off()


################################################################################
################################################################################
#####Figures 2-6 and S2-S7, S9-S10#
# 
matr<-read.csv("2geno_5patch_transient_matrixa-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)

matr_s<-read.csv("2geno_5patch_transient_matrixa-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)


N=5 #patch number
K=3000 # carrying capacity
h=1 #repetition
mo=seq(1,2,0.01)
mig=seq(0,0.05,0.0005)
###index 1: X: geno wild/robust calculation at i) migration level; ii) mortality ratio level;
X<-matrix(NA,length(mig),length(mo)) # average geno 1
Y<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality
I<-matrix(NA,length(mig),length(mo)) # average total infected host individuals

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    I1<-0
    
    for(k in 1:length(h)) # repeat
    {
    
    X1<-sum(X1,matr[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
    Y1<-sum(Y1,matr[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
    I1<-sum(I1,matr[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),3])
    }
    X[i,j]<-X1
    Y[i,j]<-Y1
    I[i,j]<-I1
    
  }
  
}
  


X_s<-matrix(NA,length(mig),length(mo)) # average geno 1
Y_s<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality
I_s<-matrix(NA,length(mig),length(mo)) # average total infected host individuals

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1_s<-0
    Y1_s<-0
    I1_s<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1_s<-sum(X1_s,matr[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1_s<-sum(Y1_s,matr[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      I1_s<-sum(I1_s,matr[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),3])
    }
    X_s[i,j]<-X1_s
    Y_s[i,j]<-Y1_s
    I_s[i,j]<-I1_s
    
  }
  
}
#####average patch densities of X and Y
X1<-X/(length(h)*N) # 5 is patch number as well as repeating number
Y1<-Y/(length(h)*N)
I1<-I/(length(h)*N)
    

#no log
geno_ratio<-(Y1/X1)
tot_pop<-((X1+Y1)/K) # 3000 is carrying capacity


XX<-melt(geno_ratio)
YY<-melt(tot_pop)
II<-melt(I1)
ZZ<-melt(X1/K)
KK<-melt(Y1/K)


######for the scenario when disease started from one local patch
X1_s<-X_s/(length(h)*N) # 5 is patch number as well as repeating number
Y1_s<-Y_s/(length(h)*N)
I1_s<-I_s/(length(h)*N)


#no log
geno_ratio_s<-(Y1_s/X1_s)
tot_pop_s<-((X1_s+Y1_s)/K) # 3000 is carrying capacity


XX_s<-melt(geno_ratio_s)
YY_s<-melt(tot_pop_s)
II_s<-melt(I1_s)
ZZ_s<-melt(X1_s/K)
KK_s<-melt(Y1_s/K)


#mig=0:0.0005:0.005;mo=1:0.15:2.5;
x_axis<-mig
y_axis<-mo


theme_set(theme_bw(20))

####FIGURE 2

tiff("heatmap_2geno5patch_transient_matrixa_12-25-19-15generations_allpatches.tiff", width=10,height=10, units='in',res=600)



v <- ggplot(XX, aes(x_axis[X1], y_axis[X2], z = value))
v + geom_contour()

v + geom_contour(aes(colour = stat(level)))
v + geom_contour(colour = "red")
v + geom_raster(aes(fill = value)) +guides(fill = guide_colorbar(title="Wild / Robust"))+theme(legend.title = element_text(size=12))+
  geom_contour(colour = "black")+xlab("Migration rate of wild type")+ylab(expression(paste('Mortality ratio ',alpha[W]/alpha[R])))+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                                       
                                                       midpoint=median(XX[,3],na.rm = TRUE), space = "rgb", na.value = "grey50", guide = 
                                                         "colourbar")+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0, 0))                                                   
                                                         
                                                                                                                          



dev.off()  





####FIGURE 3
tiff("heatmap_2geno5patch_transient_tot_pop_matrixa_12-25-19-15generations_allpatches_highacc_highdis.tiff", width=10,height=10, units='in',res=600)


v <- ggplot(YY, aes(x_axis[X1], y_axis[X2], z = value))
v + geom_contour()

v + geom_contour(aes(colour = stat(level)))
v + geom_contour(colour = "red")
v + geom_raster(aes(fill = value)) +guides(fill = guide_colorbar(title="Tot_Pop"))+theme(legend.title = element_text(size=12))+
  geom_contour(colour = "black")+xlab("Migration rate of wild type")+ylab(expression(paste('Mortality ratio ',alpha[W]/alpha[R])))+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                                                                                                                                         
                                                                                                                                                         midpoint=median(YY[,3]), space = "rgb", na.value = "grey50", guide = 
                                                                                                                                                           "colourbar")+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0, 0))                                                   



dev.off()  



#FIGURE 4

tiff("heatmap_2geno5patch_transient_matrixa_12-25-19-15generations_spread.tiff", width=10,height=10, units='in',res=600)


v <- ggplot(XX_s, aes(x_axis[X1], y_axis[X2], z = value))
v + geom_contour()

v + geom_contour(aes(colour = stat(level)))
v + geom_contour(colour = "red")
v + geom_raster(aes(fill = value)) +guides(fill = guide_colorbar(title="Wild / Robust"))+theme(legend.title = element_text(size=12))+
  geom_contour(colour = "black")+xlab("Migration rate of wild type")+ylab(expression(paste('Mortality ratio ',alpha[W]/alpha[R])))+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                                                                                                                                         
                                                                                                                                                         midpoint=median(XX_s[,3],na.rm = TRUE), space = "rgb", na.value = "grey50", guide = 
                                                                                                                                                           "colourbar")+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0, 0))                                                   



dev.off()  




###FIGURE 5
tiff("heatmap_2geno5patch_transient_tot_pop_matrixa_12-25-19-15generations_spread_highacc_highdis.tiff", width=10,height=10, units='in',res=600)

v <- ggplot(YY_s, aes(x_axis[X1], y_axis[X2], z = value))
v + geom_contour()

v + geom_contour(aes(colour = stat(level)))
v + geom_contour(colour = "red")
v + geom_raster(aes(fill = value)) +guides(fill = guide_colorbar(title="Tot_Pop"))+theme(legend.title = element_text(size=12))+
  geom_contour(colour = "black")+xlab("Migration rate of wild type")+ylab(expression(paste('Mortality ratio ',alpha[W]/alpha[R])))+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                                                                                                                                         
                                                                                                                                                         midpoint=median(YY_s[,3]), space = "rgb", na.value = "grey50", guide = 
                                                                                                                                                           "colourbar")+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0, 0))                                                   



dev.off()  




##FIGURE S2
tiff("heatmap_2geno5patch_Wild_matrixa_15generation_allpatches_12-25-19.tiff", width=10,height=10, units='in',res=200)


v <- ggplot(KK, aes(x_axis[X1], y_axis[X2], z = value/max(1)))
v + geom_contour()

v + geom_contour(aes(colour = stat(level)))
v + geom_contour(colour = "red")
v + geom_raster(aes(fill = value/max(1))) +guides(fill = guide_colorbar(title="Stand_Wild"))+theme(legend.title = element_text(size=12))+
  geom_contour(colour = "black")+xlab("Migration rate of wild type")+ylab(expression(paste('Mortality ratio ',alpha[W]/alpha[R])))+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                                                                                                                                         
                                                                                                                                                         midpoint=median(KK[,3]/max(1)), space = "rgb", na.value = "grey50", guide = 
                                                                                                                                                           
                                                                                                                                                           "colourbar")+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0, 0))                                                   

dev.off()  


###FIGURE S3


tiff("heatmap_2geno5patch_Robust_matrixa_15generation_allpatches_12-25-19.tiff", width=10,height=10, units='in',res=600)



v <- ggplot(ZZ, aes(x_axis[X1], y_axis[X2], z = value))
v + geom_contour()

v + geom_contour(aes(colour = stat(level)))
v + geom_contour(colour = "red")
v + geom_raster(aes(fill = value)) +guides(fill = guide_colorbar(title="Stand_Robust"))+theme(legend.title = element_text(size=12))+
  geom_contour(colour = "black")+xlab("Migration rate of wild type")+ylab(expression(paste('Mortality ratio ',alpha[W]/alpha[R])))+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                                                                                                                                         
                                                                                                                                                         midpoint=median(ZZ[,3],na.rm = TRUE), space = "rgb", na.value = "grey50", guide = 
                                                                                                                                                           "colourbar")+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0, 0))                                                   


dev.off()  



###FIGURE S4


tiff("heatmap_2geno5patch_tot_infected_matrixa_15generation_allpatches_12-25-19.tiff", width=10,height=10, units='in',res=200)


v <- ggplot(II, aes(x_axis[X1], y_axis[X2], z = value/K))
v + geom_contour()

v + geom_contour(aes(colour = stat(level)))
v + geom_contour(colour = "red")
v + geom_raster(aes(fill = value/K)) +guides(fill = guide_colorbar(title="Tot_Infected"))+theme(legend.title = element_text(size=12))+
  xlab("Migration rate of wild type")+ylab(expression(paste('Mortality ratio ',alpha[W]/alpha[R])))+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                                                                                                          
                                                                                                                          midpoint=median(II[,3]/K), space = "rgb", na.value = "grey50", guide = 
                                                                                                                            "colourbar")+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0, 0))                                                   

dev.off()




#FIGURE S5


tiff("heatmap_2geno5patch_Wild_matrixa_15generation_spread_12-25-19.tiff", width=10,height=10, units='in',res=200)

v <- ggplot(KK_s, aes(x_axis[X1], y_axis[X2], z = value/max(1)))
v + geom_contour()

v + geom_contour(aes(colour = stat(level)))
v + geom_contour(colour = "red")
v + geom_raster(aes(fill = value/max(1))) +guides(fill = guide_colorbar(title="Stand_Wild"))+theme(legend.title = element_text(size=12))+
  geom_contour(colour = "black")+xlab("Migration rate of wild type")+ylab(expression(paste('Mortality ratio ',alpha[W]/alpha[R])))+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                                                                                                                                         
                                                                                                                                                         midpoint=median(KK_s[,3]/max(1)), space = "rgb", na.value = "grey50", guide = 
                                                                                                                                                           
                                                                                                                                                           "colourbar")+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0, 0))                                                   

dev.off()  



###FIGURE S6
tiff("heatmap_2geno5patch_Robust_matrixa_15generation_spread_12-25-19.tiff", width=10,height=10, units='in',res=600)

v <- ggplot(ZZ_s, aes(x_axis[X1], y_axis[X2], z = value))
v + geom_contour()

v + geom_contour(aes(colour = stat(level)))
v + geom_contour(colour = "red")
v + geom_raster(aes(fill = value)) +guides(fill = guide_colorbar(title="Stand_Robust"))+theme(legend.title = element_text(size=12))+
  geom_contour(colour = "black")+xlab("Migration rate of wild type")+ylab(expression(paste('Mortality ratio ',alpha[W]/alpha[R])))+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                                                                                                                                         
                                                                                                                                                         midpoint=median(ZZ_s[,3],na.rm = TRUE), space = "rgb", na.value = "grey50", guide = 
                                                                                                                                                           "colourbar")+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0, 0))                                                   


dev.off() 




##FIGURE S7

tiff("heatmap_2geno5patch_tot_infected_matrixa_15generation_spread_12-25-19.tiff", width=10,height=10, units='in',res=200)

#

v <- ggplot(II_s, aes(x_axis[X1], y_axis[X2], z = value/K))
v + geom_contour()

v + geom_contour(aes(colour = stat(level)))
v + geom_contour(colour = "red")
v + geom_raster(aes(fill = value/K)) +guides(fill = guide_colorbar(title="Tot_Infected"))+theme(legend.title = element_text(size=12))+
  geom_contour(colour = "black")+xlab("Migration rate of wild type")+ylab(expression(paste('Mortality ratio ',alpha[W]/alpha[R])))+ scale_fill_gradient2(low = 'steelblue', mid='white',high = 'red' ,
                                                                                                                                                         
                                                                                                                                                         midpoint=median(II_s[,3]/K), space = "rgb", na.value = "grey50", guide = 
                                                                                                                                                           "colourbar")+scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0, 0))                                                   

dev.off()






############################################################################################################################
############################################################################################################################
##FIGURE_6###################################

#####plots for topology and migration


####NOTATION#######################################
#all matra to matrk collected the topological data from a to k based on Fig. S2 when disease started from all patches
#all matras to matrks collected the topological data from a to k based on Fig.S2 when disease started from one patch

##two factors: migration and mortality ratio as x and y axes
x_axis<-mig
y_axis<-mo


matras<-read.csv("2geno_5patch_transient_matrixa-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)
matrbs<-read.csv("2geno_5patch_transient_matrixb-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)
matrcs<-read.csv("2geno_5patch_transient_matrixc-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)
matrds<-read.csv("2geno_5patch_transient_matrixd-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)
matres<-read.csv("2geno_5patch_transient_matrixe-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)
matrfs<-read.csv("2geno_5patch_transient_matrixf-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)
matrgs<-read.csv("2geno_5patch_transient_matrixg-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)
matrhs<-read.csv("2geno_5patch_transient_matrixh-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)
matris<-read.csv("2geno_5patch_transient_matrixi-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)
matrjs<-read.csv("2geno_5patch_transient_matrixj-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)
matrks<-read.csv("2geno_5patch_transient_matrixk-12-25-19_15generations_spread_1_2_0_0.05_highdisease.csv",header=FALSE)



matra<-read.csv("2geno_5patch_transient_matrixa-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)
matrb<-read.csv("2geno_5patch_transient_matrixb-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)
matrc<-read.csv("2geno_5patch_transient_matrixc-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)
matrd<-read.csv("2geno_5patch_transient_matrixd-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)
matre<-read.csv("2geno_5patch_transient_matrixe-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)
matrf<-read.csv("2geno_5patch_transient_matrixf-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)
matrg<-read.csv("2geno_5patch_transient_matrixg-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)
matrh<-read.csv("2geno_5patch_transient_matrixh-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)
matri<-read.csv("2geno_5patch_transient_matrixi-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)
matrj<-read.csv("2geno_5patch_transient_matrixj-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)
matrk<-read.csv("2geno_5patch_transient_matrixk-12-25-19_15generations_allpatches_1_2_0_0.05_highdisease.csv",header=FALSE)

####for topology a
###index 1: Xa or xas: geno wild/robust calculation at i) migration level; ii) mortality ratio level;
Xa<-matrix(NA,length(mig),length(mo)) # average geno 1
Ya<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matra[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matra[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xa[i,j]<-X1
    Ya[i,j]<-Y1
    
  }
  
}



#####average patch densities of wildtype (X) and robust (Y)
X1a<-Xa/(length(h)*N) # 5 is patch number as well as repeating number
Y1a<-Ya/(length(h)*N)

####index 
geno_ratioa<-Y1a/X1a
tot_popa<-(X1a+Y1a)/3000 # 3000 is carrying capacity

XXa<-melt(geno_ratioa)
YYa<-melt(tot_popa)
ZZa<-melt(X1a/3000)
KKa<-melt(Y1a/3000)



Xb<-matrix(NA,length(mig),length(mo)) # average geno 1: wild type
Yb<-matrix(NA,length(mig),length(mo)) # average geno 2 (robust) with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrb[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrb[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xb[i,j]<-X1
    Yb[i,j]<-Y1
    
  }
  
}


###for topology b
X1b<-Xb/(length(h)*N) # 5 is patch number as well as repeating number
Y1b<-Yb/(length(h)*N)

####index 
geno_ratiob<-Y1b/X1b
tot_popb<-(X1b+Y1b)/3000 # 3000 is carrying capacity

XXb<-melt(geno_ratiob)
YYb<-melt(tot_popb)
ZZb<-melt(X1b/3000)
KKb<-melt(Y1b/3000)



#### for topology c
Xc<-matrix(NA,length(mig),length(mo)) # average geno 1
Yc<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrc[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrc[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xc[i,j]<-X1
    Yc[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1c<-Xc/(length(h)*N) # 5 is patch number as well as repeating number
Y1c<-Yc/(length(h)*N)
####index 
geno_ratioc<-Y1c/X1c
tot_popc<-(X1c+Y1c)/3000 # 3000 is carrying capacity

XXc<-melt(geno_ratioc)
YYc<-melt(tot_popc)
ZZc<-melt(X1c/3000)
KKc<-melt(Y1c/3000)



####for topology d

Xd<-matrix(NA,length(mig),length(mo)) # average geno 1
Yd<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrd[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrd[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xd[i,j]<-X1
    Yd[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1d<-Xd/(length(h)*N) # 5 is patch number as well as repeating number
Y1d<-Yd/(length(h)*N)
####index 
geno_ratiod<-Y1d/X1d
tot_popd<-(X1d+Y1d)/3000 # 3000 is carrying capacity

XXd<-melt(geno_ratiod)
YYd<-melt(tot_popd)
ZZd<-melt(X1d/3000)
KKd<-melt(Y1d/3000)



###for topology e
Xe<-matrix(NA,length(mig),length(mo)) # average geno 1
Ye<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matre[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matre[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xe[i,j]<-X1
    Ye[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1e<-Xe/(length(h)*N) # 5 is patch number as well as repeating number
Y1e<-Ye/(length(h)*N)

####index 
geno_ratioe<-Y1e/X1e
tot_pope<-(X1e+Y1e)/3000 # 3000 is carrying capacity

XXe<-melt(geno_ratioe)
YYe<-melt(tot_pope)
ZZe<-melt(X1e/3000)
KKe<-melt(Y1e/3000)


#for topology e
Xf<-matrix(NA,length(mig),length(mo)) # average geno 1
Yf<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrf[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrf[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xf[i,j]<-X1
    Yf[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1f<-Xf/(length(h)*N) # 5 is patch number as well as repeating number
Y1f<-Yf/(length(h)*N)
####index 
geno_ratiof<-Y1f/X1f
tot_popf<-(X1f+Y1f)/3000 # 3000 is carrying capacity

XXf<-melt(geno_ratiof)
YYf<-melt(tot_popf)
ZZf<-melt(X1f/3000)
KKf<-melt(Y1f/3000)


##for topology g
Xg<-matrix(NA,length(mig),length(mo)) # average geno 1
Yg<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrg[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrg[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xg[i,j]<-X1
    Yg[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1g<-Xg/(length(h)*N) # 5 is patch number as well as repeating number
Y1g<-Yg/(length(h)*N)
####index 
geno_ratiog<-Y1g/X1g
tot_popg<-(X1g+Y1g)/3000 # 3000 is carrying capacity

XXg<-melt(geno_ratiog)
YYg<-melt(tot_popg)
ZZg<-melt(X1g/3000)
KKg<-melt(Y1g/3000)



###for topology h
Xh<-matrix(NA,length(mig),length(mo)) # average geno 1
Yh<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrh[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrh[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xh[i,j]<-X1
    Yh[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1h<-Xh/(length(h)*N) # 5 is patch number as well as repeating number
Y1h<-Yh/(length(h)*N)
####index 
geno_ratioh<-Y1h/X1h
tot_poph<-(X1h+Y1h)/3000 # 3000 is carrying capacity

XXh<-melt(geno_ratioh)
YYh<-melt(tot_poph)
ZZh<-melt(X1h/3000)
KKh<-melt(Y1h/3000)



###for topology i
Xi<-matrix(NA,length(mig),length(mo)) # average geno 1
Yi<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matri[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matri[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xi[i,j]<-X1
    Yi[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1i<-Xi/(length(h)*N) # 5 is patch number as well as repeating number
Y1i<-Yi/(length(h)*N)
####index 
geno_ratioi<-Y1i/X1i
tot_popi<-(X1i+Y1i)/3000 # 3000 is carrying capacity

XXi<-melt(geno_ratioi)
YYi<-melt(tot_popi)
ZZi<-melt(X1i/3000)
KKi<-melt(Y1i/3000)



#for topology j
Xj<-matrix(NA,length(mig),length(mo)) # average geno 1
Yj<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrj[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrj[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xj[i,j]<-X1
    Yj[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1j<-Xj/(length(h)*N) # 5 is patch number as well as repeating number
Y1j<-Yj/(length(h)*N)

####index 
geno_ratioj<-Y1j/X1j
tot_popj<-(X1j+Y1j)/3000 # 3000 is carrying capacity

XXj<-melt(geno_ratioj)
YYj<-melt(tot_popj)
ZZj<-melt(X1j/3000)
KKj<-melt(Y1j/3000)


#for topology k
Xk<-matrix(NA,length(mig),length(mo)) # average geno 1
Yk<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrk[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrk[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xk[i,j]<-X1
    Yk[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1k<-Xk/(length(h)*N) # 5 is patch number as well as repeating number
Y1k<-Yk/(length(h)*N)
####index 
geno_ratiok<-Y1k/X1k
tot_popk<-(X1k+Y1k)/3000 # 3000 is carrying capacity

XXk<-melt(geno_ratiok)
YYk<-melt(tot_popk)
ZZk<-melt(X1k/3000)
KKk<-melt(Y1k/3000)



####group geno ratio by edge number
### edge 5: a
### edge 6: b, c, d
### edge 7: e,f,g
### edge 8: h,i
### edge 9: j
### edge 10:k

geno_ratio_5<-XXa  #the ratio of numbers of wild type and robust
geno1_5<-ZZa   #the number of wild_type genotype in topologies with edge number 5
geno2_5<-KKa   #the number of robust genotype in topologies with edge number 5

tot_pop_5<-YYa #total population of both wild type and robust with edge number 5

geno_ratio_6<-XXb
geno_ratio_6[,3]<-(XXb[,3]+XXc[,3]+XXd[,3])/3
geno1_6<-ZZb
geno1_6[,3]<-(ZZb[,3]+ZZc[,3]+ZZd[,3])/3
geno2_6<-KKb
geno2_6[,3]<-(KKb[,3]+KKc[,3]+KKd[,3])/3


tot_pop_6<-YYb
tot_pop_6[,3]<-(YYb[,3]+YYc[,3]+YYd[,3])/3



geno_ratio_7<-XXe
geno_ratio_7[,3]<-(XXe[,3]+XXf[,3]+XXg[,3])/3
geno1_7<-ZZe
geno1_7[,3]<-(ZZe[,3]+ZZf[,3]+ZZg[,3])/3
geno2_7<-KKe
geno2_7[,3]<-(KKe[,3]+KKf[,3]+KKg[,3])/3


tot_pop_7<-YYe
tot_pop_7[,3]<-(YYe[,3]+YYf[,3]+YYg[,3])/3

geno_ratio_8<-XXh
geno_ratio_8[,3]<-(XXh[,3]+XXi[,3])/2
geno1_8<-ZZh
geno1_8[,3]<-(ZZh[,3]+ZZi[,3])/2
geno2_8<-KKh
geno2_8[,3]<-(KKh[,3]+KKi[,3])/2

tot_pop_8<-YYh
tot_pop_8[,3]<-(YYh[,3]+YYi[,3])/2


geno_ratio_9<-XXj
geno1_9<-ZZj
geno2_9<-KKj

tot_pop_9<-YYj

geno_ratio_10<-XXk
geno1_10<-ZZk
geno2_10<-KKk

tot_pop_10<-YYk

##pick up the above results at certain level of mortality ratio across edge numbers: I picked up the row 50 in my case
geno_ratio_5_high<-geno_ratio_5[(length(mig)*50+1):(length(mig)*51),]

tot_pop_5_high<-tot_pop_5[(length(mig)*50+1):(length(mig)*51),]

geno_ratio_6_high<-geno_ratio_6[(length(mig)*50+1):(length(mig)*51),]


tot_pop_6_high<-tot_pop_6[(length(mig)*50+1):(length(mig)*51),]


geno_ratio_7_high<-geno_ratio_7[(length(mig)*50+1):(length(mig)*51),]


tot_pop_7_high<-tot_pop_7[(length(mig)*50+1):(length(mig)*51),]

geno_ratio_8_high<-geno_ratio_8[(length(mig)*50+1):(length(mig)*51),]

tot_pop_8_high<-tot_pop_8[(length(mig)*50+1):(length(mig)*51),]

geno_ratio_9_high<-geno_ratio_9[(length(mig)*50+1):(length(mig)*51),]

tot_pop_9_high<-tot_pop_9[(length(mig)*50+1):(length(mig)*51),]

geno_ratio_10_high<-geno_ratio_10[(length(mig)*50+1):(length(mig)*51),]

tot_pop_10_high<-tot_pop_10[(length(mig)*50+1):(length(mig)*51),]




#####################################################################
#########################################################################
########################################################################
##for the density and total population across different topological structures when disease started from one local patch and spreads to other places
##all Xas indicates the number of wild type at topology a
Xas<-matrix(NA,length(mig),length(mo)) # average geno 1
Yas<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matras[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matras[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xas[i,j]<-X1
    Yas[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1as<-Xas/(length(h)*N) # 5 is patch number as well as repeating number
Y1as<-Yas/(length(h)*N)

####index 
geno_ratioas<-Y1as/X1as
tot_popas<-(X1as+Y1as)/3000 # 3000 is carrying capacity

XXas<-melt(geno_ratioas)
YYas<-melt(tot_popas)
ZZas<-melt(X1as/3000)
KKas<-melt(Y1as/3000)


###topology b
Xbs<-matrix(NA,length(mig),length(mo)) # average geno 1
Ybs<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrbs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrbs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xbs[i,j]<-X1
    Ybs[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1bs<-Xbs/(length(h)*N) # 5 is patch number as well as repeating number
Y1bs<-Ybs/(length(h)*N)

####index 
geno_ratiobs<-Y1bs/X1bs
tot_popbs<-(X1bs+Y1bs)/3000 # 3000 is carrying capacity

XXbs<-melt(geno_ratiobs)
YYbs<-melt(tot_popbs)
ZZbs<-melt(X1bs/3000)
KKbs<-melt(Y1bs/3000)



#topology c
Xcs<-matrix(NA,length(mig),length(mo)) # average geno 1
Ycs<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrcs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrcs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xcs[i,j]<-X1
    Ycs[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1cs<-Xcs/(length(h)*N) # 5 is patch number as well as repeating number
Y1cs<-Ycs/(length(h)*N)
####index 
geno_ratiocs<-Y1cs/X1cs
tot_popcs<-(X1cs+Y1cs)/3000 # 3000 is carrying capacity

XXcs<-melt(geno_ratiocs)
YYcs<-melt(tot_popcs)
ZZcs<-melt(X1cs/3000)
KKcs<-melt(Y1cs/3000)


##topology d
Xds<-matrix(NA,length(mig),length(mo)) # average geno 1
Yds<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrds[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrds[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xds[i,j]<-X1
    Yds[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1ds<-Xds/(length(h)*N) # 5 is patch number as well as repeating number
Y1ds<-Yds/(length(h)*N)
####index 
geno_ratiods<-Y1ds/X1ds
tot_popds<-(X1ds+Y1ds)/3000 # 3000 is carrying capacity

XXds<-melt(geno_ratiods)
YYds<-melt(tot_popds)
ZZds<-melt(X1ds/3000)
KKds<-melt(Y1ds/3000)


##topology e
Xes<-matrix(NA,length(mig),length(mo)) # average geno 1
Yes<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matres[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matres[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xes[i,j]<-X1
    Yes[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1es<-Xes/(length(h)*N) # 5 is patch number as well as repeating number
Y1es<-Yes/(length(h)*N)

####index 
geno_ratioes<-Y1es/X1es
tot_popes<-(X1es+Y1es)/3000 # 3000 is carrying capacity

XXes<-melt(geno_ratioes)
YYes<-melt(tot_popes)
ZZes<-melt(X1es/3000)
KKes<-melt(Y1es/3000)


#topology f
Xfs<-matrix(NA,length(mig),length(mo)) # average geno 1
Yfs<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrfs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrfs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xfs[i,j]<-X1
    Yfs[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1fs<-Xfs/(length(h)*N) # 5 is patch number as well as repeating number
Y1fs<-Yfs/(length(h)*N)
####index 
geno_ratiofs<-Y1fs/X1fs
tot_popfs<-(X1fs+Y1fs)/3000 # 3000 is carrying capacity

XXfs<-melt(geno_ratiofs)
YYfs<-melt(tot_popfs)
ZZfs<-melt(X1fs/3000)
KKfs<-melt(Y1fs/3000)


#topology g
Xgs<-matrix(NA,length(mig),length(mo)) # average geno 1
Ygs<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrgs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrgs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xgs[i,j]<-X1
    Ygs[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1gs<-Xgs/(length(h)*N) # 5 is patch number as well as repeating number
Y1gs<-Ygs/(length(h)*N)
####index 
geno_ratiogs<-Y1gs/X1gs
tot_popgs<-(X1gs+Y1gs)/3000 # 3000 is carrying capacity

XXgs<-melt(geno_ratiogs)
YYgs<-melt(tot_popgs)
ZZgs<-melt(X1gs/3000)
KKgs<-melt(Y1gs/3000)


#topology h
Xhs<-matrix(NA,length(mig),length(mo)) # average geno 1
Yhs<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrhs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrhs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xhs[i,j]<-X1
    Yhs[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1hs<-Xhs/(length(h)*N) # 5 is patch number as well as repeating number
Y1hs<-Yhs/(length(h)*N)
####index 
geno_ratiohs<-Y1hs/X1hs
tot_pophs<-(X1hs+Y1hs)/3000 # 3000 is carrying capacity

XXhs<-melt(geno_ratiohs)
YYhs<-melt(tot_pophs)
ZZhs<-melt(X1hs/3000)
KKhs<-melt(Y1hs/3000)


#topology i
Xis<-matrix(NA,length(mig),length(mo)) # average geno 1
Yis<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matris[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matris[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xis[i,j]<-X1
    Yis[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1is<-Xis/(length(h)*N) # 5 is patch number as well as repeating number
Y1is<-Yis/(length(h)*N)
####index 
geno_ratiois<-Y1is/X1is
tot_popis<-(X1is+Y1is)/3000 # 3000 is carrying capacity

XXis<-melt(geno_ratiois)
YYis<-melt(tot_popis)
ZZis<-melt(X1is/3000)
KKis<-melt(Y1is/3000)



#topology j
Xjs<-matrix(NA,length(mig),length(mo)) # average geno 1
Yjs<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrjs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrjs[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xjs[i,j]<-X1
    Yjs[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1js<-Xjs/(length(h)*N) # 5 is patch number as well as repeating number
Y1js<-Yjs/(length(h)*N)

####index 
geno_ratiojs<-Y1js/X1js
tot_popjs<-(X1js+Y1js)/3000 # 3000 is carrying capacity

XXjs<-melt(geno_ratiojs)
YYjs<-melt(tot_popjs)
ZZjs<-melt(X1js/3000)
KKjs<-melt(Y1js/3000)


#topology k
Xks<-matrix(NA,length(mig),length(mo)) # average geno 1
Yks<-matrix(NA,length(mig),length(mo)) # average geno2 with larger mortality

for(i in 1:length(mig)) #ii
{
  for(j in 1:length(mo)) #i
  {
    X1<-0
    Y1<-0
    
    for(k in 1:length(h)) # repeat
    {
      
      X1<-sum(X1,matrks[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),1])
      Y1<-sum(Y1,matrks[(i+(j-1)*length(mig)+(k-1)*length(mo)*length(mig)),2])
      
    }
    Xks[i,j]<-X1
    Yks[i,j]<-Y1
    
  }
  
}


#####average patch densities of X and Y
X1ks<-Xks/(length(h)*N) # 5 is patch number as well as repeating number
Y1ks<-Yks/(length(h)*N)
####index 
geno_ratioks<-Y1ks/X1ks
tot_popks<-(X1ks+Y1ks)/3000 # 3000 is carrying capacity

XXks<-melt(geno_ratioks)
YYks<-melt(tot_popks)
ZZks<-melt(X1ks/3000)
KKks<-melt(Y1ks/3000)



####group geno ratio by edge number
### edge 5: a
### edge 6: b, c, d
### edge 7: e,f,g
### edge 8: h,i
### edge 9: j
### edge 10:k


geno_ratio_5s<-XXas #number ratio of wild type vs robust with edge number 5 when disease started from one local patch
geno1_5s<-ZZas #number of wild type with edge number 5 when disease started from one local patch
geno2_5s<-KKas #number of robust with edge number 5 when disease started from one local patch

tot_pop_5s<-YYas #total population size (wild + robust) with edge number 5 when disease started from one local patch

geno_ratio_6s<-XXbs
geno_ratio_6s[,3]<-(XXbs[,3]+XXcs[,3]+XXds[,3])/3
geno1_6s<-ZZbs
geno1_6s[,3]<-(ZZbs[,3]+ZZcs[,3]+ZZds[,3])/3
geno2_6s<-KKbs
geno2_6s[,3]<-(KKbs[,3]+KKcs[,3]+KKds[,3])/3


tot_pop_6s<-YYbs
tot_pop_6s[,3]<-(YYbs[,3]+YYcs[,3]+YYds[,3])/3 #all NA



geno_ratio_7s<-XXes
geno_ratio_7s[,3]<-(XXes[,3]+XXfs[,3]+XXgs[,3])/3 #all NA
geno1_7s<-ZZes
geno1_7s[,3]<-(ZZes[,3]+ZZfs[,3]+ZZgs[,3])/3
geno2_7s<-KKes
geno2_7s[,3]<-(KKes[,3]+KKfs[,3]+KKgs[,3])/3


tot_pop_7s<-YYes
tot_pop_7s[,3]<-(YYes[,3]+YYfs[,3]+YYgs[,3])/3 #NA

geno_ratio_8s<-XXhs
geno_ratio_8s[,3]<-(XXhs[,3]+XXis[,3])/2 #h is NA
geno1_8s<-ZZhs
geno1_8s[,3]<-(ZZhs[,3]+ZZis[,3])/2
geno2_8s<-KKhs
geno2_8s[,3]<-(KKhs[,3]+KKis[,3])/2

tot_pop_8s<-YYhs
tot_pop_8s[,3]<-(YYhs[,3]+YYis[,3])/2


geno_ratio_9s<-XXjs
geno1_9s<-ZZjs
geno2_9s<-KKjs

tot_pop_9s<-YYjs

geno_ratio_10s<-XXks
geno1_10s<-ZZks
geno2_10s<-KKks

tot_pop_10s<-YYks

#the ratio, number and total population across different edge numbers at one mortality ratio (here I picked up the mortality ratio at row 50)
geno_ratio_5_highs<-geno_ratio_5s[(length(mig)*50+1):(length(mig)*51),]

tot_pop_5_highs<-tot_pop_5s[(length(mig)*50+1):(length(mig)*51),]

geno_ratio_6_highs<-geno_ratio_6s[(length(mig)*50+1):(length(mig)*51),]


tot_pop_6_highs<-tot_pop_6s[(length(mig)*50+1):(length(mig)*51),]


geno_ratio_7_highs<-geno_ratio_7s[(length(mig)*50+1):(length(mig)*51),]


tot_pop_7_highs<-tot_pop_7s[(length(mig)*50+1):(length(mig)*51),]

geno_ratio_8_highs<-geno_ratio_8s[(length(mig)*50+1):(length(mig)*51),]

tot_pop_8_highs<-tot_pop_8s[(length(mig)*50+1):(length(mig)*51),]

geno_ratio_9_highs<-geno_ratio_9s[(length(mig)*50+1):(length(mig)*51),]

tot_pop_9_highs<-tot_pop_9s[(length(mig)*50+1):(length(mig)*51),]

geno_ratio_10_highs<-geno_ratio_10s[(length(mig)*50+1):(length(mig)*51),]

tot_pop_10_highs<-tot_pop_10s[(length(mig)*50+1):(length(mig)*51),]


theme_set(theme_bw(20))


tiff("Topology_graph_Fig6_disease_15generation_allpatches and one patch.tiff", width=10,height=10, units='in',res=200)

par(mfrow=c(1,2))
par(mar=c(6,3,5,2))

plot(x_axis,geno_ratio_10_high[,3],lwd=2,col="red",type="l",lty=1,xlab="",ylab="",cex.lab=2, cex.axis=2, cex.main=2,main="All patches")
points(x_axis,geno_ratio_9_high[,3],lwd=2,col="orange",type="l",lty=1)
points(x_axis,geno_ratio_8_high[,3],lwd=2,col="yellow",type="l",lty=1)
points(x_axis,geno_ratio_7_high[,3],lwd=2,col="light blue",type="l",lty=1)
points(x_axis,geno_ratio_6_high[,3],lwd=2,col="blue",type="l",lty=1)
points(x_axis,geno_ratio_5_high[,3],lwd=2,col="dark blue",type="l",lty=1)
legend(0.01,1.745,legend=c("edge=5", "edge=6","edge=7","edge=8","edge=9","edge=10"),cex=1.5,col=c("blue","steelblue","light blue","yellow","orange","red"),lty = rep(1,6),lwd=rep(2,6),box.lty=0)


plot(x_axis,geno_ratio_5_highs[,3],lwd=2,col="dark blue",type="l",lty=1,xlab="",ylab="",cex.lab=2, cex.axis=2, cex.main=2,ylim=c(1.82,1.95),main="One patch")
points(x_axis,geno_ratio_6_highs[,3],lwd=2,col="blue",type="l",lty=1)
points(x_axis,geno_ratio_7_highs[,3],lwd=2,col="light blue",type="l",lty=1)
points(x_axis,geno_ratio_8_highs[,3],lwd=2,col="yellow",type="l",lty=1)
points(x_axis,geno_ratio_9_highs[,3],lwd=2,col="orange",type="l",lty=1)
points(x_axis,geno_ratio_10_highs[,3],lwd=2,col="red",type="l",lty=1)
legend(0.0005,1.6707,legend=c("#5", "#6","#7","#8","#9","#10"),cex=1.5,col=c("dark blue","blue","light blue","yellow","orange","red"),lwd=rep(2,6),box.lty=0)



dev.off()  




avepatchinfect_matrixa_allpatch<-read.csv("ave_SImodel_5patch_infected_total_matrixa_allpatch_1_5_20.csv",header=FALSE)
avepatchRT_matrixa_allpatch<-read.csv("ave_SImodel_5patch_RT_matrixa_allpatch_1_5_20.csv",header=FALSE)
avepatchWT_matrixa_allpatch<-read.csv("ave_SImodel_5patch_WT_matrixa_allpatch_1_5_20.csv",header=FALSE)
avepatchinfectRT_matrixa_allpatch<-read.csv("ave_SImodel_5patch_infected_RT_matrixa_allpatch_1_5_20.csv",header=FALSE)
avepatchinfectWT_matrixa_allpatch<-read.csv("ave_SImodel_5patch_infected_WT_matrixa_allpatch_1_5_20.csv",header=FALSE)

avepatchinfect_matrixk_allpatch<-read.csv("ave_SImodel_5patch_infected_total_matrixk_allpatch_1_5_20.csv",header=FALSE)
avepatchRT_matrixk_allpatch<-read.csv("ave_SImodel_5patch_RT_matrixk_allpatch_1_5_20.csv",header=FALSE)
avepatchWT_matrixk_allpatch<-read.csv("ave_SImodel_5patch_WT_matrixk_allpatch_1_5_20.csv",header=FALSE)
avepatchinfectRT_matrixk_allpatch<-read.csv("ave_SImodel_5patch_infected_RT_matrixk_allpatch_1_5_20.csv",header=FALSE)
avepatchinfectWT_matrixk_allpatch<-read.csv("ave_SImodel_5patch_infected_WT_matrixk_allpatch_1_5_20.csv",header=FALSE)


####FIGURE S9
tiff("avepatchSImodel_matrixa_k_allpatch_1_5_20.tiff", width=10,height=10, units='in',res=600)

par(mfrow=c(3,2))
par(mar=c(2,5,1,4.5))

plot(Time,avepatchWT_matrixa_allpatch[,1],type="l",lwd=2,xlim=c(8,60),ylab=expression(bold("Wild-type Size")),cex.axis=1.5,xlab="",cex.lab=1.5,xaxt="n",ylim=c(1000,1280),col="black") #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,avepatchWT_matrixk_allpatch[,3],type="l",lwd=2,lty=2,col="red")
points(Time,avepatchWT_matrixa_allpatch[,3],type="l",lwd=2,lty=3,col="blue")
points(Time,avepatchWT_matrixa_allpatch[,5],type="l",lwd=2,lty=5,col="blue")
legend(8,1180,c("No Migration","Edge=5,Mig=0.025","Edge=5,Mig=0.05","Edge=10,Mig=0.025"),lty=c(1,3,5,2),col=c("black","blue","blue","red"),lwd=c(2,2,2,2),bty="n",cex=1.5)


plot(Time,avepatchRT_matrixa_allpatch[,1],type="l",lwd=2,xlim=c(8,60),col="black",ylab=expression(bold("Robust-type Size")),cex.axis=1.5,xlab="",cex.lab=1.5,xaxt="n",ylim=c(500,650)) #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,avepatchRT_matrixk_allpatch[,3],type="l",lwd=2,lty=2,col="red")
points(Time,avepatchRT_matrixa_allpatch[,3],type="l",lwd=2,lty=3,col="black")
points(Time,avepatchRT_matrixa_allpatch[,5],type="l",lwd=2,lty=5,col="blue")


plot(Time,avepatchWT_matrixa_allpatch[,1]/avepatchRT_matrixa_allpatch[,1],col="black",type="l",lwd=2,xlim=c(15,60),ylab=expression(bold("Host_ratio: Wild/Robust")),xlab="",cex.lab=1.5,xaxt="n",ylim=c(1.7,2.05)) #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,avepatchWT_matrixk_allpatch[,3]/avepatchRT_matrixk_allpatch[,3],type="l",lwd=2,lty=2,col="red")
points(Time,avepatchWT_matrixa_allpatch[,3]/avepatchRT_matrixa_allpatch[,3],type="l",lwd=2,lty=3,col="blue")
points(Time,avepatchWT_matrixa_allpatch[,5]/avepatchRT_matrixa_allpatch[,5],type="l",lwd=2,lty=5,col="blue")
legend(20,1.34,c("1","1.3","1.7","2"),lty=c(1,2,3,6),col=c("black","blue","brown","red"),lwd=c(2,2,2,2),title="Mortality Ratio",bty="n",cex=2)


plot(Time,avepatchWT_matrixa_allpatch[,1]+avepatchRT_matrixa_allpatch[,1],type="l",col="black",lwd=2,xlim=c(15,60),ylab=expression(bold("Host_Total: Wild+Robust")),xlab="",cex.lab=1.5,xaxt="n",ylim=c(1600,2000)) #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,avepatchWT_matrixk_allpatch[,3]+avepatchRT_matrixk_allpatch[,3],type="l",lwd=2,lty=2,col="red")
points(Time,avepatchWT_matrixa_allpatch[,3]+avepatchRT_matrixa_allpatch[,3],type="l",lwd=2,lty=3,col="blue")
points(Time,avepatchWT_matrixa_allpatch[,5]+avepatchRT_matrixa_allpatch[,5],type="l",lwd=2,lty=5,col="blue")



plot(Time,avepatchinfect_matrixa_allpatch[,1],type="l",lwd=2,xlim=c(15,60),col="black",ylab=expression(bold("Infected_total: Robust+Wild")),xlab="",cex.lab=1.5,xaxt="n",ylim=c(0,200))
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,avepatchinfect_matrixk_allpatch[,3],type="l",lwd=2,lty=2,col="red")
points(Time,avepatchinfect_matrixa_allpatch[,3],type="l",lwd=2,lty=3,col="blue")
points(Time,avepatchinfect_matrixa_allpatch[,5],type="l",lwd=2,lty=5,col="blue")


dev.off()



###FIGURE S10

avepatchinfect_matrixa_spread<-read.csv("SImodel_5patch_infected_total_matrixa_spread_1_5_20.csv",header=FALSE)
avepatchRT_matrixa_spread<-read.csv("SImodel_5patch_RT_matrixa_spread_1_5_20.csv",header=FALSE)
avepatchWT_matrixa_spread<-read.csv("SImodel_5patch_WT_matrixa_spread_1_5_20.csv",header=FALSE)
avepatchinfectRT_matrixa_spread<-read.csv("SImodel_5patch_infected_RT_matrixa_spread_1_5_20.csv",header=FALSE)
avepatchinfectWT_matrixa_spread<-read.csv("SImodel_5patch_infected_WT_matrixa_spread_1_5_20.csv",header=FALSE)


tiff("avepatchSImodel_matrixa_k_spread_1_5_20.tiff", width=10,height=10, units='in',res=600)

par(mfrow=c(3,2))
par(mar=c(2,5,1,4.5))

plot(Time,avepatchWT_matrixa_spread[,1],type="l",lwd=2,xlim=c(8,60),ylab=expression(bold("Wild-type Size")),cex.axis=1.5,xlab="",cex.lab=1.5,xaxt="n",ylim=c(1000,1280),col="black") #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,avepatchWT_matrixk_spread[,3],type="l",lwd=2,lty=2,col="red")
points(Time,avepatchWT_matrixa_spread[,3],type="l",lwd=2,lty=3,col="blue")
points(Time,avepatchWT_matrixa_spread[,5],type="l",lwd=2,lty=5,col="blue")
legend(8,1220,c("No Migration","Edge=5,Mig=0.025","Edge=5,Mig=0.05","Edge=10,Mig=0.025"),lty=c(1,3,5,2),col=c("black","blue","blue","red"),lwd=c(2,2,2,2),bty="n",cex=1.5)


plot(Time,avepatchRT_matrixa_spread[,1],type="l",lwd=2,xlim=c(8,60),col="black",ylab=expression(bold("Robust-type Size")),cex.axis=1.5,xlab="",cex.lab=1.5,xaxt="n",ylim=c(500,650)) #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,avepatchRT_matrixk_spread[,3],type="l",lwd=2,lty=2,col="red")
points(Time,avepatchRT_matrixa_spread[,3],type="l",lwd=2,lty=3,col="black")
points(Time,avepatchRT_matrixa_spread[,5],type="l",lwd=2,lty=5,col="blue")


plot(Time,avepatchWT_matrixa_spread[,1]/avepatchRT_matrixa_spread[,1],col="black",type="l",lwd=2,xlim=c(15,60),ylab=expression(bold("Host_ratio: Wild/Robust")),xlab="",cex.lab=1.5,xaxt="n",ylim=c(1.7,2.05)) #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,avepatchWT_matrixk_spread[,3]/avepatchRT_matrixk_spread[,3],type="l",lwd=2,lty=2,col="red")
points(Time,avepatchWT_matrixa_spread[,3]/avepatchRT_matrixa_spread[,3],type="l",lwd=2,lty=3,col="blue")
points(Time,avepatchWT_matrixa_spread[,5]/avepatchRT_matrixa_spread[,5],type="l",lwd=2,lty=5,col="blue")
legend(20,1.34,c("1","1.3","1.7","2"),lty=c(1,2,3,6),col=c("black","blue","brown","red"),lwd=c(2,2,2,2),title="Mortality Ratio",bty="n",cex=2)


plot(Time,avepatchWT_matrixa_spread[,1]+avepatchRT_matrixa_spread[,1],type="l",col="black",lwd=2,xlim=c(15,60),ylab=expression(bold("Host_Total: Wild+Robust")),xlab="",cex.lab=1.5,xaxt="n",ylim=c(1600,2000)) #ylim=c(1.442,1.444),
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,avepatchWT_matrixk_spread[,3]+avepatchRT_matrixk_spread[,3],type="l",lwd=2,lty=2,col="red")
points(Time,avepatchWT_matrixa_spread[,3]+avepatchRT_matrixa_spread[,3],type="l",lwd=2,lty=3,col="blue")
points(Time,avepatchWT_matrixa_spread[,5]+avepatchRT_matrixa_spread[,5],type="l",lwd=2,lty=5,col="blue")


plot(Time,avepatchinfect_matrixa_spread[,1],type="l",lwd=2,xlim=c(15,60),col="black",ylab=expression(bold("Infected_total: Robust+Wild")),xlab="",cex.lab=1.5,xaxt="n",ylim=c(0,200))
axis(1,seq(0,60,12),seq(0,15,3),cex.axis=1.5)
points(Time,avepatchinfect_matrixk_spread[,3],type="l",lwd=2,lty=2,col="red")
points(Time,avepatchinfect_matrixa_spread[,3],type="l",lwd=2,lty=3,col="blue")
points(Time,avepatchinfect_matrixa_spread[,5],type="l",lwd=2,lty=5,col="blue")

dev.off()





