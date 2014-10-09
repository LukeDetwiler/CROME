#CRUISE Runs
setwd("C:/Users/Alec/Desktop/CT_RIVER_OPTIMIZATION_MODEL/WestWeighting/STATUS_QUO_ECO_OFF")

####Import and manipulate date information####
numyrs = read.table("numyrs.txt",header=FALSE)[[1]]
StartYear = read.table("startyear.txt",header=FALSE)[[1]]
EndYear = StartYear + numyrs - 1
startdate  = paste(StartYear,"/1/1",sep="")
enddate = paste(EndYear,"/12/31",sep="")
doy = seq(as.Date(startdate),as.Date(enddate), by="day")
yy = as.numeric(format(doy,"%Y"))
mm = as.numeric(format(doy,"%m"))
dd = as.numeric(format(doy,"%d"))
no_leap = which(mm!=2 | dd!=29)
doy = doy[no_leap]
YEARS = as.numeric(format(doy,"%Y"))
MONTHS = as.numeric(format(doy,"%m"))
DAYS =  as.numeric(format(doy,"%d"))
JulianDay = rep(1:365,length(unique(YEARS)))





####Create Blank Arrays####

runs = c("STATUS_QUO_ECO_OFF",
         #"ECO_RUN_0.99",
         #"ECO_RUN_0.98",
#          "ECO_RUN_0.90",
         #"ECO_RUN_0.80", 
         #"ECO_RUN_0.70", 
#          "ECO_RUN_0.50",
#          "ECO_RUN_0.40",
#          "ECO_RUN_0.30",
#          "ECO_RUN_0.20",
#          "ECO_RUN_0.10",
         #"ECO_RUN_0.05",
         "CT_OUT_ECO_500"
)


num_runs = length(runs)

#Set Eco-Nodes of Interest
eco_node_num_list = c(33,95)


#USACE
WST_RES_BMD_ST = array(NA,c(length(doy),num_runs))
WST_RES_TWN_ST = array(NA,c(length(doy),num_runs))

WST_RES_BMD_ST_HIST = array(NA,c(length(doy),num_runs))
WST_RES_TWN_ST_HIST = array(NA,c(length(doy),num_runs))

WST_RES_BMD_R = array(NA,c(length(doy),num_runs))
WST_RES_TWN_R = array(NA,c(length(doy),num_runs))

WST_RES_BMD_R_HIST = array(NA,c(length(doy),num_runs))
WST_RES_TWN_R_HIST = array(NA,c(length(doy),num_runs))

WST_RES_BMD_ST_FLOOD_TARG_ABOVE = array(NA,c(length(doy),num_runs))
WST_RES_TWN_ST_FLOOD_TARG_ABOVE = array(NA,c(length(doy),num_runs))
WST_RES_BMD_ST_FLOOD_TARG_BELOW = array(NA,c(length(doy),num_runs))
WST_RES_TWN_ST_FLOOD_TARG_BELOW = array(NA,c(length(doy),num_runs))

WST_RES_BMD_ST_TARG_ABOVE = array(NA,c(length(doy),num_runs))
WST_RES_TWN_ST_TARG_ABOVE = array(NA,c(length(doy),num_runs))

WST_RES_BMD_ST_TARG_BELOW = array(NA,c(length(doy),num_runs))
WST_RES_TWN_ST_TARG_BELOW = array(NA,c(length(doy),num_runs))

#USACE FLOOD CHECKPOINTS
WST_G_PRXY = array (NA,c(length(doy), num_runs), dimnames = runs )
WST_RES_BMD_ST_FLOOD_TARG_ABOVE = array (NA,c(length(doy), num_runs))
WST_RES_TWN_ST_FLOOD_TARG_ABOVE = array (NA,c(length(doy), num_runs))

#ECO_NODE PENALTIES
ECO_PENALTIES = array(NA,c(length(doy),length(eco_node_num_list),num_runs))
ECO_NODE_MOD_FLOWS = array(NA,c(length(doy),length(eco_node_num_list),num_runs))
ECO_NODE_NAT_FLOWS = array(NA,c(length(doy),length(eco_node_num_list),num_runs))

ECO_PEN_HIGH1_ABOVE = array(NA,c(length(doy),length(eco_node_num_list),num_runs))
ECO_PEN_HIGH2_ABOVE = array(NA,c(length(doy),length(eco_node_num_list),num_runs))
ECO_PEN_LOW1_BELOW = array(NA,c(length(doy),length(eco_node_num_list),num_runs))
ECO_PEN_LOW2_BELOW = array(NA,c(length(doy),length(eco_node_num_list),num_runs))

ECO_PEN_HIGH1_BELOW = array(NA,c(length(doy),length(eco_node_num_list),num_runs))
ECO_PEN_HIGH2_BELOW = array(NA,c(length(doy),length(eco_node_num_list),num_runs))
ECO_PEN_LOW1_ABOVE = array(NA,c(length(doy),length(eco_node_num_list),num_runs))
ECO_PEN_LOW2_ABOVE = array(NA,c(length(doy),length(eco_node_num_list),num_runs))

#EN_TARGETS
EN_TARGETS =array(NA, c(length(doy), 4*length(eco_node_num_list),num_runs))

#DUAL PRICES
DUAL_PRICE_EN_CONSTRAINTS =array(NA, c(length(doy), 6, 1, length(eco_node_num_list)))

####Fill Blank Arrays#####
for (j in 1:num_runs) {
  #For Climate Change Runs
  #cur_dir = paste("C:/Users/janleitner/Desktop/WEST_ONLY_CC/",runs[j],sep="")
  
  cur_dir = paste("C:/Users/Alec/Desktop/CT_RIVER_OPTIMIZATION_MODEL/WestWeighting/", runs[j], sep="")
  setwd(cur_dir)
  filenames = list.files(cur_dir)  
  
  #USACE
  WST_RES_BMD_ST[,j] = scan("WST_RES_BMD_ST.txt")/43560+250
  WST_RES_TWN_ST[,j] = scan("WST_RES_TWN_ST.txt")/43560+800
  
  WST_RES_BMD_ST_HIST[,j] = scan("WST_RES_BMD_ST_HIST.txt")/43560 + 250
  WST_RES_TWN_ST_HIST[,j] = scan("WST_RES_TWN_ST_HIST.txt")/43560 + 800
  
  WST_RES_BMD_R[,j] = scan("WST_RES_BMD_R.txt")/86400
  WST_RES_TWN_R[,j] = scan("WST_RES_TWN_R.txt")/86400
  
  WST_RES_BMD_R_HIST[,j] = scan("WST_RES_BMD_R_HIST.txt")/86400
  WST_RES_TWN_R_HIST[,j] = scan("WST_RES_TWN_R_HIST.txt")/86400
  
  WST_RES_BMD_ST_FLOOD_TARG_ABOVE[,j] = scan("WST_RES_BMD_ST_FLOOD_TARG_ABOVE.txt")/43560
  WST_RES_TWN_ST_FLOOD_TARG_ABOVE[,j] = scan("WST_RES_TWN_ST_FLOOD_TARG_ABOVE.txt")/43560
  WST_RES_BMD_ST_FLOOD_TARG_BELOW[,j] = scan("WST_RES_BMD_ST_FLOOD_TARG_BELOW.txt")/43560
  WST_RES_TWN_ST_FLOOD_TARG_BELOW[,j] = scan("WST_RES_TWN_ST_FLOOD_TARG_BELOW.txt")/43560
  
  WST_RES_BMD_ST_TARG_ABOVE[,j] = scan("WST_RES_BMD_ST_TARG_ABOVE.txt")/43560
  WST_RES_TWN_ST_TARG_ABOVE[,j] = scan("WST_RES_TWN_ST_TARG_ABOVE.txt")/43560
  
  WST_RES_BMD_ST_TARG_BELOW[,j] = scan("WST_RES_BMD_ST_TARG_BELOW.txt")/43560
  WST_RES_TWN_ST_TARG_BELOW[,j] = scan("WST_RES_TWN_ST_TARG_BELOW.txt")/43560
  
  #FLOOD CHECKPOINTS
  WST_G_PRXY[,j] = scan("WST_G_PRXY.txt")/86400
  WST_RES_BMD_ST_FLOOD_TARG_ABOVE [,j]= scan("WST_RES_BMD_ST_FLOOD_TARG_ABOVE.txt")/86400
  WST_RES_TWN_ST_FLOOD_TARG_ABOVE [,j]= scan("WST_RES_TWN_ST_FLOOD_TARG_ABOVE.txt")/86400
  
  for (eco in 1:length(eco_node_num_list)) {
    pattern1 = paste("EN_",eco_node_num_list[eco],"_HIGH2_ABOVE.txt",sep="")
    pattern2 = paste("EN_",eco_node_num_list[eco],"_HIGH1_ABOVE.txt",sep="")
    pattern3 = paste("EN_",eco_node_num_list[eco],"_LOW1_BELOW.txt",sep="")
    pattern4 = paste("EN_",eco_node_num_list[eco],"_LOW2_BELOW.txt",sep="")
    pattern5 = paste("EN_",eco_node_num_list[eco],"_HIGH1_BELOW.txt",sep="")
    pattern6 = paste("EN_",eco_node_num_list[eco],"_HIGH2_BELOW.txt",sep="")
    pattern7 = paste("EN_",eco_node_num_list[eco],"_LOW1_ABOVE.txt",sep="")
    pattern8 = paste("EN_",eco_node_num_list[eco],"_LOW2_ABOVE.txt",sep="")
    ECO_PENALTIES[,eco,j] = (2*scan(filenames[grep(pattern1,filenames)]) + 
                               scan(filenames[grep(pattern2,filenames)]) + 
                               scan(filenames[grep(pattern3,filenames)]) + 
                               2*scan(filenames[grep(pattern4,filenames)]))/86400
    pattern9 = paste("EN_",eco_node_num_list[eco],"_NAT.txt",sep="")
    ECO_NODE_NAT_FLOWS[,eco,j] = scan(filenames[grep(pattern9,filenames)])/86400
    pattern10 = paste("EN_",eco_node_num_list[eco],"_MOD.txt",sep="")
    ECO_NODE_MOD_FLOWS[,eco,j] = scan(filenames[grep(pattern10,filenames)])/86400
    
    ECO_PEN_HIGH2_ABOVE[,eco,j] = scan(filenames[grep(pattern1,filenames)])/86400
    ECO_PEN_HIGH1_ABOVE[,eco,j] = scan(filenames[grep(pattern2,filenames)])/86400
    ECO_PEN_LOW1_BELOW[,eco,j] = scan(filenames[grep(pattern3,filenames)])/86400
    ECO_PEN_LOW2_BELOW[,eco,j] = scan(filenames[grep(pattern4,filenames)])/86400
    
    ECO_PEN_HIGH1_BELOW[,eco,j] = scan(filenames[grep(pattern5,filenames)])/86400
    ECO_PEN_HIGH2_BELOW[,eco,j] = scan(filenames[grep(pattern6,filenames)])/86400
    ECO_PEN_LOW1_ABOVE[,eco,j] = scan(filenames[grep(pattern7,filenames)])/86400
    ECO_PEN_LOW2_ABOVE[,eco,j] = scan(filenames[grep(pattern8,filenames)])/86400
    
    #EN Targets
    target1 = paste("WST_EN_", eco_node_num_list[eco], "_LOW2.txt", sep="")
    target2 = paste("WST_EN_", eco_node_num_list[eco], "_LOW1.txt", sep="")
    target3 = paste("WST_EN_", eco_node_num_list[eco], "_HIGH1.txt", sep="")
    target4 = paste("WST_EN_", eco_node_num_list[eco], "_HIGH2.txt", sep="")
    
    EN_TARGETS[,(4*(eco-1)+1),j] = scan(filenames[grep(target1,filenames)])
    EN_TARGETS[,(4*(eco-1)+2),j] = scan(filenames[grep(target2,filenames)])
    EN_TARGETS[,(4*(eco-1)+3),j] = scan(filenames[grep(target3,filenames)])
    EN_TARGETS[,(4*(eco-1)+4),j] = scan(filenames[grep(target4,filenames)])
  }
}

ReservoirName = c("Ball Mountain Dam", "Townshend Dam")

##CALCULATE AND NORMALIZE ECO-NODE PENALTIES JUST FOR WEST SUB-BASIN ECO-NODES####
en.name.list <- list(doy, eco_node_num_list, runs)

dimnames(ECO_PENALTIES) <- en.name.list

en.nos <- c("33", "95")

ECO_PENALTIES_WEST <- array (NA, c(num_runs))

for (i in 1:num_runs){
  ECO_PENALTIES_WEST[i] <- sum(ECO_PENALTIES[,en.nos, runs[i]])
}



####NOW LET'S MOVE ON TO Storage PENALTIES####

dimnames(WST_RES_BMD_ST_TARG_ABOVE)<-list(doy, runs)
dimnames(WST_RES_TWN_ST_TARG_ABOVE)<-list(doy, runs)
dimnames(WST_RES_BMD_ST_TARG_BELOW)<-list(doy, runs)
dimnames(WST_RES_TWN_ST_TARG_BELOW)<-list(doy, runs)

STORAGE_PENALTIES_WEST <- array(NA,c(num_runs))
STORAGE_PENALTIES_WEST_BMD <- array(NA,c(num_runs))
STORAGE_PENALTIES_WEST_TWN <- array(NA,c(num_runs))

WST_RES_BMD_ST_TARG_ABOVE_WT <- 10.008008
WST_RES_BMD_ST_TARG_BELOW_WT <- 20.016105
WST_RES_TWN_ST_TARG_ABOVE_WT <- 15.955015
WST_RES_TWN_ST_TARG_BELOW_WT <- 31.91003

for (j in 1:num_runs){
  STORAGE_PENALTIES_WEST [j] <- sum(WST_RES_BMD_ST_TARG_ABOVE[,runs[j]])*WST_RES_BMD_ST_TARG_ABOVE_WT +
    sum(WST_RES_TWN_ST_TARG_ABOVE[,runs[j]])*WST_RES_TWN_ST_TARG_ABOVE_WT +
    sum(WST_RES_BMD_ST_TARG_BELOW[,runs[j]])*WST_RES_BMD_ST_TARG_BELOW_WT +
    sum(WST_RES_TWN_ST_TARG_BELOW[,runs[j]])*WST_RES_TWN_ST_TARG_BELOW_WT
  
}
# Ball mountain Dam Penalties 
for (j in 1:num_runs){
  STORAGE_PENALTIES_WEST_BMD [j] <- sum(WST_RES_BMD_ST_TARG_ABOVE[,runs[j]])*WST_RES_BMD_ST_TARG_ABOVE_WT +
    sum(WST_RES_BMD_ST_TARG_BELOW[,runs[j]])*WST_RES_BMD_ST_TARG_BELOW_WT 
}
# Townshend Dam
for (j in 1:num_runs){
  STORAGE_PENALTIES_WEST_TWN [j] <- sum(WST_RES_TWN_ST_TARG_ABOVE[,runs[j]])*WST_RES_TWN_ST_TARG_ABOVE_WT +
    sum(WST_RES_TWN_ST_TARG_BELOW[,runs[j]])*WST_RES_TWN_ST_TARG_BELOW_WT
}

ReservoirName = c("Ball Mountain Dam", "Townshend Dam")

NumReservoirs = length(ReservoirName)

cur_run1 = which(runs=="STATUS_QUO_ECO_OFF")
cur_run2 = which(runs=="CT_OUT_ECO")
cur_run3 = which(runs=="ECO_RUN_0.30")
cur_run4 = which(runs=="ECO_RUN_0.20")
cur_run5 = which(runs=="ECO_RUN_0.10")

setwd("C:/Users/Alec/Desktop/CT_RIVER_OPTIMIZATION_MODEL/WestWeighting/Reservoir Storage Targets")

library(abind)

WST_RES_BMD_ST_TARG <- scan("WST_RES_BMD_ST_TARG.txt")/43560 + 250
WST_RES_TWN_ST_TARG <- scan("WST_RES_TWN_ST_TARG.txt")/43560 + 800

WST_RES_BMD_CAP <-54450 + 250
WST_RES_TWN_CAP <- 32900 + 800

TARGETS<-abind(WST_RES_BMD_ST_TARG, WST_RES_TWN_ST_TARG,
               along=2)

STORAGES<-abind(WST_RES_BMD_ST, WST_RES_TWN_ST,
                along=3)

HIST<-abind(WST_RES_BMD_ST_HIST, WST_RES_TWN_ST_HIST,
            along=3)

CAPS<-abind(WST_RES_BMD_CAP, WST_RES_TWN_CAP, 
            along=1)

eco.nat.select<-c(2,1)

ReservoirName = c("Ball Mountain Dam", "Townshend Dam")

NumReservoirs = length(ReservoirName)

#Original Plot 

setwd("C:/Users/Alec/Desktop/CT_RIVER_OPTIMIZATION_MODEL/WestWeighting/")

tiff("West sub-basin Storage vs Eco.tif",width = 966, height = 700, units = "px", res=100)
par(mfrow=c(1,1), mar=c(5,5,2,2), xpd=F)
plot(ECO_PENALTIES_WEST, STORAGE_PENALTIES_WEST, type="b",cex=1.5,
     
     xlim=rev(range(seq(1,max(ECO_PENALTIES_WEST)))),
     ylim=rev(range(seq(1, max(STORAGE_PENALTIES_WEST), by=1000))),
     xlab="Eco-Node Penalties", ylab="Dam Storage Penalties",
     #main="West Sub-Basin Tradeoff", 
     cex.lab=1.75, cex.axis=1.5
)
title("West Basin Tradeoffs")
points(ECO_PENALTIES_WEST[1], STORAGE_PENALTIES_WEST[1], pch=19)
points(ECO_PENALTIES_WEST[cur_run3], STORAGE_PENALTIES_WEST[cur_run3], pch=19)
points(ECO_PENALTIES_WEST[cur_run4], STORAGE_PENALTIES_WEST[cur_run4], pch=19)
points(ECO_PENALTIES_WEST[num_runs], STORAGE_PENALTIES_WEST[num_runs], pch=19)

dev.off()

resNameCounter=1
tiff("West sub-basin BMD Storage vs Eco.tif",width = 966, height = 700, units = "px", res=100)
par(mfrow=c(1,1), mar=c(5,5,2,2), xpd=F)
plot(ECO_PENALTIES_WEST, STORAGE_PENALTIES_WEST_BMD, type="b",cex=1.5,
     
     xlim=rev(range(seq(1,max(ECO_PENALTIES_WEST)))),
     ylim=rev(range(seq(1, max(STORAGE_PENALTIES_WEST), by=1000))),
     xlab="Eco-Node Penalties", ylab="Dam Storage Penalties",
     #main="West Sub-Basin Tradeoff", 
     cex.lab=1.75, cex.axis=1.5
)
title(paste0(ReservoirName[resNameCounter]," Tradeoffs"))
points(ECO_PENALTIES_WEST[1], STORAGE_PENALTIES_WEST_BMD[1], pch=19)
points(ECO_PENALTIES_WEST[cur_run3], STORAGE_PENALTIES_WEST_BMD[cur_run3], pch=19)
points(ECO_PENALTIES_WEST[cur_run4], STORAGE_PENALTIES_WEST_BMD[cur_run4], pch=19)
points(ECO_PENALTIES_WEST[num_runs], STORAGE_PENALTIES_WEST_BMD[num_runs], pch=19)

dev.off()

resNameCounter=resNameCounter+1
tiff("West sub-basin TWN Storage vs Eco.tif",width = 966, height = 700, units = "px", res=100)
par(mfrow=c(1,1), mar=c(5,5,2,2), xpd=F)
plot(ECO_PENALTIES_WEST, STORAGE_PENALTIES_WEST_TWN, type="b",cex=1.5,
     
     xlim=rev(range(seq(1,max(ECO_PENALTIES_WEST)))),
     ylim=rev(range(seq(1, max(STORAGE_PENALTIES_WEST), by=1000))),
     xlab="Eco-Node Penalties", ylab="Dam Storage Penalties",
     #main="West Sub-Basin Tradeoff", 
     cex.lab=1.75, cex.axis=1.5
)
title(paste0(ReservoirName[resNameCounter]," Tradeoffs"))
points(ECO_PENALTIES_WEST[1], STORAGE_PENALTIES_WEST_TWN[1], pch=19)
points(ECO_PENALTIES_WEST[cur_run3], STORAGE_PENALTIES_WEST_TWN[cur_run3], pch=19)
points(ECO_PENALTIES_WEST[cur_run4], STORAGE_PENALTIES_WEST_TWN[cur_run4], pch=19)
points(ECO_PENALTIES_WEST[num_runs], STORAGE_PENALTIES_WEST_TWN[num_runs], pch=19)

dev.off()

setwd("C:/Users/Alec/Desktop/CT_RIVER_OPTIMIZATION_MODEL/WestWeighting/")
#Average Storages
for (z in 1:NumReservoirs){
  
  row.select = which(HIST[,1,z] >=0)
  
  #Plot within capacity
  ymin <-min(TARGETS[,z], STORAGES[,1,z], STORAGES[,num_runs,z], max(0,HIST[,1,z]))
  ymax <-max(TARGETS[,z], STORAGES[,1,z], STORAGES[,num_runs,z], HIST[,1,z], CAPS[z])
  
  png(paste0(ReservoirName[z],"_STORAGES_CAP.png"),width = 966, height = 725, units = "px")
  par(mar=c(9,4,4,2), cex=2, xpd=FALSE)
  plot(doy[1:365], TARGETS[,z], type="l", ylim=c(ymin,ymax), xlab="", ylab="", lwd=2, col="purple4", lty=3) #target
  if (length(row.select)>0){
    lines(doy[1:365], aggregate(HIST[row.select,1,z], by=list(JulianDay[row.select]), FUN=mean)[,2], lwd=2, lty = 1 ,col="paleturquoise3")#historic
  }
  
  lines(doy[1:365], aggregate(STORAGES[,1,z], by=list(JulianDay), FUN=mean, na.rm=T)[,2], lwd=2, col="red3") #Model SQ
  lines(doy[1:365], aggregate(STORAGES[,num_runs,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="dodgerblue3") #model Eco
  abline(h =CAPS[z], lwd=2, lty=2, col="red")
  mtext(side=1,"Month",line=2.5,font=2, cex=2)
  mtext(side=2,"Storage (acre-feet)",line=2.5,font=2, cex=2)
  title(paste0(ReservoirName[z]," Storage"))
  par(xpd=TRUE)
  
  if (length(row.select)>0){
    legend("bottom",inset=c(0,-.65), c("Capacity", "Status Quo (Pt. A)", "E-Flow Target (Pt. D)", "Historic"), col=c("red", "red3", "dodgerblue3", 
                                                                                                                     "paleturquoise3"), lwd=c(2,2,2,2), lty=c(2,1,1,1), y.intersp=1.2, x.intersp=1.2, ncol=2)
  } else {
    legend("bottom",inset=c(0,-.65), c("Capacity", "Status Quo (Pt. A)", "E-Flow Target (Pt. D)"), col=c("red", "red3", "dodgerblue3"), 
           lwd=c(2,2,2), lty=c(2,1,1), y.intersp=1.2, x.intersp=1.2, ncol=2)  
    
  }
  dev.off()
  
  #Plot Extremes, Historic and Target
  ymin <-min(TARGETS[,z], aggregate(STORAGES[,1,z], by=list(JulianDay), FUN=mean)[,2], 
             aggregate(STORAGES[,num_runs,z], by=list(JulianDay), FUN=mean)[,2],
             aggregate(HIST[row.select,1,z], by=list(JulianDay[row.select]), FUN=mean)[,2])
  ymax <-max(TARGETS[,z], aggregate(STORAGES[,1,z], by=list(JulianDay), FUN=mean)[,2], 
             aggregate(STORAGES[,num_runs,z], by=list(JulianDay), FUN=mean)[,2],
             aggregate(HIST[row.select,1,z], by=list(JulianDay[row.select]), FUN=mean)[,2])
  
  png(paste0(ReservoirName[z],"_STORAGES_EXTREMES.png"),width = 966, height = 725, units = "px")
  par(mar=c(9,4,3,5), cex=2, xpd=FALSE)
  plot(doy[1:365], TARGETS[,z], type="l", ylim=c(ymin,ymax), xlab="", ylab="", lwd=2, col="purple4", lty=3) #target
  if (length(row.select)>0){
    lines(doy[1:365], aggregate(HIST[row.select,1,z], by=list(JulianDay[row.select]), FUN=mean)[,2], lwd=2, lty = 1 ,col="paleturquoise3")#historic
  }
  lines(doy[1:365], aggregate(STORAGES[,1,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="red3") #Model SQ
  lines(doy[1:365], aggregate(STORAGES[,num_runs,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="dodgerblue3") #model Eco
  
  mtext(side=1,"Month",line=2.5,font=2, cex=2)
  mtext(side=2,"Storage (acre-feet)",line=2.5,font=2, cex=2)
  title(paste0(ReservoirName[z]," Storage"))
  
  par(new=TRUE)
  plot(doy[1:365], aggregate(ECO_NODE_NAT_FLOWS[,eco.nat.select[z],1], by=list(JulianDay), FUN=mean)[,2] ,type="l",col="gray47", #natural flows
       xaxt="n",yaxt="n",xlab="",ylab="", lty=5, lwd=1)
  axis(4)
  mtext("Natural Flow (cfs)",side=4,line=3, font=2, cex=2)
  
  par(xpd=T)
  
  if (length(row.select)>0){
    legend("bottom",inset=c(2,-.625), c("Target", "Status Quo (Pt. A)", "E-Flow Target (Pt. D)", "Historic", "Natural Flows"), 
           col=c("black", "red3", "dodgerblue3", "paleturquoise3", "gray47"), 
           lwd=c(2,2,2,2,1), lty=c(3,1,1,1,5), y.intersp=1.2, x.intersp=1,
           ncol=3, cex=.9)
  } else {
    legend("bottom",inset=c(.5,-.625), c("Target", "Status Quo(Pt. A)", "E-Flow Target (Pt. D)", "Natural Flows"), 
           col=c("black", "red3", "dodgerblue3", "gray47"), 
           lwd=c(2,2,2,1), lty=c(3,1,1,5), y.intersp=1.2, x.intersp=1,
           ncol=2)
    
  }
  dev.off() 
  
  #Plot Extremes and two tradeoff points (selected under cur_run_x and y)
  png(paste0(ReservoirName[z],"_STORAGES_TRADEOFFS.png"),width = 966, height = 725, units = "px")
  par(mar=c(9,4,3,5), cex=2, xpd=FALSE, lheight=5)
  plot(doy[1:365], TARGETS[,z], type="l", ylim=c(ymin,ymax), xlab="", ylab="", lwd=2, col="purple4", lty=3) #target
  lines(doy[1:365], aggregate(STORAGES[,1,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="red3") #Model SQ
  lines(doy[1:365], aggregate(STORAGES[,num_runs,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="dodgerblue3") #model Eco
#   lines(doy[1:365], aggregate(STORAGES[,cur_run3,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="orange2", lty=1) #ECO RUN 7.0
#   lines(doy[1:365], aggregate(STORAGES[,cur_run4,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="olivedrab3", lty=1) #ECO RUN 0.50
  
  mtext(side=1,"Month",line=2.5,font=2, cex=2)
  mtext(side=2,"Storage (acre-feet)",line=2.5,font=2, cex=2)
  title(paste0(ReservoirName[z]," Storage"))
  
  par(new=TRUE)
  plot(doy[1:365], aggregate(ECO_NODE_NAT_FLOWS[,eco.nat.select[z],1], by=list(JulianDay), FUN=mean)[,2] ,type="l",col="gray47", #natural flows
       xaxt="n",yaxt="n",xlab="",ylab="", lty=5, lwd=1)
  axis(4)
  mtext("Natural Flow (cfs)",side=4,line=3, font=2, cex=2)
  par(xpd=TRUE, cex=1.55)
  legend("bottom",inset=c(-.5,-.625), c("Target","Status Quo (Pt. A)", 
                                        paste0(runs[cur_run3], " (Pt. B)"), 
                                        paste0(runs[cur_run4], " (Pt. C)"), 
                                        "E-Flow Target (Pt. D)", "Natural Flows"), 
         col=c("purple4", "red3",  "orange2", "olivedrab3", "dodgerblue3", "gray47"), 
         lwd=c(2,2,2,2,2,1), lty=c(3,1,1,1,1,5), y.intersp=1.2,x.intersp=1.2, ncol=3)
  dev.off()
  
  
  ####Make two previous plots with different limits
  
  
  ymin <-min(TARGETS[,z], aggregate(STORAGES[,1,z], by=list(JulianDay), FUN=mean)[,2], 
             aggregate(STORAGES[,num_runs,z], by=list(JulianDay), FUN=mean)[,2],
             aggregate(HIST[row.select,1,z], by=list(JulianDay[row.select]), FUN=mean)[,2]
#              ,aggregate(STORAGES[,cur_run3,z], by=list(JulianDay), FUN=mean)[,2],
#              aggregate(STORAGES[,cur_run4,z], by=list(JulianDay), FUN=mean)[,2],
#              aggregate(STORAGES[,cur_run5,z], by=list(JulianDay), FUN=mean)[,2]
             )
             
             ymax <-max(TARGETS[,z], aggregate(STORAGES[,1,z], by=list(JulianDay), FUN=mean)[,2], 
                        aggregate(STORAGES[,num_runs,z], by=list(JulianDay), FUN=mean)[,2],
                        aggregate(HIST[row.select,1,z], by=list(JulianDay[row.select]), FUN=mean)[,2]
#                         ,aggregate(STORAGES[,cur_run3,z], by=list(JulianDay), FUN=mean)[,2],
#                         aggregate(STORAGES[,cur_run4,z], by=list(JulianDay), FUN=mean)[,2],
#                         aggregate(STORAGES[,cur_run5,z], by=list(JulianDay), FUN=mean)[,2]
                        )
                        
                        
                        #Plot Extremes, Historic and Target
                        png(paste0(ReservoirName[z],"_diff_limits_STORAGES_EXTREMES_limits2.png"),width = 966, height = 725, units = "px")
                        par(mar=c(9,4,3,5), cex=2, xpd=FALSE)
                        plot(doy[1:365], TARGETS[,z], type="l", ylim=c(ymin,ymax), xlab="", ylab="", lwd=2, col="purple4", lty=3) #target
                        if (length(row.select)>0){
                          lines(doy[1:365], aggregate(HIST[row.select,1,z], by=list(JulianDay[row.select]), FUN=mean)[,2], lwd=2, lty = 1 ,col="paleturquoise3")#historic
                        }
                        lines(doy[1:365], aggregate(STORAGES[,1,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="red3") #Model SQ
                        lines(doy[1:365], aggregate(STORAGES[,num_runs,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="dodgerblue3") #model Eco
                        
                        mtext(side=1,"Month",line=2.5,font=2, cex=2)
                        mtext(side=2,"Storage (acre-feet)",line=2.5,font=2, cex=2)
                        title(paste0(ReservoirName[z]," Storage"))
                        
                        par(new=TRUE)
                        plot(doy[1:365], aggregate(ECO_NODE_NAT_FLOWS[,eco.nat.select[z],1], by=list(JulianDay), FUN=mean)[,2] ,type="l",col="gray47", #natural flows
                             xaxt="n",yaxt="n",xlab="",ylab="", lty=5, lwd=1)
                        axis(4)
                        mtext("Natural Flow (cfs)",side=4,line=3, font=2, cex=2)
                        
                        par(xpd=TRUE)
                        
                        if (length(row.select)>0){
                          legend("bottom",inset=c(0,-.625), c("Target", "Status Quo (Pt. A)", "E-Flow Target (Pt. D)", "Historic", "Natural Flows"), 
                                 col=c("black", "red3", "dodgerblue3", "paleturquoise3", "gray47"), 
                                 lwd=c(2,2,2,2,1), lty=c(3,1,1,1,5), y.intersp=1.2, x.intersp=1.2,
                                 ncol=3, adj=c(0,1))
                        } else {
                          legend("bottom",inset=c(0,-.625), c("Target", "Status Quo (Pt. A)", "E-Flow Target (Pt. D)", "Natural Flows"), 
                                 col=c("black", "red3", "dodgerblue3", "gray47"), 
                                 lwd=c(2,2,2,1), lty=c(3,1,1,5), y.intersp=1.2, x.intersp=1.2,
                                 ncol=2)
                          
                        }
                        dev.off() 
                        
                        #Plot Extremes and two tradeoff points (selected under cur_run_x and y)
                        png(paste0(ReservoirName[z],"_diff_limits_STORAGES_TRADEOFFS.png"),width = 966, height = 725, units = "px")
                        par(mar=c(9,4,3,5), cex=2, xpd=FALSE, lheight=5)
                        plot(doy[1:365], TARGETS[,z], type="l", ylim=c(ymin,ymax), xlab="", ylab="", lwd=2, col="purple4", lty=3) #target
                        lines(doy[1:365], aggregate(STORAGES[,1,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="red3") #Model SQ
                        lines(doy[1:365], aggregate(STORAGES[,num_runs,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="dodgerblue3") #model Eco
#                         lines(doy[1:365], aggregate(STORAGES[,cur_run3,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="olivedrab1", lty=1) #ECO RUN 7.0
#                         lines(doy[1:365], aggregate(STORAGES[,cur_run4,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="olivedrab3", lty=1) #ECO RUN 0.50
#                         lines(doy[1:365], aggregate(STORAGES[,cur_run5,z], by=list(JulianDay), FUN=mean)[,2], lwd=2, col="olivedrab4", lty=1)
                        
                        mtext(side=1,"Month",line=2.5,font=2, cex=2)
                        mtext(side=2,"Storage (acre-feet)",line=2.5,font=2, cex=2)
                        title(paste0(ReservoirName[z]," Storage"))
                        
                        par(new=TRUE)
                        plot(doy[1:365], aggregate(ECO_NODE_NAT_FLOWS[,eco.nat.select[z],1], by=list(JulianDay), FUN=mean)[,2] ,type="l",col="gray47", #natural flows
                             xaxt="n",yaxt="n",xlab="",ylab="", lty=5, lwd=1)
                        axis(4)
                        mtext("Natural Flow (cfs)",side=4,line=3, font=2, cex=2)
                        par(xpd=TRUE, cex=1.65)
                        legend("bottom",inset=c(-.5,-.625), c("Target","Status Quo (Pt. A)", 
                                                              paste0(runs[cur_run3], " (Pt. B)"), 
                                                              paste0(runs[cur_run4], " (Pt. C)"),
                                                              runs[cur_run5],
                                                              "E-Flow Target (Pt. D)", "Natural Flows"), 
                               col=c("purple4", "red3",  "olivedrab1", "olivedrab3", "olivedrab4", "gray47"), 
                               lwd=c(2,2,2,2,2,1), lty=c(3,1,1,1,1,5), y.intersp=1.2,x.intersp=1.2, ncol=3)
                        dev.off()
                        
                        
}