#Tremblay random walk simulation
#Note: number of random walk steps greatly reduced in this example. For full analysis following steps outlined in Barkley et al. allow ~6 days to run depending on the computer. Step length can also be increased to provide longer tracks with a faster run time.

library(glatos)
library(dplyr)
library(ggplot2)
set.seed(1)

#make the polygon for where animals will swim
mypolygon<- data.frame(x = c(1,1,400,400), y = c(1,400,400,1))

#run the random walk for both animalR and animalT
animalR<- crw_in_polygon(polyg = mypolygon, theta = c(0,8), stepLen = 1, initPos = c(NA,NA), initHeading = NA, nsteps = 10000, sp_out = FALSE)
animalT<- crw_in_polygon(polyg = mypolygon, theta = c(0,8), stepLen = 1, initPos = c(NA,NA), initHeading = NA, nsteps = 10000, sp_out = FALSE)

#Use glatos transmit_along_path function to convert positions from meters to seconds (i.e. here we are essentially getting a position for each animal for every second)
animalR.sec<- transmit_along_path(path=animalR, vel= 0.34, delayRng = c(1,1), burstDur = 0, sp_out = FALSE)
animalT.sec<- transmit_along_path(path=animalT, vel= 0.34, delayRng = c(1,1), burstDur = 0, sp_out = FALSE)

#convert to whole seconds
animalT.sec$et<- floor(animalT.sec$et)
animalR.sec$et<- floor(animalR.sec$et)

#Get animalT's transmissions

# -- modified version of transmit_along_path from glatos
transmit_along_path2<- function (path = NA, vel = 0.34, delayRng = c(120,240), burstDur = 5){
  ntrns<- max(path$et)/(delayRng[1] + burstDur)
  ints<- runif(ntrns, delayRng[1] + burstDur, delayRng[2] + burstDur)
  ints[1]<- runif(1,0, ints[1])
  etime<- cumsum(ints)
  trns<- path
  trns$trans<- path$et %in% floor(etime)
  return(trns)
}

animalT.T<- transmit_along_path2(path = animalT.sec, vel = 0.34, delayRng = c(120,240), burstDur = 5)

# Now calculate animalR detections:
#detection range function
detRngFun<- function(dm, b=c(5.5, -0.0275)) {
  p<- 1/(1+exp(-(b[1]+b[2]*dm)))
  return (p)
}

#loop for determining successful detections - pulled from glatos function 'detect_transmissions' but allows for the 'receiver' to also move - i.e. animalR
dtc<- data.frame(trns_et = NA, recv_id = NA, recv_x = NA, recv_y = NA, trns_x = NA, trns_y = NA, distance = NA)
animalT.TT<- subset(animalT.T, trans == TRUE)
for (g in 1:nrow(animalT.TT)) {
  if (g == 1)
    pb<- txtProgressBar(min = 0, max = nrow(animalT.TT), style = 3)
  posART<- match(animalT.TT$et[g], animalR.sec$et)
  distM.g<- sqrt((animalT.TT$x[g] - animalR.sec$x[posART])^2 + (animalT.TT$y[g] - animalR.sec$y[posART])^2)
  detP.g<- detRngFun(distM.g)
  succ.g<- as.logical(rbinom(length(detP.g), 1, detP.g))
  if (sum(succ.g) > 0) {
    dtc.g<- data.frame(trns_et = animalT.TT$et[g], recv_id = posART, recv_x = animalR.sec$x[posART], recv_y = animalR.sec$y[posART], trns_x = animalT.TT$x[g], trns_y = animalT.TT$y[g], distance = distM.g)
    dtc<- rbind(dtc, dtc.g)
  }
  info<- sprintf("%d%% done", round(g/nrow(animalT.TT)*100))
  setTxtProgressBar(pb,g)
  if (g == nrow(animalT.TT))
    close(pb)
}

#forget about above function for a little while, and re-calcualte the distance between animalT and animalR at each second within the simulation

animalT.T$distance<- sqrt((animalT.T$x - animalR.sec$x)^2 + (animalT.T$y - animalR.sec$y)^2)

#this can be plotted

plot(animalT.T$et, animalT.T$distance, type = 'l')

#cut the data to 'interaction events' where the animals are within 200 m of each other -- and then label each interaction event
interac<- subset(animalT.T, animalT.T$distance<200)
inter<- split(interac$et, cumsum(c(1, diff(interac$et) != 1)))
inter2<- data.frame(ID = rep(names(inter), sapply(inter, length)), Obs = unlist(inter))
interac$inter<- inter2$ID

#bring in actual detected transmissions:
interac$transDET<- interac$et %in% dtc$trns_et

#some plots to look at the interactions:

pointplot<- subset(interac, trans == TRUE)
pointplot2<- subset(interac, transDET == TRUE)

par(mfrow=c(3,3))
plot(interac$et[interac$inter==2], interac$distance[interac$inter==2], xlab = "", ylab = "")
abline(h = min(interac$distance[interac$inter == 2]), lty = 2)
points(pointplot$et[pointplot$inter == 2], pointplot$distance[pointplot$inter == 2], col = "black", bg = "white", pch = 21, cex = 2)
points(pointplot2$et[pointplot2$inter == 2], pointplot2$distance[pointplot2$inter == 2], col = "black", bg = "cyan", pch = 21, cex = 2)
plot(interac$et[interac$inter==3], interac$distance[interac$inter==3], xlab = "", ylab = "")
abline(h = min(interac$distance[interac$inter == 3]), lty = 2)
points(pointplot$et[pointplot$inter == 3], pointplot$distance[pointplot$inter == 3], col = "black", bg = "white", pch = 21, cex = 2)
points(pointplot2$et[pointplot2$inter == 3], pointplot2$distance[pointplot2$inter == 3], col = "black", bg = "cyan", pch = 21, cex = 2)
plot(interac$et[interac$inter==4], interac$distance[interac$inter==4], xlab = "", ylab = "")
abline(h = min(interac$distance[interac$inter == 4]), lty = 2)
points(pointplot$et[pointplot$inter == 4], pointplot$distance[pointplot$inter == 4], col = "black", bg = "white", pch = 21, cex = 2)
points(pointplot2$et[pointplot2$inter == 4], pointplot2$distance[pointplot2$inter == 4], col = "black", bg = "cyan", pch = 21, cex = 2)
plot(interac$et[interac$inter==5], interac$distance[interac$inter==5], xlab = "", ylab = "")
abline(h = min(interac$distance[interac$inter == 5]), lty = 2)
points(pointplot$et[pointplot$inter == 5], pointplot$distance[pointplot$inter == 5], col = "black", bg = "white", pch = 21, cex = 2)
points(pointplot2$et[pointplot2$inter == 5], pointplot2$distance[pointplot2$inter == 5], col = "black", bg = "cyan", pch = 21, cex = 2)
plot(interac$et[interac$inter==6], interac$distance[interac$inter==6], xlab = "", ylab = "")
abline(h = min(interac$distance[interac$inter == 6]), lty = 2)
points(pointplot$et[pointplot$inter == 6], pointplot$distance[pointplot$inter == 6], col = "black", bg = "white", pch = 21, cex = 2)
points(pointplot2$et[pointplot2$inter == 6], pointplot2$distance[pointplot2$inter == 6], col = "black", bg = "cyan", pch = 21, cex = 2)
plot(interac$et[interac$inter==7], interac$distance[interac$inter==7], xlab = "", ylab = "")
abline(h = min(interac$distance[interac$inter == 7]), lty = 2)
points(pointplot$et[pointplot$inter == 7], pointplot$distance[pointplot$inter == 7], col = "black", bg = "white", pch = 21, cex = 2)
points(pointplot2$et[pointplot2$inter == 7], pointplot2$distance[pointplot2$inter == 7], col = "black", bg = "cyan", pch = 21, cex = 2)
plot(interac$et[interac$inter==8], interac$distance[interac$inter==8], xlab = "", ylab = "")
abline(h = min(interac$distance[interac$inter == 8]), lty = 2)
points(pointplot$et[pointplot$inter == 8], pointplot$distance[pointplot$inter == 8], col = "black", bg = "white", pch = 21, cex = 2)
points(pointplot2$et[pointplot2$inter == 8], pointplot2$distance[pointplot2$inter == 8], col = "black", bg = "cyan", pch = 21, cex = 2)
plot(interac$et[interac$inter==9], interac$distance[interac$inter==9], xlab = "", ylab = "")
abline(h = min(interac$distance[interac$inter == 9]), lty = 2)
points(pointplot$et[pointplot$inter == 9], pointplot$distance[pointplot$inter == 9], col = "black", bg = "white", pch = 21, cex = 2)
points(pointplot2$et[pointplot2$inter == 9], pointplot2$distance[pointplot2$inter == 9], col = "black", bg = "cyan", pch = 21, cex = 2)
plot(interac$et[interac$inter==10], interac$distance[interac$inter==10], xlab = "", ylab = "")
abline(h = min(interac$distance[interac$inter == 10]), lty = 2)
points(pointplot$et[pointplot$inter == 10], pointplot$distance[pointplot$inter == 10], col = "black", bg = "white", pch = 21, cex = 2)
points(pointplot2$et[pointplot2$inter == 10], pointplot2$distance[pointplot2$inter == 10], col = "black", bg = "cyan", pch = 21, cex = 2)

#summarise the data - total number of transmissions in each approach and the min distance in that approach
detDIS<- interac %>% group_by(inter) %>% summarise(n = sum(transDET), dis = min(distance))

#plot of the raw data: min distance vs. number of detections
plot(detDIS$dis, detDIS$n, cex = 2, xlab = "min distance from animalR", ylab = "number of detected transmissions")

#now bin the minimum distance of each approach by 10m intervals
detDIS$n<- as.numeric(detDIS$n)
detDIS$bins<- cut(detDIS$dis, breaks = seq(0, 200, by = 10))
detDIS$bins<- factor(detDIS$bins)
detDIS$binNUM<- detDIS$bins
levels(detDIS$binNUM)<- seq(0,190,10)

#this line counts how many times an animal was detected n number of times in each distance bin (for example: how many times was animalT detected twice when it's minimum approach distance was 100 m)
dect<- detDIS %>% group_by(binNUM, n) %>% summarise(detectionNUM = n())
#now count how many times animalT approached animalR in each distance bin
dect1<- detDIS %>% group_by(binNUM) %>% summarise(BinCount = n())
dect2<- merge(dect, dect1, by = "binNUM")
#equation to determine %detection probability for the number of detections in each distance bin (how many times was animalT detected twice in the 100 m distance bin divided by how many times did animalT come within 100 m of animalR?)
dect2$Prob<- dect2$detectionNUM/dect2$BinCount
dect2$n<- factor(dect2$n)

#plot the raw data - add smoothing lines (geom_smooth) as appropriate
ggplot(dect2, aes(x = binNUM, y = Prob, group = n))+
geom_line(aes(group = n, colour = n), lwd = 2)


#to make a 'mock' detection dataset for animalR but with known distance:
VMTdetect<- subset(interac, transDET == TRUE)
VMTdetect<- merge(detDIS, VMTdetect, by = 'inter')
VMTdetect$computertime<- rep(as.POSIXct(Sys.time()), length.out = length(VMTdetect$et))
VMTdetect$detectime<- VMTdetect$computertime + VMTdetect$et



############# Simplified simulation from supplimentary material ########################


#Simulate an animal swimming past a stationary reciever

#modified version of 'receiver_line_det_sim' from the glatos package to output number of detections on the receiver
#and the tragectory distance between the fish and the receiver (i.e. minimum distance between the two)
VMTDetSim <- function(vel=1,delayRng=c(80, 160),burstDur=5.0,
                      recSpc=1000,maxDist=1000,rngFun,outerLim=c(0,0),nsim=1000,showPlot=FALSE)
{ 
  #check if rngFun is function
  if(any(!is.function(rngFun))) 
    stop(paste0("Error: argument 'rngFun' must be a function...\n",
                "see ?receiverLineDetSim\n check: is.function(rngFun)"))
  
  
  #Define receiver line
  
  if(any(is.na(recSpc))) recSpc <- 0 #to simulate one receiver
  xLim <- c(0,sum(recSpc)+sum(outerLim))
  recLoc <- c(outerLim[1], outerLim[1] + cumsum(recSpc))
  yLim <-  c(-maxDist, maxDist)
  
  
  #Simulate tag transmissions
  
  nTrns <- floor((diff(yLim)/vel)/delayRng[1]) #number of transmissions 
  
  #sample delays, runif is a simulation command, will generate random numbers (#of numbers to generate, min value, max value)
  del <- matrix(runif(nTrns*nsim,delayRng[1],delayRng[2]),
                nrow=nsim, ncol=nTrns) 
  del <- del + burstDur #add burst duration (for Vemco)
  trans <- t(apply(del, 1, cumsum)) #time series of signal transmissions
  #"center" the fish track over the receiver line; with some randomness
  trans <- trans - matrix(runif(nsim, trans[,nTrns/2],trans[,(nTrns/2)+1]), 
                          nrow=nsim, ncol=nTrns)
  #row = simulated fish; col = signal transmission
  fsh.x <- matrix(runif(nsim, xLim[1], xLim[2]), nrow=nsim, ncol=nTrns)
  #convert from time to distance from start
  fsh.y <- matrix(trans*vel, nrow=nsim, ncol=nTrns) 
  
  
  #Optional quick and dirty plot just to see what is happening
  if(showPlot){
    plot(NA, xlim=xLim, ylim=yLim, asp=c(1,1),
         xlab="Distance (in meters) along receiver line",
         ylab="Distance (in meters) along fish path")
    #fish tracks and transmissions
    for(i in 1:nsim){
      lines(fsh.x[i,], fsh.y[i,], col="grey") #fish tracks
      points(fsh.x[i,], fsh.y[i,], pch=20, cex=0.8) #signal transmissions
    }
    #receiver locations
    points(recLoc, rep(0,length(recLoc)), pch=21, bg='red', cex=1.2)
    legend("topleft",legend=c("receiver","sim. fish path","tag transmit"),
           pch=c(21,124,20),col=c("black","grey","black"),pt.bg=c("red",NA,NA),
           pt.cex=c(1.2,1,0.8))
  }
  
  
  #Simulate detections 
  
  #calculate distances between transmissions and receivers
  for(i in 1:length(recLoc)){ #loop through receivers
    if(i == 1) { #pre-allocate objects, if first receiver
      succ <- detP <- distM <- vector("list",length(recLoc))
      nDets <- matrix(NA, nrow=nsim, ncol=length(recLoc)) #col = receiver
    }
    #tag-receiver distances in meters
    distM[[i]] <- sqrt((fsh.x - recLoc[i])^2 + (fsh.y)^2) 
    #detection probabilities
    detP[[i]] <- matrix(rngFun(distM[[i]]), nrow=nsim) 
    #detected=1, not=0
    succ[[i]] <- matrix(rbinom(length(detP[[i]]), 1, detP[[i]]), nrow=nsim) 
    #number of times each transmitter detected on ith receiver
    nDets[,i] <- rowSums(succ[[i]]) 
  }
  
  #max detects on any one receiver for each transmitter
  maxDet <- apply(nDets, 1, max) 
  #proportion of transmitters detected more than once on any receiver
  detProb <-  mean(maxDet>1) 
  
  #Edited ANB: for VMT return a dataframe of the tag distance from the receiver and num detections
  my_list<- list(maxDet, fsh.x[,1])
  
  #return(detProb) 
  return(my_list)
}


#function for the range, dm is distance in meters and b is intercept and slope for the logistic curve
pdrf<- function(dm, b=c(5.5, -0.0275)) {
  p<- 1/(1+exp(-(b[1]+b[2]*dm)))
  return (p)
}

#preview detection range curve. Edited to have detection probability be 50% at 200m (do 5.5/200)
plot(pdrf(0:800), ylab= "Probability of detecting each coded burst", xlab= "Distance between receiver and transmitter")

#set animal and study-specific details here - see ?receiver_line_det_sim for a detailed explanation
dp<- VMTDetSim(vel=0.68,delayRng=c(120, 240),burstDur=5.0, recSpc=0, maxDist=1000,
               rngFun = pdrf,outerLim=c(200,200),nsim=1000,showPlot=FALSE)

#clean up the output
Sim.vmtDat<- data.frame('detectionNum'= dp[[1]], 'distance'= dp[[2]])

#look at the data
plot(Sim.vmtDat$distance, Sim.vmtDat$detectionNum)

#sort the data into distance bins
Sim.vmtDat$bins <- cut(Sim.vmtDat$distance, breaks= seq(0, 1600, by = 10))
Sim.vmtDat$bins<- factor(Sim.vmtDat$bins)
Sim.vmtDat$binNUM<- Sim.vmtDat$bins
levels(Sim.vmtDat$binNUM)
levels(Sim.vmtDat$binNUM)<- seq(-200, 190, 10)

#Calculate the probability of getting n number of detections within each distance bin
x= vector()
dect1= data.frame()
for (i in levels(Sim.vmtDat[ , 4])){
  
  x<- Sim.vmtDat[Sim.vmtDat$binNUM == i, 1]
  dect1<- rbind(dect1, data.frame("X1dect" = length(which(x == 1))/length(x), 
                                  "X2dect" = length(which(x == 2))/length(x),
                                  "X3dect" = length(which(x == 3))/length(x),
                                  "X4dect" = length(which(x == 4))/length(x),
                                  "X0dect" = length(which(x == 0))/length(x),
                                  "X5dect" = length(which(x == 5))/length(x),
                                  "X6dect" = length(which(x == 6))/length(x),
                                  "X7dect" = length(which(x == 7))/length(x),
                                  "X8dect" = length(which(x == 8))/length(x)
  ))
  
}
#add the distance bin
dect1$distance<- levels(Sim.vmtDat$binNUM)
dect1$distance<- as.numeric(dect1$distance)


#plot the results with smoothing (make sure to check raw data!)
plot(dect1$distance, dect1$X0dect, ylim = c(0.02,0.7), col= "white", xlab= "Distance (m)", ylab= "Quantile of overall detections")
#axis(side= 1, at= seq(0,160,10) )
lines(smooth.spline(dect1$distance, dect1$X1dect), col= "blue", lwd= 5)
lines(smooth.spline(dect1$distance, dect1$X2dect), col= "red", lwd= 5)
lines(smooth.spline(dect1$distance, dect1$X3dect), col= "green", lwd= 5)
lines(smooth.spline(dect1$distance, dect1$X4dect), col= "purple", lwd= 5)
lines(smooth.spline(dect1$distance, dect1$X0dect), col= "black", lwd= 5)
lines(smooth.spline(dect1$distance, dect1$X5dect), col= "deeppink4", lwd= 5)
lines(smooth.spline(dect1$distance, dect1$X6dect), col= "sandybrown", lwd= 5)
lines(smooth.spline(dect1$distance, dect1$X7dect), col= "turquoise4", lwd= 5)
lines(smooth.spline(dect1$distance, dect1$X8dect), col= "coral1", lwd= 5)
legend("topleft", legend = c("0 detections", "1 detection", "2 detections", "3 detections", "4 detections", "5 detections"
),
col= c("black", "blue", "red", "green", "purple", "deeppink4"), lty = 1, ncol=3, lwd= 2)


legend("topleft", legend = c("0 detections", "1 detection", "2 detections", "3 detections", "4 detections", "5 detections",
                             "6 detections", "7 detections", "8 detections"),
       col= c("black", "blue", "red", "green", "purple", "deeppink4", "sandybrown", "turquoise4",
              "coral1"), lty = 1, ncol=3, lwd= 2)
