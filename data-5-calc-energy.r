source('energy.r')

sizes=c(6,8,10,12,14,16)
frequency_k=c(200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000)
frequency_G=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6)
frequency_M_s=c(200000000,300000000,400000000,500000000,600000000,700000000,800000000,900000000,1000000000,1100000000,1200000000,1300000000,1400000000,1500000000,1600000000)
#frequency_M_s=c(200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000)

#frequency_k=c(200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600)
#frequency_G=c(200000,300000,400000,500000,600000,700000,800000,900000,1000000,1100000,1200000,1300000,1400000,1500000,1600000)
#frequency_G=c(200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600)

#frequency_k=c(200e6,300e6,400e6,500e6,600e6,700e6,800e6,900e6,1000e6,1100e6,1200e6,1300e6,1400e6,1500e6,1600e6)

# FITTING Time over Frequency ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my_benchmark_time <- read.table("data_time.txt")
my_prediction_time <- my_benchmark_time
my_prediction_time_avg <- my_benchmark_time

my_starters_time <- rbind(c(1.7500,0.072,1.27),
                  c(7.7500,0.072,1.27),
                  c(28.000,0.072,1.27),
                  c(126.00,0.072,1.27),
                  c(598.00,0.072,1.27),
                  c(2601.0,0.072,1.27))

my_parameters_time <- matrix(0,length(sizes),4)
colnames(my_parameters_time) <- c("ccb","cck","beta","assym")

x <- frequency_k/1e6
y <- log(my_benchmark_time*1e6)

par(mfrow=c(1,1))
plot(1,1,
     type="n",
     xlim=c(0,max(x)),
     ylim=c(0,max(y)),
     main="Time - Frequency",
     xlab="frequency (GHz)",
     ylab="log(Task Completion Time) (us)")

abline(v=seq(0,max(x),0.1), h=0:30, col="gray", lty="dotted")

legend("bottomleft",
       title=expression(paste("input size (2"^"N",")")),
       paste("N=",sizes,sep=""),
       col=(1:8)+1,
       lty=1,
       pch=1:8,
       ncol=4,
       bg="white",
       cex=0.75)

for(size in 1:(length(sizes)))
#for(size in 1:6)
{
  if(size < 7)
  {
    cat("Power:",sizes[size],"\n")

    #regression <- nls(y[,size] ~ log(beta1 / (x - beta2) + beta3 ),
    regression <- nls(y[,size] ~ log((beta1 / (x - beta2))  + (beta1 * beta3) ),
                      control=list(maxiter = 1000),
                      start = list(beta1=my_starters_time[size,1],
                                   beta2=my_starters_time[size,2],
                                   beta3=my_starters_time[size,3]),
                      trace=T,
                      algorithm="port",
                      lower=c(0,0,-10000))

    my_parameters_time[size,"ccb"] <- coef(regression)[1]
    my_parameters_time[size,"cck"] <- coef(regression)[2]
    my_parameters_time[size,"beta"] <- coef(regression)[3]

    errors <- (y[,size] - log(my_parameters_time[size,"ccb"] / (x - my_parameters_time[size,"cck"]) + my_parameters_time[size,"beta"])^2)

    x_fit <- seq(0,1.7,0.01)
    lines(x=x_fit,
          y=log(my_parameters_time[size,"ccb"]/(x_fit-my_parameters_time[size,"cck"])+my_parameters_time[size,"beta"]),
          type="l",
          col="blue",
          lty="dotted")

    my_parameters_time[size,"assym"] <- (my_parameters_time[size,"cck"])^(1/my_parameters_time[size,"beta"])

    my_starters_time[size,1] <- my_parameters_time[size,1]
    my_starters_time[size,2] <- my_parameters_time[size,2]
    my_starters_time[size,3] <- my_parameters_time[size,3]

    my_prediction_time[,size] <- (my_parameters_time[size,"ccb"]/(x-my_parameters_time[size,"cck"])+my_parameters_time[size,"beta"])/1e6
  }

  lines(x=x,
        y=y[,size],
        col=size+1,
        type="o",
        pch=size)

  fitted <- log(my_parameters_time[size,1]/(x-my_parameters_time[size,3])+my_parameters_time[size,2])
}

my_parameters_time_avg <- apply(my_parameters_time,2,mean)
for(size in 1:6) {
  my_prediction_time_avg[,size] <- (my_parameters_time[size,"ccb"] / (x-my_parameters_time_avg["cck"]) + my_parameters_time_avg["beta"])/1e6
}
print(((my_prediction_time-my_benchmark_time)/my_benchmark_time)*100)

cat("***********KAMESH***********:",sizes[size],"\n")
print(summary(regression))
print(my_parameters_time)


# POWER ~~~~~~~~~~~~~~~~~~~~

my_power_estimates <- read.table("data_power.txt",sep="")
colnames(my_power_estimates) <- c("freq","size","median")

my_benchmark_power <- my_benchmark_time
my_benchmark_power
for(i in 1:6) {
  my_benchmark_power[,i] <- my_power_estimates[(0:(length(frequency_M_s)-1))*6+i,"median"]
} 
my_prediction_power <- my_prediction_time

# FITTING Power over Frequency ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par(mfrow=c(1,1))

#Start Matrix
my_starters_power <- matrix(0,length(sizes), 5)
colnames(my_starters_power) <- c("freq","power","Pstatic", "Ileak", "Dyn")

#Matrix to hold the fitted parameters
my_parameters_power <- matrix(0,length(sizes),4)
colnames(my_parameters_power) <- c("power","Pstatic", "gamma", "naC")

#Matrix to hold the fitted parameters
my_parameters_power2 <- matrix(0,length(sizes),4)
colnames(my_parameters_power2) <- c("power","Istatic", "Ileak", "naC")

#Matrix to hold the fitted parameters
compute_freq <- matrix(0,length(sizes),4)
colnames(compute_freq) <- c("power","Ptotal","Pstatic","Rt")

errors <- c()
Rt <- c()
total_power <- c()

cat("PLOT TEMPERATURE TRACES ~~~~~~~~~\n")
#for(freq in 1:length(frequencies_short))
#for(freq in 1:11)
#{

for(size in 1:6)
#for(size in 3)
{    
  # INITIALIZATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #cat("\tinitializing ...\n")

  cat("Array size:",sizes[size],"\n")
  
  my_starters_power[size,"power"]    <- sizes[size]
  my_parameters_power[size,"power"]  <- sizes[size]
  my_parameters_power2[size,"power"]  <- sizes[size]
  
  f <- frequency_k[c(-14)]/1e6
  V <- get_voltage(f*1e9)
  P <- my_benchmark_power[c(-14),size]/1000
  
# P_stat <- 0.4670914*1000
#   
# compute_gamma <- function(P1,P2,f1,f2) {
#   m1 <- 0.3304
#   m2 <- 0.8077
#   V1 <- m1*f1+m2
#   V2 <- m1*f2+m2
#   gamma <- (P2*((V1)^2)*f1-P1*((V2)^2)*f2)/(P1*f2*((V2)^3)-P2*f1*((V1)^3))
#   return(gamma)
# }
#   
# my_arr <- c()
# for(k in 1:14) {
#   for(l in 1:6) {
#     my_arr <- c(my_arr,compute_gamma(P1=my_benchmark_power[k,l]-P_stat,
#                                      P2=my_benchmark_power[(k+1),l]-P_stat,
#                                      f1=frequency_G[k],
#                                      f2=frequency_G[k+1]))
#   }
#   cat(my_benchmark_power[k,l],my_benchmark_power[(k+1),l],frequency_G[k],frequency_G[(k+1)],"\n")
# }
# 
# hist(my_arr,breaks=seq(-400,100,0.5),xlim=c(-20,20))

# compute_aC <- function(P1,P2,f1,f2,gamma) {
#   m1 <- 0.3304
#   m2 <- 0.8077
#   V1 <- m1*f1+m2
#   V2 <- m1*f2+m2
#   aC <- (V1*f1-V2*f2^2+gamma*V1^2*f1-gamma*V2^2*f2)/(P1-P2)
#   return(aC)
# }
  
  
  
  # PLOTTING figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("\tplotting figures...\n")
  
  plot(x=f,
       y=P,
       col="cyan3",
       main=paste("buffer-size:",sizes[size],"at 37gC"),
       xlab="frequency (Hz)",
       ylab="power (W)",
       type="n",
       lty="dashed",
       xlim=range(frequency_k)/1e6,
       ylim=c(0,max(my_benchmark_power[,size])/1000)
       #ylim=c(0.85,1.075)
       )
  grid()
  
  # FITTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  cat("\tfitting default ...\n")
  
  r2egression <- nls(P ~ Pstat + (1+V*gamma)*1/2*naC*f*V^2,
  #r2egression <- nls(P ~ Pstat + 1/2*naC*f*V^2 + gamma * V,
                     control=list(maxiter = 1000,
                                  #minFactor=1e-25,
                                  warnOnly=TRUE),
                     start = list(Pstat = 0.5,
                                  gamma  = 1.5,
                                  naC = 0.7),
                     trace = T,
                     algorithm = "port",
                     lower = c(0,0,0),
                     upper = c(3,2000,2))  
  
  my_parameters_power[size, "Pstatic"] <- Pstatic <- coef(r2egression)[1]
  my_parameters_power[size, "gamma"]   <- gamma   <- coef(r2egression)[2]
  my_parameters_power[size, "naC"]     <- naC     <- coef(r2egression)[3]
  
  cat("\tfitting alternative ...\n")
  
  r3egression <- nls(P ~ 3.69*Istatic + Ileak * V + 1/2*naC*f*V^2,
                     control=list(maxiter = 1000,
                                  #minFactor=1e-25,
                                  warnOnly=TRUE),
                     start = list(Istatic = 0.04,
                                  Ileak  = 0.9,
                                  naC = 0.1),
                     trace = T,
                     algorithm = "port",
                     lower = c(-1,-10,-1),
                     upper = c(1,1,10))  
  
  my_parameters_power2[size, "Istatic"] <- Istatic  <- coef(r3egression)[1]
  my_parameters_power2[size, "Ileak"]   <- Ileak    <- coef(r3egression)[2]
  my_parameters_power2[size, "naC"]     <- naC2     <- coef(r3egression)[3]
  
  # PLOTTING figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cat("\tplotting fitted curve...\n")
  
  error2s <- log((P-(Pstatic+(1+V*gamma)*1/2*naC*f*V^2))^2)
  error3s <- log((P-(3.69*Istatic + Ileak * V + 1/2*naC2*f*V^2))^2)
  
  print(round(error2s,3))
  print(round(error3s,3))

  Vfit <- get_voltage(frequency_k*1e3)
  freq_fit <- frequency_k/1e6
  P <- my_benchmark_power[,size]/1000
  f <- frequency_k/1e6
  
  lines(x=freq_fit,
        y=(Pstatic+(1+Vfit*gamma)*1/2*naC*freq_fit*Vfit^2),
  #      y=Pstatic+1/2*naC*freq_fit*Vfit^2+Vfit*gamma,
        col="red",
        type="o",
        pch=3)
  
  lines(x=freq_fit,
        y=(3.69*Istatic + Vfit * Ileak + 1/2*naC2*freq_fit*Vfit^2),
        #      y=Pstatic+1/2*naC*freq_fit*Vfit^2+Vfit*gamma,
        type="o",
        col="blue",
        pch=2)
  
  #cat("Rt:")
  #print(round(Vfit*gamma,2))  
  #colnames(compute_freq) <- c("power","Pstatic","Pdynamic","Pleakage","rPstatic","rPdynamic","rPleakage")
  
  Rt       <- rbind(Rt,round(Vfit*gamma,2))
  Pleakage <- (1/2*naC*freq_fit*Vfit^2)*Vfit*gamma + Pstatic
  Pdynamic <- (1/2*naC*freq_fit*Vfit^2) + Pleakage
  
  lines(x=f,
        y=P,
        col="cyan3",
        type="l",
        lty="dashed")
  
  #lines(x=freq_fit,
  #      y=rep(Pstatic,length(freq_fit)),
  #      col=3,
  #      type="l")
  #lines(x=freq_fit,
  #      y=Pleakage,
  #      col=4,
  #      type="l")
  
  legend("topleft",
         c("measured power", "fitted power", "static power","leakage power"),
         col=c("cyan3","red","green","blue"),
         lty=c("dashed",rep("solid",3)))
  
#   Pleakage <- (1/2*naC*freq_fit*Vfit^2)*Vfit*gamma
#   Pdynamic <- (1/2*naC*freq_fit*Vfit^2)
#   Ptotal <- Pdynamic + Pleakage + Pstatic
#   
#   compute_freq[table_index,"power"]     <- table_index
#   compute_freq[table_index,"Pstatic"]   <- Pstatic
#   compute_freq[table_index,"Rt"]        <- gamma
  
  f <- frequency_k/1e6
  V <- get_voltage(f*1e9)
  my_prediction_power[,size] <- (Pstatic + (1+V*gamma) * 1/2*naC*f*V^2)*1000
}

# my_parameters_time_avg <- apply(my_parameters_time,2,mean)
# for(size in 1:6) {
#   my_prediction_time_avg[,size] <- (my_parameters_time[size,"ccb"] / (x-my_parameters_time_avg["cck"]) + my_parameters_time_avg["beta"])/1e6
# }

cat("***********KAMESH***********:",sizes[size],"\n")
print(summary(r2egression))
print(my_parameters_power)

# ENERGY modeling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  meanPstatic <- mean(my_parameters_power[,"Pstatic"])
  
  weigths_os <- 1-(0.115/frequency_G)

  my_prediction_time_avg_os <- my_prediction_time#*weigths_os
  my_benchmark_time_os <- my_benchmark_time#*weigths_os

  my_theoretical_energy  <- (my_prediction_power/1000-meanPstatic)*my_prediction_time_avg_os
  my_measured_energy     <- (my_benchmark_power/1000-meanPstatic)*my_benchmark_time_os

  my_theoretical_energy_norm  <- my_theoretical_energy
  my_measured_energy_norm     <- my_measured_energy

  for(i in 1:6) {
    my_theoretical_energy_norm[,i] <- my_theoretical_energy_norm[,i]/(2^sizes[i])*1e9
    my_measured_energy_norm[,i]         <- my_measured_energy[,i]/(2^sizes[i])*1e9
  }

  plot(x=1,
       y=1,
       type="n",
       ylim=range(my_measured_energy_norm[,1:6],my_theoretical_energy_norm[,1:6]),
       xlim=range(frequency_G))

  for(i in 1:6) {
    lines(x=frequency_G,
          y=my_theoretical_energy_norm[,i],
          col="black",
          type="o",
          lty="dashed",
          pch=i)
    
    lines(x=frequency_G,
          y=my_measured_energy_norm[,i],
          col=i+1,
          type="o",
          lty="solid",
          pch=i)
    
  }

legend("top",
       paste(sizes),
       col=2:7,
       lty="solid",
       pch=1:6,
       ncol=2)


abs(round(apply((((my_prediction_time[,1:6]-my_benchmark_time[,1:6]))/my_benchmark_time[,1:6])*100, 1, mean),2))
abs(round(apply((((my_prediction_power[,1:6]-my_benchmark_power[,1:6]))/my_benchmark_power[,1:6])*100, 1, mean),2))
abs(round(apply((((my_theoretical_energy_norm[,1:6]-my_measured_energy_norm[,1:6]))/my_measured_energy_norm[,1:6])*100, 1, mean),2))

