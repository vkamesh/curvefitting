library(TTR)
library(MASS)
library(NORMT3)
library(minpack.lm)

#DATA
temperature_table <- cbind(c(173,180,188,196,204,212,220,228,236,244,252,259,266,273,289,304,314,325,337,347,361,376,391,406,417,431,447,474,491,499,511,519,547,568,585,597,614,629,647,672,690,720,735,755,775,795,818,841,864,887,909,932,954,976,999,1021,1051,1077,1103,1129,1155,1177,1199,1220,1242,1263,1284,1306,1326,1349,1369,1390,1411,1433,1454,1474,1486,1499,1512,1531,1548,1570),c(800,790,780,770,760,750,740,730,720,710,700,690,680,670,660,650,640,630,620,610,600,590,580,570,560,550,540,530,520,510,500,490,480,470,460,450,440,430,420,410,400,390,380,370,360,350,340,330,320,310,300,290,280,270,260,250,240,230,220,210,200,190,180,170,160,150,140,130,120,110,100,90,80,70,60,50,40,30,20,10,0,-10))
temperature_table_inv <- cbind(c(800,790,780,770,760,750,740,730,720,710,700,690,680,670,660,650,640,630,620,610,600,590,580,570,560,550,540,530,520,510,500,490,480,470,460,450,440,430,420,410,400,390,380,370,360,350,340,330,320,310,300,290,280,270,260,250,240,230,220,210,200,190,180,170,160,150,140,130,120,110,100,90,80,70,60,50,40,30,20,10,0,-10),c(173,180,188,196,204,212,220,228,236,244,252,259,266,273,289,304,314,325,337,347,361,376,391,406,417,431,447,474,491,499,511,519,547,568,585,597,614,629,647,672,690,720,735,755,775,795,818,841,864,887,909,932,954,976,999,1021,1051,1077,1103,1129,1155,1177,1199,1220,1242,1263,1284,1306,1326,1349,1369,1390,1411,1433,1454,1474,1486,1499,1512,1531,1548,1570))
frequency_levels <- as.integer(c(25,50,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600)*1e6)
voltage_levels   <- c(0.850,0.875,0.900,0.925,0.950,0.950,0.950,0.975,1.000,1.025,1.075,1.125,1.175,1.225,1.250,1.275,1.325,1.350)

# FUNCTIONS

fconvertTemp <- function(temp) {
  
  apply(as.matrix(my_data[,2]), 1,
        function(x) {
          index <- which.min(abs(x-temperature_table[,1]))
          #kl <- temperature_table[index,2]/10+0.5
          up <- temperature_table[index,1]
          bo <- temperature_table[index-1,1]
          diff <- 
          
          kl <- temperature_table[index,2]/10+0.5
          return(kl)
  })

}

f_netwonstefan <- function(x,alpha,beta,gamma) {return((alpha*exp(-beta*x)/(1+gamma*alpha*(1-exp(-beta*x)))))}
f_netwonstefan_inv <- function(x,alpha,beta,gamma) {return(alpha-(alpha*exp(-beta*x)/(1+gamma*alpha*(1-exp(-beta*x)))))}
f_netwonstefan_sim <- function(x,alpha,beta,gamma) {return(alpha-(alpha*exp(-beta*x)/(1+gamma*(1-exp(-beta*x)))))}
fi_netwonstefan <- function(x,alpha,beta,gamma) {return(-(1/beta)*log(1-x/(alpha*(1+gamma*alpha-gamma*x))))}

f_newton <- function(x,alpha,beta,gamma) {return((beta - beta*(1-exp(-x/alpha)))+gamma)}


fit_newton_stefan_inv <- function(x,y, s_alpha=99, s_beta=0.01, s_gamma=1) {
  reg <- nls(y ~ f_netwonstefan_inv(x,alpha,beta,gamma),
             control=list(maxiter = 1000),
             start = list(alpha=s_alpha,
                          beta=s_beta,
                          gamma=s_gamma),
             trace=T,
             algorithm="port",
             lower=c(0,0,0),
             upper=c(1000,1,1000))
  
  my_data <- NULL
  my_data$coefs  <- c(coef(reg)[1],coef(reg)[2],coef(reg)[3])
  my_data$errors <- sum((y-f_netwonstefan(x,coef(reg)[1],coef(reg)[2],coef(reg)[3]))^2)
  return(my_data)
}

fit_newton_stefan <- function(x,y, s_alpha=99, s_beta=0.01, s_gamma=1) {
  reg <- nls(y ~ f_netwonstefan(x,alpha,beta,gamma),
             control=list(maxiter = 1000),
             start = list(alpha=s_alpha,
                          beta=s_beta,
                          gamma=s_gamma),
             trace=T,
             algorithm="port",
             lower=c(0,0,0),
             upper=c(1000,1,1000))
  
  my_data <- NULL
  my_data$coefs  <- c(coef(reg)[1],coef(reg)[2],coef(reg)[3])
  my_data$errors <- sum((y-f_netwonstefan(x,coef(reg)[1],coef(reg)[2],coef(reg)[3]))^2)
  return(my_data)
}



fit_stefan <- function(x,y) {
  reg <- nls(y ~ (gamma*exp(-beta*x)-alpha^4)^(1/4),
             control=list(maxiter = 1000),
             start = list(alpha=99,
                          beta=0.01,
                          gamma=1),
             trace=T,
             algorithm="port",
             lower=c(0,0,0),
             upper=c(1000,1,1000))
  
  my_data <- NULL
  my_data$coefs  <- c(coef(reg)[1],coef(reg)[2],coef(reg)[3])
  my_data$errors <- sum(((gamma*exp(-beta*x)-alpha^4)^(1/4))^2)
  return(my_data)
}



fit_netwon_inv <- function(x,y) {
  regression <- nls(y ~ beta*(1-exp(-((x^alpha)/gamma))),
                    #regression <- nls(y ~ beta1*(1-(beta3/t)^beta2),
                    control=list(maxiter = 10000),
                    start = list(beta=94,
                                 gamma=20,
                                 alpha=1),
                    trace=T,
                    algorithm="port",
                    lower=c(0,0,0),
                    upper=c(200,100,100))
  
  my_data <- c(coef(regression)[3],coef(regression)[1], coef(regression)[2])
  #colnames(my_data) <- c("beta","gamma","alpha")
  return(my_data)
}

fit_newton <- function(x,y, s_alpha=1, s_beta=2, s_gamma=10) {
  regression <- nls(y ~ (beta - beta*(1-exp(-(x/alpha))))+gamma,
                    #regression <- nls(y ~ beta1*(1-(beta3/t)^beta2),
                    control=list(maxiter = 10000),
                    start = list(alpha=s_alpha,
                                 beta=s_beta,
                                 gamma=s_gamma),
                    trace=T,
                    algorithm="port",
                    lower=c(0,0,0),
                    upper=c(1000,1000,1000))
  
  #colnames(my_data) <- c("beta","gamma","alpha")
  my_data <- NULL
  my_data$coefs <- c(coef(regression)[1], coef(regression)[2], coef(regression)[3])
  my_data$errors <- sum( (y - ((coef(regression)[2] - coef(regression)[2]*(1-exp(-(x/coef(regression)[1]))))+coef(regression)[3])    )^2)
  return(my_data)
}


extrapolate_exp <- function(c1, c2, c3, offset, x) {
  return((-c2*log(1-(x-offset)/c1))^(1/c3))
}

fEb <- function(f,V,naC,gamma,ccb,cck,beta){
  c1 <- ccb*naC/2
  c2 <- ccb*naC/2*gamma
  #return((c2*f*V^3+c1*f*V^2+ccb*Pstatic)/(f^beta-cck))
  return((c2*f*V^3+c1*f*V^2)/(f^beta-cck))
}


fptfunct <- function(x,alpha,beta,gamma)
  {return((-gamma*log(1-x/beta))^(1/alpha))}

fiptfunct  <- function(x,alpha,beta,gamma)
  {return(beta*(1-exp(-(x^alpha)/gamma)))}

fconvert <- function(x,gammaa,gammab,betaa,betab)
  {return(betaa*(1-(1-x/betab)^(gammab/gammaa)))}


fconvert2 <- function(x,gammaa,gammab,betaa,betab,alphaa,alphab)
{return( betaa* ( 1-(exp((-log((1-x/betab)^gammab))^(alphaa/alphab)))^(-1/gammaa)) )}

fconvert3 <- function(x,gammaa,gammab,betaa,betab,alphaa,alphab)
{return( betaa* ( 1-(exp((-gammab*log((1-x/betab)))^(alphaa/alphab)))^(-1/gammaa)) )}

fconvert4 <- function(x,gammaa,gammab,betaa,betab,alphaa,alphab)
{return( betaa* ( 1-(1/(1-x/betab))^(-gammab/gammaa)) )}


# naC   <- 5000
# beta  <- 1.7
# m1    <- 0.3304
# m2    <- 0.8077
# f     <- seq(0.1,1.6,0.0001)
# V     <- m1*f+m2
# ccb   <- 1000
# 
# cck   <- 0.07
# gamma <- 3.14
# 
# #my_Eb <- fEb(f,V,naC,gamma,ccb,cck,beta)
# #plot(x=f, y=my_Eb, col="red", type="l")
# 
# top <- ((1+gamma*V))*f*V^2
# bottom <- 1/(f^beta-cck)
# plot(x=f, y=bottom, type="l", col="red", log="x")
# 
# 
# 
# > cat("cck:",range(my_data[my_index,1]),"\n")
# cck: 0 0.71 
# 
# 
# my_data <- read.table("data.txt");
# colnames(my_data) <- c("cck","gamma", "beta", "convex")
# 
# my_index <- which(my_data[,"convex"] == 1)
# 
# cat("cck:",range(my_data[my_index,1]),"\n")
# cat("gamma:",range(my_data[my_index,2]),"\n")
# cat("beta:",range(my_data[my_index,4]),"\n")
# 
# plot3d(my_data[,"cck"],
#        my_data[,"gamma"],
#        my_data[,"beta"],
#        col=c("red","green")[my_data[,"convex"]+1],
#        size=3)
# 
# 
# 
# plot(x=f, y=1/(f^1.28-0.128), type="l", col="red")


fdEb <- function(f,V,naC,gamma,ccb,cck,beta,m1,m2){
  
  t1 <- (-(1+gamma*V)*V^2*f^beta*beta)/(f^beta-cck)^2
  t2 <- gamma*m1*f*V^2/(f^beta-cck)
  t3 <- ((1+gamma*V)*V^2)/(f^beta-cck)
  t4 <- (2*(1+gamma*V)*f*V*m1)/(f^beta-cck)

  return(t1+t2+t3+t4)
}
  
fEb_cross <- function(f,V,naC,gamma,Pstatic,ccb,cck,beta, m1, m2){
  c1 <- ccb*naC/2
  c2 <- ccb*naC/2*gamma
  #return((c2*f*V^3+c1*f*V^2+ccb*Pstatic)/(f^beta-cck))
  return(((c2*f*V^3+c1*f*V^2)/(f^beta-cck)-m2)/(m1*f))
}


fexp <- function(x, Ts, Istat, expo, shift) {
  return(Ts+exp(Istat*((x-shift)^expo)))
}

fiexp <- function(x, Ts, Istat, expo, shift) {
  return((log(x-Ts,base=exp(1))/Istat)^(1/expo)+shift)
}

fnewton <- function(t, Tambient, Tstart, beta1, beta2) {
  return(Tambient-(Tambient-Tstart)*exp(-beta1*(t)^beta2))
}

finewton <- function(t, Tambient, Tstart, beta1, beta2) {
  return(((-1/beta1)*log((t-Tambient)/(-(Tambient-Tstart)),base=exp(1)))^(1/beta2))
}

ftime <- function(f, cck, ccb, beta) {
  return(ccb/(f^beta-cck))
}

fpower <- function(f, V, naC, gamma, Pstatic) {
  #return(Pstatic + (1+V*gamma)*1/2*naC*f*V^2)
  return((1+V*gamma)*1/2*naC*f*V^2)
}

interpolate <- function(value, conversion_table) {
  table_source <- conversion_table[,1]
  table_destination <- conversion_table[,2]
  
  index <- which(table_source > value)
  top <- index[1]
  base <- top-1
  scale <- table_destination[top]-table_destination[base]
  return((value-table_source[base])/(table_source[top]-table_source[base])*scale+table_destination[base])
}

find_max_ccf<- function(a,b)
{
  d <- ccf(a, b, plot = FALSE, lag.max=500)
  cor = d$acf[,,1]
  lag = d$lag[,,1]
  res = data.frame(cor,lag)
  res_max = res[which.max(res$cor),]
  return(res_max)
}

myF <- function (x,t,T1,T2) {
  return(T1+(T2-T1)*x/t)
}

temperature_correction <- function(T, P, Tt) {
  A <- P/(T+337)
  return(P+(Tt-T)*A)
}

interpolate_single <- function(value, conversion_table) {
  table_source <- conversion_table[,1]
  table_destination <- conversion_table[,2]
  
  #cat(table_source[1] > table_source[2], table_source[1], table_source[2], "\n")
  if(table_source[1] > table_source[2]) {
    temp <- table_source[1]
    table_source[1] <- table_source[2]
    table_source[2] <- temp
    temp <- table_destination[1]
    table_destination[1] <- table_destination[2]
    table_destination[2] <- temp
  }
  
  scale <- table_destination[2]-table_destination[1]
  return((value-table_source[1])/(table_source[2]-table_source[1])*scale+table_destination[1])
}


translate <- function(value, table_source, table_destination){
  return(table_destination[which(table_source == value)])
}

get_voltage <- function(hertz) {
  return(apply(as.matrix(hertz),1,function(x) translate(x, frequency_levels, voltage_levels)))
}

# Returns every nth sample in data
resample <- function(data, n, middle=FALSE) {
  down_sample_count <- floor(length(data)/n)
  data_index <- seq(1,length(data),n)[1:down_sample_count]
  down_data <- data[data_index]
  
  if(middle)
    down_data <- down_data+(down_data[2]-down_data[1])/2
    
  return(down_data)
}


# Returns each average of consecutive groups of n samples
downsample <- function(data, n) {
  down_sample_count <- floor(length(data)/n)
  data_index <- seq(1,length(data),n)[1:down_sample_count]
  down_data <- apply(as.matrix(data_index),1,function(x) {median(data[x:(x+n-1)])})
  return(down_data)
}


get_max <- function(data) {
  return(which(data == max(data)));
}

# Returns the borders of a time series defined by a threshold
#   threshold: the threshold
#   data: the data to analyse, must be a ts!
get_borders <- function(threshold=1300, my_data, start=1.8) {
  
  if(class(my_data) != "data.frame") {
    cat("Error: the data is not a Data Frame!\n")
    silent();
  }
  
  time <- 1
  power <- 2
  borders <- c()
  
  borders$threshold <- threshold
  borders$start     <- start
  
  my_data <- my_data[which(my_data[,time] > start),]
  
  #threshold <- 1300
  #x <- seq(0,length(temp[,power]),10)
  #plot(x=temp[x,time], y=temp[x,power], type="l")
  #abline(h=threshold, col="blue")
  
  index <- which(my_data[,power] > threshold)
  #Must add 1 to get the first rising flank
  index <- c(1,index)
  my_borders <- which(diff(index) > 20)
  flanks_up   <- my_borders+1
  #flanks_down   <- borders-1
  #length(flanks_up)
  #length(flanks_down)
  
  borders$index <- array(index[flanks_up])
  
  flanks_up_time <- my_data[index[flanks_up],time]
  flanks_up_time <- flanks_up_time[which(diff(flanks_up_time) > 1.5)]
  
  #flanks_down_time <- temp[index[flanks_down],time]
  #flanks_down_time <- flanks_down_time[which(diff(flanks_down_time) > 0.5)]
  #abline(v=flanks_up_time, col="red")
  #abline(v=flanks_down_time, col="green")
  
  borders$time <- array(flanks_up_time)
  
  
  return(borders)
}







# my_linear <- function(x, x1, x2, y1, y2) {
#   m <- (y2-y1)/(x2-x1)
#   b <- y2-m*x2
#   return(m*x+b)
# }
# 
# my_table <- rbind(c(200,29.4,568.83,37.2,582.38),
#                   c(400,30.4,730.32,41.0,754.95),
#                   c(600,31.8,909.75,45.5,948.90),
#                   c(700,32.9,1024.60,47.8,1070.75),
#                   c(800,34.5,1304.0,53.4,1391.25))
# 
# 
# x <- -1000:1000
# plot(x=0, y=0, xlim=c(-400,60), ylim=c(0,1000), type="n")
# 
# for (i in 1:5) {
#   lines(x=my_table[i,c(2,4)], y=my_table[i,c(3,5)], type="p")
#   lines(x=x, y=my_linear(x,my_table[i,x1],my_table[i,x2],my_table[i,y1],my_table[i,y2]), type="l", lty="dashed")
# }
# abline(v=-273.15, col="blue")



P.a15 <- function(T,freq,core) {
  f.slope <- 0.220 - 0.315 * freq + 0.467 * freq^2
  f.offset <- f.slope/2.202
  a0 <- f.slope * core + f.offset
  a1 <- -56.652 * freq + 165.896 + (5-core) * 8.430
  a2 <- 33.105
  return(a0+exp((T-a1)/a2))
}

P.a7 <- function(T,freq,core) {
  f.slope <- 0.028 - 0.093 * freq + 0.371 * freq^2
  f.offset <- f.slope/2.202
  a0 <- f.slope * core + f.offset
  a1 <- -38.242 * freq + 187.668 + (5-core) * 8.430
  a2 <- 33.105
  return(a0+exp((T-a1)/a2))
}



