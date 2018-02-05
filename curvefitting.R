#### Set working directory ####
#setwd("~/Documents/Telecom-Paristech/paper2/curvefitting/kamesh")

#### Initialize parameters ####

benchmarks = c("Gold-rader","Blowfish","SHA")
frequency_M = c(100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500)
frequency_G = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5)
Voltage = c(0.95,0.97,0.99,1.01,1.03,1.05,1.07,1.09,1.11,1.13,1.15,1.17,1.19,1.21,1.23)

#### Fitting execution time over frequency ####

# Read data from the file and assign coloumn names
benchmark_time <- read.table("data_time_kamesh.txt")
colnames(benchmark_time) <- c("Gold-rader","Blowfish","SHA")

my_benchmark_time <- benchmark_time

prediction_time <- benchmark_time

prediction_time_avg <- benchmark_time

starting_values <- rbind(c(3.0,0.072,1.27),
                         c(10.0,0.072,1.27),
                         c(10.0,0.072,1.27))

# Build a matrix for the constants and initilize them to '0'
# and assign coloumn names to the matrix. assym = asymptote
constant_parameters <- matrix(0,length(benchmarks),4)
colnames(constant_parameters) <- c("ccb","cck","beta","assym")

x <- frequency_M
y <- benchmark_time

# Plot the graph (without data curves)
par(mfrow=c(1,1))
plot(1,1,
     type="n",
     xaxp  = c(0, 1600, 16),
     xlim=c(0,max(x)),
     ylim=c(0,max(y)),
     main="Execution time vs Frequency",
     xlab="Frequency (MHz)",
     ylab="Execution time (s)")

#axis(side=1, at=seq(0, 1600, by=100))
abline(v=seq(0,max(x),100), h=seq(0,max(y),100), col="gray", lty="dotted")

legend("topright",
       #title=expression(paste("Benchmarks")),
       paste("",benchmarks,sep=""),
       col=(1:8)+1,
       lty=1,
       pch=1:8,
       ncol=1,
       bg="white",
       cex=0.75)

# Run the regression analysis (fitting) for each benchmark
for(size in 1:length(benchmarks)) {
  if(size < (length(benchmarks)+1)) {
    
    cat("Benchmark:",benchmarks[size],"\n")
    
    # regression analysis
    regression <- nls(y[,size] ~ (constant1 * ((1 / (x - constant2))  + constant3)),
                         control=list(maxiter = 1000),
                         start = list(constant1=starting_values[size,1],
                                      constant2=starting_values[size,2],
                                      constant3=starting_values[size,3]),
                         trace = T,
                         algorithm = "port",
                         lower = c(0,0,-10000))

    constant_parameters[size,"ccb"] <- coef(regression)[1]
    constant_parameters[size,"cck"] <- coef(regression)[2]
    constant_parameters[size,"beta"] <- coef(regression)[3]
    
    errors <- (y[,size] - (constant_parameters[size,"ccb"] * ((1 / (x - constant_parameters[size,"cck"])) + constant_parameters[size,"beta"])))
    
    # Fit the predicted data
    x_fit <- seq(100,1500,10)
    lines(x=x_fit,
          y=(constant_parameters[size,"ccb"] * ((1 / (x_fit - constant_parameters[size,"cck"])) + constant_parameters[size,"beta"])),
          type="l",
          col="blue",
          lty="dotted")
    
    # Vertical asymptote is at infinity
    constant_parameters[size,"assym"] <- (constant_parameters[size,"cck"])^(1/constant_parameters[size,"beta"])
     
    starting_values[size,1] <- constant_parameters[size,1]
    starting_values[size,2] <- constant_parameters[size,2]
    starting_values[size,3] <- constant_parameters[size,3]

    prediction_time[,size] <- (constant_parameters[size,"ccb"] / (x-constant_parameters[size,"cck"]) + constant_parameters[size,"beta"])
    
  }
  
  # Plot the data on the graph
  lines(x=x,
        y=y[,size],
        col=size+1,
        type="o",
        pch=size)
  
  fitted <- (constant_parameters[size,1] / (x-constant_parameters[size,3]) + constant_parameters[size,2])
  
}

# Calculating average prediction time
constant_parameters_avg <- apply(constant_parameters,2,mean)
for(size in 1:length(benchmarks)) {
  prediction_time_avg[,size] <- (constant_parameters[size,"ccb"] / (x - constant_parameters_avg["cck"]) + constant_parameters_avg["beta"])
  }

# Print percentage increase/decrease in the prediction times
print(((prediction_time - benchmark_time) / benchmark_time) * 100)


#### Fitting power over frequency ####

my_power_estimates <- read.table("av-power_vs_frequency.txt",sep="")
colnames(my_power_estimates) <- c("Gold-rader","Blowfish","SHA")

my_benchmark_power <- my_power_estimates

prediction_power <- my_power_estimates

par(mfrow=c(1,1))

#Matrix to hold the fitted parameters
my_parameters_power <- matrix(0,length(benchmarks),4)
colnames(my_parameters_power) <- c("power","Pstatic", "gamma", "naC")

#Matrix to hold the fitted parameters
compute_freq <- matrix(0,length(benchmarks),4)
colnames(compute_freq) <- c("power","Ptotal","Pstatic","Rt")

errors <- c()
Rt <- c()
total_power <- c()

for (size in 1:length(benchmarks)) {
  
  cat("Benchmark:",benchmarks[size],"\n")
  my_parameters_power[size,"power"]  <- benchmarks[size]
  
  f <- frequency_M
  V <- Voltage
  P <- my_benchmark_power[,size]
  
  cat("\tplotting figures...\n")
  
  plot(x=f,
       y=P,
       col="red",
       main=paste("Average power vs Frequency:", benchmarks[size]),
       xlab="Frequency (MHz)",
       ylab="Average power (mW)",
       type="n",
       lty="dashed",
       xlim=range(frequency_M),
       ylim=c(0,max(my_benchmark_power[,size]+5))
       #ylim=c(0,60)
  )
  
  abline(v=(seq(0,1600,100)), h=(seq(0,60,5)), col="lightgray", lty="dotted")
 
  # FITTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  cat("\tfitting default ...\n")
  
  regression_power <- nls(P ~ Pstatic + (1+V*gamma)*1/2*naC*f*V^2,
                     #r2egression <- nls(P ~ Pstat + 1/2*naC*f*V^2 + gamma * V,
                     control=list(maxiter = 1000,
                                  #minFactor=1e-5,
                                  warnOnly=T),
                     start = list(Pstatic = 0.1,
                                  gamma = 0.1,
                                  naC   = 0.1),
                     trace = T,
                     algorithm = "port",
                     lower = c(0,0,0),
                     upper = c(2,30,2))
  
  my_parameters_power[size, "Pstatic"] <- Pstatic <- coef(regression_power)[1]
  my_parameters_power[size, "gamma"]   <- gamma   <- coef(regression_power)[2]
  my_parameters_power[size, "naC"]     <- naC     <- coef(regression_power)[3]
  
  print(my_parameters_power)

  Vfit <- Voltage
  freq_fit <- frequency_M
  P <- my_benchmark_power[,size]
  f <- frequency_M

  lines(x=f,
        y=P,
        col="blue",
        type="l")

  lines(x=freq_fit,
        y=(Pstatic+(1+Vfit*gamma)*1/2*naC*freq_fit*Vfit^2),
        #y=Pstatic+1/2*naC*freq_fit*Vfit^2+Vfit*gamma,
        col="red",
        type="o",
        pch=3,
        lty="dashed")
  
  prediction_power[,size] <- (Pstatic+(1+Vfit*gamma)*1/2*naC*freq_fit*Vfit^2)

}

#### Energy modeling and fitting ####

experimental_energy <- my_benchmark_power * my_benchmark_time

theoretical_energy <- prediction_power * prediction_time

plot(x=1,
     y=1,
     main=paste("Total energy consumption vs Frequency"),
     xaxp  = c(0, 1600, 16),
     yaxp  = c(0, 2000, 10),
     xlab= "Frequency (MHz)",
     ylab= "Total energy consumption (mJ)",
     type="n",
     ylim=range(experimental_energy[,1:3],theoretical_energy[,1:3]),
     xlim=range(frequency_M))

abline(v=(seq(0,1600,100)), h=(seq(0,2000,200)), col="lightgray", lty="dotted")

for (size in 1:length(benchmarks)) {
  lines(x=frequency_M,
        y=theoretical_energy[,size],
        col="black",
        type="o",
        lty="dashed",
        pch=size)
  
  lines(x=frequency_M,
        y=experimental_energy[,size],
        col=size+1,
        type="o",
        lty="solid",
        pch=size)

}

#seq.int(0.95, 1.25, length.out = 10)


